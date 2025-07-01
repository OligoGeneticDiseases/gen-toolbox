import datetime
import gc
import os
import traceback
from random import random

import hail as hl
from pathlib import Path

import math
import pandas as pd
from hail.utils import info

from src.utils.file.file_meta import get_metadata
from src.utils.general.data_manipulation import batcher
from src.data_processing.file_io.readers import find_filetype
from src.data_processing.file_io.writers import write_filelist, trim_prefix
from src.data_processing.vcf.read import import_and_annotate_vcf_batch, load_db_batch, load_mt, import_and_annotate_vcf
from src.data_processing.vcf.transform import reduce_to_2d_table, create_frequency_bins, extract_genotype_dataframe
from src.data_processing.hail.genomic_operations import merge_matrix_tables_rows, merge_matrix_tables_cols, \
    create_related_samples_table, get_variants_counts_per_gene
from src.data_processing.pca.analysis import pca_graphing


N_BATCH = 10

def set_reference_file():
    project_root = Path(__file__).resolve().parents[2]
    reference_path = project_root / "GRCh37_MT.json"

    hl.default_reference = hl.ReferenceGenome.read(str(reference_path))

def handle_quit():
    # TODO: log some extra output for Nextlow / ensure correct error codes are sent
    hl.stop()


def find_elements(dictionary, x, pos=0, sep=None):
    """
    Find elements matching a dict with structure key: tuple(list(), list()) structure.
    :param dictionary: The dict to be searched.
    :param x: Lookup value
    :param pos: Position of list that is searched (i.e. tuple[0] or tuple[1]. Default=0.
    :param sep: Default None. Split elements according to this separator.
    :return: Matching key: item pairs where x in items() sublist.
    """
    result = dict()
    for key, item in dictionary.items():
        it = item[pos]
        if sep is not None:
            it = item[pos].split(sep)
        if x in it:
            result[key] = item
    return result


def write_frequency_table(result, args, name, extra_tag=""):
    table_path = hl.utils.timestamp_path(
        os.path.join(args.dest, f"{name}_{extra_tag}.csv")
    )
    if not Path(table_path).exists() or args.overwrite:
        result.to_csv(table_path)
        hl.utils.info(f"OLIGO: wrote output to {table_path}")
    else:
        hl.utils.info(f"OLIGO: skipped writing output:\n{result}.")


def create_null_dist_match_anti(vcfs, commandhandler, sample_size=None):
    """
    Create a match set of equal size to phenotype samples by randomly sampling from non-phenotype samples,
    and use the rest (including phenotype samples) as anti-match.

    Args:
        vcfs (List[Path]): List of VCF paths
        commandhandler: Object with .metadata and .args.phenotype
        sample_size (int, optional): Number of random controls to sample as match.
                                     If None, defaults to number of phenotype cases.

    Returns:
        List[List[Path]]: [matched_vcfs, anti_matched_vcfs]
    """
    assert hasattr(commandhandler, "metadata"), "Missing metadata in commandhandler"
    assert hasattr(commandhandler.args, "phenotype"), "Missing phenotype argument"

    phenotype = commandhandler.args.phenotype
    sep = ","
    pos = 0

    # Identify phenotype and non-phenotype keys
    keys_with_pheno = [k for k, v in commandhandler.metadata.items()
                       if len(v) > pos and phenotype in v[pos].split(sep)]
    keys_without_pheno = [k for k in commandhandler.metadata if k not in keys_with_pheno]

    if not keys_with_pheno:
        raise ValueError(f"No metadata entries found with phenotype key: {phenotype}")
    if not keys_without_pheno:
        raise ValueError("No metadata entries found without phenotype")

    if sample_size is None:
        sample_size = len(keys_with_pheno)

    if len(keys_without_pheno) < sample_size:
        raise ValueError(f"Not enough controls to sample: requested {sample_size}, available {len(keys_without_pheno)}")

    sampled_keys = set(random.sample(keys_without_pheno, sample_size))

    def trim(vcf):
        return trim_prefix(vcf.stem)

    matched_vcfs = [vcf for vcf in vcfs if trim(vcf) in sampled_keys]
    all_other_keys = set(keys_without_pheno + keys_with_pheno) - sampled_keys
    anti_matched_vcfs = [vcf for vcf in vcfs if trim(vcf) in all_other_keys]

    if not matched_vcfs or not anti_matched_vcfs:
        raise ValueError("One of the returned VCF groups is empty.")

    return [matched_vcfs, anti_matched_vcfs]



def collect_vcfs_by_phenotype(commandhandler, allow_null_prefix=False, symlink=True):
    """
    Load VCF paths and split them by phenotype (cases and controls).
    Creates symlinks to temp folders for case/control.

    Returns:
        List[List[Path]]: [case_vcfs, control_vcfs]
    """
    set_reference_file()
    # --- Check for required arguments ---
    assert hasattr(commandhandler, "args"), "commandhandler is missing .args"
    assert getattr(commandhandler.args, "file", None), "commandhandler.args.file is required"
    assert getattr(commandhandler.args, "phenotype", None), "commandhandler.args.phenotype is required"
    assert getattr(commandhandler.args, "globals", None), "commandhandler.args.globals is required"
    assert getattr(commandhandler.args, "temp", None), "commandhandler.args.temp is required"

    # Gather all VCFs
    full_paths = [Path(p) for p in commandhandler.args.file]
    vcfs = []
    for path in full_paths:
        if path.is_file() and path.suffix == ".vcf":
            vcfs.append(path)
        else:
            vcfs.extend(path.rglob("*.vcf"))

    metadata = get_metadata(commandhandler.args.globals)
    assert metadata is not None

    match_and_anti_match = []

    if commandhandler.args.phenotype and not commandhandler.args.phenotype.startswith("*"):
        matches = find_elements(metadata, commandhandler.args.phenotype, sep=",")
        match_keys = set(matches.keys())

        filtered_vcfs = [path for path in vcfs if trim_prefix(path.stem) in match_keys]
        assert filtered_vcfs, f"No matches found with phenotype {commandhandler.args.phenotype}, quitting."
        anti_match = list(set(vcfs) - set(filtered_vcfs))

        match_and_anti_match.append(filtered_vcfs)
        match_and_anti_match.append(anti_match)

        info(
            f"Found {len(filtered_vcfs)} matches, {len(anti_match)} anti-matches for phenotype: {commandhandler.args.phenotype}")
        if symlink:
            for label, group in zip(['positive', 'negative'], [filtered_vcfs, anti_match]):
                symlink_root = Path(f"{commandhandler.args.temp}/{label}_{commandhandler.args.phenotype}_{len(group)}")
                for link in group:
                    symlink_path = symlink_root / link.name
                    symlink_path.parent.mkdir(parents=True, exist_ok=True)
                    if not symlink_path.exists():
                        symlink_path.symlink_to(link)

    elif allow_null_prefix and commandhandler.args.phenotype and commandhandler.args.phenotype.startswith("*"):
        match_and_anti_match.append(create_null_dist_match_anti(vcfs, commandhandler))  # returns list of [cases, controls]

    else:
        raise ValueError("Phenotype must be provided or start with '*' to create null distribution.")

    return match_and_anti_match


class CommandHandler:

    def __init__(self, args):
        self.args = args
        self.metadata = None

    def handle_find_type_command(self):
        files = find_filetype(self.args.source, self.args.type, verbose=False)
        write_filelist(
            self.args.directory,
            "{0}.{1}.txt".format(
                os.path.basename(os.path.normpath(self.args.source)), self.args.type
            ),
            files,
            regex=self.args.regex,
        )

    def handle_read_vcfs_command(self):
        """
        Handle read_vcfs command with the command input. This will create 2 output tables (create_frequency_bins)
        for a match and anti_match set of a given phenotype (matches from globals file).
        Creates batches of VCF files so that Hail would not crash.
        """
        match_and_anti_match = collect_vcfs_by_phenotype(self,
            allow_null_prefix=True
        )

        # Create frequency tables for both match and anti-match
        for k, positive_or_negative_set in enumerate(match_and_anti_match):
            n_batches = int(math.ceil(len(positive_or_negative_set) / N_BATCH))
            batch_matrix_tables = []
            set_tag = "positive" if k == 0 else "negative"
            set_tag = f"{self.args.phenotype}_" + set_tag
            info(f"Starting with set: {set_tag} ({len(positive_or_negative_set)}). "
                 f"Dividing work into batches of {N_BATCH}, with {n_batches} batches.")
            if k == 0:
                info(f"Skipping {set_tag} due to debugging tag active.")
                stack = traceback.extract_stack()
                # Reports this location in the logfile.
                caller = stack[-1]
                info(f"File: {caller.filename}, Function: {caller.name}, Line: {caller.lineno}")
                continue #uncommenting this will skip reading positive set


            for i, batch in enumerate(batcher(positive_or_negative_set, N_BATCH)):
                combined_mt_path = Path(self.args.dest).joinpath(
                    "multi_batch_{0}_{1}_dataset_{2}_{3}.mt".format(
                        i+1, len(batch), len(positive_or_negative_set), set_tag
                    )).as_posix()
                if not Path(combined_mt_path).exists():
                    matrix_tables = import_and_annotate_vcf_batch(
                        batch, metadata=self.metadata, annotate=self.args.annotate, location=self.args.dest,
                        interval=self.args.interval
                    )  # hail.import_vcf() and hail.VEP() within this function
                    batch_combined_mt = merge_matrix_tables_rows(
                        matrix_tables, self.args.phenotype
                    )  # Merge batch of matrix tables into one
                    if self.args.write:
                        batch_combined_mt.write(combined_mt_path)
                        # Write MatrixTable to disk TODO: log values don't make sense after a certain number of batches i.e. batch 9/8
                        info(
                            "Wrote batch {0}/{1} onto disk with {2} subelements.".format(
                                i + 1, n_batches, len(batch)
                            )
                        )
                        info(f"reading{combined_mt_path}")
                        #batch_combined_mt = hail.read_matrix_table(combined_mt_path)
                        batch_combined_mt.unpersist()
                else:
                    info(
                        "Skipping batch {0}/{1}. Read from {2} disk with {3} subelements.".format(
                            i + 1, n_batches, combined_mt_path, len(batch)
                        ))
                    #batch_combined_mt = hail.read_matrix_table(combined_mt_path)
                #  None in case of small batch with phenotype negative, another batch might contain results
                if Path(combined_mt_path).exists() is not None:
                    batch_matrix_tables.append(combined_mt_path)
            if len(batch_matrix_tables) > int(N_BATCH / 2):
                super_batch = []
                for j, batch in enumerate(batcher(batch_matrix_tables, int(N_BATCH / 10))):
                    super_batch.append(merge_matrix_tables_rows(load_db_batch(batch)))
                final_combined_mt = merge_matrix_tables_rows(super_batch)
            else:
                final_combined_mt = merge_matrix_tables_rows(load_db_batch(batch_matrix_tables))
            reduced_mt = reduce_to_2d_table(
                final_combined_mt
            )  # Reduce the matrix table to a 2D matrix table with gene and frequency as keys

            result = create_frequency_bins(
                reduced_mt
            )  # Create a final output frequency table with 16 bins
            # Export the final frequency table to a TSV file
            write_frequency_table(
                result,
                self.args,
                name="frequency_table_{0}".format(
                    len(positive_or_negative_set)
                ),
                extra_tag=set_tag,
            )
        handle_quit()

    def handle_read_vcf_to_pickle_command(self):
        match_and_anti_match = collect_vcfs_by_phenotype(self,
            allow_null_prefix=True
        )

        final_dfs = []
        batch_size = N_BATCH

        for k, vcf_group in enumerate(match_and_anti_match):
            all_batches = []
            for i in range(0, len(vcf_group), batch_size):
                batch = vcf_group[i:i + batch_size]
                batch_df = pd.DataFrame()

                for vcf in batch:
                    mt = import_and_annotate_vcf(vcf, metadata=self.metadata, annotate=self.args.annotate,
                                                 interval=self.args.interval)
                    bin_dict = get_variants_counts_per_gene(mt)
                    freq_bins = create_frequency_bins(bin_dict, _raw_out=False)
                    batch_df = batch_df.add(freq_bins, fill_value=0)
                    mt.unpersist()
                    del mt
                    gc.collect()

                batch_df = batch_df.astype("UInt64")
                batch_path = Path(self.args.dest).joinpath(f"partial_N_batch_{N_BATCH}_{self.args.phenotype}_{k}_batch{i}.pkl")
                batch_df.to_pickle(batch_path)
                info(f"Saved batch {i} to {batch_path}")
                all_batches.append(batch_path)

            # Combine batches
            combined_df = pd.DataFrame()
            for batch_path in all_batches:
                df = pd.read_pickle(batch_path)
                combined_df = combined_df.add(df, fill_value=0)

            combined_df = combined_df.astype("UInt64")
            final_dfs.append(combined_df)

        assert len(final_dfs) == 2, "Expected both case and control sets."
        case_df, control_df = final_dfs
        case_df.columns = [f"{col}_case" for col in case_df.columns]
        control_df.columns = [f"{col}_control" for col in control_df.columns]

        joined_df = case_df.join(control_df, how="outer").fillna(0)
        out_path = Path(self.args.dest).joinpath(f"gene_bin_aggregate_counts_case_control_{self.args.phenotype}.pkl")
        joined_df.to_pickle(out_path)
        info(f"✅ Saved merged case/control frequency data to {out_path}")

    def handle_check_relatedness(self):
        """
        Handles the check relatedness command, using KING inference.
        Outputs a table of probably related samples and phi scores in the specified output folder
        """
        mt_paths = list()
        full_paths = [Path(path) for path in self.args.file]
        for path in full_paths:
            if path.is_file() and path.suffix == ".vcf":
                mt_paths.append(path)
            else:
                mt_paths.extend(path.glob("*.vcf"))

        batch_matrix_tables = []
        # don't annotate, downfilter to chr2
        for i, batch in enumerate(batcher(mt_paths, N_BATCH)):
            mts = import_and_annotate_vcf_batch(batch, annotate=False, interval=["2"])

            mt_combined = merge_matrix_tables_cols(mts)
            info(
                "Combined batch {0} of {1}".format(i, int(len(mt_paths) / N_BATCH))
            )
            batch_matrix_tables.append(mt_combined)
        mt_combined = merge_matrix_tables_cols(batch_matrix_tables)
        king = create_related_samples_table(mt_combined)
        king.filter_entries(king.phi > 0.4, keep=True).show()
        dest = Path(self.args.dest).joinpath("relatedness.tsv")
        if not dest.exists():
            dest.parent.mkdir()
        king.entries().to_pandas().to_csv(dest.__str__())
        info("Wrote relatedness file to: {0}".format(dest))

    def handle_pca(self):
        """
        This function will graph PCA relatedness from a relatedness table. Unused.
        """
        pca_graphing(self.args.pca, self.args.pca_tsv)

    def handle_run_skat_command(self):
        """
        Reads VCFs per phenotype group (cases and controls), processes in batches,
        extracts genotype matrices, joins them across samples, and saves final matrices
        as pickled pandas DataFrames for SKAT-compatible input.
        | gene\:pos | Sample\_1 | Sample\_2 | Sample\_3 | ... |
        | --------- | --------- | --------- | --------- | --- |
        | VAR\_1    | 0         | 1         | 2         |     |
        | VAR\_2    | 1         | 0         | 1         |     |
        | ...       | ...       | ...       | ...       |     |

        """
        match_and_anti_match = collect_vcfs_by_phenotype(self, allow_null_prefix=False)
        phenotype_vector = []
        sample_ids = []
        for k, vcf_list in enumerate(match_and_anti_match):
            cohort_tag = "positive" if k == 0 else "negative"
            cohort_tag = f"{self.args.phenotype}_" + cohort_tag + "_SKAT"

            out_dir = Path(self.args.dest) / cohort_tag
            out_dir.mkdir(parents=True, exist_ok=True)

            info(f"⚙️  Starting SKAT preprocessing for: {cohort_tag} ({len(vcf_list)} VCFs)")

            batch_idx = 0
            batch_df_list = []

            for i, batch in enumerate(batcher(vcf_list, N_BATCH)):
                info(f"🧪 Processing batch {i + 1} / {math.ceil(len(vcf_list) / N_BATCH)}")
                if i > 3:
                    break
                batch_combined_df = []
                for vcf in batch:
                    sample_id = trim_prefix(vcf.stem)
                    sample_ids.append(sample_id)
                    phenotype_vector.append(1 if k == 0 else 0) # one phenotype per column

                    mt = import_and_annotate_vcf(
                        vcf,
                        metadata=self.metadata,
                        annotate=self.args.annotate,
                        interval=self.args.interval,
                    )
                    batch_combined_df.append(mt)

                batch_combined_df = merge_matrix_tables_cols(batch_combined_df)
                out_path = out_dir / f"genotype_batch_{batch_idx+1}.mt"
                batch_combined_df.checkpoint(out_path, overwrite=True)

                batch_combined_df.to_pickle(out_path)
                info(f"📦 Saved batch {batch_idx+1} to: {out_path}")
                batch_df_list.append(out_path)
                batch_idx += 1

            # Merge all batch dataframes into final cohort-wide matrix
            final_mt = merge_matrix_tables_cols()

            final_df = final_df.fillna(0).astype("int32")
            final_out_path = Path(self.args.dest) / f"genotype_matrix_{cohort_tag}.pkl"
            final_df.to_pickle(final_out_path)

            phenotype_df = pd.DataFrame({
                "sample_id": sample_ids,
                "phenotype": phenotype_vector
            })
            phenotype_df.to_csv(Path(self.args.dest) / "phenotype_vector.csv", index=False)
            phenotype_df.to_pickle(Path(self.args.dest) / "phenotype_vector.pkl")
            info(f"✅ Saved final SKAT genotype matrix for {cohort_tag} to {final_out_path}")
            info(f"Saved phenotype-vector pickle to {Path(self.args.dest).joinpath('phenotype_vector.pkl')}")



