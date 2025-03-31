import datetime
import os
import traceback

import hail
from pathlib import Path

import math
from hail.utils import info

from src.utils.file.file_meta import get_metadata
from src.utils.general.data_manipulation import batcher
from src.data_processing.file_io.readers import find_filetype
from src.data_processing.file_io.writers import write_filelist, trim_prefix
from src.data_processing.vcf.read import import_and_annotate_vcf_batch, load_db_batch, load_mt
from src.data_processing.vcf.transform import reduce_to_2d_table, create_frequency_bins
from src.data_processing.hail.genomic_operations import merge_matrix_tables_rows, merge_matrix_tables_cols, create_related_samples_table
from src.data_processing.pca.analysis import pca_graphing


N_BATCH = 500


def handle_quit():
    # TODO: log some extra output for Nextlow / ensure correct error codes are sent
    hail.stop()


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
    table_path = hail.utils.timestamp_path(
        os.path.join(args.dest, "{0}_{1}.csv".format(name, extra_tag))
    )
    if not Path(table_path).exists() or args.overwrite:
        result.to_csv(table_path)
        hail.utils.info("OLIGO: wrote output to {0}".format(table_path))
    else:
        hail.utils.info("OLIGO: skipped writing output:\n{0}.".format(result))


class CommandHandler:

    def __init__(self, args):
        self.args = args

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
        hail.default_reference = hail.ReferenceGenome.read("/home/villem/branches/stats-burden-analysis/GRCh37_MT.json")

        full_paths = [Path(path) for path in self.args.file]
        vcfs = []
        #  If phenotype is set, write two tables of both phenotype match and the opposite set
        match_and_anti_match = list()
        for path in full_paths:
            if path.is_file() and path.suffix == ".vcf":
                vcfs.append(path)
            else:
                vcfs.extend(path.rglob("*.vcf"))
        metadata = get_metadata(self.args.globals)
        assert metadata is not None
        if self.args.phenotype is not None:
            # Downfilter according to the info in the globals file
            matches = find_elements(metadata, self.args.phenotype, sep=",")
            filtered_vcfs = list()
            for path in vcfs:
                if (
                    trim_prefix(path.stem) in matches.keys()
                ):  # TODO: make it into dict calling
                    filtered_vcfs.append(path)
            assert (
                len(filtered_vcfs) > 0
            ), "No matches found with given phenotype {0}, quitting.".format(self.args.phenotype)  # Quit the pipeline if no matches are found, otherwise the antimatch will be all
            anti_match = list(set(vcfs) - set(filtered_vcfs))
            match_and_anti_match.append(filtered_vcfs)
            match_and_anti_match.append(anti_match)
            hail.utils.info("Found {0} matches for phenotype {1}. Anti-matches: {2}".format(len(filtered_vcfs),
                                                                                            self.args.phenotype,
                                                                                            len(anti_match)))

            # Create symlink folders from base data (i.e. VCFs)
            for link in filtered_vcfs:
                symlink_path = Path(f"{self.args.temp}/positive_{self.args.phenotype}_{len(filtered_vcfs)}/{link.name}")
                symlink_path.parent.mkdir(parents=True, exist_ok=True)
                if not symlink_path.exists():
                    symlink_path.symlink_to(link)
            for link in anti_match:
                symlink_path = Path(f"{self.args.temp}/negative_{self.args.phenotype}_{len(anti_match)}/{link.name}")
                symlink_path.parent.mkdir(parents=True, exist_ok=True)
                if not symlink_path.exists():
                    symlink_path.symlink_to(link)
        else:
            match_and_anti_match.append(vcfs)  # Append all

        # Create frequency tables for both match and anti-match
        for k, positive_or_negative_set in enumerate(match_and_anti_match):
            n_batches = int(math.ceil(len(positive_or_negative_set) / N_BATCH))
            batch_matrix_tables = []
            set_tag = "positive" if k == 0 else "negative"
            set_tag = "{0}_".format(self.args.phenotype) + set_tag
            info(f"Starting with set: {set_tag} ({len(positive_or_negative_set)}). "
                 f"Dividing work into batches of {N_BATCH}, with {n_batches} batches.")
            if k == 0:
                info(f"Skipping {set_tag} due to debugging tag active.")
                stack = traceback.extract_stack()
                # Reports this location in the logfile.
                caller = stack[-1]
                info(f"File: {caller.filename}, Function: {caller.name}, Line: {caller.lineno}")
                continue #uncommenting this will skip reading positive set
                pass

            for i, batch in enumerate(batcher(positive_or_negative_set, N_BATCH)):
                combined_mt_path = Path(self.args.dest).joinpath(
                    "multi_batch_{0}_{1}_dataset_{2}_{3}.mt".format(
                        i+1, len(batch), len(positive_or_negative_set), set_tag
                    )).as_posix()
                if not Path(combined_mt_path).exists():
                    matrix_tables = import_and_annotate_vcf_batch(
                        batch, metadata=metadata, annotate=self.args.annotate, location=self.args.dest
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

    def handle_load_db_command(self):
        """
        Loads Hail MatrixTables as a database, i.e. load several intermediary batch files into a single output.
        Outputs the combined frequency bins table (create_frequency_bins). Uses batching like in read_vcf
        This command is not required in normal analysis workflow.
        """
        metadata_dict = None
        mt_paths = []
        if self.args.globals is not None:
            metadata_dict = get_metadata(self.args.globals)
        full_paths = [Path(path) for path in self.args.file]
        for path in full_paths:
            if path.is_dir():
                if path.suffix == "negative.mt":
                    mt_paths.append(path)
                else:
                    mt_paths.extend(path.glob("*.mt"))
        mt_paths.sort()
        batch_matrix_tables = []
        steps = 0
        if len(mt_paths) > 1:
            # Reduce that batch size by a factor of 10 when loading large subsets
            # TODO: Batching doesn't actually do much
            #  unless write to disk is active due to Pyspark only handling lazily references until actual data is called
            for i, batch in enumerate(batcher(mt_paths, int(N_BATCH / 10))):
                batch_combined_mt = merge_matrix_tables_rows(load_db_batch(batch))
                batch_matrix_tables.append(batch_combined_mt)
                if self.args.write:
                    batch_combined_mt.write(
                        Path(self.args.dest)
                        .joinpath(
                            "multi_batch_dataset_{0}_{1}_{2}.mt_combined".format(
                                i+1, len(batch), len(mt_paths)+1
                            )
                        )
                        .__str__()
                    )
                hail.utils.info(
                    "OLIGO: multi_batch_dataset_{0}_{1}_{2}.mt_combined complete".format(
                        i+1, len(batch), len(mt_paths)+1
                    )
                )
                steps += 1
            super_batch = []
            if len(batch_matrix_tables) > 5:
                # TODO: don't save into memory if not required, Hail will persist files anyway
                for j, batch in enumerate(
                    batcher(batch_matrix_tables, int(N_BATCH / 5))
                ):
                    super_batch.append(merge_matrix_tables_rows(batch))
                    hail.utils.info(
                        "OLIGO: super_multi_batch_dataset_{0}_{1} complete".format(
                            j, len(batch)
                        )
                    )
                    steps += 1
            else:
                super_batch = batch_matrix_tables
            final_combined_mt = merge_matrix_tables_rows(super_batch)
        else:
            steps = 1
            if len(mt_paths) > 1:
                batch_matrix_tables = load_db_batch(mt_paths)
                final_combined_mt = merge_matrix_tables_rows(batch_matrix_tables)
            else:
                final_combined_mt = load_mt(mt_paths[0].__str__())
        dest = Path(self.args.dest).joinpath(
            "multi_batch_dataset_merged_{0}.mt_combined".format(steps)
        )
        if self.args.write and (not dest.exists() or self.args.overwrite):
            # Can't write and load recursively from the same file, therefore skip writing if file exists
            final_combined_mt.write(dest.__str__())
        hail.utils.info("OLIGO: Merge load_db command in {0} steps".format(steps))
        result = create_frequency_bins(reduce_to_2d_table(final_combined_mt))
        write_frequency_table(
            result,
            self.args,
            name="frequency_table_{0}".format(steps),
            extra_tag="negative",
        )
        hail.utils.info("OLIGO: Finished LoadDB command.")
        handle_quit()

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
            hail.utils.info(
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
        hail.utils.info("Wrote relatedness file to: {0}".format(dest))

    def handle_pca(self):
        """
        This function will graph PCA relatedness from a relatedness table. Unused.
        """
        pca_graphing(self.args.pca, self.args.pca_tsv)
