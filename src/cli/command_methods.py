import datetime
import os
import hail
from pathlib import Path
from hail.utils import info

from src.utils.file.file_meta import get_metadata
from src.utils.general.data_manipulation import batcher
from src.data_processing.file_io.readers import find_filetype
from src.data_processing.file_io.writers import write_filelist, trim_prefix
from src.data_processing.vcf.read import (
    import_and_annotate_vcf_batch_read,
    import_and_annotate_vcf_batch,
    load_db_batch,
    load_mt, import_and_annotate_gvcf_batch,
)
from src.data_processing.vcf.transform import reduce_to_2d_table, create_frequency_bins
from src.data_processing.hail.genomic_operations import (
    merge_matrix_tables_rows,
    merge_matrix_tables_cols,
    create_related_samples_table,
)
from src.data_processing.pca.analysis import pca_graphing


unique = hash(datetime.datetime.utcnow())  # TODO: delete, variable not used
N_BATCH = 75


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
        self.handle_pca = None
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
        print("Handle_read_vcf_command starting...")
        full_paths = [Path(path) for path in self.args.file]
        print("full_paths = [Path(path) for path in self.args.file] done now")

        vcfs = []
        gvcfs = []

        print(" vcfs, gvcfs = [] done now")


        #  If phenotype is set, write two tables of both phenotype match and the opposite set
        match_and_anti_match = list()

        # Collect VCF and GVCF files from specified paths
        for path in full_paths:
            if path.is_file() and path.suffix in {".vcf", ".gvcf"}:
                (gvcfs if path.suffix == ".gvcf" else vcfs).append(path)
            else:
                vcfs.extend(path.rglob("*.vcf"))
                gvcfs.extend(path.rglob("*.gvcf"))

        metadata = get_metadata(self.args.globals)
        assert metadata is not None

        match_and_anti_match = []
        if self.args.phenotype:
            # Filter according to the phenotype in metadata
            matches = find_elements(metadata, self.args.phenotype, sep=",")
            filtered_vcfs = [p for p in vcfs if trim_prefix(p.stem) in matches]
            filtered_gvcfs = [p for p in gvcfs if trim_prefix(p.stem) in matches]

            assert filtered_vcfs or filtered_gvcfs, f"No matches for phenotype {self.args.phenotype}, quitting."
            anti_match_vcfs = list(set(vcfs) - set(filtered_vcfs))
            anti_match_gvcfs = list(set(gvcfs) - set(filtered_gvcfs))
            match_and_anti_match = [(filtered_vcfs, filtered_gvcfs), (anti_match_vcfs, anti_match_gvcfs)]

            # Log matched and anti-matched counts
            hail.utils.info(
                f"Found {len(filtered_vcfs)} VCF and {len(filtered_gvcfs)} GVCF matches for phenotype {self.args.phenotype}. "
                f"Anti-matches: {len(anti_match_vcfs)} VCF and {len(anti_match_gvcfs)} GVCF."
            )
        else:
            match_and_anti_match = [(vcfs, gvcfs)] # append all

        # Create frequency tables for both match and anti-match
        for k, (vcf_set, gvcf_set) in enumerate(match_and_anti_match):
            n_batches = (len(vcf_set) + len(gvcf_set)) // N_BATCH
            batch_matrix_tables = []
            set_tag = f"{self.args.phenotype}_{'positive' if k == 0 else 'negative'}"

            for i, batch in enumerate(batcher(vcf_set + gvcf_set, N_BATCH)):
                # Initialize as None for each batch
                vcf_combined_mt = None
                gvcf_combined_mt = None

                vcf_batch = [file for file in batch if file.suffix == ".vcf"]
                gvcf_batch = [file for file in batch if file.suffix == ".gvcf"]

                # Process VCFs and GVCFs in separate batches
                if vcf_batch:
                    vcf_combined_mt = import_and_annotate_vcf_batch_read(
                        vcf_batch, metadata=metadata, annotate=self.args.annotate
                    )
                    if vcf_combined_mt:
                        batch_matrix_tables.append(vcf_combined_mt)

                if gvcf_batch:
                    gvcf_combined_mt = import_and_annotate_gvcf_batch(
                        gvcf_batch, metadata=metadata, annotate=self.args.annotate
                    )
                    if gvcf_combined_mt:
                        batch_matrix_tables.append(gvcf_combined_mt)

                if self.args.write:
                    for mt, file_type in zip([vcf_combined_mt, gvcf_combined_mt], ["vcf", "gvcf"]):
                        if mt is not None:  # Check if mt is defined
                            mt.write(
                                Path(self.args.dest)
                                .joinpath(f"batch_{i}_{len(batch)}_{file_type}_{set_tag}.mt")
                                .__str__()
                            )
                            info(f"Wrote batch {i}/{n_batches} for {file_type.upper()} files.")

                final_combined_mt = merge_matrix_tables_rows(batch_matrix_tables)
                reduced_mt = reduce_to_2d_table(final_combined_mt)
                result = create_frequency_bins(reduced_mt)
                write_frequency_table(result, self.args, name="frequency_table", extra_tag=set_tag)
            handle_quit()


def handle_load_db_command(self):
    """
    Loads Hail MatrixTables as a database, i.e. loads VCF and GVCF batches separately into outputs.
    Outputs the combined frequency bins table (create_frequency_bins).
    Uses batching similar to handle_read_vcfs_command.
    """
    metadata_dict = get_metadata(self.args.globals) if self.args.globals else None
    full_paths = [Path(path) for path in self.args.file]
    mt_paths, gvcf_paths = [], []

    # Separate VCF MatrixTables and GVCF files
    for path in full_paths:
        if path.suffix == ".mt":
            mt_paths.append(path)
        elif path.suffix == ".gvcf":
            gvcf_paths.append(path)

    # Process VCF MatrixTables
    batch_matrix_tables = []
    for i, batch in enumerate(batcher(mt_paths, N_BATCH // 10)):
        batch_combined_mt = merge_matrix_tables_rows(load_db_batch(batch))
        batch_matrix_tables.append(batch_combined_mt)
        if self.args.write:
            batch_combined_mt.write(
                Path(self.args.dest)
                .joinpath(f"multi_batch_dataset_{i}_{len(batch)}.mt_combined")
                .__str__()
            )
        hail.utils.info(
            f"Processed VCF batch {i + 1}/{len(mt_paths) // (N_BATCH // 10)} with {len(batch)} files."
        )

    # Process GVCFs if any
    gvcf_matrix_tables = []
    if gvcf_paths:
        hail.utils.info(f"Found {len(gvcf_paths)} GVCF files to process.")
        for i, batch in enumerate(batcher(gvcf_paths, N_BATCH // 10)):
            batch_combined_gvcf_mt = import_and_annotate_gvcf_batch(
                batch, metadata=metadata_dict, annotate=self.args.annotate
            )
            gvcf_matrix_tables.append(batch_combined_gvcf_mt)
            if self.args.write:
                batch_combined_gvcf_mt.write(
                    Path(self.args.dest)
                    .joinpath(f"multi_batch_gvcf_{i}_{len(batch)}.mt_combined")
                    .__str__()
                )
            hail.utils.info(
                f"Processed GVCF batch {i + 1}/{len(gvcf_paths) // (N_BATCH // 10)} with {len(batch)} files."
            )

    # Merge VCF tables only (do not combine with GVCFs to avoid format inconsistencies)
    if batch_matrix_tables:
        hail.utils.info("Merging all VCF MatrixTables into final combined table.")
        final_combined_mt = merge_matrix_tables_rows(batch_matrix_tables)
        dest = Path(self.args.dest).joinpath("multi_batch_dataset_merged.mt")
        if self.args.write and (not dest.exists() or self.args.overwrite):
            final_combined_mt.write(dest.__str__())

        # Create frequency bins and output the final frequency table for VCFs
        result = create_frequency_bins(reduce_to_2d_table(final_combined_mt))
        write_frequency_table(result, self.args, name="frequency_table", extra_tag="db_vcf")
        hail.utils.info("Completed processing and merging for VCF tables.")

    # Process and save GVCF outputs separately if any
    if gvcf_matrix_tables:
        hail.utils.info("Merging all GVCF MatrixTables into final combined table.")
        final_combined_gvcf_mt = merge_matrix_tables_rows(gvcf_matrix_tables)
        dest_gvcf = Path(self.args.dest).joinpath("multi_batch_dataset_gvcf_merged.mt")
        if self.args.write and (not dest_gvcf.exists() or self.args.overwrite):
            final_combined_gvcf_mt.write(dest_gvcf.__str__())

        # Create frequency bins and output the final frequency table for GVCFs
        result_gvcf = create_frequency_bins(reduce_to_2d_table(final_combined_gvcf_mt))
        write_frequency_table(result_gvcf, self.args, name="frequency_table_gvcf", extra_tag="db_gvcf")
        hail.utils.info("Completed processing and merging for GVCF tables.")

    # Finish command and exit
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
            if path.is_file() and (
                path.suffix in {".vcf", ".gvcf"}
            ):  # gvcf files will have .g.vcf, .genome.vcf which wil be captured by .vcf only
                mt_paths.append(path)
            else:
                mt_paths.extend(path.glob("*.vcf"))
                mt_paths.extend(path.glob("*.gvcf"))

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
