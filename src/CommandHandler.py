import os
import sys
from pathlib import Path
import datetime

from hail.utils import info

from .hail_methods import import_and_annotate_vcf, merge_matrix_tables, multi_way_union_mts, reduce_to_2d_table, create_frequency_bins, \
    import_and_annotate_vcf_batch, load_db
from .utils import find_filetype, batcher, get_metadata
from .file_utility import write_filelist

unique = hash(datetime.datetime.utcnow())
N_BATCH = 200

class CommandHandler:
    def __init__(self, args):
        self.args = args

    def handle_find_type_command(self):
        files = find_filetype(self.args.source, self.args.type, verbose=False)
        write_filelist(self.args.directory, "{0}.{1}.txt".format(os.path.basename(os.path.normpath(self.args.source)), self.args.type), files, regex=self.args.regex)

    def handle_read_vcfs_command(self):
        full_paths = [Path(path) for path in self.args.file]
        vcfs = []

        for path in full_paths:
            if path.is_file() and path.suffix == ".vcf":
                vcfs.append(path)
            else:
                vcfs.extend(path.glob("*.vcf"))
        n_batches = int(len(vcfs)/N_BATCH)
        batch_matrix_tables = []
        for i, batch in enumerate(batcher(vcfs, N_BATCH)):
            matrix_tables = import_and_annotate_vcf_batch(batch, self.args.annotate)  # hail.import_vcf() and hail.VEP() within this function

            batch_combined_mt = merge_matrix_tables(matrix_tables)  # Merge batch of matrix tables into one
            if self.args.write:
                batch_combined_mt.write(os.path.join(self.args.dest, "multi_batch_{0}_dataset_{1}.mt".format(i, len(vcfs)))) # Write MatrixTable to disk
                info("Wrote batch {0}/{1} onto disk with {2} subelements.".format(i, n_batches, N_BATCH))

            batch_matrix_tables.append(batch_combined_mt)
        if len(batch_matrix_tables) > 10:
            super_batch = []
            for j, batch in enumerate(batch_matrix_tables, N_BATCH % 5):
                super_batch.append(merge_matrix_tables(batch))
            final_combined_mt = merge_matrix_tables(super_batch)
        else:
            final_combined_mt = merge_matrix_tables(batch_matrix_tables)
        reduced_mt = reduce_to_2d_table(
            final_combined_mt)  # Reduce the matrix table to a 2D matrix table with gene and frequency as keys

        final_frequency_table = create_frequency_bins(reduced_mt,
                                                      num_bins=16)  # Create a final output frequency table with 16 bins

        # Export the final frequency table to a TSV file
        #final_frequency_table.flatten().export(os.path.join(self.args.dest, "final_frequency_table.tsv"))

    def handle_load_db_command(self):
        metadata_dict = None
        mt_paths = []
        if self.args.globals is not None:
            metadata_dict = get_metadata(self.args.globals)
        full_paths = [Path(path) for path in self.args.file]
        for path in full_paths:
            if path.is_dir() and path.suffix == ".mt":
                mt_paths.append(path.__str__())
            else:
                mt_paths.extend(path.glob("*.mt"))

        batch_matrix_tables = load_db(mt_paths)
        if len(batch_matrix_tables) > 10:
            super_batch = []
            for j, batch in enumerate(batch_matrix_tables, N_BATCH % 5):
                super_batch.append(merge_matrix_tables(batch))
            final_combined_mt = merge_matrix_tables(super_batch)
        else:
            final_combined_mt = merge_matrix_tables(batch_matrix_tables)