import os
import sys
from pathlib import Path
import datetime

import hail.utils
from hail.utils import info

from .hail_methods import import_and_annotate_vcf, merge_matrix_tables, multi_way_union_mts, reduce_to_2d_table, \
    create_frequency_bins, \
    import_and_annotate_vcf_batch, load_db_batch, load_mt
from .utils import find_filetype, batcher, get_metadata
from .file_utility import write_filelist

unique = hash(datetime.datetime.utcnow())
N_BATCH = 50

def handle_quit():
    #TODO: log some extra output for Nextlow / ensure correct error codes are sent
    hail.stop()


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
                batch_combined_mt.write(Path(self.args.dest).joinpath(
                                                     "multi_batch_{0}_{1}_dataset_{2}.mt".format(i, len(batch),
                                                                                                 len(vcfs))).__str__())  # Write MatrixTable to disk
                info("Wrote batch {0}/{1} onto disk with {2} subelements.".format(i, n_batches, N_BATCH))

            batch_matrix_tables.append(batch_combined_mt)
        if len(batch_matrix_tables) > 10:
            super_batch = []
            for j, batch in enumerate(batch_matrix_tables, int(N_BATCH / 5)):
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
        handle_quit()

    def handle_load_db_command(self):
        metadata_dict = None
        mt_paths = []
        if self.args.globals is not None:
            metadata_dict = get_metadata(self.args.globals)
        full_paths = [Path(path) for path in self.args.file]
        for path in full_paths:
            if path.is_dir():
                if path.suffix == ".mt":
                    mt_paths.append(path)
                else:
                    mt_paths.extend(path.glob("*.mt"))

        batch_matrix_tables = []
        steps = 0
        if len(mt_paths) > 5:
            # Reduce that batch size by a factor of 10 when loading large subsets
            # TODO: Batching doesn't actually do much
            #  unless write to disk is active due Pyspark only handling references until actual data is called
            for i, batch in enumerate(batcher(mt_paths, int(N_BATCH / 10))):

                batch_combined_mt = merge_matrix_tables(load_db_batch(batch))
                batch_matrix_tables.append(batch_combined_mt)
                if self.args.write:
                    batch_combined_mt.write(Path(self.args.dest).joinpath("multi_batch_dataset_{0}_{1}.mt_combined".format(i, len(batch))).__str__())
                hail.utils.info("OLIGO: multi_batch_dataset_{0}_{1}.mt_combined complete".format(i, len(batch)))
                steps += 1
            super_batch = []
            if len(batch_matrix_tables) > 5:
                # TODO: don't save into memory if not required, Hail will persist files anyway
                for j, batch in enumerate(batcher(batch_matrix_tables, int(N_BATCH / 5))):
                    super_batch.append(merge_matrix_tables(batch))
                    hail.utils.info("OLIGO: super_multi_batch_dataset_{0}_{1} complete".format(j, len(batch)))
                    steps += 1
            else:
                super_batch = batch_matrix_tables
            final_combined_mt = merge_matrix_tables(super_batch)
        else:
            steps = 1
            if len(mt_paths) > 1:
                batch_matrix_tables = load_db_batch(mt_paths)
                final_combined_mt = merge_matrix_tables(batch_matrix_tables)
            else:
                final_combined_mt = load_mt(mt_paths[0].__str__())
        dest = Path(self.args.dest).joinpath("multi_batch_dataset_merged_{0}.mt".format(steps))
        if self.args.write and (not dest.exists() or self.args.overwrite):
            # Can't write and load recursively from the same file, therefore skip writing if file exists
            final_combined_mt.write(dest.__str__())
        hail.utils.info("OLIGO: Merge load_db command in {0} steps".format(steps))
        result = create_frequency_bins(reduce_to_2d_table(final_combined_mt))
        if self.args.write and (not dest.exists() or self.args.overwrite):
            result.write(hail.utils.timestamp_path(os.path.join(self.args.dest, "out_{0}.tsv")))
        hail.utils.info("OLIGO: Finished LoadDB command.")
