import os
from pathlib import Path
import hail as hl
from hail_methods import import_and_annotate_vcf, merge_matrix_tables, reduce_to_2d_table, create_frequency_bins

unique = hash(datetime.datetime.utcnow())


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

        matrix_tables = []
        for vcf_path in vcfs:
            annotated_vcf = import_and_annotate_vcf(vcf_path)  # hail.import_vcf() and hail.VEP() within this function
            matrix_tables.append(annotated_vcf)

        combined_mt = merge_matrix_tables(matrix_tables)  # Merge all matrix tables into one

        reduced_mt = reduce_to_2d_table(
            combined_mt)  # Reduce the matrix table to a 2D matrix table with gene and frequency as keys

        final_frequency_table = create_frequency_bins(reduced_mt,
                                                      num_bins=16)  # Create a final output frequency table with 16 bins

        # Export the final frequency table to a TSV file
        final_frequency_table.flatten().export(os.path.join(self.args.dest, "final_frequency_table.tsv"))

    def handle_load_db_command(self):
        metadata_dict = None
        if self.args.globals is not None:
            metadata_dict = get_metadata(self.args.globals)

        dirpath = Path(self.args.directory)
        gnomad_tb = load_hailtables(dirpath, self.args.number, self.args.out, metadata_dict, self.args.overwrite, self.args.phenotype)
        gnomad_tb.describe()
        gnomad_tb.flatten().export(Path(self.args.out).parent.joinpath("gnomad_tb{0}.tsv".format(unique)).__str__())