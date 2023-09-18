import os


class CommandFactory:

    def __init__(self, parser):
        self.find_type = None
        self.parser = parser
        self.subparsers = parser.add_subparsers(title="commands", dest="command")

    def create_find_type_command(self):
        self.find_type = self.subparsers.add_parser(
            "findtype", help="Find all specific files of a given filetype."
        )
        self.find_type.add_argument(
            "-s", "--source", help="Directory to be searched.", type=str, required=True
        )

        # TODO: Default value returns None type
        self.find_type.add_argument(
            "-d",
            "--directory",
            help="Directory to be saved to.",
            nargs="?",
            const=".",
            type=str,
            required=True,
        )

        self.find_type.add_argument(
            "-t",
            "--type",
            help="Filetype to be used.",
            action="store",
            type=str,
            required=True,
        )
        self.find_type.add_argument(
            "-r", "--regex", help="Filter strings by including parsed regex.", type=str
        )

    def create_read_vcfs_command(self):
        read_vcfs = self.subparsers.add_parser(
            "readvcfs", help="Turn VCF file(s) into a Hail MatrixTable."
        )
        read_vcfs.add_argument(
            "-f",
            "--file",
            help="The VCF file(s) [comma seperated], .txt/.list of VCF paths to be parsed or folder containing VCF "
            "files.",
            nargs="+",
            required=True,
        )
        read_vcfs.add_argument(
            "-d",
            "--dest",
            help="Destination folder to write the frequency tables and Hail MatrixTable files (if --write is active).",
            nargs="?",
            const=os.path.abspath(".."),
            required=True,
        )
        read_vcfs.add_argument(
            "-r",
            "--overwrite",
            help="Overwrites any existing output MatrixTables, HailTables.",
            action="store_true",
        )
        read_vcfs.add_argument(
            "-g",
            "--globals",
            help="Tab delimited input file containing globals string  for a given unique sample (e.g. "
            "Identifier\\t.Phenotype\\tMutations",
            required=True,
            type=str,
        )
        read_vcfs.add_argument(
            "-p",
            "--phenotype",
            help="Filter according to this specific phenotype.",
            required=False,
            default=None,
            type=str,
        )
        read_vcfs.add_argument(
            help="Annotate flag will annotate input VCFs with Ensembl VEP. Set to false to skip annotations."
            "Annotated VCFs must be contain VEP annotations (CSQ string must be present with gene names, "
            "allele frequencies).",
            dest="annotate",
            default=True,  # sets the default flag for "annotate" as active implicitly
            action="store_true",
        )
        read_vcfs.add_argument(
            "--no-annotate",
            help="Settings this flag will skip annotations using VEP.",
            dest="annotate",
            action="store_false",
        )
        read_vcfs.add_argument(
            "--write",
            help="Write the unioned VCF MatrixTable to disk for later loading and handling.",
            default=False,
            action="store_true",
        )
        read_vcfs.add_argument(
            "-t",
            "--temp",
            help="Destination folder for all Hail temp files",
            nargs="?",
            const=os.path.abspath(".."),
            default="/tmp/",
            required=False,
        )

    def create_load_db_command(self):
        loaddb = self.subparsers.add_parser(
            "loaddb", help="Load a folder containing HailTables."
        )
        loaddb.add_argument(
            "-f",
            "--file",
            help="The VCF MatrixTable folder(s) [comma seperated], .txt/.list of VCF paths to be parsed or folder containing VCF "
            "files.",
            nargs="+",
            required=True,
        )
        loaddb.add_argument(
            "-r",
            "--overwrite",
            help="Overwrites any existing output MatrixTables, HailTables.",
            action="store_true",
            default=False,
            required=False,
        )
        loaddb.add_argument(
            "-d",
            "--dest",
            help="Destination folder to write the frequency tables and Hail MatrixTable files (if --write is active).",
            nargs="?",
            const=os.path.abspath(".."),
            required=True,
        )
        # Deprecated? Load the entire folder, might be useful if Nextflow scripts are used to add massive parallelization
        loaddb.add_argument(
            "-n",
            "--number",
            help="Number of tables to be collated.",
            nargs="?",
            required=False,
            type=int,
            default=-1,
        )
        loaddb.add_argument(
            "-g",
            "--globals",
            help="Tab delimited input file containing globals string  for a given unique sample (e.g. "
            "Identifier\\t.Phenotype\\tMutations",
            required=True,
            type=str,
        )
        loaddb.add_argument(
            "--phenotype",
            help="Filter a subset of samples with given phenotype. Regex strings accepted e.g. r'NA\d+",
            required=False,
            type=str,
        )
        loaddb.add_argument(
            "--write",
            help="Write the unioned VCF MatrixTable to disk for later loading and handling.",
            default=False,
            action="store_true",
        )

    def create_check_relatedness_command(self):
        check_relatedness = self.subparsers.add_parser(
            "relatedness2", help="Create a table of KING phi relatedness of input VCFs"
        )
        check_relatedness.add_argument(
            "-f",
            "--file",
            help="The VCF file(s) [comma seperated], .txt/.list of VCF paths to be parsed or folder containing VCF "
            "files.",
            nargs="+",
            required=True,
        )
        check_relatedness.add_argument(
            "-d",
            "--dest",
            help="Destination folder to write the frequency tables and Hail MatrixTable files (if --write is active).",
            nargs="?",
            const=os.path.abspath(".."),
            required=True,
        )
        check_relatedness.add_argument(
            "-t",
            "--temp",
            help="Destination folder for all Hail temp files",
            nargs="?",
            const=os.path.abspath(".."),
            default="/tmp/",
            required=False,
        )

    def create_pca_command(self):
        pca_cmd = self.subparsers.add_parser(
            "pca", help="Output a graph from PCA data."
        )
        pca_cmd.add_argument(
            "--pca",
            help="Input PCA file.",
            nargs="?",
            const=os.path.abspath(".."),
            required=True,
        )

        pca_cmd.add_argument(
            "--pca_tsv",
            help="Input PCA_tsv file.",
            nargs="?",
            const=os.path.abspath(".."),
            required=True,
        )
        pca_cmd.add_argument(
            "-d",
            "--dest",
            help="Destination folder to write the frequency tables and Hail MatrixTable files (if --write is active).",
            nargs="?",
            const=os.path.abspath(".."),
            required=True,
        )
