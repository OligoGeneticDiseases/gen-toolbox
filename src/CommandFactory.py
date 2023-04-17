import os


class CommandFactory:
    def __init__(self, parser):
        self.find_type = None
        self.parser = parser
        self.subparsers = parser.add_subparsers(title="commands", dest="command")

    def create_find_type_command(self):
        self.find_type = self.subparsers.add_parser("findtype", help="Find all specific files of a given filetype.")
        self.find_type.add_argument("-s", "--source", help="Directory to be searched.", type=str, required=True)

        # TODO: Default value returns None type
        self.find_type.add_argument("-d", "--directory", help="Directory to be saved to.", nargs='?', const=".",
                                    type=str, required=True)

        self.find_type.add_argument("-t", "--type", help="Filetype to be used.", action="store", type=str,
                                    required=True)
        self.find_type.add_argument("-r", "--regex", help="Filter strings by including parsed regex.", type=str)

    def create_read_vcfs_command(self):
        read_vcfs = self.subparsers.add_parser("readvcfs", help="Turn VCF file(s) into a Hail MatrixTable.")
        read_vcfs.add_argument(
            "-f",
            "--file",
            help="The VCF file(s) [comma seperated], .txt/.list of VCF paths to be parsed or folder containing VCF "
                 "files.",
            nargs='+',
            required=True
        )
        read_vcfs.add_argument(
            "-d",
            "--dest",
            help="Destination folder to write the Hail MatrixTable files.",
            nargs='?',
            const=os.path.abspath(".."),
            required=True
        )
        read_vcfs.add_argument(
            "-r",
            "--overwrite",
            help="Overwrites any existing output MatrixTables, HailTables.",
            action="store_true"
        )
        read_vcfs.add_argument(
            "-g",
            "--globals",
            help="Tab delimited input file containing globals string  for a given unique sample (e.g. "
                 "Identifier\\t.Phenotype\\tMutations",
            required=True,
            type=str
        )

    def create_load_db_command(self):
        loaddb = self.subparsers.add_parser("loaddb", help="Load a folder containing HailTables.")
        loaddb.add_argument("-d",
                            "--directory",
                            help="Folder to load the Hail MatrixTable files from.",
                            nargs='?',
                            const=os.path.abspath(".."),
                            required=True
                            )
        loaddb.add_argument("-r",
                            "--overwrite",
                            help="Overwrites any existing output MatrixTables, HailTables.",
                            action="store_true",
                            required=False
                            )
        loaddb.add_argument("-o",
                            "--out",
                            help="Output destination.",
                            required=True,
                            type=str
                            )
        loaddb.add_argument("-n",
                            "--number",
                            help="Number of tables to be collated.",
                            nargs="?",
                            required=False,
                            type=int,
                            default=-1,
                            )
        loaddb.add_argument("-g",
                            "--globals",
                            help="Tab delimited input file containing globals string for a given unique sample.",
                            type=str
                            )
        loaddb.add_argument("--phenotype",
                            help="Filter a subset of samples with given phenotype. Regex strings accepted e.g. r'NA\d+",
                            required=False,
                            type=str
                            )
