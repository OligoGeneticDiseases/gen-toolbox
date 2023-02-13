import os
from pathlib import Path
import datetime
from utils import find_filetype, write_filelist, get_metadata, load_hailtables

from utils import write_gnomad_table

unique = hash(datetime.datetime.utcnow())


class CommandHandler:
    def __init__(self, args):
        self.args = args

    def handle_find_type_command(self):
        files = find_filetype(self.args.source, self.args.type, verbose=False)
        write_filelist(self.args.directory, "{0}.{1}.txt".format(os.path.basename(os.path.normpath(self.args.source)), self.args.type), files, regex=self.args.regex)

    def handle_read_vcfs_command(self):
        full_paths = [Path(path) for path in self.args.file]
        files = set()
        for path in full_paths:
            if path.is_file():
                if path.suffix == ".vcf":  # VCF files are parsed.
                    files.add(path)
                else:  # might be a list of VCFs
                    with open(path, "r") as filelist:
                        for line in filelist:
                            # coerce lines into path
                            p = Path(line.strip())
                            if p.suffix == ".vcf":
                                files.add(p)
            else:  # Glob folder for *.VCF
                files |= set(path.glob("*.vcf"))
        gnomad_path = Path(self.args.dest).joinpath(Path("gnomad_tb"))
        if gnomad_path.exists():
            if self.args.overwrite:
                metadata = self.args.globals
            else:
                FileExistsError("The combined gnomad_tb exists and --overwrite is not active! "
                                "Rename or move the folder {0}".format(gnomad_path.__str__()))
        else:
            gnomad_tb = write_gnomad_table(files, self.args.dest, overwrite=self.args.overwrite,
                                           metadata=self.args.globals)
        # gnomad_tb.describe()
        gnomad_tb.flatten().export(Path(self.args.dest).parent.joinpath("gnomad.tsv").__str__())

    def handle_load_db_command(self):
        metadata_dict = None
        if self.args.globals is not None:
            metadata_dict = get_metadata(self.args.globals)

        dirpath = Path(self.args.directory)
        gnomad_tb = load_hailtables(dirpath, self.args.number, self.args.out, metadata_dict, self.args.overwrite, self.args.phenotype)
        gnomad_tb.describe()
        gnomad_tb.flatten().export(Path(self.args.out).parent.joinpath("gnomad_tb{0}.tsv".format(unique)).__str__())