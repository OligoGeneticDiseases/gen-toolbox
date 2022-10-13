import argparse
import os.path
import shutil
from pathlib import Path
from sys import stderr

import hail as hl
from hailtop.utils import time

import file_utility

mt_path = "/mnt/c/Users/ville/Documents/hail_vep"
vep_json = "/mnt/c/Users/ville/Documents/VEP_cfg.json"


def vcfs_to_matrixtable(f, destination=None, write=True):
    files = list()
    if type(f) is list:
        for vcf in f:
            files.append(vcf)

    elif not f.endswith(".vcf") and not f.endswith(".gz"):
        with open(f) as vcflist:
            for vcfpath in vcflist:
                stripped = vcfpath.strip()
                assert os.path.exists(stripped)
                files.append(stripped)
    else:
        assert os.path.exists(f), "Path {0} does not exist.".format(f)
        files.append(f)  # Only one file

    # recode = {f"chr{i}":f"{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    # Can import only samples of the same key (matrixtable join), if input is list of vcfs
    table = hl.import_vcf(files, force=True, reference_genome='GRCh37', contig_recoding={"chr1": "1",
                                                                                         "chr2": "2",
                                                                                         "chr3": "3",
                                                                                         "chr4": "4",
                                                                                         "chr5": "5",
                                                                                         "chr6": "6",
                                                                                         "chr7": "7",
                                                                                         "chr8": "8",
                                                                                         "chr9": "9",
                                                                                         "chr10": "10",
                                                                                         "chr11": "11",
                                                                                         "chr12": "12",
                                                                                         "chr13": "13",
                                                                                         "chr14": "14",
                                                                                         "chr15": "15",
                                                                                         "chr16": "16",
                                                                                         "chr17": "17",
                                                                                         "chr18": "18",
                                                                                         "chr19": "19",
                                                                                         "chr20": "20",
                                                                                         "chr21": "21",
                                                                                         "chr22": "22",
                                                                                         "chrX": "X",
                                                                                         "chrY": "Y"})
    if write:
        if not os.path.exists(destination):
            table.write(destination)
        else:
            raise FileExistsError(destination)
    return table


def parse_empty(text):
    return hl.if_else(text == "", hl.missing(hl.tint32), hl.float(text))


def append_table(table, out=None, write=False):
    mt_a = table.annotate_rows(CSQ=table.info.CSQ.first().split("\\|"))
    # mt_a = mt_a.drop(mt_a.info) # Drop the already split string
    mt_a = mt_a.annotate_rows(impact=mt_a.CSQ[2])
    mt_a = mt_a.annotate_rows(gene=mt_a.CSQ[3])
    mt_a = mt_a.annotate_rows(Entrez_ID=hl.int(parse_empty(mt_a.CSQ[4])))
    mt_a = mt_a.annotate_rows(AC=mt_a.info.AC)
    mt_a = mt_a.annotate_rows(CADD_phred=hl.float(parse_empty(mt_a.CSQ[33])))
    mt_a = mt_a.filter_entries((hl.len(mt_a.filters) == 0), keep=True)  # Remove all not PASS
    mt_a = mt_a.annotate_rows(gnomAD_exomes_AF=hl.float(parse_empty(mt_a.CSQ[36])))
    mt_a = mt_a.annotate_rows(gnomAD_genomes_AF=hl.float(parse_empty(mt_a.CSQ[40])))
    if write and out is not None:
        mt_a.write(out)
    return mt_a


def parse_tables(tables):
    mt_tables = []
    # Positional arguments from VEP annotated CSQ string. TODO: Query from VCF header
    for i, table in enumerate(tables):
        mt_tables.append(append_table(table))

    return mt_tables


def mts_to_table(tables):
    for i, tb in enumerate(tables):
        tb = tb.entries()  # Convert from MatrixTable to Table
        tables[i] = tb.key_by(tb.gene, tb.Entrez_ID)  # Key by gene
    return tables


def mt_join(mt_list):
    mt_final = None
    for i, mt in enumerate(mt_list):
        if i == 0:
            mt_final = mt
        mt_final = hl.experimental.full_outer_join_mt(mt_final, mt)  # An outer join of MatrixTables
        # mt_final.write(mt_path)
    return mt_final


def table_join(tables_list):
    # Join Tables into one Table.
    if tables_list is not None:
        unioned = tables_list[0]  # Initialize with a single table
    if len(tables_list) > 1:
        unioned = unioned.union(*tables_list[1:])
    return unioned


def gnomad_table(unioned):
    gnomad_tb = unioned.group_by(unioned.gene).aggregate(
        modifier=hl.struct(
            gnom0001=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.001) & (unioned.impact.contains(hl.literal("MODIFIER"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0005=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.001) & (unioned.gnomAD_genomes_AF < 0.005) & (
                unioned.impact.contains(hl.literal("MODIFIER"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0010=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.005) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("MODIFIER")),
                                      hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("MODIFIER")), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.10) & (unioned.impact.contains("MODIFIER")),
                                     hl.agg.array_sum(unioned.AC)[0])),
        low=hl.struct(
            gnom0001=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.001) & (unioned.impact.contains(hl.literal("LOW"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0005=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.001) & (unioned.gnomAD_genomes_AF < 0.005) & (
                unioned.impact.contains(hl.literal("LOW"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0010=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.005) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("LOW")),
                                      hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("LOW")), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.10) & (unioned.impact.contains("LOW")),
                                     hl.agg.array_sum(unioned.AC)[0])),
        moderate=hl.struct(
            gnom0001=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.001) & (unioned.impact.contains(hl.literal("MODERATE"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0005=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.001) & (unioned.gnomAD_genomes_AF < 0.005) & (
                unioned.impact.contains(hl.literal("MODERATE"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0010=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.005) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("MODERATE")),
                                      hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("MODERATE")), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.10) & (unioned.impact.contains("MODERATE")),
                                     hl.agg.array_sum(unioned.AC)[0])),
        high=hl.struct(
            gnom0001=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.001) & (unioned.impact.contains(hl.literal("HIGH"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0005=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.001) & (unioned.gnomAD_genomes_AF < 0.005) & (
                unioned.impact.contains(hl.literal("HIGH"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0010=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.005) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("HIGH")),
                                      hl.agg.array_sum(unioned.AC)[0]),
            gnomad_0100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.01) & (
                unioned.impact.contains("HIGH")), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.10) & (unioned.impact.contains("HIGH")),
                                     hl.agg.array_sum(unioned.AC)[0])))
    return gnomad_tb


def write_gnomad_table(vcfs, dest, overwrite=False):
    files = dict()
    with open(vcfs) as vcflist:
        for vcfpath in vcflist:
            stripped = vcfpath.strip()
            vcfpath = Path(stripped)
            assert vcfpath.exists()
            prefix = file_utility.trim_prefix(vcfpath.stem)
            destination = Path(dest).joinpath(Path(stripped).stem)
            if not destination.exists():
                # Read all vcfs and make a dict, keeps in memory!
                files[prefix] = append_table(vcfs_to_matrixtable(stripped, destination.__str__(), False),
                                         out=destination.__str__(), write=True)
            else:
                files[prefix] = hl.read_matrix_table(destination.__str__())

    # Turn MatrixTables into HailTables, keyed by gene, join
    unioned_table = table_join(mts_to_table(list(files.values())))
    gnomad_tb = gnomad_table(unioned_table)
    gnomadpath = Path(dest).joinpath(Path("gnomad_tb"))
    if gnomadpath.exists():
        if not overwrite:
            raise FileExistsError(gnomadpath)
        else:
            stderr.write("WARNING: Overwrite is active. Deleting pre-existing filetree {0}\n".format(gnomadpath))
            gnomadpath.rmdir()
    else:
        gnomad_tb.write(gnomadpath.__str__())
    return gnomad_tb


if __name__ == '__main__':
    try:

        parser = argparse.ArgumentParser(prog="Annotation pipeline command-line file tool.")
        subparsers = parser.add_subparsers(title="commands", dest="command")
        findtype = subparsers.add_parser("Findtype", help="Find all specific files of a given filetype.")
        findtype.add_argument("-s", "--source", help="Directory to be searched.", action="store", type=str)
        findtype.add_argument("-d", "--directory", help="Directory to be saved to.", nargs='?', const=".",
                              action="store", type=str)  # TODO: Default value returns None type
        findtype.add_argument("-t", "--type", help="Filetype to be used.", action="store", type=str)
        findtype.add_argument("-r", "--regex", help="Filter strings by including parsed regex.", action="store",
                              type=str)
        readvcfs = subparsers.add_parser("Readvcfs", help="Turn VCF file(s) into a Hail MatrixTable.")
        readvcfs.add_argument("-f", "--file",
                              help="The VCF file to be parsed or file containing VCF pathnames on each line, LF (Unix)",
                              action="store", type=str)  # TODO: Default value returns None type
        readvcfs.add_argument("-d", "--dest", help="Destination folder to write the Hail MatrixTable files.",
                              nargs='?', const=os.path.abspath("."))
        readvcfs.add_argument("-r", "--overwrite", help="Overwrites any existing output MatrixTables, HailTables.",
                              action="store_true")

        args = parser.parse_args()
        if args.command is not None:

            if str.lower(args.command) == "findtype":
                # print(findtypes(args.directory, args.type))
                files = file_utility.find_filetype(args.source, args.type)

                file_utility.write_filelist(args.directory, "{0}.{1}.txt".format(os.path.basename(
                    os.path.normpath(args.source)), args.type), files, regex=args.regex)
                # Convert the directory into a name for the file, passing found files with
                # Regex in files matching only with a matching regex (e.g. *.vep.vcf wildcard).
                # Unique files only, duplicates written to duplicates_*.txt
            elif str.lower(args.command) == "readvcfs":
                print("Turning {0} into HailTable in directory {1}".format(args.file, args.dest))
                if Path(args.dest).joinpath(Path("gnomad_tb")).exists():
                    if not args.overwrite:
                        gnomad_tb = hl.read_table(Path(args.dest).joinpath(Path("gnomad_tb")).__str__())
                    else:
                        gnomad_tb = write_gnomad_table(args.file, args.dest, overwrite=args.overwrite)
                else:
                    gnomad_tb = write_gnomad_table(args.file, args.dest, overwrite=args.overwrite)
                gnomad_tb.describe()
                gnomad_tb.flatten().export(Path(args.dest).parent.joinpath("gnomad.tsv").__str__())

        #  Empty stdin string / command
        else:
            parser.print_usage()
            # print("Invalid input, quitting.")
            assert os.path.exists(mt_path)
            mt = hl.read_matrix_table(mt_path)

    except KeyboardInterrupt:
        print("Quitting.")
        raise

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
