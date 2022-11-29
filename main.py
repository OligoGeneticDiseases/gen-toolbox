import argparse
import os.path
import shutil
import sys
from pathlib import Path
from sys import stderr

import hail as hl
from pyspark import *

import file_utility

hail_home = Path(hl.__file__).parent.__str__()


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


def append_table(table, prefix, out=None, write=False, metadata=None):
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
    if metadata is not None:
        mt_a = mt_a.annotate_globals(phenotype=metadata.get(prefix, "NA"))

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
    if tables_list is not None and len(tables_list) > 0:
        unioned = tables_list[0]  # Initialize with a single table
    else:
        raise Exception("No tables to be joined based on current configuration.")
    if len(tables_list) > 1:
        unioned = unioned.union(*tables_list[1:])
    return unioned


def get_metadata(metadata_path):
    p = Path(metadata_path)
    metadata_dict = dict()
    assert p.exists()
    with p.open() as f:
        for line in f.readlines():
            s = line.strip().split("\t")
            ecode = file_utility.trim_prefix(s[0])
            if ecode not in metadata_dict:
                if len(s) >= 2:
                    metadata_dict[ecode] = s[1]
                else:
                    metadata_dict[ecode] = "NA"
            else:
                sys.stderr.write("Found duplicate key {0} for line {1}. Existing object {2}\n"
                                 .format(ecode, s, (ecode, metadata_dict[ecode])))
    return metadata_dict


def gnomad_table(unioned):
    gnomad_tb = unioned.group_by(unioned.gene).aggregate(
        modifier=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.01) & (unioned.impact.contains(hl.literal("MODIFIER"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_1_5=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.05) & (
                unioned.impact.contains(hl.literal("MODIFIER"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_5_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.05) & (
                unioned.impact.contains("MODIFIER")), hl.agg.array_sum(unioned.AC)[0])),
        low=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.01) & (unioned.impact.contains(hl.literal("LOW"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_1_5=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.05) & (
                unioned.impact.contains(hl.literal("LOW"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_5_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.05) & (
                unioned.impact.contains("LOW")), hl.agg.array_sum(unioned.AC)[0])),
        moderate=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.01) & (unioned.impact.contains(hl.literal("MODERATE"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_1_5=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.05) & (
                unioned.impact.contains(hl.literal("MODERATE"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_5_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.05) & (
                unioned.impact.contains("MODERATE")), hl.agg.array_sum(unioned.AC)[0])),
        high=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.gnomAD_genomes_AF < 0.01) & (unioned.impact.contains(hl.literal("HIGH"))),
                hl.agg.array_sum(unioned.AC)[0]),
            gnomad_1_5=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.01) & (unioned.gnomAD_genomes_AF < 0.05) & (
                unioned.impact.contains(hl.literal("HIGH"))), hl.agg.array_sum(unioned.AC)[0]),
            gnomad_5_100=hl.agg.filter((unioned.gnomAD_genomes_AF > 0.05) & (
                unioned.impact.contains("HIGH")), hl.agg.array_sum(unioned.AC)[0])))
    return gnomad_tb


def write_gnomad_table(vcfs, dest, overwrite=False, metadata=None):
    hailtables = dict()
    metadata_dict = get_metadata(metadata)
    with open(vcfs) as vcflist:
        for vcfpath in vcflist:
            stripped = vcfpath.strip()
            vcfpath = Path(stripped)
            assert vcfpath.exists()
            prefix = file_utility.trim_prefix(vcfpath.stem)
            destination = Path(dest).joinpath(Path(stripped).stem)
            if overwrite or not destination.exists():
                # Read all vcfs and make a dict, keeps in memory!
                hailtables[prefix] = append_table(vcfs_to_matrixtable(stripped, destination.__str__(), False),
                                                  prefix, destination.__str__(), True, metadata_dict)
            elif destination.exists():
                hailtables[prefix] = hl.read_matrix_table(destination.__str__())

    # Turn MatrixTables into HailTables, keyed by gene, join
    unioned_table = table_join(mts_to_table(list(hailtables.values())))
    gnomad_tb = gnomad_table(unioned_table)
    gnomadpath = Path(dest).joinpath(Path("gnomad_tb"))
    if gnomadpath.exists():
        if not overwrite:
            raise FileExistsError(gnomadpath)
        else:
            stderr.write("WARNING: Overwrite is active. Deleting pre-existing directory {0}\n".format(gnomadpath))
            shutil.rmtree(gnomadpath)
    else:
        gnomad_tb.write(gnomadpath.__str__())
    return gnomad_tb


def load_hailtables(dest, number, overwrite=False, phenotype=None):
    hailtables = dict()
    gnomadpath = Path(dest).joinpath(Path("gnomad_tb"))

    for folder in dest.iterdir():
        if folder.is_dir():
            vcfname = folder.name
            print(vcfname)
            if vcfname != "gnomad_tb":  # Skip the folder containing the end product
                prefix = file_utility.trim_prefix(vcfname)
                ht = hl.read_matrix_table(folder.__str__())
                hailtables[prefix] = ht
    # TODO: Slicing
    if number == -1:
        number = len(hailtables)
    if phenotype is not None:
        # Union HailTables with a given phenotype, thereby filtering
        print("Filtering tables based on phenotype {0}".format(phenotype))
        matched_tables = list(filter(lambda t: hl.eval(t.phenotype.matches(phenotype, True)),
                                     hailtables.values()))
        if len(matched_tables) > 0:
            sys.stderr.write("Found {0} matching table(s) with given phenotype key.\n".format(len(matched_tables)))
            unioned_table = table_join(mts_to_table(matched_tables))
        else:
            sys.stderr.write("NO tables matched to phenotype \"{0}\".\n".format(phenotype))
            phens = list(hl.eval(t.phenotype) for t in hailtables.values())
            raise KeyError("Phenotype keys available: {0}".format(phens))
    else:
        # Else union all tables
        unioned_table = table_join(mts_to_table(list(hailtables.values())))

    gnomad_tb = gnomad_table(unioned_table)
    if gnomadpath.exists():
        if not overwrite:
            raise FileExistsError(gnomadpath)
        else:
            stderr.write("WARNING: Overwrite is active. Deleting pre-existing filetree {0}\n".format(gnomadpath))
            shutil.rmtree(gnomadpath)
            gnomad_tb.write(gnomadpath.__str__())
    else:
        gnomad_tb.write(gnomadpath.__str__())
    return gnomad_tb


if __name__ == '__main__':
    try:

        parser = argparse.ArgumentParser(prog="Gnomad frequency table burden analysis pipeline command-line tool using "
                                              "Hail\n{0}".format(hl.cite_hail()))
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
        readvcfs.add_argument("-g", "--globals", help="Tab delimited input file containing globals string "
                                                      "for a given unique sample.", action="store", type=str)
        loaddb = subparsers.add_parser("Loaddb", help="Load a folder containing HailTables.")
        loaddb.add_argument("-d", "--directory", help="Folder to load the Hail MatrixTable files from.",
                            nargs='?', const=os.path.abspath("."))
        loaddb.add_argument("-r", "--overwrite", help="Overwrites any existing output MatrixTables, HailTables.",
                            action="store_true")
        loaddb.add_argument("-o", "--out", help="Output destination.", action="store", type=str)
        loaddb.add_argument("-n", "--number", help="Number of tables to be collated.", nargs="?",
                            type=int,
                            default=-1)
        loaddb.add_argument("-g", "--globals", help="Tab delimited input file containing globals string "
                                                    "for a given unique sample.", action="store", type=str)
        loaddb.add_argument("--phenotype", help="Filter a subset of samples with given phenotype.", action="store",
                            type=str)

        args = parser.parse_args()
        if args.command is not None:

            if str.lower(args.command) == "findtype":
                # print(findtypes(args.directory, args.type))
                files = file_utility.find_filetype(args.source, args.type, verbose=False)

                file_utility.write_filelist(args.directory, "{0}.{1}.txt".format(os.path.basename(
                    os.path.normpath(args.source)), args.type), files, regex=args.regex)
                # Convert the directory into a name for the file, passing found files with
                # Regex in files matching only with a matching regex (e.g. *.vep.vcf wildcard).
                # Unique files only, duplicates written to duplicates_*.txt

            else:
                conf = SparkConf()
                conf.set('spark.sql.files.maxPartitionBytes', '60000000000')
                conf.set('spark.sql.files.openCostInBytes', '60000000000')
                conf.set('spark.submit.deployMode', u'client')
                conf.set('spark.app.name', u'HailTools-TSHC')
                conf.set('spark.executor.memory', "4g")
                conf.set("spark.jars", "{0}/backend/hail-all-spark.jar".format(hail_home))
                conf.set("spark.executor.extraClassPath", "./hail-all-spark.jar")
                conf.set("spark.driver.extraClassPath", "{0}/backend/hail-all-spark.jar".format(hail_home))
                conf.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
                conf.set("spark.kryo.registrator", "is.hail.kryo.HailKryoRegistrator")
                sc = SparkContext(conf=conf)
                hl.init(backend="spark", sc=sc)
                if str.lower(args.command) == "readvcfs":

                    if Path(args.dest).joinpath(Path("gnomad_tb")).exists():
                        if not args.overwrite:
                            gnomad_tb = hl.read_table(Path(args.dest).joinpath(Path("gnomad_tb")).__str__())
                        else:
                            gnomad_tb = write_gnomad_table(args.file, args.dest, overwrite=args.overwrite,
                                                           metadata=args.globals)
                    else:
                        gnomad_tb = write_gnomad_table(args.file, args.dest, overwrite=args.overwrite,
                                                       metadata=args.globals)
                    gnomad_tb.describe()
                    gnomad_tb.flatten().export(Path(args.dest).parent.joinpath("gnomad.tsv").__str__())
                elif str.lower(args.command) == "loaddb":
                    if args.globals is not None:
                        metadata_dict = get_metadata(args.globals)

                    dirpath = Path(args.directory)
                    gnomad_tb = load_hailtables(dirpath, args.number, args.overwrite, args.phenotype)
                    gnomad_tb.describe()
                    gnomad_tb.flatten().export(Path(args.out).parent.joinpath("{0}_gnomad.tsv".format(args.phenotype)).__str__())
                    input("Waiting to exit. Press any key.")

        #  Empty stdin string / command
        else:
            parser.print_usage()
            # print("Invalid input, quitting.")
            # assert os.path.exists(mt_path)
            # mt = hl.read_matrix_table(mt_path)

    except KeyboardInterrupt:
        print("Quitting.")
        raise

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
