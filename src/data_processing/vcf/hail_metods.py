import os
import sys
import datetime
import shutil
import hail as hl
from pathlib import Path
from src.data_processing.file_io.readers import trim_prefix

unique = hash(datetime.datetime.utcnow())

def write_gnomad_table(vcfs, dest, overwrite=False, metadata=None):
    gnomad_tb = None #TODO: define as global variable?
    hailtables = dict()
    metadata_dict = get_metadata(metadata)
    for vcfpath in vcfs:
        assert vcfpath.exists()
        prefix = trim_prefix(vcfpath.stem)
        destination = Path(dest).joinpath(vcfpath.stem)
        if overwrite or not destination.exists():
            # Read all vcfs and make a dict, keeps in memory!
            hailtables[prefix] = append_table(vcfs_to_matrixtable(vcfpath.__str__(), destination.__str__(), False),
                                              prefix, destination.__str__(), True, metadata_dict)
        elif destination.exists():
            hailtables[prefix] = hl.read_matrix_table(destination.__str__())
            sys.stderr.write("Overwrite is not active, opening existing file instead: {0}\n"
                             .format(destination.__str__()))
        else:
            FileExistsError("The output HailTable exists and --overwrite is not active in destination {0}"
                            .format(destination.__str__()))

    # Turn MatrixTables into HailTables, keyed by gene, join
    unioned_table = table_join(mts_to_table(list(hailtables.values())))

    gnomad_tb = gnomad_table(unioned_table)
    gnomadpath = Path(dest).joinpath(Path("gnomad_tb"))
    if gnomadpath.exists():
        if not overwrite:
            raise FileExistsError(gnomadpath)
        else:
            sys.stderr.write("WARNING: Overwrite is active. Deleting pre-existing directory {0}\n".format(gnomadpath))
            shutil.rmtree(gnomadpath)
    else:
        gnomad_tb.write(gnomadpath.__str__())
        pass
    return gnomad_tb


def get_metadata(metadata_path):
    p = Path(metadata_path)
    metadata_dict = dict()
    assert p.exists()
    with p.open(encoding="latin-1") as f:
        for line in f.readlines():
            s = line.strip().split("\t")
            ecode = trim_prefix(s[0])
            if ecode not in metadata_dict:
                if len(s) >= 2:
                    metadata_dict[ecode] = [s[1], s[2]]
                else:
                    metadata_dict[ecode] = ["NA", "NA"]
            else:
                sys.stderr.write("Found duplicate key {0} for line {1}. Existing object {2}.\n"
                                 .format(ecode, s, (ecode, metadata_dict[ecode])))
    return metadata_dict


def append_table(table, prefix, out=None, write=False, metadata=None):
    mt_a = table
    mt_a = mt_a.annotate_rows(VEP_str=mt_a.vep.first().split("\\|"))
    mt_a = mt_a.annotate_entries(AC=mt_a.GT.n_alt_alleles(),
                                 VF=hl.float(mt_a.AD[1] / mt_a.DP))
    mt_a = mt_a.annotate_rows(impact=mt_a.VEP_str[0],
                              gene=mt_a.VEP_str[1],
                              HGNC_ID=hl.int(parse_empty(mt_a.VEP_str[2])),
                              MAX_AF=hl.float(parse_empty(mt_a.VEP_str[3])))
    mt_a = mt_a.drop(mt_a.info)
    mt_a = mt_a.filter_entries(mt_a.VF >= 0.3, keep=True)  # Remove all not ALT_pos < 0.3 / DP > 20
    if metadata is not None:
        phen, mut = metadata.get(prefix, ["NA", "NA"])
        if len(phen) == 0: phen = "NA"
        if len(mut) == 0: mut = "NA"
        mt_a = mt_a.annotate_globals(metadata=hl.struct(phenotype=phen, mutation=mut))

    if write and out is not None:
        mt_a.write(out)
    return mt_a


def vcfs_to_matrixtable(f, destination=None, write=True, annotate=True):
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

    contig_prefix = "chr"
    contig_recoding = {f"{contig_prefix}{i}": str(i) for i in range(1, 23)}
    contig_recoding.update({"chrX": "X", "chrY": "Y"})
    table = hl.import_vcf(files, force=True, reference_genome='GRCh37', contig_recoding=contig_recoding)

    if annotate:
        table = table.filter_rows(table.alleles[1] != '*')  # These alleles break VEP, filter out star alleles.
        table = hl.methods.vep(table, config="vep_settings.json", csq=True)
    if write:
        if not os.path.exists(destination):
            table.write(destination)
        else:
            raise FileExistsError(destination)
    return table


def table_join(tables_list):
    # Join Tables into one Table.
    if tables_list is not None and len(tables_list) > 0:
        unioned = tables_list[0]  # Initialize with a single table
    else:
        raise Exception("No tables to be joined based on current configuration.")
    if len(tables_list) > 1:
        unioned = unioned.union(*tables_list[1:], unify=True)
    return unioned.cache()


def mts_to_table(tables):
    for i, tb in enumerate(tables):
        tb = tb.key_cols_by()
        tb = tb.entries()  # Convert from MatrixTable to Table
        tables[i] = tb.key_by(tb.gene)  # Key by gene
    return tables


def gnomad_table(unioned, text="modifier"):
    sys.stderr.write("Creating MAX_AF_frequency table\n")
    gnomad_tb = unioned.group_by(unioned.gene).aggregate(
        modifier=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.MAX_AF < 0.01) & (unioned.impact.contains(hl.literal("MODIFIER"))),
                hl.agg.sum(unioned.AC)),
            gnomad_1_5=hl.agg.filter((unioned.MAX_AF > 0.01) & (unioned.MAX_AF < 0.05) & (
                unioned.impact.contains(hl.literal("MODIFIER"))), hl.agg.sum(unioned.AC)),
            gnomad_5_100=hl.agg.filter((unioned.MAX_AF > 0.05) & (
                unioned.impact.contains(hl.literal("MODIFIER"))), hl.agg.sum(unioned.AC))),
        low=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.MAX_AF < 0.01) & (unioned.impact.contains(hl.literal("LOW"))),
                hl.agg.sum(unioned.AC)),
            gnomad_1_5=hl.agg.filter((unioned.MAX_AF > 0.01) & (unioned.MAX_AF < 0.05) & (
                unioned.impact.contains(hl.literal("LOW"))), hl.agg.sum(unioned.AC)),
            gnomad_5_100=hl.agg.filter((unioned.MAX_AF > 0.05) & (
                unioned.impact.contains(hl.literal("LOW"))), hl.agg.sum(unioned.AC))),
        moderate=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.MAX_AF < 0.01) & (unioned.impact.contains(hl.literal("MODERATE"))),
                hl.agg.sum(unioned.AC)),
            gnomad_1_5=hl.agg.filter((unioned.MAX_AF > 0.01) & (unioned.MAX_AF < 0.05) & (
                unioned.impact.contains(hl.literal("MODERATE"))), hl.agg.sum(unioned.AC)),
            gnomad_5_100=hl.agg.filter((unioned.MAX_AF > 0.05) & (
                unioned.impact.contains(hl.literal("MODERATE"))), hl.agg.sum(unioned.AC))),
        high=hl.struct(
            gnomad_1=hl.agg.filter(
                (unioned.MAX_AF < 0.01) & (unioned.impact.contains(hl.literal("HIGH"))),
                hl.agg.sum(unioned.AC)),
            gnomad_1_5=hl.agg.filter((unioned.MAX_AF > 0.01) & (unioned.MAX_AF < 0.05) & (
                unioned.impact.contains(hl.literal("HIGH"))), hl.agg.sum(unioned.AC)),
            gnomad_5_100=hl.agg.filter((unioned.MAX_AF > 0.05) & (
                unioned.impact.contains(hl.literal("HIGH"))), hl.agg.sum(unioned.AC)))
    )
    return gnomad_tb


def parse_empty(text):
    return hl.if_else(text == "", hl.missing(hl.tint32), hl.float(text))





def load_hailtables(dest, number, out=None, metadata=None, overwrite=False, phenotype=None):
    hailtables = dict()
    gnomadpath = Path(dest).joinpath(Path("gnomad_tb", str(unique)))
    ### TODO: Remove temporary fix
    gnomad_tb = None
    ecode_phenotype = dict()
    inverse_matches = dict()
    ###
    count = sum(1 for t in dest.iterdir())
    sys.stderr.write("{0} items in folder {1}\n".format(count, str(dest)))
    i = 0
    toolbar_width = 1 if count//10 == 0 else count//10
    # setup toolbar
    sys.stderr.write("Loading MatrixTables\n Progress\n")

    for idx, folder in enumerate(dest.iterdir(), 1):
        if folder.is_dir():
            vcfname = folder.name
            outpath = Path(out).joinpath(vcfname)
            #print(vcfname)
            if vcfname.rfind("gnomad_tb") == -1:  # Skip the folders containing the end product
                prefix = trim_prefix(vcfname)
                mt_a = hl.read_matrix_table(folder.__str__())
                #mt_a.cache()
                if metadata is not None:
                    # mt_a.write(outpath.__str__()) #Hail scripts in here fix loaded MatrixTables
                    # and outputs into new args.out
                   #mt_a.drop('phenotype')
                    #phen, mut = metadata.get(prefix, ["NA", "NA"])
                    #mt_a = mt_a.annotate_globals(metadata=hl.struct(phenotype=phen, mutation=mut))
                    #mt_a.write(outpath.__str__())
                    #mt_a.describe()

                    pass
                hailtables[prefix] = mt_a.filter_entries(mt_a.DP > 30)
                if idx//toolbar_width >= i:
                    sys.stderr.write("[{0}] Done {1}%\n".format("x"*(toolbar_width//10)*(i)+"-"*(toolbar_width//10)*(10-i), idx//toolbar_width*10))
                    i+=1

    # TODO: Slicing
    print("Read {0} HailTables".format(len(hailtables.values())))
    if number == -1:
        number = len(hailtables)
    if phenotype is not None:
        # Union HailTables with a given phenotype, thereby filtering
        sys.stderr.write("Filtering tables based on phenotype \"{0}\"\n".format(phenotype))
        for key, ht in hailtables.items():
            #ecode_phenotype[key] = hl.eval(ht.metadata.phenotype)
            #TODO: Load phenotypes once to increase speed of dict searches
            pass

        #matched_tables = list(filter(lambda k, t: t.startswith(phenotype), ecode_phenotype.items()))
        matched_tables = list(filter(lambda t: hl.eval(t.metadata.phenotype.matches(phenotype)),
                    hailtables.values()))
        if len(matched_tables) > 0:
            sys.stderr.write("Found {0} matching table(s) with given phenotype key\n".format(len(matched_tables)))
            unioned_table = table_join(mts_to_table(matched_tables))
        else:
            sys.stderr.write("NO tables matched to phenotype \"{0}\"\n".format(phenotype))
            phens = list(hl.eval(t.metadata.phenotype) for t in hailtables.values())
            raise KeyError("Phenotype keys available: {0}".format(phens))
    else:
        # Else union all tables
        unioned_table = table_join(mts_to_table(list(hailtables.values())))
        #sys.stderr.write("Writing intermediary unioned table to {0}\n".format(gnomadpath.parent.__str__() + "\gnomad_tb_unioned" + str(unique)))
        #unioned_table.write(gnomadpath.parent.__str__() + "\gnomad_tb_unioned" + str(unique))

    gnomad_tb = gnomad_table(unioned_table)
    if gnomadpath.exists():
        if not overwrite:
            raise FileExistsError(gnomadpath)
        else:
            sys.stderr.write("WARNING: Overwrite is active. Deleting pre-existing filetree {0}\n".format(gnomadpath))
            shutil.rmtree(gnomadpath)
            gnomad_tb.write(gnomadpath.__str__())
    else:
        gnomad_tb.write(gnomadpath.__str__())
        pass
    return gnomad_tb