import pathlib

import hail as hl

from src.data_processing.vcf.hail_metods import parse_empty
from src.utils.general.string_operations import trim_prefix


# Content: Functions related to reading various file formats.


def load_db_batch(mts):
    batch = []
    for mt_path in mts:
        batch.append(load_mt(mt_path.__str__()))
    return batch


def load_mt(mt_path):
    mt = hl.read_matrix_table(mt_path)
    hl.utils.info(f"Read MatrixTable {mt_path}")
    return mt


def import_and_annotate_vcf_batch(vcfs, metadata=None, annotate=True, interval=None, location=None):
    batch = []
    for vcf in vcfs:
        mt = import_and_annotate_vcf(
                vcf, metadata=metadata, annotate=annotate, interval=interval)
        mt_path = pathlib.Path(location).joinpath(f"{hl.eval(mt.prefix)}.mtx").as_posix()
        #mt.write(mt_path)
        #batch.append(load_mt(mt_path))
        batch.append(mt)
    return batch


def import_and_annotate_vcf(vcf_path, metadata=None, annotate=True, interval=None):
    """
    Import a VCF file and annotate it using hail.import_vcf() and hail.VEP().

    :param vcf_path: The path to the VCF file.
    :param annotate: Annotates the input VCF file using VEP (default=True). Set to false to skip annotation
    (if already annotated with VEP)
    :param interval hl.LocusInterval to import only certain genomic interval
    :return: Annotated MatrixTable.
    """
    prefix = trim_prefix(vcf_path.stem)
    contig_prefix = "chr"
    contig_recoding = {f"{contig_prefix}{i}": str(i) for i in range(1, 23)}
    contig_recoding.update({"chrX": "X", "chrY": "Y", "chrM": "MT"})

    mt = hl.import_vcf(
        vcf_path.__str__(), reference_genome=hl.default_reference, contig_recoding=contig_recoding
    )
    # mt = hl.filter_alleles(mt, lambda allele, i: hl.is_star(mt.alleles[0], allele))
    # updated_info = mt.info.annotate(AC=mt.new_to_old.map(lambda i: mt.info.AC[i - 1]))
    # ds_result = mt.annotate_rows(info=updated_info)
    # mt = ds_result
    if interval is not None:
        mt = hl.filter_intervals(
            mt,
            [
                hl.parse_locus_interval(
                    x,
                )
                for x in interval
            ],
        )
    mt = mt.filter_rows(mt.locus.contig != "MT") # Filter mitochondrial
    mt = mt.filter_rows(mt.alleles[1] != "*")  # Filter star alleles as these break VEP

    if annotate:
        mt = hl.vep(mt, "./src/config/vep_settings.json")
        mt = mt.annotate_rows(
            impact=mt.vep.IMPACT,
            gene=mt.vep.SYMBOL,
            HGNC_ID=mt.vep.HGNC_ID,
            MAX_AF=mt.vep.MAX_AF,
        )
    else:  # get the data from the CSQ string
        mt = mt.annotate_rows(
            vep=mt.info.CSQ.first().split("\\|")
        )  # Convert CSQ string into the expected VEP output

        mt = mt.annotate_rows(
            # For self-annotated the indexes might be very different in a CSQ string, i.e. for DRAGEN output it might be
            # 2, 3, 22, 23
            impact=mt.vep[0],
            gene=mt.vep[1],
            HGNC_ID=hl.int(parse_empty(mt.vep[2])),
            MAX_AF=hl.float(parse_empty(mt.vep[3])),
        )
    mt = mt.filter_rows(mt.MAX_AF < 0.01) # Filter out common variants (>1%)
    mt = mt.annotate_entries(AC=mt.GT.n_alt_alleles(), VF=hl.float(mt.AD[1] / mt.DP))
    if metadata is not None:
        phen, mut = metadata.get(prefix, ["NA", "NA"])
        if len(phen) == 0:
            phen = "NA"
        if len(mut) == 0:
            mut = "NA"
        mt.annotate_globals(
            metadata=hl.struct(phenotype=phen, mutation=mut)
        )  # Annotate all rows with corresponding meta
    mt = mt.drop(mt.vep)  # Drop now duplicated field
    mt = mt.drop(mt.info)
    mt = mt.filter_entries(
        mt.VF >= 0.3, keep=True
    )  # Remove all not ALT_pos < 0.3 / DP > 30
    mt = mt.filter_entries(mt.DP > 30, keep=True)
    mt = mt.annotate_globals(prefix=prefix)

    return mt
