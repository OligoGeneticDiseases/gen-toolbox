import hail as hl

from src.data_processing.vcf.hail_metods import parse_empty
from src.utils.general.string_operations import trim_prefix

# Content: Functions related to reading various file formats.

print("read.py is starting......")

def load_db_batch(mts):
    print("load_db_batch is starting....")
    batch = []
    for mt_path in mts:
        batch.append(load_mt(mt_path.__str__()))
    return batch


def load_mt(mt_path):
    print("mt_path starting...")
    return hl.read_matrix_table(mt_path)


def import_and_annotate_vcf_batch_read(
    vcfs, metadata=None, annotate=True, interval=None
):
    print("import_and_annotate_vcf_batch_read() is starting..")
    # function for avoiding conflicts with columns merging in check command_factory.check_relatedness command
    batch = []
    for vcf in vcfs:
        batch.append(
            import_and_annotate_vcf(
                vcf, metadata=metadata, annotate=annotate, interval=interval
            )
        )
    row_keys = [set(i.row.info.keys()) for i in batch]
    entry_keys = [set(i.entry._fields.keys()) for i in batch]
    row_common = list(set.intersection(*row_keys))
    entry_common = list(set.intersection(*entry_keys))
    fields_to_convert = {
        "GQ": hl.tfloat64,
        "PL": hl.tarray(hl.tint32),
        "AD": hl.tarray(hl.tint32),
        "GT": hl.tcall,
        "DP": hl.tint32,
    }
    batch = [mt.annotate_rows(info=mt.row.info.select(*row_common)) for mt in batch]
    batch = [mt.select_entries(*entry_common) for mt in batch]
    batch = [
        mt.annotate_entries(
            **{
                field: hl.or_else(mt[field], hl.missing(fields_to_convert[field]))
                for field in fields_to_convert
                if field in entry_common
            }
        )
        for mt in batch
    ]
    mt = hl.MatrixTable.union_rows(*batch, _check_cols=False)
    mt = annotation(mt, metadata=metadata, annotate=annotate)
    return mt


def annotation(mt, metadata=None, annotate=True):
    print("annotation is starting .....")
    if annotate:
        mt = hl.vep(mt, "./src/config/vep_settings.json")
        print("prindime skeemi siin välja:")
        mt.describe()
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
            impact=mt.vep[0],
            gene=mt.vep[1],
            HGNC_ID=hl.int(parse_empty(mt.vep[2])),
            MAX_AF=hl.float(parse_empty(mt.vep[3])),
        )
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
    mt = mt.filter_entries(
        mt.VF >= 0.3, keep=True
    )  # Remove all not ALT_pos < 0.3 / DP > 30
    mt.filter_entries(mt.DP > 30, keep=True)
    # mt = mt.select_entries(mt.AD,mt.DP,mt.GT,mt.VF,mt.AC,)
    # mt = mt.select_rows(mt.impact, mt.gene, mt.HGNC_ID, mt.MAX_AF)
    return mt


def import_vcf(vcf_path, metadata=None, annotate=True, interval=None):
    print("import vcf is starting....")
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
    contig_recoding.update({"chrX": "X", "chrY": "Y"})
    mt = hl.import_vcf(
        vcf_path.__str__(), reference_genome="GRCh37", contig_recoding=contig_recoding
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
    # mt.alleles[1] != "*" have been replaced with (hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*")
    # hl.len(mt.alleles) > 1) will filter out nucleotide bases without mutation
    # No need to handle string characters other than "*" like "." as hl.import_vcf does not read them in records
    mt = mt.filter_rows(
        (hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*")
    )  # Filter star alleles as these break VEP

    return mt


def import_and_annotate_vcf_batch(vcfs, metadata=None, annotate=True, interval=None):
    print("import_and_annotate_vcf_batch is starting ....")

    batch = []
    for vcf in vcfs:
        batch.append(
            import_and_annotate_vcf(
                vcf, metadata=metadata, annotate=annotate, interval=interval
            )
        )
    return batch


def import_and_annotate_vcf(vcf_path, metadata=None, prefix=None, annotate=True, interval=None):
    print("import_and_annotate_vcf is starting ....")
    """
    Import a VCF file and annotate it using hail.import_vcf() and hail.VEP().

    :param vcf_path: The path to the VCF file.
    :param annotate: Annotates the input VCF file using VEP (default=True). Set to false to skip annotation.
    :param interval hl.LocusInterval to import only certain genomic interval.
    :return: Annotated MatrixTable.
    """
    # Use filename stem as prefix for metadata, if necessary
    prefix = trim_prefix(vcf_path.stem) if 'stem' in dir(vcf_path) else None
    contig_prefix = "chr"
    contig_recoding = {f"{contig_prefix}{i}": str(i) for i in range(1, 23)}
    contig_recoding.update({"chrX": "X", "chrY": "Y"})

    # Import VCF with contig recoding
    mt = hl.import_vcf(
        vcf_path.__str__(),
        reference_genome="GRCh37", 
        contig_recoding=contig_recoding
    )
    
    
    # filter out non-standard contigs
    standard_contigs = {f"{contig_prefix}{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    mt = mt.filter_rows(hl.literal(standard_contigs).contains(mt.locus.contig))

    # Handle intervals, if specified
    if interval is not None:
        mt = hl.filter_intervals(
            mt,
            [hl.parse_locus_interval(x) for x in interval],
        )
    
    # Filter alleles and apply annotation as needed
    mt = mt.filter_rows(mt.alleles[1] != "*")  # Filter star alleles as these break VEP
    if annotate:
        mt = hl.vep(mt, "./src/config/vep_settings.json")
        mt.describe()
        mt = mt.annotate_rows(
            impact=mt.vep.IMPACT,
            gene=mt.vep.SYMBOL,
            HGNC_ID=mt.vep.HGNC_ID,
            MAX_AF=mt.vep.MAX_AF,
        )
    else:
        mt = mt.annotate_rows(
            vep=mt.info.CSQ.first().split("\\|")
        )
        mt = mt.annotate_rows(
            impact=mt.vep[0],
            gene=mt.vep[1],
            HGNC_ID=hl.int(parse_empty(mt.vep[2])),
            MAX_AF=hl.float(parse_empty(mt.vep[3])),
        )
    
    mt = mt.annotate_entries(AC=mt.GT.n_alt_alleles(), VF=hl.float(mt.AD[1] / mt.DP))
    
    # Handle metadata if present
    if metadata is not None and prefix:
        phen, mut = metadata.get(prefix, ["NA", "NA"])
        phen = phen if phen else "NA"
        mut = mut if mut else "NA"
        mt = mt.annotate_globals(metadata=hl.struct(phenotype=phen, mutation=mut))
    
    mt = mt.drop(mt.vep)  # Drop duplicated VEP fields
    mt = mt.filter_entries(mt.VF >= 0.3)  # Filter by VF threshold
    mt = mt.filter_entries(mt.DP > 30, keep=True)

    return mt


def import_and_annotate_gvcf_batch(gvcfs, metadata=None, annotate=True, interval=None):
    print("import_and_annotate_gvcf_batch is starting....")
    batch = []
    for gvcf in gvcfs:
        prefix = trim_prefix(gvcf.stem) if hasattr(gvcf, 'stem') else None
        batch.append(
            import_and_annotate_gvcf(
                gvcf, metadata=metadata, prefix=prefix, annotate=annotate, interval=interval
            )
        )

    merged_mt = hl.MatrixTable.union_rows(*batch, _check_cols=False)
    return merged_mt


def filter_empty_rows(mt):
    """
    Filter out empty rows in MatrixTable if no relevant annotations are found.
    """
    print("Filtering empty rows...")
    mt = mt.filter_rows(
        hl.is_defined(mt.locus) & (hl.len(mt.alleles) > 1)
    )  # Ensure at least one allele and defined locus
    return mt


def annotate_gvcf(mt, annotate=True):
    """
    Annotate the GVCF MatrixTable. Use VEP if available; otherwise, use fallback CSQ-based annotation.
    """
    print("Annotating GVCF...")
    if annotate:
        try:
            mt = hl.vep(mt, "./src/config/vep_settings.json")
            mt = mt.annotate_rows(
                impact=hl.or_missing(hl.is_defined(mt.vep), mt.vep.IMPACT),
                gene=hl.or_missing(hl.is_defined(mt.vep), mt.vep.SYMBOL),
                HGNC_ID=hl.or_missing(hl.is_defined(mt.vep), mt.vep.HGNC_ID),
                MAX_AF=hl.or_missing(hl.is_defined(mt.vep), mt.vep.MAX_AF),
            )
            print("VEP annotation completed.")
        except Exception as e:
            print(f"VEP annotation failed: {e}. Attempting fallback annotation.")
            if 'CSQ' in mt.info:
                mt = mt.annotate_rows(
                    vep=hl.or_missing(hl.is_defined(mt.info.CSQ), mt.info.CSQ.first().split("\\|"))
                )
                mt = mt.annotate_rows(
                    impact=hl.or_missing(hl.len(mt.vep) > 0, mt.vep[0]),
                    gene=hl.or_missing(hl.len(mt.vep) > 1, mt.vep[1]),
                    HGNC_ID=hl.int(hl.or_missing(hl.len(mt.vep) > 2, mt.vep[2])),
                    MAX_AF=hl.float(hl.or_missing(hl.len(mt.vep) > 3, mt.vep[3])),
                )
            else:
                print("No valid CSQ field found for fallback annotation.")

    # Apply entry annotations
    mt = mt.annotate_entries(
        AC=hl.or_missing(hl.is_defined(mt.GT), mt.GT.n_alt_alleles()),
        VF=hl.if_else(
            hl.is_defined(mt.AD) & (mt.DP > 0),
            hl.float(mt.AD[1] / mt.DP),
            hl.missing(hl.tfloat)
        )
    )
    return mt


def import_and_annotate_gvcf(gvcf_path, metadata=None, prefix=None, annotate=True, interval=None):
    """
    Import a GVCF file, filter non-standard contigs and alleles, and annotate using Hail's VEP or CSQ-based annotation.

    :param gvcf_path: The path to the GVCF file.
    :param metadata: Dictionary with metadata to annotate globally.
    :param annotate: Annotate using VEP (default=True).
    :param interval: List of intervals to filter the data within specific genomic regions.
    :return: Annotated MatrixTable.
    """
    prefix = trim_prefix(gvcf_path.stem)

    # Define contig recoding and allowed contigs
    contig_prefix = "chr"
    contig_recoding = {f"{contig_prefix}{i}": str(i) for i in range(1, 23)}
    contig_recoding.update({"chrX": "X", "chrY": "Y", "chrM": "MT"})

    # Import the GVCF file with contig recoding to match GRCh37 format
    mt = hl.import_vcf(
        gvcf_path.__str__(),
        reference_genome="GRCh37",
        contig_recoding=contig_recoding,
    )
    print("MatrixTable imported with contig recoding.")

    # Debug: Describe imported MatrixTable before filtering
    mt.describe()

    # Define allowed contigs as standard chromosomes
    primary_contigs = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                       "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"}

    # Filter out non-standard contigs based on mt.locus.contig
    mt = mt.filter_rows(
        (hl.literal(primary_contigs).contains(mt.locus.contig)) &
        (hl.len(mt.alleles) > 1) &
        (mt.alleles[1] != "*")
    )
    print("Filtered for primary chromosomes and non-star alleles.")

    # Debug: Describe MatrixTable after filtering
    mt.describe()

    # Annotate VF if not present
    if 'VF' not in mt.entry:
        mt = mt.annotate_entries(
            VF=hl.if_else(
                hl.is_defined(mt.AD) & (mt.DP > 0),
                hl.float(mt.AD[1] / mt.DP),
                hl.missing(hl.tfloat)
            )
        )
        print("Annotated VF entry field based on AD and DP.")

    # Proceed with annotation if enabled
    if annotate:
        try:
            mt = hl.vep(mt, "./src/config/vep_settings.json")
            mt = mt.annotate_rows(
                impact=hl.or_missing(hl.is_defined(mt.vep), mt.vep.IMPACT),
                gene=hl.or_missing(hl.is_defined(mt.vep), mt.vep.SYMBOL),
                HGNC_ID=hl.or_missing(hl.is_defined(mt.vep), mt.vep.HGNC_ID),
                MAX_AF=hl.or_missing(hl.is_defined(mt.vep), mt.vep.MAX_AF),
            )
        except Exception as e:
            print(f"VEP annotation failed: {e}. Attempting fallback annotation.")
            if 'CSQ' in mt.row:
                mt = mt.annotate_rows(
                    vep=hl.or_missing(hl.is_defined(mt.info.CSQ), mt.info.CSQ.first().split("\\|"))
                )
                mt = mt.annotate_rows(
                    impact=hl.or_missing(hl.len(mt.vep) > 0, mt.vep[0]),
                    gene=hl.or_missing(hl.len(mt.vep) > 1, mt.vep[1]),
                    HGNC_ID=hl.int(hl.or_missing(hl.len(mt.vep) > 2, mt.vep[2])),
                    MAX_AF=hl.float(hl.or_missing(hl.len(mt.vep) > 3, mt.vep[3])),
                )
            else:
                print("No valid CSQ field found for fallback annotation.")

    # Apply VF-based filtering if VF exists
    if 'VF' in mt.entry:
        mt = mt.filter_entries((mt.VF >= 0.3) & (mt.DP > 30), keep=True)
    else:
        print("Warning: VF field not found in MatrixTable. Skipping VF-based filtering.")

    return mt