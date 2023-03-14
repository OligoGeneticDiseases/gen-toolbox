// Set input and output directories
params.input_dir = './tests/'
params.output_dir = './tests/output'

// Validate VCF process
process validateVCF {
    input:
    // Take input files ending in '.vcf.gz' from the input directory
    file '*.vcf.gz' from params.input_dir

    output:
    // Output filtered files into filtered_ch
    file "${name}.vcf.gz" into filtered_ch

    script:
    """
    // Use vcftools to filter out indels and keep only variants with AC=1
    vcftools --gzvcf ${input} --remove-indels --recode --stdout | \
    grep -E "^#|AC=1" > ${name}.vcf
    // Use bgzip to compress the output
    bgzip ${name}.vcf
    // Use tabix to index the compressed output
    tabix -p vcf ${name}.vcf.gz
    """

    errorStrategy 'ignore' // Ignore errors
}

// Annotate variants process
process annotateVariants {
    input:
    // Take filtered input files from filtered_ch
    file '*.vcf.gz' from filtered_ch

    output:
    // Output annotated files into annotated_ch
    file "${name}.annotated.vcf.gz" into annotated_ch

    script:
    """
    // Start the Hail dataproc cluster
    hailctl dataproc start my_cluster
    // Submit annotation.py to the cluster with the input and output files and the GRCh38 reference genome
    hailctl dataproc submit my_cluster annotation.py \
    --input ${input} \
    --output ${name}.annotated.vcf.gz \
    --vep GRCh38
    // Stop the Hail dataproc cluster
    hailctl dataproc stop my_cluster
    """

    errorStrategy 'ignore' // Ignore errors
}

// Combine annotations process
process combineAnnotations {
    input:
    // Take annotated input files from annotated_ch
    file '*.annotated.vcf.gz' from annotated_ch.collect()

    output:
    // Output the combined annotated file to the output directory
    file "${params.output_dir}/annotated_combined.vcf.gz"

    script:
    """
    // Use bcftools to concatenate the annotated files into a single file
    bcftools concat ${input} --threads 4 -O z -o ${output}
    // Use tabix to index the compressed output
    tabix -p vcf ${output}
    """

    errorStrategy 'ignore' // Ignore errors
}

// Define workflow
workflow {
    // Create a channel of input files
    Channel.fromPath(params.input_dir)
           // Filter to keep only files ending in '.vcf.gz'
           .filter{ it.name.endsWith('.vcf.gz') }
           // Split each file into its individual lines
           .set{ vcf_ch }
           .splitCsv()
           // Filter each line to keep only those starting with '#' or containing 'AC=1'
           .map{ it.filter { line -> line.startsWith('#') || line.contains("AC=1") } }
           // Annotate the variants using the Hail dataproc cluster
           .annotate(sorted: true)
           // Set the output to filtered_ch
           .set{ filtered_ch }

    // Retry processes if they fail
    validateVCF.filterMeta { isRetry() }
    annotateVariants.filterMeta { isRetry() }
    combineAnnotations.filterMeta { isRetry() }

    // View the output of the combineAnnotations process
    combineAnnotations.view()
}