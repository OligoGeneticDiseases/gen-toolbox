params.vcf_file = "/Users/markus/gen-toolbox/tests/in_silico_sorted.vcf"
params.destination = "output"
params.type = "vcf"

if (!params.vcf_file || !params.destination) {
    exit 1, "Error: VCF file and destination directory must be provided. Use --vcf_file {path to VCF file} --destination {destination dir}"
}

process findType {
    input:
    val vcf_file from params.vcf_file
    val destination from params.destination
    val type from params.type

    output:
    path("${destination}") into result

    script:
    """
    python3 /Users/markus/gen-toolbox/main.py findtype -s ${vcf_file} -d ${destination} -t ${type}
    """
}

result.view { outputDir ->
    println "Pipeline completed successfully. Output saved to ${outputDir}"
}
