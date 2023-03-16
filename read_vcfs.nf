params.vcf_files = ""
params.destination = "output"
params.globals = ""

if (!params.vcf_files || !params.destination || !params.globals) {
    exit 1, "Error: VCF files, destination directory, and globals file must be provided. Use --vcf_files {path to VCF files} --destination {destination dir} --globals {path to globals file}"
}

process readVCFs {
    input:
    val vcf_files from params.vcf_files
    val destination from params.destination
    val globals from params.globals

    output:
    path("${destination}") into result

    script:
    """
    python3 /Users/markus/gen-toolbox/main.py readvcfs -f ${vcf_files} -d ${destination} -g ${globals}
    """
}

result.view { outputDir ->
    println "Pipeline completed successfully. Output saved to ${outputDir}"
}
