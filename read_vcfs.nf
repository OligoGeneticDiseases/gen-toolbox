nextflow.enable.dsl=2

params.vcf_files = "*"
params.destination = "output"
params.globals = ""
params.python_script_path = "$PWD/main.py"

vcfsCh = Channel.fromPath("$params.vcf_files/*", checkIfExists: true )
                      | map { file -> return tuple(file.baseName, file) }

//vcfsCh = Channel.fromPath(params.vcf_files, checkIfExists: true ).splitText(){ it -> [basename(it), it] }

process readVCFs {
    conda "$PWD/environment.yml"
    publishDir("$params.destination/$name", mode: 'copy')
    tag { batch }
    input:
	tuple val(name), path(vcf_path)
	
	output:
    path("*")

    script:
    """
	python ${params.python_script_path} annotatevcfs --file $vcf_path --dest ${params.destination} -g ${params.globals}
    """
}

workflow {
    vcfsCh | view()
    readVCFs(vcfsCh)
}
