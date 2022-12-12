# Oligogenicity
The concept of this project is to collate large numbers of VCF files of unique samples using [Hail](https://hail.is/) from Broad Institute [[github]](https://github.com/hail-is/hail), 
annotate variants using [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and create variant frequency tables in a pipeline fashion. 
The tool is set up to run locally on a large server (e.g., 64 core CPU with 2GB of RAM per core). Input VCF formats may differ between generations, but the pipeline is dependent on correct INFO fields.
Input VCF files can be validated using the VCF-validator from the [vcftools](https://vcftools.github.io/index.html) toolset.
An additional statistical module enables analysis of the frequency tables to reveal potential relations between the variant burden within gene sets and correlating phenotypes – a form of variant enrichment analysis using Monte-Cristo permutation analysis. 
The project is designed to be dockerized. 
The project is directed by Tartu University Hospital Centre of Medical Genetics / Tartu University Institute of Clinical Medicine. 

## Setup
1.	Clone the project into a local folder
2.	Insert the VEP annotator perl script location to the config file
3.	Docker image build -t oligogenicity:latest .
4.	Docker run -v {path to VEP}:{path in config} –help
Alternatively 
1.	clone the project
2.	Insert the VEP annotator perl script location to the config file
3.	Activate a local python 3 venv env
4.	Pip install -t requirements.txt
5.	Activate the venv environment
6.	Python3 main.py --help

