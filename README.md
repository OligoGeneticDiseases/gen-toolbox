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

## Alternatively non-docker setup on local machine
1.	clone the project
2. install perl from https://www.perl.org/get.html
3. install VEP from https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html
4. update vep_settings.json, set VEP executable, VEP cache dir 
5. create local virtual environment by executing command python3 -m venv {path-to-your-environment}
6. activate virtual environment by using command {path-to-your-environment}/bin/activate
7. install requirements by command pip install -r requirements.txt
8.	Python3 main.py --help

# files information
1. main.py is root file to execute with your command ex. python main.py findtype -s {source dir} -d {destination dir} -t vcf 
2. CommandFactory.py this is a command creator class, you will find all respective command creation methods in this class
3. CommandHaldler.py this is command handler, you will find all commands invoke methods in this class
4. utils.py, utility file holds utility functions which gets called from commandHandler.py