# Oligogenicity
The gen-toolbox-dev project is a comprehensive tool for collating large numbers of VCF files of unique samples, annotating variants, and creating variant frequency tables in a pipeline fashion. The project is designed to run locally on a large server and is directed by Tartu University Hospital Centre of Medical Genetics / Tartu University Institute of Clinical Medicine.

This repository contains tools for genomic data processing and analysis. It provides a command-line interface (CLI) for executing various commands related to genomic data.

## Concept

The concept of this project is to collate large numbers of VCF files of unique samples using [Hail](https://hail.is/) from Broad Institute [GitHub](https://github.com/hail-is/hail), annotate variants using [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), and create variant frequency tables in a pipeline fashion. The tool is set up to run locally on a large server (e.g., 64 core CPU with 2GB of RAM per core). Input VCF formats may differ between generations, but the pipeline is dependent on correct INFO fields. Input VCF files can be validated using the VCF-validator from the [vcftools](https://vcftools.github.io/index.html) toolset. An additional statistical module enables analysis of the frequency tables to reveal potential relations between the variant burden within gene sets and correlating phenotypes – a form of variant enrichment analysis using Monte Carlo permutation analysis. The project is designed to be dockerized.

## Libraries Used

- Hail: A Python library for scalable genomic data analysis.
- Variant Effect Predictor (VEP): A tool for annotating and predicting the effects of genetic variants.


## Libraries Used

- Hail: A Python library for scalable genomic data analysis.
- Variant Effect Predictor (VEP): A tool for annotating and predicting the effects of genetic variants.

## File Tree

```
├── main.py
└── src
    ├── cli
    │   ├── command_factory.py
    │   └── command_methods.py
    ├── config
    │   ├── gene_config.json
    │   ├── settings.py
    │   └── vep_settings.json
    ├── data_processing
    │   ├── file_io
    │   │   ├── readers.py
    │   │   └── writers.py
    │   ├── hail
    │   │   └── genomic_operations.py
    │   ├── pca
    │   │   └── analysis.py
    │   └── vcf
    │       ├── hail_metods.py
    │       ├── read.py
    │       └── transform.py
    ├── stats.py
    └── utils
        ├── file
        │   ├── file_meta.py
        │   └── file_search.py
        └── general
            ├── data_manipulation.py
            └── string_operations.py
```

## File Descriptions

- `main.py`: Main entry point of the application. Initializes the CLI and executes the appropriate command based on user input.
- `src/cli/command_factory.py`: Defines the `CommandFactory` class for creating and executing CLI commands.
- `src/cli/command_methods.py`: Defines methods for available CLI commands.
- `src/config/gene_config.json`: Contains configuration settings for genes.
- `src/config/settings.py`: Defines the `Settings` class for managing configuration settings.
- `src/config/vep_settings.json`: Contains configuration settings for VEP.
- `src/data_processing/file_io/readers.py`: Defines classes for reading data from various file formats.
- `src/data_processing/file_io/writers.py`: Defines classes for writing data to various file formats.
- `src/data_processing/hail/genomic_operations.py`: Defines functions for genomic operations using Hail.
- `src/data_processing/pca/analysis.py`: Defines functions for PCA on genomic data.
- `src/data_processing/vcf/hail_metods.py`: Defines functions for working with VCF files using Hail.
- `src/data_processing/vcf/read.py`: Defines functions for reading VCF files.
- `src/data_processing/vcf/transform.py`: Defines functions for transforming VCF files.
- `src/stats.py`: Defines functions for statistical analysis on genomic data.
- `src/utils/file/file_meta.py`: Defines functions for working with file metadata.
- `src/utils/file/file_search.py`: Defines functions for searching for files.
- `src/utils/general/data_manipulation.py`: Defines functions for data manipulation.
- `src/utils/general/string_operations.py`: Defines functions for string operations.

## Setup
1.	Clone the project into a local folder
2.	Insert the VEP annotator perl script location to the config file
3. `docker image build -t oligogenicity:latest`
4. `docker run -v {path to VEP}:{path in config} –help`

## Alternatively non-docker setup on local machine
1. Clone the project
2. Install perl from https://www.perl.org/get.html
3. Install VEP from https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html
4. Update vep_settings.json, set VEP executable, VEP cache dir 
5. Create local virtual environment by executing command `python3 -m venv {path-to-your-environment} `
6. Activate virtual environment by using command `{path-to-your-environment}/bin/activate `
7. Install requirements by command `pip install -r requirements.txt`
8. Run: `python3 main.py --help`