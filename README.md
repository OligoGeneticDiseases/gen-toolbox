# Oligogenicity Analysis Toolbox

This pipeline provides a comprehensive suite of tools for genomic variant burden analysis focused on oligogenic interactions. The project uses the [Hail](https://hail.is/) framework from Broad Institute for scalable genomic data processing, [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) for variant annotation, and Monte Carlo simulation methods for statistical analysis of rare variant enrichment in gene sets. 

This repository contains tools for variant data processing and analysis. It provides a command-line interface (CLI) for executing various commands related to variant data.

> **Table of Contents**  
> - [Concept](#concept)
> - [Installation and Requirements](#installation-and-requirements)  
> - [Usage Overview](#usage-overview)  
> - [Commands](#commands)  
>   - [1. `findtype`](#1-findtype)  
>   - [2. `readvcfs`](#2-readvcfs)  
>   - [3. `runskat`](#3-runskat)  
> - [Statistical Analysis](#statistical-analysis)
> - [Data Files](#data-files)
> - [Example Workflows](#example-workflows)

---    

## Concept

The concept of this project is to collate large numbers of VCF files of unique samples using [Hail](https://hail.is/) from Broad Institute [GitHub](https://github.com/hail-is/hail), annotate variants using [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), and create variant frequency tables in a pipeline fashion. The tool is set up to run locally on a large server (e.g., 64 core CPU with 2GB of RAM per core). 

Input VCF formats may differ between generations, but the pipeline is dependent on correct INFO fields. Input VCF files can be validated using the VCF-validator from the [vcftools](https://vcftools.github.io/index.html) toolset. 

An additional statistical module enables analysis of the frequency tables to reveal potential relations between the variant burden within gene sets and correlating phenotypes – a form of variant enrichment analysis using Monte Carlo permutation analysis. The project is designed to be dockerized.

**Key Features:**
- **Variant Burden Analysis:** Focuses on oligogenic variant interactions
- **Monte Carlo Simulation:** Statistical analysis for rare variant enrichment
- **Scalable Processing:** Uses Hail for large-scale genomic data processing
- **VEP Integration:** Comprehensive variant annotation
- **Frequency Tables:** Creates detailed variant frequency bins for statistical analysis

---    

### Libraries Used

- **Hail:** A Python library for scalable genomic data analysis
- **Variant Effect Predictor (VEP):** A tool for annotating and predicting the effects of genetic variants
- **Pandas:** Data manipulation and analysis for statistical computations
- **NumPy/SciPy:** Scientific computing for Monte Carlo simulations

## Setup

### Docker Setup
1. Clone the project into a local folder
2. Insert the VEP annotator perl script location to the config file
3. `docker image build -t oligogenicity:latest`
4. `docker run -v {path to VEP}:{path in config} --help`

### Local Machine Setup
1. Clone the project
2. Install perl from https://www.perl.org/get.html
3. Install VEP from https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html
4. Update `vep_settings.json`, set VEP executable, VEP cache dir 
5. Create local virtual environment by executing command `python3 -m venv {path-to-your-environment}`
6. Activate virtual environment by using command `{path-to-your-environment}/bin/activate`
7. Install requirements by command `pip install -r requirements.txt`
8. Run: `python3 main.py --help`

---    

## Installation and Requirements

1. **Python 3.9+** (or newer)
2. **[Hail 0.2+](https://hail.is/docs/0.2/install/)**  
3. A **Spark** environment or config that Hail can use (the pipeline sets some default Spark conf parameters in `spark_conf.json`)
4. Any additional environment setups for [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) if you plan to annotate

Ensure you have installed Hail and VEP as appropriate for your environment. The pipeline expects a working `hail` library and that your `spark.jars` references the `hail-all-spark.jar`. If you see messages about "pip-installed Hail requires additional configuration," set the Spark properties as recommended or adjust them in `src/config/spark_conf.json`.

---

## Usage Overview

Call `main.py` with one of the subcommands (`findtype`, `readvcfs`, `runskat`) and the relevant arguments. For instance:

```bash
python3 main.py <command> [options...]
```
Each command has unique flags. Below are the details.
          

## Commands

### 1. `findtype`

**Description:** Searches through a directory (recursively) to find files with a specific extension. Can optionally filter by a regex, and writes a list of the found files to disk.

**CLI Args:**
* `-s / --source` : Directory to be searched (required)
* `-d / --directory` : Directory to save the output file listing (required)
* `-t / --type` : Extension/filetype to look for (required)
* `-r / --regex` : Optional regex string to further filter results

**Example Usage:**

```bash
python3 main.py findtype \
  --source /path/to/input/folder \
  --directory /path/to/output/folder \
  --type .vcf \
  --regex "S1" 
```
This searches `/path/to/input/folder` for files ending in `.vcf, then writes a text file listing results into `/path/to/output/folder`. This text file can be used to collate files using a number of shell tools or edited to create a 

### 2. `readvcfs`

**Description:** Imports VCF(s) into Hail MatrixTables, annotates them (optionally via VEP), merges them in memory, and writes frequency tables out. Also requires a `globals` file to attach metadata (e.g., phenotype) to the input samples. Creates frequency bins for statistical analysis.

**CLI Args:**
* `-f / --file` : VCF file(s), or a text file containing VCF paths, or a folder of VCFs (required)
* `-d / --dest` : Output directory for final tables/MatrixTables (required)
* `-r / --overwrite` : Overwrite existing output data if set
* `-g / --globals` : Tab-delimited file for sample-level metadata (required)
  * The first column is typically the sample ID; additional columns (e.g., phenotype) can be used
* `-p / --phenotype` : If set, the code filters samples that match this phenotype (optional)
* `--annotate / --no-annotate` : Whether to run VEP annotation (defaults to --annotate)
* `--write` : If set, writes the final merged MatrixTable to disk in addition to frequency tables
* `-t / --temp` : Temporary directory for Hail scratch data
* `--interval` : Interval to downfilter i.e. chrX or 1:12000-20000

**Example Usage:**

```bash
python3 main.py readvcfs \
  --file /data/my_vcfs/*.vcf \
  --dest /output/folder \
  --globals /data/phenotype-metadata.txt \
  --phenotype "Case" \
  --annotate \
  --write \
  --temp /tmp/hail/
```

This finds all .vcf files under `/data/my_vcfs/`, annotates them with VEP, merges them, writes frequency tables and MatrixTables in `/output/folder`, and references metadata from `/data/phenotype-metadata.txt`.

### 3. *Optional `runskat`*

**Description:** Reads VCFs per phenotype group (cases and controls), processes in batches, extracts genotype matrices, joins them across samples, and saves final matrices as pickled pandas DataFrames for SKAT-compatible input.

**CLI Args:**
* `-f / --file` : VCF file(s), or a text file containing VCF paths, or a folder of VCFs (required)
* `-d / --dest` : Output directory for final matrices (required)
* `-r / --overwrite` : Overwrites existing output data if set
* `-g / --globals` : Tab-delimited file for sample-level metadata (required)
* `-p / --phenotype` : Filter according to this specific phenotype (required)
* `--annotate / --no-annotate` : Whether to run VEP annotation (defaults to --annotate, will fail if VEP is not installed and configured)
* `--write` : Write intermediate MatrixTables to disk for later loading
* `-t / --temp` : Temporary directory for Hail scratch data
* `--interval` : Interval to downfilter i.e. chrX or 1:12000-20000

**Example Usage:**

```bash
python3 main.py runskat \
  --file /data/my_vcfs/*.vcf \
  --dest /output/skat \
  --globals /data/phenotype-metadata.txt \
  --phenotype "Case" \
  --annotate \
  --temp /tmp/hail/
```

This processes VCFs for cases and controls, creates genotype matrices suitable for SKAT analysis, and saves them as pickled DataFrames in `/output/skat`.

---

## Statistical Analysis

The statistical analysis component of this pipeline is implemented in a comprehensive Jupyter notebook located at `extras/stats-notebook.ipynb`. This notebook contains:

- **Frequency Table Processing:** Analysis of variant burden summary stats analysis
- **Monte Carlo Permutation Analysis:** Simulating random gene pathways for variant burden analysis
- **Oligogenic Interaction Models:** Logistic regression for detecting burden patterns in the simulations
- **Visualization Tools:** Plots and charts for interpreting results
- **Result Interpretation:** Methods for evaluating statistical significance

The notebook processes the frequency tables generated by the pipeline commands and performs sophisticated statistical analyses to identify potential oligogenic interactions and variant burden correlations with phenotypes.

---

## Data Files

The `extras/data/` directory contains example frequency table datasets for analysis:

### Frequency Tables
- **`frequency_table_639_LIHAS_positive_rarevariants.csv`** (191KB, 6,378 lines)
  - Neuromuscular disease dataset positive samples (N=639), also including samples with monogenic rare variants
- **`frequency_table_639_NULLDIST_positive.csv`** (190KB, 6,358 lines)  
  - Expected null distribution of random samples as comparison (N=639)
- **`frequency_table_9059_LIHAS_negative_rarevariants.csv`** (411KB, 12,944 lines)
  - Dataset of negative samples for NMD (N=9,059) with rare variants for other phenotypes 
- **`frequency_table_9059_NULLDIST_negative.csv`** (415KB, 13,073 lines)
  - Expected null distribution for negative samples (N=9,059), i.e. random phenotypes

These files represent processed variant frequency data that can be used for statistical analysis in the notebook. The summary variant burden files were generated from a dataset of 9698 samples sequenced for a variety of diseases using Illumina TruSight One and TruSight Expanded targeted gene sequencing panels and represents the data that was analysed for the primary publication. You may generate your own frequency data with targeted data/WES/WGS VCF inputs.

---

## Example Workflows

Below is a typical multi-step usage:

### 1. Discovery and Listing
```bash
python3 main.py findtype \
  -s /path/to/inputs \
  -d /path/to/listing \
  -t .vcf
```
This writes a file in `/path/to/listing` enumerating .vcf files found in `/path/to/inputs`. A samples file containing the phenotypes of available samples should be generated from this output:
i.e. metadata.txt, *tab* seperated 3 columns of sample \t phenotype \t comments. One sample per line. The phenotype column may contain several phenotypes per sample, defaults to comma seperated.
*sample1.vcf*\t*NMD,DMD*\t*NC_000023.10:g.33229398C>A*
This enables grouping of one or more phenotypes for analysis. This sample-phenotype file is used as input in *readvcfs,runskat* (--globals) and further filtered (--phenotype DMD) to generate a positive set. 

### 2. VCF Processing and Frequency Table Generation
```bash
python3 main.py readvcfs \
  -f /path/to/my_vcfs/*.vcf \
  -d /path/to/readvcfs_out \
  -g /path/to/metadata.txt \
  --phenotype "Case" \
  --annotate \
  --write
```
This merges VCFs, annotates via VEP, creates frequency tables, and writes out MatrixTables in `/path/to/readvcfs_out`.

### 3. SKAT Matrix Generation
```bash
python3 main.py runskat \
  -f /path/to/my_vcfs/*.vcf \
  -d /path/to/skat_output \
  -g /path/to/metadata.txt \
  --phenotype "Case" \
  --annotate
```
Optionally, this creates genotype matrices suitable for SKAT analysis and saves them as pickled DataFrames. SKAT analysis might be useful for a comparison of findings from Monte Carlo simulation analysis.

### 4. Statistical Analysis
Open and run the Jupyter notebook:
```bash
jupyter notebook extras/stats-notebook.ipynb
```
Follow the notebook to perform Monte Carlo simulations and variant burden analysis on the generated frequency tables.

---

### 5. Limitations
- **Single Workstation VEP:** Using Hail to activate VEP annotations on a single workstation can be very slow in comparison to modern clinical annotation pipelines. One might therefore be more inclined to use this tool with a single generation VCF input such as Illumina DRAGEN output VCF files that already contain variant annotations.
- **Limited annotations:** The pipeline is set to work on 4 minimum annotations from VEP: impact, gene label, HGNC ID, MAX_AF
- **Different VCF generations:** Using different generations of input VCF files may pose unexpected faults if formatting changes throughout, Hail expects all VCFs to have the 4 annotations to be in the same format and location within the VCF (i.e. CSQ string formatting). More annotations are allowed, but then the code snippet should be modified accordingly for correct CSQ string handling.
- **Unsupported variants:** Star-alleles are filtered due to the way Hail processes variants.
- **Hail exceptions and logging:** The pipeline is setup to fail gracefully after catching exceptions, which are contained in Hail logs. Hail logfiles can grow to a large size if logging levels are inapproriately low for large datasets. As the pipeline is sequential, the pipeline will skip intermediary output files and continue from the last known position if the command is restarted with the same CLI inputs.

## Contributing

Ethics approval was granted by the University of Tartu Research Ethics Committee for reanalysis of variant files. In accordance with the research permit 374/M-6. This project is a collaboration between the Institute of Clinical Medicine, Tartu University and the Genetics and Personalized Medicine Clinic, Tartu University Hospital. The work is financed by the Estonian Research Council grant PSG774. All issues in repository should be added to the issues section in code repository.

## License

See LICENSE.md for license information.

