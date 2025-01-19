# Oligogenicity
This pipeline provides a suite of commands (via `main.py`) to manage, annotate, and analyze VCFs with the Hail framework. The pipeline includes commands for locating specific files, importing/annotating VCFs, combining MatrixTables into a “database,” computing relatedness (via KING), and performing PCA plotting.

The gen-toolbox project is a comprehensive tool for collating large numbers of VCF files of unique samples, annotating variants, and creating variant frequency tables in a pipeline fashion. The project is designed to run locally on a large server and is directed by Tartu University Hospital Centre of Medical Genetics / Tartu University Institute of Clinical Medicine. This work was supported by the Estonian Research Council grant PSG774.

This repository contains tools for genomic data processing and analysis. It provides a command-line interface (CLI) for executing various commands related to genomic data.

> **Table of Contents**  
> - [Concept](#Concept)
> - [Installation and Requirements](#installation-and-requirements)  
> - [Usage Overview](#usage-overview)  
> - [Commands](#commands)  
>   - [1. `findtype`](#1-findtype)  
>   - [2. `readvcfs`](#2-readvcfs)  
>   - [3. `loaddb`](#3-loaddb)  
>   - [4. `relatedness2`](#4-relatedness2)  
>   - [5. `pca`](#5-pca)  
> - [Example Workflows](#example-workflows)

---    

## Concept

The concept of this project is to collate large numbers of VCF files of unique samples using [Hail](https://hail.is/) from Broad Institute [GitHub](https://github.com/hail-is/hail), annotate variants using [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), and create variant frequency tables in a pipeline fashion. The tool is set up to run locally on a large server (e.g., 64 core CPU with 2GB of RAM per core). Input VCF formats may differ between generations, but the pipeline is dependent on correct INFO fields. Input VCF files can be validated using the VCF-validator from the [vcftools](https://vcftools.github.io/index.html) toolset. An additional statistical module enables analysis of the frequency tables to reveal potential relations between the variant burden within gene sets and correlating phenotypes – a form of variant enrichment analysis using Monte Carlo permutation analysis. The project is designed to be dockerized.

---    

### Libraries Used

- Hail: A Python library for scalable genomic data analysis.
- Variant Effect Predictor (VEP): A tool for annotating and predicting the effects of genetic variants.

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

---    

## Installation and Requirements

1. **Python 3.9+** (or newer).
2. **[Hail 0.2+](https://hail.is/docs/0.2/install/).**  
3. A **Spark** environment or config that Hail can use (the pipeline sets some default Spark conf parameters in `spark_conf.json`).
4. Any additional environment setups for [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) if you plan to annotate.

Ensure you have installed Hail and VEP as appropriate for your environment. The pipeline expects a working `hail` library and that your `spark.jars` references the `hail-all-spark.jar`. If you see messages about “pip-installed Hail requires additional configuration,” set the Spark properties as recommended or adjust them in `src/config/spark_conf.json`.

---

## Usage Overview

Call `main.py` with one of the subcommands (`findtype`, `readvcfs`, `loaddb`, `relatedness2`, `pca`) and the relevant arguments. For instance:

```bash
python3 main.py <command> [options...]
```
Each command has unique flags. Below are the details.
          

## Commands

1) ```findtype```

Description: Searches through a directory (recursively) to find files with a specific extension. Can optionally filter by a regex, and writes a list of the found files to disk.

CLI Args:
* -s / --source : Directory to be searched (required).
* -d / --directory : Directory to save the output file listing (required).
* -t / --type : Extension/filetype to look for (required).
* -r / --regex : Optional regex string to further filter results.

Example Usage:

```bash
python3 main.py findtype \
  --source /path/to/input/folder \
  --directory /path/to/output/folder \
  --type .vcf \
  --regex "S1" 
```
This searches /path/to/input/folder for files ending in .vcf whose paths also match the regex "S1", then writes a text file listing results into /path/to/output/folder.

2) ```readvcfs```

Description: Imports VCF(s) into Hail MatrixTables, annotates them (optionally via VEP), merges them in memory, and writes frequency tables or the combined MT out. Also reads a globals file to attach metadata (e.g., phenotype).

CLI Args:
*	-f / --file : VCF file(s), or a text file containing VCF paths, or a folder of VCFs (required).
*	-d / --dest : Output directory for final tables/MatrixTables (required).
*	-r / --overwrite : Overwrite existing output data if set.
*	-g / --globals : Tab-delimited file for sample-level metadata (required).
*	The first column is typically the sample ID; additional columns (e.g., phenotype) can be used.
*	-p / --phenotype : If set, the code filters samples that match this phenotype (optional).
*	--annotate / --no-annotate : Whether to run VEP annotation (defaults to --annotate).
*	--write : If set, writes the final merged MatrixTable to disk in addition to a summary table.
*	-t / --temp : Temporary directory for Hail scratch data.

Example Usage:

```bash
python3 main.py readvcfs \
  --file /data/my_vcfs/*.vcf \
  --dest /output/folder \
  --globals /data/phenotype-pseudon.txt \
  --annotate \
  --write \
  --temp /tmp/hail/
```

This finds all .vcf files under /data/my_vcfs/, annotates them with VEP, merges them, writes a MatrixTable in /output/folder, and references metadata from /data/phenotype-pseudon.txt.

3) ```loaddb```

Description: Loads Hail MatrixTables from disk (previously written by readvcfs or other steps), merges them, and writes out a combined “database” (also merges frequency bins). Optionally filters by phenotype.

CLI Args:
*	-f / --file : One or more .mt folder paths, or a text list of them, or a directory of MatrixTables.
*	-r / --overwrite : Overwrite existing output data if set.
*	-d / --dest : Output directory for final frequency tables / MT.
*	-n / --number : Number of tables to load/merge (default -1 loads all).
*	-g / --globals : Same style metadata file for sample IDs, phenotypes, etc.
*	--phenotype : Regex-based sample filtering.
*	--write : If set, writes the final merged MatrixTable to disk.

Example Usage:

```bash
python3 main.py loaddb \
  --file /output/folder/*.mt \
  --dest /merged/db \
  --globals /data/phenotype-pseudon.txt \
  --write
```
This merges all .mt files in /output/folder/*.mt, applies sample metadata from /data/phenotype-pseudon.txt, and writes out a combined dataset to /merged/db.

4) ```relatedness2```

Description: Reads one or more VCF(s), merges them by columns, then runs Hail’s KING method to compute pairwise relatedness (phi). Writes a table of relatedness or “related” pairs to disk. By default, it only processes chromosome 2 (the code specifically filters to interval=["2"] as an example).

CLI Args:
*	-f / --file : VCF file(s) or directory of VCFs (required).
*	-d / --dest : Destination for the resulting relatedness.tsv file (required).
*	-t / --temp : Temporary folder for Hail’s intermediate data.

Example Usage:

```bash
python3 main.py relatedness2 \
  --file /data/my_vcfs/sample1.vcf /data/my_vcfs/sample2.vcf \
  --dest /results/relatedness \
  --temp /tmp/hail/
```

This merges sample1.vcf and sample2.vcf, filters to chromosome 2, runs KING, and writes /results/relatedness/relatedness.tsv with phi scores for each sample pair.

5) ```pca```

Description: Plots PCA results from external tab-delimited data. The pipeline’s usage expects you to have a PCA output and an associated .tsv with population or sample metadata.

CLI Args:
*	--pca : Input PCA file (required).
*	--pca_tsv : Second PCA metadata file (required).
*	-d / --dest : Output directory for any resulting images or logs.

Example Usage:

```bash
python3 main.py pca \
  --pca /data/pca_out/pca_scores.txt \
  --pca_tsv /data/pca_out/pca_meta.tsv \
  --dest /figures/
```

This reads the PCA components from pca_scores.txt, merges them with metadata from pca_meta.tsv, and produces a 2D scatterplot (PC1 vs. PC2) stored in /figures/.

## Example Workflows

Below is a typical multi-step usage:
1.	List VCFs:

```bash
python3 main.py findtype \
  -s /path/to/inputs \
  -d /path/to/listing \
  -t .vcf
```

This writes a file in /path/to/listing enumerating .vcf files found in /path/to/inputs.

2.	Read & Annotate:

```bash
python3 main.py readvcfs \
  -f /path/to/my_vcfs/*.vcf \
  -d /path/to/readvcfs_out \
  -g /path/to/metadata.txt \
  --annotate \
  --write
```

This merges VCFs, annotates via VEP, and writes out one or more MatrixTables in /path/to/readvcfs_out.

3.	Load & Merge:

```bash
python3 main.py loaddb \
  -f /path/to/readvcfs_out/*.mt \
  -d /path/to/final_db \
  -g /path/to/metadata.txt \
  --phenotype "Case.*" \
  --write
```

Loads each .mt from step 2, merges them with metadata, filters by a phenotype regex "Case.*", and writes the combined data in /path/to/final_db.

4.	KING Relatedness:

```bash
python3 main.py relatedness2 \
  -f /path/to/my_vcfs/*.vcf \
  -d /results/ \
  -t /tmp/hailtemp/
```

Merges the listed VCFs by columns, filters to chromosome 2, then computes pairwise relatedness using KING. Outputs relatedness.tsv in /results.

5.	PCA:

```bash
python3 main.py pca \
  --pca /results/pca_scores.txt \
  --pca_tsv /results/pca_metadata.tsv \
  -d /results/plots
```

Loads the PCA components (pca_scores.txt), merges with sample data in pca_metadata.tsv, and plots PC1 vs. PC2 in /results/plots.

