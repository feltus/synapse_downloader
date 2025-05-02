# synapse_downloader

# Background

This script will download files from synapse.org.  This specific example will download MCF10A RNAseq datasets.

# Preparing for data download
#Get a list of Synapse IDs from https://www.synapse.org/Synapse:syn18486863/tables/ and place in a file.
#Example MCF10A RNAseq dataset list is included as 'mcf10a-lincs_rnaseq-datasets.txt'

#Set  up a token at synapse.org and place in ~/.synapseConfig
```
[authentication]
authtoken = <TOKEN>
```

#Install Synapse Data Transfer client.  You might want to install in a conda or venv environment.
```
pip install synapseclient
```
# Download the datasets 

```
python DL_SYNAPSE.py -f mcf10a-lincs_rnaseq-datasets.txt
```

# Manipulate HDFH5 files
## Installing HDF5 Tools
conda install -c conda-forge hdf5 h5py pytables

## Viewing File Structure
### List the contents/structure
h5ls syn15574159_abundance.h5

### List with more details
h5ls -v syn15574159_abundance.h5

### Recursive listing (show all groups/datasets)
h5ls -r syn15574159_abundance.h5

## Examining File Contents
### Dump entire contents
h5dump syn15574159_abundance.h5

### Export to text format
h5dump -o syn15574159_abundance.output.txt syn15574159_abundance.h5

# Convert HDF5 files to Gene Expression Matrices

## H5 to GEM Converter
A Python tool for converting Kallisto HDF5 RNA-seq outputs into gene expression matrices (GEM), with support for count data and TPM calculations.

## Overview
This utility processes Kallisto RNA-seq output files (.h5 format) and generates standardized gene expression matrices for downstream analysis. It can process a single file or batch process all files in a directory, combining them into a unified matrix.

## Features
Extracts integer read counts from Kallisto HDF5 files
Calculates Transcripts Per Million (TPM) using gene lengths from a GTF file
Generates combined expression matrices from multiple samples
Creates library size statistics for quality control
Handles ENSEMBL transcript IDs with robust parsing

## Requirements
Python 3.6+
Required packages:
numpy
pandas
h5py

## Install dependencies
```
pip install numpy pandas h5py

## Usage
Basic Command Structure
```
python h5_to_gem.py [-f FILE] [-d DIRECTORY] -g GTF_FILE [-o OUTPUT] [-l]
```

Command-line Arguments
Argument	Description
-f, --file	Path to a single HDF5 file to process
-d, --directory	Directory containing multiple HDF5 files to batch process
-g, --gtf	Path to GTF annotation file for gene lengths (required)
-o, --output	Base name for output files (default: derived from input filename)
-l, --list	List structure of HDF5 file and exit
Note: Either -f or -d must be specified.

### Examples
```
#Process a Single File
python h5_to_gem.py -f syn15574159_abundance.h5 -g Homo_sapiens.GRCh38.113.chr.gtf

#Process All Files in a Directory
python h5_to_gem.py -d ./kallisto_results -g Homo_sapiens.GRCh38.113.chr.gtf -o breast_cancer_rnaseq

# Examine H5 File Structure
python h5_to_gem.py -f syn15574159_abundance.h5 -l

# Specify Custom Output Name
python h5_to_gem.py -f syn15574159_abundance.h5 -g Homo_sapiens.GRCh38.113.chr.gtf -o mcf10a_rnaseq
```

##Output Files
The script generates three files:

{base_name}.count.gem: Tab-separated matrix of integer read counts
Rows: Transcript IDs
Columns: Sample names
{base_name}.tpm.gem: Tab-separated matrix of TPM values
Rows: Transcript IDs
Columns: Sample names
{base_name}.library.specs: CSV file with library size information
Columns: filename, library_size
TPM Calculation
Transcripts Per Million (TPM) normalizes for both sequencing depth and transcript length:

Divide read count by transcript length (in kb) for each transcript
Sum all these values to get the normalization factor
Divide each count/length value by this sum and multiply by 1,000,000
Transcript ID Handling
The script extracts ENSEMBL transcript IDs (e.g., ENST00000488147.1) from complex identifiers such as:

ENSTR0000302805.7|ENSGR0000168939.11|OTTHUMG00000022675.3|OTTHUMT00000058823.3|SPRY3-001|SPRY3|9019|protein_coding|
It handles IDs with or without version numbers for maximum compatibility with the GTF file.

Troubleshooting
Missing Transcript Lengths
If transcripts have zero length in the TPM calculation:

Check that the GTF file is compatible with your transcript IDs
Ensure the GTF file contains exon features with correct transcript_id attributes
Run with -l option to examine the H5 file structure
Sample ID Mismatch
If processing multiple files and seeing "Gene IDs do not match reference gene set" warnings:

Verify that all H5 files were generated using the same reference transcriptome
Check for sample contamination or mixed species data

## Example Workflow
```
#Download GTF annotation
wget ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz

#Run the conversion
python h5_to_gem.py -d ./kallisto_results -g Homo_sapiens.GRCh38.113.gtf -o breast_cancer_study
```

This will generate three files:

breast_cancer_study.count.gem
breast_cancer_study.tpm.gem
breast_cancer_study.library.specs
