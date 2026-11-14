# WGS Analysis Pipeline - Variant Calling & Annotation
### An End-to-end Whole Genome Sequencing (WGS) Analysis Workflow for Human Genomic Data  


> ### **Project Overview**
### This project provides a full pipeline analysis of human whole genome sequencing data from raw fastq files to functional annotations with integration of bioinformatics tools, this workflow includes: 
* Retrieval of sample files from **1000G Project**
* Quality control of raw reads using **Fastqc**
* Trimming & filtration of lower quality reads based on QC **(optional)**
* Mapping of reads to the reference genome **GRCh38**
* Manipulation & Processing of **BAM** files using **Samtools** 
* Variants calling, genotyping and filtration of variants **GATK4**
* Manipulation & Processing of **VCF** files **Bcftools** 
* Functional annotation of variants **SnpEff**
##

## Main Workflow
```
[A: Raw FASTQs Reads]
        |
        ↓
[B: QC of FASTQs Reads (FastQC)]
        ├────────────────────────┐
        |                        |
        |               [Trimming (Optional)]
        |                        |
        ├────────────────────────┘ 
        ↓
[C.1: Mapping to Reference Genome (BWA-MEM)]
        |
        ↓
[C.2: Sorting & Indexing of BAM File (Samtools)]
        |
        ↓
[D: Variants Analysis]───────────────────────────────────────┐
        |                                                    |
        ↓                                                    |
[D.1: Marking Duplicates]                                    |
        |                                                    |
        ↓                                                    |
[D.2: Base Quality Score Recalibrator (BQSR)]                |
        |                                                    |
        ↓                                                    |
[D.3: Post-Alignment QC]                                  [GATK4]
        |                                                    |
        ↓                                                    |
[D.4: HaplotypeCaller]                                       |
        |                                                    |
        ↓                                                    |
[D.5: GenotypeGVCFs]                                         |
        |                                                    |
        ↓                                                    |
[D.6: Filtering of Variants]─────────────────────────────────┘
        |
        ↓
  --Filtered VCF--
        |
        ↓
[E: Annotation of Variants (SnpEff)]
        |
        ↓
  --Annotated VCF--
```
##

## Requirements
### Input Files & Reference
|   Resource    |      Input Files      |
|-----------|-------------------|
| Sample Reads | `R1.fastq` `R2.fastq` |
| Reference Genome | `Homo_sapiens_assembly38.fasta` |
| Known Variants VCFs| `dbsnp_146.hg38.vcf.gz` `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` |
| Annotation| `GRCh38.115`|

### Tools Dependencies
|   Tool    |      Purpose      |           Version          |
|-----------|-------------------|----------------------------|
| Fastqc    | Quality Check & Control of Reads | `≥ v0.13.0` |
| Fastp     | Trimming of Lower Quality Reads  | `≥ v1.4.0`  |
| BWA-MEM   | Mapping of Reads to Reference    | `≥ v0.7.19` |
| Samtools  | Processing of SAM & BAM Files    | `≥ v1.24`   |
| GATK4     | Calling & Filtering of Variants  | `≥ v4.6.2.0`|
| Bcftools  | Processing of VCF Files          | `≥ v1.24`   |
| SnpEff    | Functional Annotation of Variants| `≥ v5.4c`    |
##

## Setup Environment & Required Tools

### Clone Repository
```
git clone https://github.com/Fathallah808/WGS-Analysis-Pipeline.git
cd WGS-Analysis-Pipeline
```

### Create & Activate Conda Environment
```
# Creating conda environment
conda create -n WGS_Analysis

# Activation of conda environment
conda activate WGS_Analysis
```

### Install Dependencies
```
conda install -c bioconda -c conda-forge \
fastqc fastp bwa samtools gatk4 bcftools snpEff -y
```

### Shell Script Permission
```
chmod +x Analysis.sh
```
##

## Project Structure
```
WGS-Analysis-Pipeline/
├──FASTQs/
|  └──SRR077311 Raw FASTQ Files
|
├──Fastqc_Results/
|  └──Fastqc Html Reports & Zip Files
|
├──Trimmed_Reads/
|  └──Trimmed FASTQ Reads
|
├──Reference/
|  └──Ref Genome, Dictionary, Fai, Index Files
|
├──Aligned_Reads/
|  ├─Bam Files
|  ├─BQSR Table
|  └─Collection_Metrics
|
├──Known_VCFs/
|  └─Known Variants VCFs
|
├──Variants_Files/
|  └─Variant Calling VCFs Files
|
├──Annotation_Results
|  └─Annotated VCF File
|
├──Analysis.sh
└──README.md
```
##

## References & Resources

### Seq data were derived from public 1000G project data: 
```
https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00107/sequence_read/
```

### GATK Resources:
```
https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels
https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
```
##

## Disclaimer
This pipeline is intended for educational and portfolio project purposes, the analyzed data were obtained from public 
**1000 genome project** available data and it's intended for practicing and it doesn't provide any original findings
##