#!/bin/bash

#—————————————————————————————————————————————————————————————————————————————————————————————————#
# Whole Genome Sequencing (WGS) Variant Calling & Annotation Pipeline
# Tools: FastQC | Fastp | BWA-MEM | Samtools | GATK | SnpEff
# Data were derived from: https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00107/sequence_read/
#—————————————————————————————————————————————————————————————————————————————————————————————————#

#———————————————————————————————————————————————————————————————————#
# Define Paths & Set Variables 
#———————————————————————————————————————————————————————————————————#

# Dir Paths
FASTQs="FASTQs"
Fastqc_Results="Fastqc_Results"
Trimmed_Reads="Trimmed_Reads"
Reference="Reference"
Aligned_Reads="Aligned_Reads"
Known_VCFs="Known_VCFs"
Variants_Files="Variants_Files"
Annotation_Results="Annotation_Results"

mkdir -p $FASTQs $Fastqc_Results $Trimmed_Reads $Reference \
 $Aligned_Reads $Known_VCFs $Variants_Files $Annotation_Results

# Sample Configurations 
Sample_ID="SRR077311"
Sample_Name="HG00107"


#———————————————————————————————————————————————————————————————————#
# Retrieving of Data
#———————————————————————————————————————————————————————————————————#

echo "====== Download Sample Files ======"

wget -P ${FASTQs} https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00107/sequence_read/SRR077311_1.filt.fastq.gz
wget -P ${FASTQs} https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00107/sequence_read/SRR077311_2.filt.fastq.gz

echo "====== Download Reference Genome ======"

wget -P ${Reference} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta

echo "====== Download GATK Bundle & Known Sites ======"

# Hg38 dbSNP
wget -P ${Known_VCFs} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/dbsnp_146.hg38.vcf.gz
wget -P ${Known_VCFs} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/dbsnp_146.hg38.vcf.gz.tbi

# Mills and 1000G Gold Standard Indels
wget -P ${Known_VCFs} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -P ${Known_VCFs} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

echo "——————————————————Files Retrieving [Done]—————————————————————"


#———————————————————————————————————————————————————————————————————# 
# 1: Quality Control of Raw Reads
#———————————————————————————————————————————————————————————————————#
echo "----- 1: Quality Control of Raw Reads [Running] -----"

fastqc \
 ${FASTQs}/SRR077311_1.filt.fastq.gz \
 ${FASTQs}/SRR077311_2.filt.fastq.gz \
 -o ${Fastqc_Results}
 
echo "----- Fastqc Run [Done] -----"     


#———————————————————————————————————————————————————————————————————#
# 2: Trimming of Lower Quality Reads (Optional)
#———————————————————————————————————————————————————————————————————#
#echo "----- 2: Trimming of Lower Quality Reads [Running] -----"

#fastp \
# -i ${FASTQs}/SRR077311_1.filt.fastq.gz \
# -I ${FASTQs}/SRR077311_2.filt.fastq.gz \
# -o ${Trimmed_Reads}/${Sample_Name}_R1_trimmed.fastq \
# -O ${Trimmed_Reads}/${Sample_Name}_R2_trimmed.fastq \
# --detect-adapter-for-pe \
# -q 30 \
# -l 50 \
# -h ${Trimmed_Reads}/${Sample_Name}_fastp_report.html

# QC of Trimmed Reads
#fastqc \
# ${Trimmed_Reads}/${Sample_Name}_R1_trimmed.fastq \
# ${Trimmed_Reads}/${Sample_Name}_R2_trimmed.fastq \
# -o ${Trimmed_Reads}

#echo "----- Trimming of Lower Quality Reads [Done] -----"     


#———————————————————————————————————————————————————————————————————#
# 3: Prepare the Reference Genome
#———————————————————————————————————————————————————————————————————#
echo "----- 3: Prepare the Reference Genome [Running] -----"

# Indexing Reference Genome 
bwa index ${Reference}/Homo_sapiens_assembly38.fasta

# Creating Fai File
samtools faidx ${Reference}/Homo_sapiens_assembly38.fasta

# Creating Sequence Dictionary
gatk createSequenceDictionary \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -O ${Reference}/Homo_sapiens_assembly38.dict

echo "----- Reference Genome Indexing & Processing [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 4: Mapping of Reads to Reference Genome (GRCh38)
#———————————————————————————————————————————————————————————————————#
echo "----- 4: Mapping of Reads to Reference Genome [Running] -----"

# Alignment to Reference Genome
bwa mem -R "@RG\tSM:${Sample_Name}\tLB:LIB1\tID:${Sample_ID}\tPL:ILLUMINA" \
 ${Reference}/Homo_sapiens_assembly38.fasta \
 ${FASTQs}/SRR077311_1.filt.fastq.gz \
 ${FASTQs}/SRR077311_2.filt.fastq.gz \
 | samtools sort -o ${Aligned_Reads}/${Sample_Name}_Sorted.bam

# Alignment Statistics
samtools flagstat \
 ${Aligned_Reads}/${Sample_Name}_Sorted.bam \
 > ${Aligned_Reads}/${Sample_Name}_Alignment_Flagstat.txt

# Indexing of BAM File
samtools index ${Aligned_Reads}/${Sample_Name}_Sorted.bam

# Marking Duplicates
gatk MarkDuplicates \
 -I ${Aligned_Reads}/${Sample_Name}_Sorted.bam \
 -O ${Aligned_Reads}/${Sample_Name}_Marked_Duplicates.bam \
 -M ${Aligned_Reads}/Metrics.txt

# Indexing
samtools index ${Aligned_Reads}/${Sample_Name}_Marked_Duplicates.bam

echo "----- Mapping to Ref. & Indexing [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 5: Base Quality Score Recalibration (BQSR)
#———————————————————————————————————————————————————————————————————#
echo "----- 5: Base Quality Score Recalibration [Running] -----"

# Generating Recalibration Table
gatk BaseRecalibrator \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -I ${Aligned_Reads}/${Sample_Name}_Marked_Duplicates.bam \
 --known-sites ${Known_VCFs}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
 --known-sites ${Known_VCFs}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
 -O ${Aligned_Reads}/BQSR.Recal.table

# Apply BQSR
gatk ApplyBQSR \
 -I ${Aligned_Reads}/${Sample_Name}_Marked_Duplicates.bam \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 --bqsr-recal-file ${Aligned_Reads}/BQSR.Recal.table \
 -O ${Aligned_Reads}/${Sample_Name}.Recalibrate.bam

# Indexing of Recalibrated BAM File
samtools index ${Aligned_Reads}/${Sample_Name}.Recalibrate.bam

echo "----- Base Quality Score Recalibration [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 6: Post-Alignment QC
#———————————————————————————————————————————————————————————————————#
echo "----- 6: Post-Alignment QC [Running] -----"

# Alignment Summary
gatk CollectAlignmentSummaryMetrics \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -I ${Aligned_Reads}/${Sample_Name}.Recalibrate.bam \
 -O ${Aligned_Reads}/Alignment.Summary.txt

# Insert Size Metrics
gatk CollectInsertSizeMetrics \
 -I ${Aligned_Reads}/${Sample_Name}.Recalibrate.bam \
 -O ${Aligned_Reads}/Insert.Size.txt \
 -H ${Aligned_Reads}/Insert.Size.Histogram.pdf

echo "----- Post-Alignment QC [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 7: HaplotypeCaller - Calling of Variants
#———————————————————————————————————————————————————————————————————#
echo "----- 7: HaplotypeCaller - Calling of Variants [Running] -----"

# HaplotypeCaller
gatk HaplotypeCaller \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -I ${Aligned_Reads}/${Sample_Name}.Recalibrate.bam \
 -O ${Variants_Files}/${Sample_Name}.g.vcf \
 -ERC GVCF

# Genotyping 
gatk GenotypeGVCFs \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -V ${Variants_Files}/${Sample_Name}.g.vcf \
 -O ${Variants_Files}/${Sample_Name}.raw.vcf

echo "----- Calling of Variants [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 8: Selection & Filtration of Variants (SNPs/Indels)
#———————————————————————————————————————————————————————————————————#
echo "-----8: Selection of Variants SNPs/Indels [Running] -----"

SelectVariants_Function() {

gatk SelectVariants \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -V ${Variants_Files}/${Sample_Name}.raw.vcf \
 --select-type $1 \
 -O $2
}

# SNPs
SelectVariants_Function \
 SNP \ 
 ${Variants_Files}/${Sample_Name}_Raw.Snps.vcf

# Indels
SelectVariants_Function \
 INDEL \ 
 ${Variants_Files}/${Sample_Name}_Raw.Indels.vcf

echo "----- Selection of Variants SNPs/Indels [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 9: Variants Filtering (Hard Filtering)
#———————————————————————————————————————————————————————————————————#
echo "----- 9: Variants Filtering [Running] -----"

# SNPs Filtering
gatk VariantFiltration \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -V ${Variants_Files}/${Sample_Name}_Raw.Snps.vcf \
 --filter-name "QD_Filter" -filter "QD < 2.0" \
 --filter-name "FS_Filter" -filter "FS > 60.0" \
 --filter-name "MQ_Filter" -filter "MQ < 40.0" \
 --filter-name "SOR_Filter" -filter "SOR > 3.0" \
 --filter-name "MQRankSum_Filter" -filter "MQRankSum < -12.0" \
 --filter-name "ReadPosRankSum_Filter" -filter "ReadPosRankSum < -8.0" \
 -O ${Variants_Files}/${Sample_Name}_Filtered.Snps.vcf

# Selection of PASS Filters
gatk SelectVariants \
 --exclude-filtered true \
 -V ${Variants_Files}/${Sample_Name}_Filtered.Snps.vcf \
 -O ${Variants_Files}/${Sample_Name}_Pass_Filtered.Snps.vcf

# Indels Filtering
gatk VariantFiltration \
 -R ${Reference}/Homo_sapiens_assembly38.fasta \
 -V ${Variants_Files}/${Sample_Name}_Raw.Indels.vcf \
 --filter-name "QD_Filter" -filter "QD < 2.0" \
 --filter-name "FS_Filter" -filter "FS > 200.0" \
 --filter-name "SOR_Filter" -filter "SOR > 10.0" \
 --filter-name "ReadPosRankSum_Filter" -filter "ReadPosRankSum < -20.0" \
 -O ${Variants_Files}/${Sample_Name}_Filtered.Indels.vcf 

# Selection of PASS Filters
gatk SelectVariants \
 --exclude-filtered true \
 -V ${Variants_Files}/${Sample_Name}_Filtered.Indels.vcf \
 -O ${Variants_Files}/${Sample_Name}_Pass_Filtered.Indels.vcf

echo "----- Variant Filtration [Done] -----"


#———————————————————————————————————————————————————————————————————#
# 10: Annotation of Variants (SnpEff)
#———————————————————————————————————————————————————————————————————#
echo "----- 10: Annotation of Variants (SnpEff) -----"

# Annotation Function
Annotation_Function() {
    snpEff ann \
     GRCh38.115 \
     $1 \
     > $2
}

# SNPs Annotation
Annotation_Function \
 ${Variants_Files}/${Sample_Name}_Pass_Filtered.Snps.vcf \
 ${Annotation_Results}/${Sample_Name}_Annotated.Snps.vcf

# Indels Annotation
Annotation_Function \
 ${Variants_Files}/${Sample_Name}_Pass_Filtered.Indels.vcf \
 ${Annotation_Results}/${Sample_Name}_Annotated.Indels.vcf

echo "----- Annotation of Variants (SnpEff) [Done] -----"

#=================================================================================================#
echo "============================================================================================"
echo "=====================[Variant Calling & Annotation Processing Done]========================="
echo "============================================================================================"
#=================================================================================================#