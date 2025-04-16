#!/bin/bash

set -euo pipefail

#########################################
# Script: 02_alignment.sh
# Purpose: Align demultiplexed paired-end reads with Bowtie2 and process BAM
#
# Dependencies: bowtie2, samtools, sambamba
# Input: clean FASTQ from demultiplexing step
# Output: sorted BAM ready for deduplication
#########################################

### User-defined variables
PROJECT="ITscCT"
TISSUE="MG"
TIMEPOINT="10W"
REPLICATE="Rep1"
MARK="H3K4me3"
SAMPLE="${PROJECT}-${TISSUE}-${TIMEPOINT}-${REPLICATE}-${MARK}"
THREADS=30
BASE_DIR="/path/to/project"
FASTQ_DIR="${BASE_DIR}/demultiplex/clean"
ALIGN_DIR="${BASE_DIR}/align"
REF_GENOME="/path/to/reference/bowtie2/index_prefix"  # <-- customize this

mkdir -p ${ALIGN_DIR}/summary
cd ${ALIGN_DIR}

### Alignment with Bowtie2
bowtie2 \
  --very-sensitive-local \
  --soft-clipped-unmapped-tlen \
  --no-mixed --no-discordant --dovetail \
  --phred33 -I 10 -X 1000 \
  -x ${REF_GENOME} \
  -1 ${FASTQ_DIR}/${SAMPLE}_R1_clean.fq.gz \
  -2 ${FASTQ_DIR}/${SAMPLE}_R2_clean.fq.gz \
  -S ${ALIGN_DIR}/${SAMPLE}.bowtie2.sam \
  -p ${THREADS} \
  &> ${ALIGN_DIR}/summary/${SAMPLE}_bowtie2.log

### Convert to BAM, filter low quality, sort
samtools view -hbSq 10 -@ ${THREADS} ${SAMPLE}.bowtie2.sam > ${SAMPLE}.bowtie2.bam
sambamba sort -t ${THREADS} -o ${SAMPLE}.bowtie2.sort.bam ${SAMPLE}.bowtie2.bam

### Clean up intermediate files
rm ${SAMPLE}.bowtie2.sam
rm ${SAMPLE}.bowtie2.bam

echo "Alignment and sorting complete for ${SAMPLE}."
