#!/bin/bash

set -euo pipefail

#########################################
# Script: 03_deduplication.sh
# Purpose: Mark duplicates using Picard, using CB tag for single-cell barcodes
#
# Dependencies: picard >=2.0, java, samtools
# Input: Sorted BAM from alignment step
# Output: Deduplicated BAM and duplication metrics
#########################################

### User-defined variables
PROJECT="ITscCT"
TISSUE="MG"
TIMEPOINT="10W"
REPLICATE="Rep1"
MARK="H3K4me3"
SAMPLE="${PROJECT}-${TISSUE}-${TIMEPOINT}-${REPLICATE}-${MARK}"
THREADS=20
BASE_DIR="/path/to/project"
ALIGN_DIR="${BASE_DIR}/align"
OUT_DIR="${ALIGN_DIR}/dedup"
TEMP_DIR="${OUT_DIR}/temp"

mkdir -p ${OUT_DIR} ${TEMP_DIR}

### Deduplication with Picard
picard MarkDuplicates \
  I=${ALIGN_DIR}/${SAMPLE}.bowtie2.sort.bam \
  O=${OUT_DIR}/${SAMPLE}.dedup.bam \
  METRICS_FILE=${OUT_DIR}/${SAMPLE}_dedup_metrics.txt \
  BARCODE_TAG=CB \
  REMOVE_DUPLICATES=true \
  ASSUME_SORT_ORDER=coordinate \
  TMP_DIR=${TEMP_DIR} \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

### Index deduplicated BAM
samtools index -@ ${THREADS} ${OUT_DIR}/${SAMPLE}.dedup.bam

echo "Deduplication complete for ${SAMPLE}."
