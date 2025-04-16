#!/bin/bash

set -euo pipefail

#########################################
# Script: 04_filter_split.sh
# Purpose: Filter BAM by barcode and split by sample/mark using Sinto
#
# Dependencies: sinto, samtools
# Input:
#   - Deduplicated BAM file (with CB tags)
#   - Barcode-to-sample mapping file (TSV with barcode \t sample)
# Output:
#   - Split BAM files in a dedicated output folder
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
ALIGN_DIR="${BASE_DIR}/align"
DEDUP_BAM="${ALIGN_DIR}/dedup/${SAMPLE}.dedup.bam"
INDEX_MAP_FILE="${BASE_DIR}/index/"${SAMPLE}.barcode_to_sample.txt"  # format: CB\tsample_name
SPLIT_DIR="${ALIGN_DIR}/dedup/split/${REPLICATE}"

mkdir -p ${SPLIT_DIR}

### Run sinto filterbarcodes to split by sample identity
sinto filterbarcodes \
  -b ${DEDUP_BAM} \
  -c ${INDEX_MAP_FILE} \
  -p ${THREADS} \
  --outdir ${SPLIT_DIR}

### Index the split BAM files
cd ${SPLIT_DIR}
for bam in *.bam; do
  samtools index -@ ${THREADS} $bam
done

echo "Barcode filtering and sample splitting complete for ${SAMPLE}."
