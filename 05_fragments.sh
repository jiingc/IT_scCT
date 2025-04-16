#!/bin/bash

set -euo pipefail

#########################################
# Script: 05_fragments.sh
# Purpose: Generate fragment BED files per sample using Sinto and prepare for indexing
#
# Dependencies: sinto, samtools, bedtools, bgzip, tabix
# Input: Split BAMs (from filter step)
# Output: Sorted, compressed, and indexed BED files
#########################################

### User-defined variables
PROJECT="ITscCT"
TISSUE="MG"
TIMEPOINT="10W"
REPLICATE="Rep1"
THREADS=30
BASE_DIR="/path/to/project"
SPLIT_DIR="${ALIGN_DIR}/dedup/split/${REPLICATE}"
FRAG_DIR="${BASE_DIR}/fragments/${REPLICATE}"

mkdir -p ${FRAG_DIR}

### Generate fragments from split BAMs
cd ${SPLIT_DIR}

for bam in *.bam; do
  sample_name=${bam%%.bam}  # get file prefix without extension
  output_bed="${FRAG_DIR}/${sample_name}.fragments.bed"

  sinto fragments \
    -b ${bam} \
    -f ${output_bed} \
    -p ${THREADS} \
    --collapse_within

  sort -k1,1 -k2,2n ${output_bed} > ${output_bed%.bed}.sort.bed
  bgzip -@ ${THREADS} ${output_bed%.bed}.sort.bed
  tabix -p bed ${output_bed%.bed}.sort.bed.gz

done

echo "Fragment BED generation and indexing complete for ${REPLICATE}."
