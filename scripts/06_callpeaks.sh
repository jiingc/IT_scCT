#!/bin/bash

set -euo pipefail

#########################################
# Script: 06_callpeaks.sh
# Purpose: Call peaks from BAM using MACS2 (supports histone marks & TFs)
#
# Dependencies: macs2
# Input:
#   - Genome size (e.g., mm, hs)
#   - Mark type: histone/tf/lads
# Output:
#   - Peak files in broadPeak or narrowPeak format
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

GENOME=$1       # mm / hs
MARK_TYPE=$2    # histone / tf / lads

SPLIT_DIR="${BASE_DIR}/align/dedup/split/${REPLICATE}"
BAM_FILE="${SPLIT_DIR}/${SAMPLE}.bam"
OUT_DIR="${BASE_DIR}/peaks"
mkdir -p ${OUT_DIR}

### Mode selection
if [[ "$MARK_TYPE" == "histone" || "$MARK_TYPE" == "lads" ]]; then
  PEAK_MODE="--broad"
  FORMAT="BAMPE"
else
  PEAK_MODE=""
  FORMAT="BAM"
fi

### Call peaks with no lambda and no summit
macs2 callpeak \
  -t ${BAM_FILE} \
  -f ${FORMAT} \
  -g ${GENOME} \
  -n ${SAMPLE} \
  ${PEAK_MODE} \
  -q 0.05 \
  --keep-dup all \
  --nolambda \
  --outdir ${OUT_DIR} \
  2> ${OUT_DIR}/${SAMPLE}_macs2.log

echo "MACS2 peak calling complete for ${SAMPLE}."
