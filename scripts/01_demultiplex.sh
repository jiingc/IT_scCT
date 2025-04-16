#!/bin/bash

set -euo pipefail

#########################################
# Script: 01_demultiplex.sh
# Purpose: Demultiplex raw FASTQ files using Cutadapt for IT-scC&T-seq
#
# Dependencies: cutadapt >=4.0
# Input:
#   - Raw FASTQ files: ${RAW_DIR}/${SAMPLE}_R1.fq.gz and _R2.fq.gz
#   - Index fasta: H5/H7 + Q5/Q7 index files
#   - Linker sequences
# Output:
#   - Demultiplexed FASTQ files in ${OUT_DIR}
#########################################

### User-defined variables
PROJECT="ITscCT"
TISSUE="MG"
TIMEPOINT="10W"
REPLICATE="Rep1"
MARK="H3K4me3"
SAMPLE="${PROJECT}-${TISSUE}-${TIMEPOINT}-${REPLICATE}-${MARK}"
THREADS=4
BASE_DIR="/path/to/project"
RAW_DIR="${BASE_DIR}/raw"
INDEX_DIR="${BASE_DIR}/index"
OUT_DIR="${BASE_DIR}/demultiplex"

mkdir -p ${OUT_DIR}/{round1_H5H7,round2_linkertrim,round3_Q5Q7,clean,unassigned,logs}

### Step 1: Trim H5/H7 index (round 1)
cutadapt \
  -j ${THREADS} \
  -e 0.1 \
  --action=trim \
  --rename='{id} CB:Z:{r1.adapter_name}{r2.adapter_name}' \
  -g ^file:${INDEX_DIR}/H501-H548.fasta \
  -G ^file:${INDEX_DIR}/H701-H724.fasta \
  -o ${OUT_DIR}/round1_H5H7/${SAMPLE}_round1_H5H7_R1.fq.gz \
  -p ${OUT_DIR}/round1_H5H7/${SAMPLE}_round1_H5H7_R2.fq.gz \
  ${RAW_DIR}/${SAMPLE}_R1.fq.gz ${RAW_DIR}/${SAMPLE}_R2.fq.gz \
  --untrimmed-o ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_round1_R1.fq.gz \
  --untrimmed-p ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_round1_R2.fq.gz \
  &>> ${OUT_DIR}/logs/${SAMPLE}_cutadapt_round1_H5H7.log

### Step 2: Trim linker sequences (round 2)
cutadapt \
  -j ${THREADS} \
  -e 0.25 \
  --action=trim \
  -g GGCGGTAGGCGTGCTC \
  -G CCAACACCCGTGCGCTG \
  -o ${OUT_DIR}/round2_linkertrim/${SAMPLE}_round2_linkertrim_R1.fq.gz \
  -p ${OUT_DIR}/round2_linkertrim/${SAMPLE}_round2_linkertrim_R2.fq.gz \
  ${OUT_DIR}/round1_H5H7/${SAMPLE}_round1_H5H7_R1.fq.gz ${OUT_DIR}/round1_H5H7/${SAMPLE}_round1_H5H7_R2.fq.gz \
  --untrimmed-o ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_round2_R1.fq.gz \
  --untrimmed-p ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_round2_R2.fq.gz \
  &>> ${OUT_DIR}/logs/${SAMPLE}_cutadapt_round2_linker.log

### Step 3: Trim Q5/Q7 index (round 3)
cutadapt \
  -j ${THREADS} \
  -e 0.1 \
  --action=trim \
  --rename='{header}{r1.adapter_name}{r2.adapter_name}' \ ## final CB structure CB:Z:H5XXH7XXQ5XXQ7XX
  -g ^file:${INDEX_DIR}/Q501-Q508.fasta \
  -G ^file:${INDEX_DIR}/Q701-Q703.fasta \
  -o ${OUT_DIR}/round3_Q5Q7/${SAMPLE}_round3_Q5Q7_R1.fq.gz \
  -p ${OUT_DIR}/round3_Q5Q7/${SAMPLE}_round3_Q5Q7_R2.fq.gz \
  ${OUT_DIR}/round2_linkertrim/${SAMPLE}_round2_linkertrim_R1.fq.gz ${OUT_DIR}/round2_linkertrim/${SAMPLE}_round2_linkertrim_R2.fq.gz \
  --untrimmed-o ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_round3_R1.fq.gz \
  --untrimmed-p ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_round3_R2.fq.gz \
  &>> ${OUT_DIR}/logs/${SAMPLE}_cutadapt_round3_Q5Q7.log

### Step 4: Trim ME-Tn5 sequences (final clean step)
cutadapt \
  -j ${THREADS} \
  -e 0.25 \
  --action=trim \
  -g AGATGTGTATAAGAGACAG \
  -G AGATGTGTATAAGAGACAG \
  -o ${OUT_DIR}/clean/${SAMPLE}_R1_clean.fq.gz \
  -p ${OUT_DIR}/clean/${SAMPLE}_R2_clean.fq.gz \
  ${OUT_DIR}/round3_Q5Q7/${SAMPLE}_round3_Q5Q7_R1.fq.gz ${OUT_DIR}/round3_Q5Q7/${SAMPLE}_round3_Q5Q7_R2.fq.gz \
  --untrimmed-o ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_clean_R1.fq.gz \
  --untrimmed-p ${OUT_DIR}/unassigned/${SAMPLE}_unassigned_clean_R2.fq.gz \
  &>> ${OUT_DIR}/logs/${SAMPLE}_cutadapt_clean_ME.log


### Final output FASTQs are: ${OUT_DIR}/clean/${SAMPLE}_R1_clean.fq.gz and _R2_clean.fq.gz

echo "Demultiplexing complete for ${SAMPLE}."
