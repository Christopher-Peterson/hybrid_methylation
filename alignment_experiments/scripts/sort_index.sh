#!/bin/bash
# Sort and index a bam file

FILE=${1}
IN_DIR="bams/dedup"
OUT_DIR="bams/sorted"
N_CORES=128

# Sort
SORT_ARGS="-o $OUT_DIR/$FILE -@$N_CORES $IN_DIR/$FILE"
samtools sort $SORT_ARGS

# Index
samtools index -b -@$N_CORES $OUT_DIR/$FILE
