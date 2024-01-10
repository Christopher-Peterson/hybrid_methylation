#!/bin/bash
ID=${1}
# BASE_FILE=`ls offspring/bams/lane?-$ID-*bam`
# DD_OUT_DIR=${2:-"offspring/bams/dedup"}
N_CORES=16
# 
# # Remove the "e" from "pe" or "se" and add the dash to the front
# ENDS_ARG=-p
# 
# ARGS="--output_dir $DD_OUT_DIR --outfile $ID $ENDS_ARG --bam  $BASE_FILE"
# deduplicate_bismark $ARGS
# 
# # Now sort
DEDUP_FILE=offspring/bams/dedup/lane*${ID}*.deduplicated.bam
OUT_FILE=offspring/sorted_bams/${ID}_sorted.bam
# mkdir -p $IN_DIR $OUT_DIR

# Sort
SORT_ARGS="-o $OUT_FILE -@$N_CORES $DEDUP_FILE"   #$DD_OUT_DIR/${ID}.deduplicated.bam"
samtools sort $SORT_ARGS

# Index
samtools index -b -@$N_CORES $OUT_FILE


# Exit w/ error if fail
if [ ! -f $OUT_FILE ]; then 
  exit 1
fi
# # Move the report
# REPORT=$OUT_DIR/${FILE/.bam/.deduplication_report.txt}
# mv $REPORT reports/dedup
