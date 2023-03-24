#!/bin/bash
# De-duplicate
# 
# For some reason, the alias isn't working

FILE=${1}
OUT_DIR="offspring/bams/dedup"
mkdir -p $OUT_DIR
# Handle paired vs. single ends

ENDS_ARG="-p" # Should be -p or -s

ARGS="--output_dir $OUT_DIR $ENDS_ARG --bam  $FILE"
deduplicate_bismark $ARGS

OUT_FILE=`ls $OUT_DIR/*$FILE*deduplication_report.txt`

# Exit w/ error if fail
if [ ! -f $OUT_FILE ]; then 
  exit 1
fi
# Move the report
# REPORT=$OUT_DIR/${FILE/.bam/.deduplication_report.txt}
# mv $REPORT reports/dedup

# This line works:
# deduplicate_bismark --output_dir bams/dedup -p --bam bams/lane1-A1-A_S1_L001_1.trim_bismark_bt2_pe.bam
# deduplicate_bismark --output_dir bams/dedup -p --bam bams/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.bam
# 
# 
# 
# TO DO NEXT:
# 1. Dedup all offspring
# 2. Sort/index all offspring
# 3. calmeth all offspring
# 4. MethHaplo 