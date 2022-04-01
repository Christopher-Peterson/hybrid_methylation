#!/bin/bash
# De-duplicate
# 
# For some reason, the alias isn't working

FILE=${1}
OUT_DIR=${2:-"bams/dedup"}

# Handle paired vs. single ends
TYPE=$(cat config/ends.type)
# Remove the "e" from "pe" or "se" and add the dash to the front
ENDS_ARG="-${TYPE/e/}" # Should be -p or -s

ARGS="--output_dir $OUT_DIR $ENDS_ARG --bam  bams/$FILE"
deduplicate_bismark $ARGS

# Move the report
REPORT=$OUT_DIR/${FILE/.bam/.deduplication_report.txt}
mv $REPORT reports/dedup

# This line works:
# deduplicate_bismark --output_dir bams/dedup -p --bam bams/lane1-A1-A_S1_L001_1.trim_bismark_bt2_pe.bam