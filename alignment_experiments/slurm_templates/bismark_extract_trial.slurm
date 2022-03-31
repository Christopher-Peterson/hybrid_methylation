#!/bin/bash

#-------------------------------------------------------
#SBATCH -J methyl_extract_trial  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/extract_trial.o
#SBATCH -e logs/extract_trial.e
#SBATCH -t 2:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# This is an attempt to make a generic slurm script that can handle both single and paired ends

## Figure out arguments ####

OUT_DIR="methyl_extract"

FILE="bams/dedup/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam"

scripts/bismark_extract.sh $FILE $OUT_DIR

