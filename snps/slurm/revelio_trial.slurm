#!/bin/bash

#-------------------------------------------------------
#SBATCH -J sort_trial  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/sort_trial.o
#SBATCH -e logs/sort_trial.e
#SBATCH -t 5:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# Sort and index a bam file
module unload xalt

mkdir -p temp

source genome.sh
BAM_IN=bams/lane1-A1-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam
BAM_OUT=masked_bams/A1_masked.bam
pysam scripts/revelio.py -f $GENOME_PATH -T 128 -t temp $BAM_IN $BAM_OUT 

samtools index -b -@$32 $BAM_OUT