#!/bin/bash

#-------------------------------------------------------
#SBATCH -J sort_trial  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/sort_trial.o
#SBATCH -e logs/sort_trial.e
#SBATCH -t 2:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# Sort and index a bam file

FILE="lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam"

scripts/sort_index.sh $FILE