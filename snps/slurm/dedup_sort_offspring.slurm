#!/bin/bash

#-------------------------------------------------------
# SBATCH -J dedup_sortoffspring  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 2
#SBATCH -n 8
#SBATCH -p normal
#SBATCH -o logs/dedup_sort.o
#SBATCH -e logs/dedup_sort.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=$PWD/jobs/dedup_sort

mkdir -p offspring/bams/dedup offspring/sorted_bams 
function add_job { 
  local ID=${1}
  local RAW_BAM=`ls offspring/bams/lane*${ID}*bam`
  local DUP_OUT=offspring/bams/dedup/lane?-${ID}*report.txt
  local SORT_OUT=offspring/sorted_bams/${ID}_sorted.bam.bai
  
  if [ ! -s $DUP_OUT ]; then
    echo "scripts/dedup_offspring.sh $RAW_BAM && scripts/sort_offspring.sh $ID"  >> $LAUNCHER_JOB_FILE  
  elif [ ! -s $SORT_OUT ]; then # deduped, but not sorted
    echo "scripts/sort_offspring.sh $ID" >> $LAUNCHER_JOB_FILE  #&& " 
  fi
}

for i in {1..8}; do 
  add_job X$i 
done


export LAUNCHER_WORKDIR=$PWD
$LAUNCHER_DIR/paramrun
# # 540906
# sbatch --dependency=afterok:540918 slurm/calmeth.slurm
