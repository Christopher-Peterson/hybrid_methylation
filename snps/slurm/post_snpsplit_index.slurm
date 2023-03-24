#!/bin/bash

#-------------------------------------------------------
# SBATCH -J post_snpsplit_index  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p normal
#SBATCH -o logs/post_snpsplit_index.o
#SBATCH -e logs/post_snpsplit_index.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=$PWD/jobs/post_snpsplit_index
N_CORES=8

mkdir -p offspring/split_sort
function add_job { 
  local in_bam=${1}
  local out_bam=${in_bam/snpsplit/split_sort}
  local out_bam=${out_bam/.genome/_genome} # Replace dot w/ underscore, as required for InformMe
  if [ ! -s ${out_bam}.bai ]; then
    echo "samtools sort -o $out_bam -@$N_CORES $in_bam && samtools index -b -@$N_CORES $out_bam"  >> $LAUNCHER_JOB_FILE  
  fi
}

for bam in offspring/snpsplit/X?.genome?.bam; do 
  add_job $bam
done


export LAUNCHER_WORKDIR=$PWD
$LAUNCHER_DIR/paramrun
# # 540906
# sbatch --dependency=afterok:540918 slurm/calmeth.slurm
