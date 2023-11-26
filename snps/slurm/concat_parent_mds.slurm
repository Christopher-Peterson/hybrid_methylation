#!/bin/bash

#-------------------------------------------------------
#SBATCH -J concat_parent_mds  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p normal
#SBATCH -o logs/concat_parent_mds.o
#SBATCH -e logs/concat_parent_mds.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Concat the paired md bam reads of each parent; also, rename them to follow the 
# Axy format (or should I not...)

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/concat_md_bams
source genome.sh

# Pull in the regions as an array
N_CORES=32

# What should I do w/ the two pairs of reads? 
# Concat them?

# Parent_offspring lists
mkdir -p md_bam_merge

add_job() {
  local md1=$1
  local md2=${md1/-md1_/-md2_}
  local main_name=`echo $md1 | grep -Po '[s,m][0-9]+' `
  local out="md_bam_merge/${main_name}.bam"
  echo "samtools merge -@${N_CORES} $out $md1 $md1 && samtools index $out" >> $LAUNCHER_JOB_FILE
}

# Assemble the job file
> $LAUNCHER_JOB_FILE
for bam in md_bams/A-*-md1_*.bam; do
add_job $bam
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

