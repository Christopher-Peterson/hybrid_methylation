#!/bin/bash

#-------------------------------------------------------
#SBATCH -J vcf_filter_join_mask  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p development
#SBATCH -o logs/vcf_filter_join_mask.o
#SBATCH -e logs/vcf_filter_join_mask.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Filter the bcf files, join the mdRAD and WGBS BCF's, and run the genome mask script.

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/vcf_filter_join_mask


function add_job {
  local md_BCF=$1
  local wgbs_BCF="${md_BCF/concat/wgbs}"
  local bn=`basename $md_BCF`
  local name=${bn/.bcf/}
  if [ ! -d masked_genomes/$name ]; then
    echo "scripts/merge_md_wgbs_vcf.sh $name $md_BCF $wgbs_BCF" >> $LAUNCHER_JOB_FILE
  fi
}

# Assemble the job file
> $LAUNCHER_JOB_FILE
for bcf in vcf/concat/A*.bcf; do
  add_job $bcf
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

