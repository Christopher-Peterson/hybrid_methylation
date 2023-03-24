#!/bin/bash

#-------------------------------------------------------
#SBATCH -J vcf_filter_parent  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p development
#SBATCH -o logs/vcf_filter_parent.o
#SBATCH -e logs/vcf_filter_parent.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/vcf_filter

function add_job {
  local BCF=$1
  local OUT="${BCF/.bcf/_filtered.vcf}"
  local FILTER='FMT/GP>0.95 && (FMT/GT="1/0" || FMT/GT="0/1)"'
  # 95% Genotype prob heterozygotes
  # Use SED to rename the samples in the files
  local SED="sed 's%masked_bams/%%g' - | sed 's/_masked.bam//g' - "
  echo "bcftools filter --include '${FILTER}' ${BCF} | $SED > ${OUT}" >> $LAUNCHER_JOB_FILE
}

# Assemble the job file
> $LAUNCHER_JOB_FILE
for bcf in vcf/concat/parents/A?.bcf; do
  add_job $bcf
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

