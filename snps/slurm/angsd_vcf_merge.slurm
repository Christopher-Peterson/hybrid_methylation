#!/bin/bash

#-------------------------------------------------------
#SBATCH -J angsd_vcf_merge  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p normal
#SBATCH -o logs/angsd_vcf_merge.o
#SBATCH -e logs/angsd_vcf_merge.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/angsd_vcf_merge
# source genome.sh

mkdir -p vcf/concat

DIRECS=`ls vcf | grep -Po "A[1-4][1-4]"`
REGIONS=( `cat genome_regions.txt` )

# Assemble the files one at a time

function add_job {
  local DIR=$1
  local OUTPUT='vcf/concat/'`basename $DIR`".bcf"
  if [ ! -f $OUTPUT ]; then
    local CONCAT_CMD="bcftools concat -n -o $OUTPUT"
    for reg in "${REGIONS[@]}"; do
      local CONCAT_CMD="${CONCAT_CMD} vcf/${DIR}/vcf_mod_${reg}.bcf"
    done
    echo $CONCAT_CMD >> $LAUNCHER_JOB_FILE
  fi
}

# Assemble the job file
> $LAUNCHER_JOB_FILE
for dir in $DIRECS; do
    echo OUTPUT='vcf/concat/'`basename $dir`".bcf"
  # add_job $dir
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

