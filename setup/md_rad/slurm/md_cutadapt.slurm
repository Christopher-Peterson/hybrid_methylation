#!/bin/bash

#-------------------------------------------------------
#SBATCH -J md_cutadapt
#SBATCH -N 1 
#SBATCH -n 76
#SBATCH -p normal
#SBATCH -o logs/md_cutadapt.o
#SBATCH -e logs/md_cutadapt.e
#SBATCH -t 00:15:00
#### #SBATCH -A tagmap

#------------------------------------------------------
module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`
export LAUNCHER_JOB_FILE=jobs/trim

mkdir -p trimmed

function add_job {
  local file=$1
  local out="trimmed/`basename ${file/_filt0/_trim}`"
  if [ ! -z $out ]; then 
    echo "cutadapt \
      -u 6 \
      -j 16 \
      -a GATCGGAAGAGCA \
      -a AGATCGGAAGAGC \
      --minimum-length 20 \
      -o $out $file ">> $LAUNCHER_JOB_FILE
  fi
}

# Populate the job file
>$LAUNCHER_JOB_FILE
for file in filt0/*fastq; do 
  add_job $file
done

# Run the launcher
$LAUNCHER_DIR/paramrun


