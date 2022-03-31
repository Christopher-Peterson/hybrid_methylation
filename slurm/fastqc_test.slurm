#!/bin/bash

#-------------------------------------------------------
#SBATCH -J fastqc
#SBATCH -N 1 
#SBATCH -n 124
#SBATCH -p development
#SBATCH -o logs/fastqc.o
#SBATCH -e logs/fastqc.e
#SBATCH -t 01:00:00
#### #SBATCH -A tagmap

#------------------------------------------------------
module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`

OUT_DIR=fastqc_results_raw

# Populate the job file
JOB_FILE=jobs/fastqc
: > $JOB_FILE

for file in fastqs/*.fastq
do echo "fastqc -o $OUT_DIR/ -f fastq $file" >> $JOB_FILE
done

# Run the launcher
export LAUNCHER_JOB_FILE=$JOB_FILE
$LAUNCHER_DIR/paramrun

