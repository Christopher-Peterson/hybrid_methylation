#!/bin/bash

#-------------------------------------------------------
#SBATCH -J fastqc_test
#SBATCH -N 1 
#SBATCH -n 48
#SBATCH -p development
#SBATCH -o logs/fastqc_test.o
#SBATCH -e logs/fastqc_test.e
#SBATCH -t 01:00:00
#### #SBATCH -A tagmap

#------------------------------------------------------
module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`


# Populate the job file
export LAUNCHER_JOB_FILE=jobs/fastqc_test
: > $LAUNCHER_JOB_FILE

# For the raw fastqs
IN_DIR=trim_test/fastqs
OUT_DIR=trim_test/fastqc_results_raw
for file in $IN_DIR/*.fastq
do echo "fastqc -o $OUT_DIR/ -f fastq $file" >> $LAUNCHER_JOB_FILE
done
# And the trimmed ones
IN_DIR=trim_test/fastqs_trimmed
OUT_DIR=trim_test/fastqc_results_trimmed
for file in $IN_DIR/*.fastq
do echo "fastqc -o $OUT_DIR/ -f fastq $file" >> $LAUNCHER_JOB_FILE
done


# Run the launcher
$LAUNCHER_DIR/paramrun

