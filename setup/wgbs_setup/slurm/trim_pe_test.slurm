#!/bin/bash

#-------------------------------------------------------
#SBATCH -J trim_pe
#SBATCH -N 1 
#SBATCH -n 12
#SBATCH -p development
#SBATCH -o logs/trim_pe_test.o
#SBATCH -e logs/trim_pe_test.e
#SBATCH -t 00:20:00

#------------------------------------------------------
# TRIM PAIRED-END READS

module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`


# Populate the job file
export LAUNCHER_JOB_FILE=jobs/trim_pe_test
> $LAUNCHER_JOB_FILE

# Create the job script
for file in trim_test/fastqs/*_R2_001.fastq
do R1IN=${file/_R2_001.fastq/}_R1_001.fastq
   R2IN=${file}
   R1OUT=${file/_R2_001.fastq/}_1.trim
   R2OUT=${file/_R2_001.fastq/}_2.trim
  echo "cutadapt \
  -a GATCGGAAGAGCA \
  -A GATCGGAAGAGCA \
  -a AGATCGGAAGAGC \
  -A AGATCGGAAGAGC \
  --minimum-length 20 \
  -q 20 \
  -o $R1OUT \
  -p $R2OUT \
  $R1IN \
  $R2IN" >> $LAUNCHER_JOB_FILE
done

# Run the launcher
$LAUNCHER_DIR/paramrun

echo "Done with run"
cd trim_test
mkdir fastqs_trimmed
mv fastqs/*trim fastqs_trimmed
