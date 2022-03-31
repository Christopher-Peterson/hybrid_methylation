#!/bin/bash

#-------------------------------------------------------
#SBATCH -J trim_pe
#SBATCH -N 4 
#SBATCH -n 12
#SBATCH -p development
#SBATCH -o logs/trim_pe.o
#SBATCH -e logs/trim_pe.e
#SBATCH -t 02:00:00
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------
# TRIM PAIRED-END READS

module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`

CORES_PER_JOB=42
# Populate the job file
JOB_FILE=jobs/trim_pe
: > $JOB_FILE

# Create the job script
for file in fastqs/*_R2_001.fastq
do R1IN=${file/_R2_001.fastq/}_R1_001.fastq
R2IN=${file}
R1OUT=${file/_R2_001.fastq/}_1.trim
R2OUT=${file/_R2_001.fastq/}_2.trim
echo "cutadapt -j $CORES_PER_JOB \
  -a GATCGGAAGAGCA \
  -A GATCGGAAGAGCA \
  -a AGATCGGAAGAGC \
  -A AGATCGGAAGAGC \
  --minimum-length 20 \
  -q 20 \
  -o $R1OUT \
  -p $R2OUT \
  $R1IN \
  $R2IN" >> $JOB_FILE
done

# Run the launcher
export LAUNCHER_JOB_FILE=$JOB_FILE
$LAUNCHER_DIR/paramrun

echo "Done with run"
mkdir fastqs_trimmed
mv fastqs/*trim fastqs_trimmed