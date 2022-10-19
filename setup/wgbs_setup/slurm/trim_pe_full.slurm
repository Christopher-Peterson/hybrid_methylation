#!/bin/bash

#-------------------------------------------------------
#SBATCH -J trim_pe
#SBATCH -N 4 
#SBATCH -n 12
#SBATCH -p normal
#SBATCH -o logs/trim_pe.o
#SBATCH -e logs/trim_pe.e
#SBATCH -t 00:40:00
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------
# TRIM PAIRED-END READS

module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`

CORES_PER_JOB=42
# Populate the job file
export LAUNCHER_JOB_FILE=jobs/trim_pe_full
: > $LAUNCHER_JOB_FILE

# Create the job script
for file in $SCRATCH/reads/raw/*_R2_001.fastq
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
    $R2IN" >> $LAUNCHER_JOB_FILE
  done

# Run the launcher
$LAUNCHER_DIR/paramrun

echo "Done with run"
cds hybrid_methylation/reads
mkdir trimmed
mv raw/*trim trimmed
rm trimmed/*trim.temp.*
