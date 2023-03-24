#!/bin/bash

#-------------------------------------------------------
#SBATCH -J mask_genomes  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p development
#SBATCH -o logs/mask_genomes.o
#SBATCH -e logs/mask_genomes.e
#SBATCH -t 0:30:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt
source genome.sh
mkdir -p masked_genomes

GENOME_NAME=Amil # Genericize if necessary
# Setup the job file
export LAUNCHER_JOB_FILE=jobs/mask_genomes

# Assemble the job file
> $LAUNCHER_JOB_FILE
for vcf in vcf/concat/A??_filtered.vcf; do
  echo "scripts/run_SNPsplit_genome_prep.sh $vcf $GENOME_NAME" >> $LAUNCHER_JOB_FILE
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

