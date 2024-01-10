#!/bin/bash

#-------------------------------------------------------
# SBATCH -J dss_independent  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p development
#SBATCH -o logs/dss_independent.o
#SBATCH -e logs/dss_independent.e
#SBATCH -t 1:00:00  

#------------------------------------------------------
export LAUNCHER_JOB_FILE=$PWD/jobs/independent

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt

# Setup the job file

RSCRIPT=r-methylkit
# Maternal/Paternal information; in the repo
parental_roles=$SCRATCH/hybrid_methylation/metadata/dam_sire.csv
# Which allele-specific genomes are attached to which parents
genome_order=offspring_parents.tsv

# Setup runs
$RSCRIPT scripts/dss_independent_setup.r "$parental_roles" "$genome_order"

> $LAUNCHER_JOB_FILE
for i in {1..8}; do
  echo "$RSCRIPT scripts/dss_independent_run.r X${i}" >> $LAUNCHER_JOB_FILE
done

export LAUNCHER_WORKDIR=$PWD
$LAUNCHER_DIR/paramrun
