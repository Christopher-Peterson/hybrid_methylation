#!/bin/bash

#-------------------------------------------------------
# SBATCH -J dss_pairwise  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 10
#SBATCH -n 10
#SBATCH -p normal
#SBATCH -o logs/dss_pairwise.o
#SBATCH -e logs/dss_pairwise.e
#SBATCH -t 4:05:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=$PWD/jobs/dss_pairwise

r-methylkit scripts/dss_setup_pairwise.r $LAUNCHER_JOB_FILE pairwise_out

export LAUNCHER_WORKDIR=$PWD
$LAUNCHER_DIR/paramrun
# # 540906
# sbatch --dependency=afterok:540918 slurm/calmeth.slurm
