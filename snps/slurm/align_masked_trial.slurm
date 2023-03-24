#!/bin/bash

#-------------------------------------------------------
#SBATCH -J align_offspring_trial  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/align_masked_trial.o
#SBATCH -e logs/align_masked_trial.e
#SBATCH -t 5:00:00
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# This is an attempt to make a generic slurm script that can handle both single and paired ends

N_CORE=${N_CORE:-24} # Adjust via passing env variables

PARMS=$(cat config/bismark.parms) # Read parm files
module unload xalt
ID=X5
N_CORE=$N_CORE PARMS="$PARMS" scripts/align_masked.sh $ID


