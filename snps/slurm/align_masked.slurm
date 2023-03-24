#!/bin/bash

#-------------------------------------------------------
#SBATCH -J align_offspring  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 8
#SBATCH -n 8
#SBATCH -p normal
#SBATCH -o logs/align_offspring.o
#SBATCH -e logs/align_offspring.e
#SBATCH -t 5:30:00
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# This is an attempt to make a generic slurm script that can handle both single and paired ends

N_CORE=${N_CORE:-24} # Adjust via passing env variables

PARMS="$(cat config/bismark.parms)" # Read parm files
ml launcher
module unload xalt
mkdir -p offspring/bams/dedup offspring/reports offspring/sorted_bams

### Functions
# Define differences of pe and se; each gets a different setup_job function
SCRIPT="scripts/align_masked.sh"
function setup_job {
    local ID=$1
    # Using a temp directory in scratch
    BAM=offspring/bams/lane?-${ID}-?_S?_L00?_1.trim_bismark_bt2_pe.bam
    if [ ! -f ${BAM} ]; then # IF $BAM does not exist
     # Append the job for file to the job script
     echo "N_CORE=$N_CORE PARMS=\"$PARMS\" $SCRIPT $ID" >> $LAUNCHER_JOB_FILE
  fi
  }

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/align_masked.sh
> $LAUNCHER_JOB_FILE
for id in {1..8}
  do setup_job X${id}
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

