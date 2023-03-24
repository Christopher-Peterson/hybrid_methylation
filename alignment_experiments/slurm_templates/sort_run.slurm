#!/bin/bash

#-------------------------------------------------------
#SBATCH -J sort_run
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -p normal
#SBATCH -o logs/sort_run.o
#SBATCH -e logs/sort_run.e
#SBATCH -t 1:00:00  
#------------------------------------------------------

# Sort and Index the bam files

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/sort_run

## Figure out arguments ####

SCRIPT="scripts/sort_index.sh"
mkdir -p bams/sorted
# Check to see if the full output files exist & have non-zero size
# If not, add the bam file to the job list
function add_job {
  local FULL_FILE=$1
  local FILE=$(basename $FULL_FILE) # Remove directories & suffix for pattern matching
  # echo $FILE
  # There should be a bai file here if the script has already run
  local BAI_FILE=bams/sorted/${FILE}.bai
  # echo $BAI_FILE
#  if [ ! -s $BAI_FILE ]; then # If it doesn't exist/is empty, add the job
    # Append the job for file to the job script
#   echo "$SCRIPT $FILE"   >> $LAUNCHER_JOB_FILE
  fi
}


# Assemble the job file
> $LAUNCHER_JOB_FILE
for file in bams/dedup/lane*deduplicated.bam
do add_job $file
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

