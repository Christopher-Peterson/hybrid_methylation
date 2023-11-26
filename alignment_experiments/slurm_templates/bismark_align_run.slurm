#!/bin/bash

#-------------------------------------------------------
#SBATCH -J sed_NAME  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N sed_NODES
#SBATCH -n sed_NODES
#SBATCH -p normal
#SBATCH -o logs/bismark_run.o
#SBATCH -e logs/bismark_run.e
#SBATCH -t sed_TIME
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# This is an attempt to make a generic slurm script that can handle both single and paired ends

TYPE=$(cat config/ends.type)
N_CORE=${N_CORE:-24} # Adjust via passing env variables

PARMS=$(cat config/bismark.parms) # Read parm files
ml launcher

### Functions
# Checks to see if a BAM file exists; if it doesn't adds the script to create it (from $FILES) to the job file
function add_job_if_no_bam {
  # This requires BAM, FILES, SCRIPT, TEMP_DIR, N_CORE, and LAUNCHER_JOB_FILE to be defined in calling environment
  if [ ! -f ${BAM} ]; then # IF $BAM does not exist
     # Append the job for file to the job script
     echo "TEMP_DIR=$TEMP_DIR N_CORE=$N_CORE PARMS=\"$PARMS\" $SCRIPT $FILES" >> $LAUNCHER_JOB_FILE
  fi
}

# Define differences of pe and se; each gets a different setup_job function
if [ "$TYPE" == "pe" ]; then
  SCRIPT="scripts/bismark_align_pe.sh"
  function setup_job {
    FILE_1=$(basename $1)
    FILE_2=${FILE_1/_1.trim/}_2.trim
    # Using a temp directory in scratch
    TEMP_DIR=temp/${FILE_1}_wd
    BAM=bams/${FILE_1}_bismark_bt2_pe.bamd
    FILES="$FILE_1 $FILE_2"
    add_job_if_no_bam
  }
else
  SCRIPT="scripts/bismark_align_se.sh"
  TEMP_DIR=/tmp/bismark
  function setup_job {
    FILE=$(basename $1)
    # Using a temp directory in scratch
    BAM=bams/${FILE}_bismark_bt2.bam
    FILES=$FILE
    add_job_if_no_bam
  }
fi

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/bismark_run
> $LAUNCHER_JOB_FILE



# Create the job file
# Assemble the commands
for file in reads/trimmed/lane1-A*_1.trim
  do setup_job $file
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

