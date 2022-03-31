#!/bin/bash

#-------------------------------------------------------
#SBATCH -J methyl_extract_run
#SBATCH -N 11
#SBATCH -n 11
#SBATCH -p normal
#SBATCH -o logs/extract_run.o
#SBATCH -e logs/extract_run.e
#SBATCH -t 2:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# This is an attempt to make a generic slurm script that can handle both single and paired ends

ml launcher
# Setup the job file
export LAUNCHER_JOB_FILE=jobs/extract_run

## Figure out arguments ####

OUT_DIR="methyl_extract"
SCRIPT="scripts/bismark_extract.sh"

# Check to see if the full output files exist & have non-zero size
# If not, add the bam file to the job list
function add_job {
  local FILE=$1
  local NAME=$(basename ${FILE/.bam/}) # Remove directories & suffix for pattern matching
  
  # There should be eight files that match $NAME in $OUT_DIR
  # Check that there are 8 of them and that each has a size > 0
  # If any of them are, set ADD to TRUE
  if [ $(ls $OUT_DIR/*$NAME* 2> /dev/null | wc -l) == 8 ]; then # The 2> redirects errors to /dev/null so that zeroes can be handled
    ADD="FALSE"
    for out_file in $OUT_DIR/*$NAME* 
    do if [ -s out_file ]; then # IF out file is size 0, add the job
      ADD="TRUE"
      fi
    done
  else # Not 8 matching files
    ADD="TRUE"
  fi

  if [ $ADD == "TRUE" ]; then # IF the report file has a size > 0
  # Append the job for file to the job script
  echo "$SCRIPT $FILE $OUT_DIR" >> $LAUNCHER_JOB_FILE
  fi
}


# Assemble the job file
> $LAUNCHER_JOB_FILE
for file in bams/dedup/*deduplicated.bam
do add_job $file
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

