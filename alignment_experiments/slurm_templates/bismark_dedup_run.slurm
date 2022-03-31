#!/bin/bash

#-------------------------------------------------------
#SBATCH -J dedup_run
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p normal
#SBATCH -o logs/dedup_run.o
#SBATCH -e logs/dedup_run.e
#SBATCH -t 2:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# This is an attempt to make a generic slurm script that can handle both single and paired ends

ml launcher

OUT_DIR="bams/dedup"
REPORT_DIR="reports/dedup"
mkdir -p $OUT_DIR $REPORT_DIR

# read the type to determine -p or -s
SCRIPT="scripts/dedup.sh"

# Check to see if the output text file has non-zero size
# If so, add the bam file to the job list
function add_job {
  FILE=$(basename $1)
  REPORT=$REPORT_DIR/${FILE/.bam/.deduplication_report.txt}
  if [ ! -s ${REPORT} ]; then # IF the report file has a size > 0
     # Append the job for file to the job script
     echo "$SCRIPT $FILE" $OUT_DIR >> $LAUNCHER_JOB_FILE
  fi
}

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/dedup_run
> $LAUNCHER_JOB_FILE

# Create the job file
# Assemble the commands
for file in bams/*.bam
  do add_job $file
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

