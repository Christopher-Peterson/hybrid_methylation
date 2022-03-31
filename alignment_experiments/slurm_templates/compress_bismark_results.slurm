#!/bin/bash

#-------------------------------------------------------
#SBATCH -J compress  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -p normal
#SBATCH -o logs/compress.o
#SBATCH -e logs/compress.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Compress the methyl_extract files in parallel

ml launcher

IN_DIR="methyl_extract"

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/compress_bismark
> $LAUNCHER_JOB_FILE

add_job = function {
  local FILE=$1
  local MESSAGE="$FILE Complete"
  echo "pigz -k -9 $FILE && echo $MESSAGE" >> LAUNCHER_JOB_FILE
}
# Create the job file
# Assemble the commands
for file in $IN_DIR/CpG_context*
  do add_job $file
done
for file in $IN_DIR/Non_CpG_context*
  do add_job $file
done
for file in $IN_DIR/*CpG_report.txt
  do add_job $file
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

