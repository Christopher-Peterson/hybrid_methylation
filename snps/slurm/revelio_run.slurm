#!/bin/bash

#-------------------------------------------------------
#SBATCH -J revelio_run  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 3
#SBATCH -n 3
#SBATCH -p normal
#SBATCH -o logs/revelio_run.o
#SBATCH -e logs/revelio_run.e
#SBATCH -t 1:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/revelio

# Sort and index a bam file
module unload xalt

mkdir -p temp
N_THREAD=128
source genome.sh
# BAM_IN=bams/lane1-A1-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam
# BAM_OUT=masked_bams/A1_masked.bam
# pysam scripts/revelio.py -f $GENOME_PATH -T 128 -t temp $BAM_IN $BAM_OUT &
# 
# 

## Figure out arguments ####

SCRIPT="pysam scripts/revelio.py -f $GENOME_PATH -T $N_THREAD -t temp"
# Check to see if the full output files exist & have non-zero size
# If not, add the bam file to the job list
function add_job {
  local BAM_IN=$1
  local SAMPLE=`echo $BAM_IN | grep -o 'A[1-4]'` # Remove directories & suffix for pattern matching
  # There should be a bai file here if the script has already run
  local BAM_OUT=masked_bams/"${SAMPLE}_masked.bam"
  # echo $BAI_FILE
  if [ ! -s $BAM_OUT ]; then # If it doesn't exist/is empty, add the job
    # Append the job for file to the job script
    echo "$SCRIPT $BAM_IN $BAM_OUT && samtools index -b -@$32 $BAM_OUT" >> $LAUNCHER_JOB_FILE
  fi
}

# Assemble the job file
> $LAUNCHER_JOB_FILE
for file in bams/lane?-A*bam
do add_job $file
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

