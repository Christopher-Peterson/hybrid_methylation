#!/bin/bash
# Align a trimmed read file w/ bismark; 

# Pass $PARMS as env variable to adjust parameters

# First argument:
FILE1=${1:-lane1-A1-A_S1_L001_1.trim}
FILE2=${2:-lane1-A1-A_S1_L001_2.trim}
GENOME=${3:-genome.sh} # Could also be genome2.sh

N_CORE=${N_CORE:-24} # Adjust via passing env variables
# Loads GENOME_FILE
PARMS=${PARMS:-""}

BASE_DIR=$SCRATCH/hybrid_methylation
LOCAL_DIR=$(pwd)

source $GENOME

# Temporary directory should either be on /tmp or /scratch
mkdir -p ${TEMP_DIR:=/tmp/bismark} # Pass TEMP_DIR as an env variable to change
cd $TEMP_DIR

# Link to temp directory
ln -s $LOCAL_DIR/reads/trimmed/$FILE1 $FILE1
ln -s $LOCAL_DIR/reads/trimmed/$FILE2 $FILE2
# Would this be faster copying it instead?

bismark \
 --parallel $N_CORE \
--bowtie2  \
--non_directional $PARMS \
$GENOME_FOLDER -1 $FILE1 -2 $FILE2


mv *.bam $LOCAL_DIR/bams
mv *.txt $LOCAL_DIR/reports

