#!/bin/bash
# Align a trimmed read file w/ bismark; 

# Pass $PARMS as env variable to adjust parameters
module unload xalt
# First argument:
DIREC=${1:-lane2-X2-E_S7_L002_1.trim_wd}
# FILE1=${1:-lane1-A1-A_S1_L001_1.trim}
# FILE2=${2:-lane1-A1-A_S1_L001_2.trim}
GENOME=${2:-genome2.sh} # Could also be genome2.sh
BASE_FILE=`echo $DIREC | sed 's/_1.trim_wd//' -`
FILE1="$DIREC/${BASE_FILE}_1.trim_unmapped_reads_1.fq.gz"
FILE2="$DIREC/${BASE_FILE}_2.trim_unmapped_reads_2.fq.gz"



N_CORE=${N_CORE:-24} # Adjust via passing env variables
# Loads GENOME_FILE
PARMS=${PARMS:-""}

BASE_DIR=$SCRATCH/hybrid_methylation
LOCAL_DIR=$(pwd)

source $LOCAL_DIR/$GENOME

# Temporary directory should either be on /tmp or /scratch
# mkdir -p ${TEMP_DIR:=/tmp/bismark} # Pass TEMP_DIR as an env variable to change
# cd $TEMP_DIR
cd temp

# Link to temp directory
# ln -s $LOCAL_DIR/reads/trimmed/$FILE1 $FILE1
# ln -s $LOCAL_DIR/reads/trimmed/$FILE2 $FILE2
# Would this be faster copying it instead?

bismark \
 --parallel $N_CORE \
--bowtie2  \
--non_directional $PARMS \
$GENOME_FOLDER -1 $FILE1 -2 $FILE2


mv *.bam $LOCAL_DIR/bams/align2
mv *.txt $LOCAL_DIR/reports/align2

