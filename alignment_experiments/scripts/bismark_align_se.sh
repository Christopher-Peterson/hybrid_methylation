#!/bin/bash
# Align a trimmed read file w/ bismark, using the /tmp directory for faster processing
# Pass $PARMS as env variable to adjust parameters

# First argument:
FILE=${1:-lane1-A1-A_S1_L001_2.trim}
FILE2=${2:-lane1-A1-A_S1_L001_2.trim}
GENOME=${3:-genome.sh} # Could also be genome2.sh


module unload xalt
# module load launcher
# export LAUNCHER_WORKDIR=`pwd`

PARMS=${PARMS:-""}

N_CORE=${N_CORE:-32} # Adjust via passing env variables
# Loads GENOME_FILE

BASE_DIR=$SCRATCH/hybrid_methylation
SAVE_DIR=`pwd`

source $SAVE_DIR/$GENOME
# source $BASE_DIR/singularity_aliases.sh

# We're running this on /tmp to increase speed
mkdir /tmp/bismark
cd /tmp/bismark

# Link to temp directory
ln -s $BASE_DIR/reads/trimmed/$FILE $FILE
# Would this be faster copying it instead?

bismark \
--parallel $N_CORE \
--non_directional $PARMS\
--bowtie2  \
$GENOME_FOLDER $FILE


mv *.bam $SAVE_DIR/bams
mv *.txt $SAVE_DIR/reports

