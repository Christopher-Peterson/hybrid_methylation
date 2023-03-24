#!/bin/bash

ID=${1:-X5}
PARMS=${PARMS:-"-N 1 --score_min L,0,-0.8 -np 0"}
PARENTS=`cat parents.txt | grep $ID | grep -o A..`
echo $PARENTS

# GENOME=masked_genomes/registry/${PARENTS}.sh
# We'll be changin directories, so let's save this command w/ an absolute path
BISMARK2=$PWD/scripts/bismark2


# Reads
BASE_DIR=$SCRATCH/hybrid_methylation
FILE1="${BASE_DIR}/reads/trimmed/lane?-${ID}-*_1.trim"
FILE2="${BASE_DIR}/reads/trimmed/lane?-${ID}-*_2.trim"

N_CORE=${N_CORE:-24} # Adjust via passing env variables
LOCAL_DIR=$PWD
GENOME_FOLDER=`ls -d $PWD/masked_genomes/$PARENTS/A?_A?_dual_hybrid*N-masked`
# source $LOCAL_DIR/$GENOME

# Temporary directory should either be on /tmp or /scratch
TEMP_DIR=$LOCAL_DIR/temp/bismark/$ID
mkdir -p $TEMP_DIR

cd $TEMP_DIR
# Link to temp directory
ln -sf $FILE1 "${ID}_1.fq"
ln -sf $FILE2 "${ID}_2.fq"

$BISMARK2 \
--parallel $N_CORE \
--bowtie2  \
--non_directional $PARMS \
$GENOME_FOLDER -1 $FILE1 -2 $FILE2


mv *.bam $LOCAL_DIR/offspring/bams
mv *.txt $LOCAL_DIR/offspring/reports

