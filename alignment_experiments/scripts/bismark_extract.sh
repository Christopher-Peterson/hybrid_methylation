#!/bin/bash
# Extract methylation reports via bismark

FILE=$1
OUT_DIR=${2:-"methyl_extract"}
GENOME=${3:-genome.sh} # Could also be genome2.sh

BASE_DIR=$SCRATCH/hybrid_methylation
LOCAL_DIR=$(pwd)
source $LOCAL_DIR/$GENOME

## Figure out arguments ####
# Number of parallel cores
PROCESSES="--multicore 42" # Uses 3 cores per process, so this is the max on LS6

FOLDER_ARGS="-o $OUT_DIR --genome_folder $GENOME_FOLDER"

# Handle paired vs. single ends
TYPE=$(cat config/ends.type)
# Remove the "e" from "pe" or "se" and add the dash to the front
ENDS_ARG="-${TYPE/e/}" # Should be -p or -s

# read the type to determine -p or -s
SCAFFOLD_ARGS=""
# Are scaffolds really necessary? It looks like a lot of files can be open at once
# SCAFFOLD_ARGS="--scaffolds --buffer_size 89%"

OUTPUT_ARGS="--merge_non_CpG  --comprehensive --cytosine_report"

module load tacc-singularity 
module unload xalt 

# I'm having some alias issues
# Irksome
# singularity exec -H /home \
#   /work/04386/crpeters/ls6/nextflow/work/singularity/nfcore-methylseq-1.6.1.img \
#   bismark_methylation_extractor \
#   $ENDS_ARG $FOLDER_ARGS $SCAFFOLD_ARGS $OUTPUT_ARGS $PROCESSES $FILE

bismark_methylation_extractor $ENDS_ARG $FOLDER_ARGS $SCAFFOLD_ARGS $OUTPUT_ARGS $PROCESSES $FILE

# Should probably add some cleanup to this afterward