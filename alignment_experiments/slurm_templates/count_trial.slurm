#!/bin/bash

#-------------------------------------------------------
#SBATCH -J count_trial  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/count_trial.o
#SBATCH -e logs/count_trial.e
#SBATCH -t 2:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# Sort and index a bam file
module unload xalt

FILE="lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam"

mkdir -r bams/counts

IGVTOOLS=$WORK/IGV_Linux_2.12.3/igvtools
IN_FILE=bams/sorted/$FILE
OUT_FILE=bams/counts/${FILE/bam/wig}
source ../../genome.sh
$IGVTOOLS count $IN_FILE $OUT_FILE $GENOME_PATH

# ml tacc-singularity
# IMAGE=/work/04386/crpeters/ls6/singularity/r-tidyverse-optparse_4.1.2.sif
# HDIR=/home
# singularity shell -H $HDIR $IMAGE
# R
