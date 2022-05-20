#!/bin/bash

# Visually inspect the alignment counts over the genome

# Things that should be done first:
# mkdir -r bams/counts

# Stagess:
# 1. Sort and Index the bam file
# 2. Generate the count file (.wig) w/ IGVtools
# 3. Convert wig file into a data frame
# 4. Create visualization and link it to the summary directory

# Input (should be deduplicated)
FILE_DEFAULT="lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam"
if [ $(cat config/ends.type) == "se" ]; then
  FILE_DEFAULT=${FILE_DEFAULT/_pe/}
fi
FILE=${1:-$FILE_DEFAULT}
STAGE=${2:-"1"}
# Derivatives
SORT_FILE=bams/sorted/$FILE # IGV Input
WIG_FILE=bams/counts/${FILE/bam/wig} # IGV Output
RDS_FILE=${WIG_FILE/wig/rds}
FIG_FILE=${RDS_FILE/rds/png}
# Some useful paths
IGVTOOLS=$WORK/IGV_Linux_2.12.3/igvtools
source $SCRATCH/hybrid_methylation/genome.sh


# Sort and Index the Bam
if [ $STAGE -le "1" ]; then
  scripts/sort_index.sh $FILE
fi
## Run IGVtools to create the WIG file

# Make WIG file
if [ $STAGE -le "2" ]; then
  mkdir -p bams/counts
  $IGVTOOLS count $SORT_FILE $WIG_FILE $GENOME_PATH
fi
# Turn the wig file into a data frame
if [ $STAGE -le "3" ]; then
  TAPESTRY_ROWS="3400"
  r-tidy-opt scripts/read_counts.r $WIG_FILE $RDS_FILE  $TAPESTRY_ROWS # Fgure out other parameters
fi
# Turn the data frame into a figure
if [ $STAGE -le "4" ]; then
  TAPESTRY_THRESHOLDS="30 200"
  r-tidy-opt scripts/plot_counts.r $RDS_FILE $FIG_FILE  $TAPESTRY_THRESHOLDS # Fgure out other parameters
  # Link the png file to alignment_experiments/summary
  WD=$(pwd)
  cd ../summary
  NAME=$(basename $WD)
  ln -s -f $WD/$FIG_FILE ${NAME}.png
  cd $WD
fi
