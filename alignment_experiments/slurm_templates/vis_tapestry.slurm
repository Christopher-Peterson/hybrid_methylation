#!/bin/bash

#-------------------------------------------------------
#SBATCH -J vis_tapestry 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/tapestry.o
#SBATCH -e logs/tapestry.e
#SBATCH -t 1:00:00  
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

# Sort and index a bam file
module unload xalt

FILE="lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam"
if [ $(cat config/ends.type) == "se" ]; then
  FILE=${FILE/_pe/}
fi
STAGE=${STAGE:-"1"}
scripts/visualize_alignment_tapestry.sh $FILE $STAGE