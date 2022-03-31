#!/bin/bash

#-------------------------------------------------------
#SBATCH -J sed_NAME
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/bismark_trial.o
#SBATCH -e logs/bismark_trial.e
#SBATCH -t 05:00:00
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------

TYPE=$(cat config/ends.type)
PARMS=$(cat config/bismark.parms)
N_CORE=${N_CORE:-24}

FILE1=lane1-A2-B_S2_L001_1.trim

if [ "$TYPE" == "pe" ]; then
  FILE2=lane1-A2-B_S2_L001_2.trim
  SCRIPT="scripts/bismark_align_pe.sh"
  TEMP_DIR="temp/${FILE1}_wd"
else
  FILE2="" # not passed for se
  SCRIPT="scripts/bismark_align_se.sh"
  TEMP_DIR="/tmp/bismark"
fi
# Execute the script
TEMP_DIR=$TEMP_DIR N_CORE=$N_CORE PARMS="$PARMS" $SCRIPT $FILE1 $FILE2
# lane1-A2-B_S2_L001_1.trim
