# A Set of functions to be used for setting up experiments
# source these into bash before running experiments
# Save as alignment_experiments/functions.sh

function sbatch_job {
  # Like sbatch, but only returns the job id; good for queuing
  local job=`sbatch "${@}" | grep "Submitted batch job" | grep -Po "[0-9]+$" `
  echo $job
}

function new_alignment_experiment {
  # Args: Name, pe/se, parameters (in quotes), genome (without the '.sh')
  local NAME=$1
  local TYPE=${2:-"pe"} # Should be pe or se for paired or single ends
  local BISMARK_PARMS=${3:-"-q -N 1 --score-min L,0,-0.6"} # Parms to pass to bismark; should be quoted
  local GENOME=${4:-"Amil"}.sh
  
  # Setup the folder structure
  local WD=$(pwd)
  local BASE_DIR=$SCRATCH/hybrid_methylation
  mkdir -p $BASE_DIR/alignment_experiments/$NAME
  cd $BASE_DIR/alignment_experiments/$NAME
  mkdir slurm logs jobs temp bams reports config methyl_extract 
  mkdir bams/dedup reports/dedup bams/sorted bams/align reports/align
  # Save some basic configuration details
  echo $TYPE > config/ends.type
  echo $BISMARK_PARMS > config/bismark.parms
  
  # Copy the slurm template files
  cp $BASE_DIR/alignment_experiments/slurm_templates/* slurm
  
  # Replace some of the templated contents of the slurm files w/ their appropriate values
  local RUN_NAME="align_run_$NAME"
  local TRIAL_NAME="align_trial_$NAME"
  if [ $TYPE == "pe" ]; then 
    local RUN_NODES="4"
    local RUN_TIME='05:00:00' # Must be in single quotes to avoid escapes
  else
    local RUN_NODES="4"
    local RUN_TIME='07:00:00' # Must be in single quotes to avoid escapes
  fi  
  sed -e "s/sed_TIME/$RUN_TIME/g" \
      -e "s/sed_NODES/$RUN_NODES/g" \
      -e "s/sed_NAME/$RUN_NAME/g" \
    -i slurm/bismark_align_run.slurm slurm/bismark_align_run.slurm
  sed "s/sed_NAME/$TRIAL_NAME/g" \
    -i slurm/bismark_align_trial.slurm slurm/bismark_align_trial.slurm
  
  # symlink the important names
  ln -s $BASE_DIR/alignment_experiments/scripts scripts
  ln -s $BASE_DIR/reads reads
  
  # copy the genome(s) in 
  ln -s $BASE_DIR/genomes/$GENOME genome.sh
  
  # Done
  echo "$NAME set up with TYPE = $TYPE, parameters $BISMARK_PARMS"
  cd $WD
}

function run_alignment {
  local OPTION=${1:-"trial"} # OR run
  local NAME=$2
  local afterok=$3
  local BASE_DIR=$SCRATCH/hybrid_methylation
  local DIREC=$BASE_DIR/alignment_experiments/$NAME
  
  if [ -z $afterok ]; then
    afterok="-d afterok:${afterok}"
  fi
  if [ -z $NAME ]; then
    echo "Directory name must be suplied as second argument"
  elif [ ! -d $DIREC ]; then
    echo "$NAME is not a directory"
  else
    local WD=$(pwd)
    cd $DIREC
    sbatch slurm/bismark_align_${OPTION}.slurm $afterok
    cd $WD
  fi
}

function dedup_alignment {
  local NAME=$1
  local BASE_DIR=$SCRATCH/hybrid_methylation
  local DIREC=$BASE_DIR/alignment_experiments/$NAME
  local SLURM=bismark_dedup_run.slurm
  # If the slurm file doesn't exist, copy it from the templates
  if [ ! -f $DIREC/slurm/$SLURM  ]; then
    cp $BASE_DIR/alignment_experiments/slurm_templates/$SLURM $DIREC/slurm/$SLURM  
  fi
  
  if [ -z $NAME ]; then
    echo "Directory name must be suplied as first argument"
  elif [ ! -d $DIREC ]; then
    echo "$NAME is not a directory"
  else
    local WD=$(pwd)
    cd $DIREC
    sbatch slurm/$SLURM
    cd $WD
  fi
}
# Add a sort run version of this
function vis_tapestry {
  local NAME=$1
  local STAGE=${2:-"1"}
  local BASE_DIR=$SCRATCH/hybrid_methylation
  local DIREC=$BASE_DIR/alignment_experiments/$NAME
  local SLURM=vis_tapestry.slurm
  # If the slurm file doesn't exist, copy it from the templates
  if [ ! -f $DIREC/slurm/$SLURM  ]; then
    cp $BASE_DIR/alignment_experiments/slurm_templates/$SLURM $DIREC/slurm/$SLURM  
  fi
  
  if [ -z $NAME ]; then
    echo "Directory name must be suplied as first argument"
  elif [ ! -d $DIREC ]; then
    echo "$NAME is not a directory"
  else
    local WD=$(pwd)
    cd $DIREC
    STAGE=$STAGE sbatch slurm/$SLURM
    cd $WD
  fi
}