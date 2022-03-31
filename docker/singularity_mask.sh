# Create binary wrappers for singularity images
# This should be saved in $WORK/singularity
# Source this in the install slurm script

# Be sure to add $WORK/singularity/bin to $PATH-

function singularity_mask {
  # Args: bin_name, image file, command = bin_name
  local BIN_NAME=$1
  local SIF=$2 # Should contain whole path
  local CMD=${3:-BIN_NAME}
  local BIN=$WORK/singularity/bin/$BIN_NAME
  # Setup the base of the file
  echo "#!/bin/bash
    module load tacc-singularity 
    HDIR=/home
    IMAGE=$SIF
    CMD=$CMD" > $BIN
  # Add the final line using single quotes, to avoid escaping
  echo 'singularity exec -H $HDIR $IMAGE $CMD "$@" ' >> $BIN
  # Make it executable
  chmod +x $BIN
}

