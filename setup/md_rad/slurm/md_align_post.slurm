#!/bin/bash

#-------------------------------------------------------
#SBATCH -J md_align_post
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p normal
#SBATCH -o logs/md_align_post.o
#SBATCH -e logs/md_align_post.e
#SBATCH -t 01:00:00
#### #SBATCH -A tagmap

#------------------------------------------------------
module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`
export LAUNCHER_JOB_FILE=jobs/trim_post

# source ../../genomes/Amil.sh
mkdir -p bams
N_CORES=16
GENOME_PATH=genome/Amil

function add_job {
  local file=$1
  local out="bams/`basename ${file/.sam/.bam}`"
  if [ ! -z $out ]; then 
    local args="-u -@${N_CORES}"
    # Sort, deduplicate, and index
    echo "samtools collate -O $args $file | \
    samtools fixmate -m $args - - | \
    samtools sort $args - | \
    samtools markdup -@${N_CORES} - $out && \
    samtools index $out" >> $LAUNCHER_JOB_FILE
  fi
}

# Populate the job file
>$LAUNCHER_JOB_FILE
for file in sams/*sam; do 
  add_job $file
done

# Run the launcher
$LAUNCHER_DIR/paramrun


