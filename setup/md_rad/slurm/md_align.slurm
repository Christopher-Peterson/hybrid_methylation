#!/bin/bash

#-------------------------------------------------------
#SBATCH -J md_align
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p normal
#SBATCH -o logs/md_align.o
#SBATCH -e logs/md_align.e
#SBATCH -t 01:00:00
#### #SBATCH -A tagmap

#------------------------------------------------------
module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`
export LAUNCHER_JOB_FILE=jobs/trim

# source ../../genomes/Amil.sh
mkdir -p sams
N_CORES=16
mkdir -p trimmed
GENOME_PATH=genome/Amil

function add_job {
  local file=$1
  local out="sams/`basename ${file/_trim.fastq/.sam}`"
  if [ ! -z $out ]; then 
    echo "bowtie2 -x $GENOME_PATH -U ${file} -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p $N_CORES -S $out" >> $LAUNCHER_JOB_FILE
    # echo "bowtie2 -x $GENOME_PATH -U ${file} --local  -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p $N_CORES -S $out" >> $LAUNCHER_JOB_FILE
  fi
}

# Populate the job file
>$LAUNCHER_JOB_FILE
for file in trimmed/*fastq; do 
  add_job $file
done

# Run the launcher
$LAUNCHER_DIR/paramrun


