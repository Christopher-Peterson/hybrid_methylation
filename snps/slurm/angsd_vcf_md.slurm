#!/bin/bash

#-------------------------------------------------------
#SBATCH -J angsd_vcf_md  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 4
#SBATCH -n 64
#SBATCH -p normal
#SBATCH -o logs/angsd_vcf_md.o
#SBATCH -e logs/angsd_vcf_md.e
#SBATCH -t 2:00:00  

#------------------------------------------------------

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=jobs/angsd_vcf_md
source genome.sh

# Pull in the regions as an array
REGIONS=( `cat $GENOME_PATH | grep -oP ">N._[0-9]+.1" | grep -oP "N.+"` )
# write regions out
echo "${REGIONS[@]}" > genome_regions.txt
# REGIONS=( `cat genome_regions.txt` )


# Create the list of parent-combination BAMs

parent_lst=("parents" "s4" "s6" "m8" "m11") # "parents" is there to enforce 1-indexing

# First, alias the parents to the A-naming system
for i in {1..4}; do
ln -sfL ${parent_lst[$i]}.bam md_bam_merge/A${i}.bam
ln -sfL ${parent_lst[$i]}.bam.bai md_bam_merge/A${i}.bam.bai

done

# Parent_offspring lists
mkdir -p bam_lists
for i in {1..3}; do
  for (( j=$i+1; j<=4; j++ )); do
    ls md_bam_merge/A{$i,$j}*bam > "bam_lists/A$i$j"
  done
done

################### Plan:
## Make a list of paired bams for each offspring's parents
## Function: takes Bamlist name and chrom name as args
  ## Adds job to the job queue if resulting files don't exist
  ## output should be in vcfs/bamlist/...
## Iterate through all of it
## 
## Save the chromosome list elsewhere
## Use the chrom list to append each vcf group in order using bcftools cat (or similar?)

function add_job {
  local BAM_LIST=$1
  local NAME=$2
  local REGION=$3
  local THREADS=$4
  local OUTFILE="vcf/$NAME/vcf_mod_${REGION}"
  # Here's where you check the output
  
  if [ ! -s $OUTFILE.bcf ]; then
    echo "angsd -bam $BAM_LIST -r $REGION -doBcf 1 -GL 1  \
     -doPost 1  -doMajorMinor 1 -doMaf 1\
     -doGeno 1 -doCounts 1 \
    --ignore-RG 0 -geno_minDepth -1 -minInd 0 \
    -nThreads $THREADS  -out $OUTFILE" >> $LAUNCHER_JOB_FILE
  fi
}

function run_bamlist {
  local BAM_LIST=$1
  local NAME=`basename $BAM_LIST`
  # regions should be externally defined
  mkdir -p vcf/$NAME
  # Figure out cores
  for reg in "${REGIONS[@]}"; do
    # For the real chromosomes (Start w/ NC), use 8 cores; for everything else, use 1
    local HEADER=`echo $reg | grep -oP "^N."`
    if [ $HEADER = "NC" ]; then
      local CORES=8
    else 
      local CORES=1
    fi
    add_job $BAM_LIST $NAME $reg $CORES
  done; 
}

# Assemble the job file
> $LAUNCHER_JOB_FILE
for bamlist in bam_lists/*; do
  run_bamlist $bamlist
done

export LAUNCHER_WORKDIR=`pwd`
$LAUNCHER_DIR/paramrun

