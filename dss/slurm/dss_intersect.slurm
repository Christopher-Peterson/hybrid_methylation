#!/bin/bash

#-------------------------------------------------------
# SBATCH -J dss_intersect  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p vm-small
#SBATCH -o logs/dss_intersect.o
#SBATCH -e logs/dss_intersect.e
#SBATCH -t 1:00:00  

#------------------------------------------------------

# Merge the chrom-level rds files for each individula & re-export as BED

ml launcher
module unload xalt

# Setup the job file
export LAUNCHER_JOB_FILE=$PWD/jobs/dss_intersect_r
export LAUNCHER_WORKDIR=$PWD

ind_dir=pairwise_out/pairwise_indiv
joint_dir=pairwise_out/pairwise_joint
mkdir -p $ind_dir $joint_dir
# Run the first intersections
# source scripts/intersect_bed_files.sh

# R Script: concats all chroms for a pair & converts to bedgraph file
RSCRIPT="r-tidy-opt"
add_r_job() {
  local p1=$1
  local p2=$2
  local out_file="$ind_dir/${p1}_${p2}_all_dml.bed"
  echo $out_file
  if [ ! -s $out_file  ]; then
    echo "$RSCRIPT scripts/dss_pairwise_to_bed.r $p1 $p2" # >> $LAUNCHER_JOB_FILE
  fi
}

# Get list of pairs, separated by "_"
pairs=`ls pairwise_out/dml_test/*chr1.rds | sed -e 's%pairwise_out/dml_test/%%' -e 's/_chr1.rds//'`
> $LAUNCHER_JOB_FILE
# Split pair into 2 args inside the loop
for p in $pairs; do add_r_job ${p/_/ }; done
# Run launcher if there are jobs in it
if [ ! -z $LAUNCHER_JOB_FILE ]; then $LAUNCHER_DIR/paramrun ; fi

# Use bedtools to intersect parent and offspring pairs
export LAUNCHER_JOB_FILE=$PWD/jobs/dss_intersect_bed
add_bedtools_job() {
  local offspring=$1
  local parents=( `grep $offspring offspring_parents.tsv` )
  # for 7 and 7, the 1 & 2 are turned around
  # I'm not sure why putting this in double quotes doesn't work, but it doesn't.
  local in_offspring=`ls ${ind_dir}/${offspring}*.bed`
  # Finding the parents is a bit of a pain, because the order may be swapped
  local parent_pat="(${parents[1]}|${parents[2]})_(${parents[1]}|${parents[2]})"
  local in_parents=`ls ${ind_dir}/*bed | grep -P $parent_pat`
  
  local out_file="${joint_dir}/${offspring}_joint.bed"
  local cmd="bedtools intersect -wa -wb -a ${in_parents} -b ${in_offspring} > ${out_file}"
  if [[ ! -s $out_file || -z $out_file ]]; then
    echo $cmd >> $LAUNCHER_JOB_FILE
  fi
}
> $LAUNCHER_JOB_FILE
for i in {1..8}; do  add_bedtools_job X${i}; done
$LAUNCHER_DIR/paramrun 

# r-tidy-opt scripts/dss_intersect_combine.r
