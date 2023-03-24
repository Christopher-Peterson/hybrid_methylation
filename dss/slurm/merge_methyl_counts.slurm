#!/bin/bash

#-------------------------------------------------------
# SBATCH -J merge_methyl_counts  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p vm-small
#SBATCH -o logs/merge_methyl_counts.o
#SBATCH -e logs/merge_methyl_counts.e
#SBATCH -t 2:00:00  

# Add N reads per loci to DSS data (and filter down)

#------------------------------------------------------
# Setup the job file
export LAUNCHER_JOB_FILE=$PWD/jobs/merge_methyl_counts
export LAUNCHER_WORKDIR=$PWD


ml launcher
module unload xalt

dss_dir=pairwise_out/trimmed_joint
n_dir=count_beds/indiv
combo_n_dir=count_beds/joint
out_dir=pairwise_out/pairwise_joint_n
mkdir -p $out_dir $n_dir $combo_n_dir

# First, run the setup, converting the data in rds form to bed
RSCRIPT="r-tidy-opt"
$RSCRIPT scripts/methyl_counts_to_bed.r $n_dir

# Use bedtools to intersect parent and offspring pairs
add_job() {
  local offspring=$1
  
  # Ordering matters
  # The offspring pair we'll get by looking at the pairwise_out/pairwise_indiv file names
  local offspring_pair=( `ls pairwise_out/pairwise_indiv/${offspring}*bed | \
        grep -oP "${offspring}.[1|2]_" | \
        sed -e 's/_/.bed/g' -e "s|X|${n_dir}/X|g"` )
  local parents=( `grep $offspring offspring_parents.tsv | sed -e "s|A|${n_dir}/A|g"` )
  
  local dss_file=${dss_dir}/${offspring}_joint.bed
  local combo_file=${combo_n_dir}/${offspring}.bed
  local out_file=${out_dir}/${offspring}.bed
  
  # awk script to find min N
  local awk_min='{min=$4;for(i=9;i<=19;i=i+5){if ($i<min)min=$i;} print $1, $2, $3, min; }'
  
  # Intersect the four count files together and report the minimum count
  # per locus
  local combo_cmd="bedtools intersect -wa -wb -a ${offspring_pair[0]} -b ${offspring_pair[1]} |\
  bedtools intersect -wa -wb -a ${parents[1]}.bed -b - | \
  bedtools intersect -wa -wb -a ${parents[2]}.bed -b - | \
  awk -v OFS='\t' '${awk_min}' > ${combo_file}"
  
  
  # Define column headers for the final output
  local out_header="$(head ${dss_file} -n1) min_N"
  # Awk command to remove duplicate chrom/start/stop columns
  local awk_post='{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $15}'
  
  # Command that writes out the header, tab-separated
  local dss_header="echo '${out_header}' | tr ' ' '\\t' > ${out_file}"
  # Merge the dss data w/ the sample size file
  local dss_cmd="${dss_header} && \
  bedtools intersect -wa -wb -a ${dss_file} -b ${combo_file} | \
  awk -v OFS='\t' '${awk_post}'  >> ${out_file}"

  if [[ ! -s $out_file || -z $out_file ]]; then # No output file
    if [[ ! -s $combo_file || -z $combo_file ]]; then # no combo file either
      echo "$combo_cmd && $dss_cmd" >> $LAUNCHER_JOB_FILE
    else
      echo "$dss_cmd" >> $LAUNCHER_JOB_FILE
    fi
  fi 
}

> $LAUNCHER_JOB_FILE
for i in {1..8}; do  add_job X${i}; done
$LAUNCHER_DIR/paramrun 
# Combine & filter outputs
r-tidy-opt scripts/filter_delta_n.r
