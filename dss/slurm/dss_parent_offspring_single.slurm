#!/bin/bash

#-------------------------------------------------------
# SBATCH -J dss_po_single  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 4
#SBATCH -n 16
#SBATCH -p development
#SBATCH -o logs/dss_po_single.o
#SBATCH -e logs/dss_po_single.e
#SBATCH -t 2:00:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

ml launcher
module unload xalt

export LAUNCHER_JOB_FILE=$PWD/jobs/dss_po_single
export LAUNCHER_WORKDIR=$PWD

# override via commandline
delta=${delta:-0.5}
N=${N:-5}
# Second version
 
out_dir="par_off_single_delta_${delta}_N${N}"
base_data="out/joint_with_N_delta_${delta}_N${N}.rds"
main_data="out/joint_with_N_delta_${delta}_N${N}.bed" # produced by merge_methyl_counts.slurm
delta_data="${out_dir}/dml_with_po_delta_${delta}_N${N}.bed"
ordered_data="out/dss_filtered_data_delta_${delta}_N${N}.bed"
mkdir -p $out_dir
# Setup the job file
N_CORE=${N_CORE:-32}
po_dir=${out_dir}/dml_test
merge_dir=${out_dir}/dml_merge
mkdir -p $merge_dir $po_dir

 
# setup launcher

# Set up the data file and the launcher script for estimating parent-offspring differences
r-tidy-opt scripts/dml_parent_offspring_single_setup.r $delta $N $LAUNCHER_JOB_FILE $N_CORE $base_data $out_dir
# Estimate parent-offspring differences
# if [-s $LAUNCHER_JOB_FILE ];then 
$LAUNCHER_DIR/paramrun; #fi

# after that's done, merge the results:

export LAUNCHER_JOB_FILE=$PWD/jobs/dss_po_single_p2
po_groups=`ls $po_dir | sed -e 's/_chr.*rds//g' | uniq`
# Merge all files for po pair, then intersect with the main data
> $LAUNCHER_JOB_FILE
for po in $po_groups; do
  echo "r-tidy-opt scripts/dml_parent_offspring_single_merge.r  ${po_dir}/${po}*rds | \
  bedtools intersect -wa -wb -a 'stdin' -b $main_data >  ${merge_dir}/${po}.bed" >> $LAUNCHER_JOB_FILE
done
$LAUNCHER_DIR/paramrun
# Combine the po groups
merged_files=`echo ${po_groups}' ' | sed -e 's/ /.bed /g' -e "s|A|${merge_dir}/A|g"`
# Define new header for read_tsv
new_header="#chrom start stop po_delta po_delta_se parent offspring chrom2 start2 stop2 delta_p delta.se_p delta_o delta.se_o offspring2 dam sire prob_5 min_N"
# Combine header & all merged files
echo $new_header | tr ' ' '\t' > $delta_data && cat  $merged_files >> $delta_data

### Merge and order

# Maternal/Paternal information
parental_roles=$SCRATCH/hybrid_methylation/metadata/dam_sire.csv
# Which allele-specific genomes are attached to which parents
genome_order=offspring_parents.tsv
# Determine which order dss put everything in
samples_delta_order=dss_delta_order.tsv
ls pairwise_out/dml_test/*chr1.rds | sed -e 's|.*/||g' -e 's|_chr.*||g' -e 's|_|\t|g' > $samples_delta_order

RSCRIPT=r-tidy-opt

$RSCRIPT scripts/merge_and_order_data.r "$delta_data" "$ordered_data" "$parental_roles" "$genome_order" "$samples_delta_order"
