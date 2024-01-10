#!/bin/bash

# Use bedgraph to attach gene names to everything...

window_dir=$STOCKYARD/tagmap-share/genomes/Amil_v2.1_ForBismark
db="${window_dir}/geneBoundaries.bed ${window_dir}/promoterBoundaries.bed ${window_dir}/tssBoundaries.bed"
dss_dir=$SCRATCH/hybrid_methylation/dss


mkdir -p annotated_loci
add_gene_id() {
  local in_file=$1
  local out_file=$2
  # Add them to it 
  
  # Append to the out file
  bedtools intersect -a $in_file -b $db -wa -wb >> $out_file
}
gene_id_for_group() {
  local group=$1
  local dir=`dirname $group`
  local base=`basename ${group/_chr*bed/}`
  local out_file="annotated_loci/${base}.bed"
  > $out_file
  for i in {1..14}; do
    add_gene_id "$dir/${base}_chr${i}.bed" $out_file
  done
}

for group in $dss_dir/pairwise_out/dml_bed/*chr1.bed; do
  gene_id_for_group $group & 
done

# intersect_bed() {
#   local dat=${1:-bed_filter/parents_0.5_0.05.bed}
#   local window=${2:-$window_dir/promoterBoundaries.bed}
#   local bn_dat=`basename $dat`
#   local bn_win=`basename $window`
#   local out_file="intersecting_bed/${bn_dat/.bed/}_${bn_win}"
#   
#   bedtools intersect -wa -a $dat  -b $window > $out_file
# }

