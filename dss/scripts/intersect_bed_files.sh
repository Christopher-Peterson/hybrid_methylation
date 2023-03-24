#!/bin/bash

# the goal here: Filter the bed files in bed_dir
# with the ones in bed_filter
window_dir=$STOCKYARD/tagmap-share/genomes/Amil_v2.1_ForBismark

mkdir -p intersecting_bed
intersect_bed() {
  local dat=${1:-bed_filter/parents_0.5_0.05.bed}
  local window=${2:-$window_dir/promoterBoundaries.bed}
  local bn_dat=`basename $dat`
  local bn_win=`basename $window`
  local out_file="intersecting_bed/${bn_dat/.bed/}_${bn_win}"
  
  bedtools intersect -wa -a $dat  -b $window > $out_file
}

for dat in bed_filter/paren*; do
  for wind in $window_dir/*bed; do
    intersect_bed $dat $wind
  done
done

# Some of these have illegal characters in them
# I think NA's?
# Go back and filter them
# 