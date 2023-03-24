# Take the raw methylation counts & convert to bed files
# Then create a script that will concat the bed files
suppressPackageStartupMessages({
  library(tidyverse)
  # library(rlang)
})
# if(!exists('argv')) argv = commandArgs(TRUE)
# out_dir = argv[1] %|% 'count_beds/indiv'
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
in_files = dir('pairwise_out/pairwise_joint_n', pattern = "bed", full.names = TRUE)
out_full = 'out/joint_with_N_all.rds'
# out_filt = 'out/joint_with_N_delta_0.5_N1.rds'
# out_filtN = 'out/joint_with_N_delta_0.5_N5.rds'
# out_filtN_bed = str_replace(out_filtN, '.rds', '.bed')
# out_filtN2 = 'out/joint_with_N_delta_0.3_N8.rds'
# out_filtN2_bed = str_replace(out_filtN2, '.rds', '.bed')

combined_data = read_tsv(in_files, id = 'file') |> select(-file)
# combined_data = read_rds(out_full)
write_rds(combined_data, out_full)
# MOST OF THIS FUNCTIONALITY HAS BEEN MOVED TO DML_PARENT_OFFSPRING_SINGLE_SETUP.R
# filt_p = 0.995
# 
# filt_delta = 0.5
# filt_N = 5
# 
# filt_delta = combined_data |> filter(prob_5 >= filt_p)
# filt_delta |> select(-start) |>   write_rds(out_filt)
# 
# # Filter for N
# filt_delta_N = filt_delta |> filter(min_N >= filt_N)
# filt_delta_N |> select(-start) |> write_rds(out_filtN)
# filt_delta_N |> write_tsv(out_filtN_bed)
# 
# # Output filt_delta_N as a bed file
# get_dml = \(delta, delta_se, threshold = 0.5) {
#   # taken from DSS::callDML
#   p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
#   p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
#   p1 + p2
# }
# 
# #alternate filter
# filt_alt = combined_data |> filter(min_N >= 7, get_dml(delta_p, delta.se_p, threshold = 0.3) >= filt_p)
# filt_alt |> select(-start) |> write_rds(out_filtN2)
# filt_alt |> write_tsv(out_filtN2_bed)
