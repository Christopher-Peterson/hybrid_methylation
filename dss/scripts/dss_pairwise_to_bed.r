# What to do now:

# Start w/ the pairwise dml files
# Read in all of a single individual's,
# write them out as a properly indexed bed
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(rlang)
  library(glue)
})
if(!exists('argv'))argv = commandArgs(TRUE)
parent1 = argv[1]  %|% "A1"
parent2 = argv[2]  %|% "A2"

out_dir = 'pairwise_indiv'
out_file = glue('pairwise_out/{out_dir}/{parent1}_{parent2}_all_dml.bed')
in_files = glue("pairwise_out/dml_test/{parent1}_{parent2}_chr{chr}.rds", chr = 1:14)

in_data = map_dfr(in_files, \(x) read_rds(x) |> as_tibble() )

bed_data = in_data |> 
  select(`#chrom`=chr, start = pos, end = pos, diff, diff.se) |> 
    mutate(start = start - 1L)
  
write_tsv(bed_data, out_file)