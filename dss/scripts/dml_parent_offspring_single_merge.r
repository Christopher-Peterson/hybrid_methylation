# merge par-off-single file w/ same parent/offsprin gcombo
# Write result to stdout


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
#  library(rlang)
})

if(!exists('argv')) argv = commandArgs(TRUE)

# In this case, ALL of the args should be input beds
input_files = argv

parent_offspring = basename(input_files[1]) |> str_remove("_chr.+") |> str_split('_') |> unlist()
  

input_data = map_dfr(input_files, read_rds) |> 
  select(`#chrom` = chr, start = pos, stop = pos, po_delta = diff, po_delta_se = diff.se) |>
  mutate(start = start -1L, parent = parent_offspring[1], offspring = parent_offspring[2])
  
write_tsv(input_data, stdout())


# Next:
# Pipe into bedtools
# Filter against the focal set
# Combine output with 'cat' (Add a header too)
