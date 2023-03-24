suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})


in_files = dir("pairwise_out/dml_test", full.names = TRUE)

out_dir = 'pairwise_out/dml_bed'
dir.create(out_dir, showWarnings = FALSE)

# tst = read_rds(in_files[1])
# tst |> as_tibble() |> head()

offspring_info = read_tsv('offspring_parents.tsv', col_names = c('off', 'p1', 'p2')) |> 
  mutate(pair = 1:n()) |> 
  rowwise() |> 
  mutate(type = case_when(
    all(c_across(p1:p2) %in% c("A1", "A2")) ~ "A. selago F1",
    all(c_across(p1:p2) %in% c("A3", "A4")) ~ "A. milepora F1",
    TRUE ~ "Hybrid F1"
  )) |> ungroup()

file_info = local({
  dat = offspring_info |> 
    mutate(off1 = paste0(off, ".1"), off2 = paste0(off, ".2")) |> 
    select(parent_a = p1, parent_b = p2, offspring_a = off1, offspring_b = off2, pair) |> 
    pivot_longer(-pair, names_to = c('generation', '.value'), names_sep = '_')
    
  bind_rows(dat |> rename(id1 = a, id2 = b),
            dat |> rename(id1 = b, id2 = a))
})
  # mutate(offspring = glue("{off}.1_{off}.2"),
  #        parents = glue("{p2}_{p1}")) |> 
  # select(parents, offspring, pair) |> 
  # pivot_longer(-pair, names_to = 'generation', values_to = 'prefix')


# Ah crap, got to sort...
# damn
# needs simplifying
# 
# X7 & 8 are backwards too ...
# Figure it out

run_data = tibble(in_file = in_files) |> 
  mutate(out_file = file.path(out_dir, str_replace(basename(in_file), 'rds', 'bed')),
         ids = in_files |> basename() |> str_remove("_chr.+rds")
         ) |> 
  separate(ids, into = c('id1', 'id2'), sep = '_') |> 
  left_join(file_info, by = c('id1', 'id2')) |> 
  filter(!is.na(pair))

convert_to_bed = \(in_file, out_file, id1, id2, pair, generation) {
  dat = read_rds(in_file) |> as_tibble() |> 
    select(`#chrom` = chr, start = pos, end = pos, mu1, mu2,
           diff, diff.se, pval, fdr, phi1, phi2) |> 
    mutate(id1 = id1, id2 = id2, pair = pair, generation = generation, start = start - 1)
  write_tsv(dat, out_file)
  
}

pwalk(run_data, convert_to_bed)
