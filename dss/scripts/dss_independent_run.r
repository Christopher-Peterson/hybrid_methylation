if(!exists('argv')) argv = commandArgs(TRUE)
suppressPackageStartupMessages({
  library(DSS)
  library(tidyverse)
  library(rlang)
  library(glue)
})
cross = argv[1] %|%  'X1'

out_dir = 'independent_loci/dml'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
in_file = 'independent_loci/input_data/cross_{cross}.tsv' |> glue()
out_file = glue('{out_dir}/independent_{cross}.rds')
in_data = read_tsv(in_file)

# Minimum size of N
min_n = in_data |> 
  mutate(N_min = pmin(N_mat_o, N_pat_o, N_mat_p, N_pat_p)) |> 
  select(chr = `#chrom`, pos, N_min)

# Pull in the data & re-order
crossing_data = read_tsv('crossing_metadata.tsv') |> 
  filter(cross_name == cross) |> 
  rename(cross = cross_name) |> 
  pivot_longer(-cross, names_to = 'id', values_to = 'name')
data_long = in_data |> 
  pivot_longer(-c(1, 2, cross), 
               names_to = c(".value", "lineage", "generation"), names_sep = '_', ) |> 
  mutate(id = paste(lineage, generation, sep = '_')) |> 
  left_join(crossing_data) |> 
  select(-lineage, -generation)

# Get a list of chromosomes
chrom_lists = local({
  chroms = data_long |> distinct(`#chrom`) |> unlist() |> str_subset("NC")
  # Split scaffolds into 6 groups that will all be run together
  scaffolds = data_long |> distinct(`#chrom`) |> unlist() |> str_subset("NC",   negate = TRUE) |> 
    tibble(scaf = _) |> mutate(group = rep_along(scaf, 1:6)) |> group_split(group) |> map('scaf')
  chrom_lists = c(scaffolds, chroms |> as.list() )
  chrom_lists
})

rm(offspring_parents, in_data) # clean up memory

# Quickly run DML
quick_dml = function(dat, names) {
  # browser()
  main_df = dat |> 
    select(-cross) |> 
    filter(id %in% names) |> 
    rename(real_pos = pos, real_chr = `#chrom`) |> 
    chop(c(N, X, id, name)) |> 
    mutate(chr = paste0('pseudo_chrom_', 1:n()),
           pos = 1) |> 
    unchop(c(N, X, id, name)) |> 
    relocate(chr, pos, N, X, everything()) 
  
  data_list = main_df |>  group_split(id) 
  bsseq = data_list |> makeBSseqData(sampleNames = names) |> suppressWarnings()
  dml = DMLtest(bsseq,  group1 = names[1], group2 = names[2],  equal.disp = TRUE, smoothing = FALSE, ncores = 1L)
  # browser()
  as_tibble(dml) |> left_join(main_df |> distinct(chr, pos, real_chr, real_pos), by = c("chr", "pos")) |> 
    select(chr = real_chr, pos = real_pos, delta = diff, delta.se = diff.se, mu1, mu2)
}
joint_dml = \(chroms, ...) {
  # browser()
  par = data_long |>  filter(`#chrom` %in% chroms) |> quick_dml(c('mat_p', 'pat_p')) 
  off = data_long |>  filter(`#chrom` %in% chroms) |> quick_dml(c('mat_o', 'pat_o'))
  mat = data_long |>  filter(`#chrom` %in% chroms) |> quick_dml(c('mat_o', 'mat_p'))
  pat = data_long |>  filter(`#chrom` %in% chroms) |> quick_dml(c('pat_o', 'pat_p'))
  # browser()
  # Inner joining parent and offspring, because I need them both together
  # full joining the rest because I can filter out missing parts and they don't require parity
  # Seriously consider some alternate methods, because you're missing too many parts with this
  # Maybe look at the beta model again? 
  full_join(
    inner_join(par, off, c('chr', 'pos'), suffix = c('_p', '_o')),
    full_join(mat, pat, c('chr', 'pos'), suffix = c('_mat', '_pat')),
    c('chr', 'pos'))
}
all_chrs = lapply(chrom_lists , joint_dml) |> bind_rows() |> mutate(cross = cross)
all_chrs |> left_join(min_n) |> write_rds(out_file)
