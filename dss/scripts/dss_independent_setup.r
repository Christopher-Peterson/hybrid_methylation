if(!exists('argv')) argv = commandArgs(TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rlang)
  library(glue)
})
out_dir = 'independent_loci/input_data'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_count = argv[1] |> as.integer() %|% 5L
max_count = argv[2] |> as.integer() %|% 400L
parent_role_data =  argv[3] %|% (Sys.getenv('SCRATCH') |> file.path('hybrid_methylation/metadata/dam_sire.csv'))
genome_order_data = argv[4] %|% 'offspring_parents.tsv'

### Determine Sample Metadata ####

# For the final version of this analysis, we want $\Delta$ to be Maternal - Paternal.
# For some cases, this has been lost; however, the data trail exist to identify 
# which samples need to have their signs flipped. 

parent_role = read_csv(parent_role_data, show_col_types = FALSE) |> 
  pivot_longer(DAM:SIRE, names_to = 'role', values_to = 'parent') |> 
  rename(offspring = ID) |> 
  mutate(role = recode(role, DAM='maternal', SIRE='paternal'))

genome_order = read_tsv(genome_order_data, show_col_types = FALSE,
                        col_names = c('offspring', 'g2', 'g1')) |> 
  pivot_longer(g1:g2, names_to = 'genome', names_prefix = 'g', values_to = 'parent') |> 
  mutate(offspring_asm = glue("{offspring}.{genome}"))

# Combine ordering info

sample_metadata = genome_order |> 
  left_join(parent_role, by = c('parent', 'offspring')) |> 
  select(-genome) |> 
  pivot_longer(c(parent, offspring_asm), names_to = 'generation', values_to = 'id') |> 
  pivot_wider(names_from = role, values_from = id) |> 
  mutate(concat = glue("{maternal}_{paternal}"))

crossing_data = sample_metadata |> select(offspring:paternal) |> 
  mutate(generation = substr(generation, 1,1)) |> 
  rename(mat = maternal, pat = paternal) |> 
  pivot_wider(names_from = 'generation', values_from = c('mat', 'pat')) |> 
  rename(cross_name = offspring)
write_tsv(crossing_data, 'crossing_metadata.tsv')

# Prepare Input Files
adult_list = local({
  adult_names = c("A1", "A2", "A3", "A4")
  dir("parents", pattern = "*A[1-4].*cov.gz", full.names = TRUE) |> 
    set_names(adult_names)
})
offspring_list = local({
  files = dir("offspring", pattern = "X[1-8].genome[1|2].bismark.cov.gz", full.names = TRUE) 
  nms = files |> basename() |> str_remove('.bismark.cov.gz') |> str_remove("genome")
  files |> set_names(nms)
})

read_bsseq_cov = \(file) {
  dat = read_tsv(file, col_names = c('chr', 'pos','end', 'meth_prec', 'N_meth', 'N_unmeth'), 
                 col_select = -c(end, meth_prec))
  # browser()
  dat |> dplyr::mutate(N = N_meth + N_unmeth) |> 
    dplyr::select(chr, pos, N, X = N_meth)
}

# 
get_intersecting_data = \(cross_name, mat_p, pat_p, mat_o, pat_o, min_reads = min_count, max_reads = max_count) {
  p1 = read_bsseq_cov(adult_list[mat_p])
  o1 = read_bsseq_cov(offspring_list[mat_o])
  join1 = inner_join(p1, o1, by = c('chr', 'pos'), suffix = c('_mat_p', '_mat_o'))
  rm(p1, o1)
  p2 = read_bsseq_cov(adult_list[pat_p])
  o2 = read_bsseq_cov(offspring_list[pat_o])
  join2 = inner_join(p2, o2, by = c('chr', 'pos'), suffix = c('_pat_p', '_pat_o'))
  rm(p2, o2)
  # browser()
  full_join = 
    inner_join(join1, join2, by = c('chr', 'pos'))
  filtered_join = full_join |> 
    filter(if_all(c(N_mat_p, N_pat_p, N_mat_o, N_pat_o), 
                  list(\(x) x >= min_reads, \(x) x <= max_reads)))
  out_file = file.path(out_dir, paste0('cross_', cross_name, '.tsv'))
  filtered_join |> rename(`#chrom` = chr) |> 
    mutate(cross = cross_name) |> write_tsv(out_file)
  invisible()
}

pwalk(crossing_data, get_intersecting_data) 
