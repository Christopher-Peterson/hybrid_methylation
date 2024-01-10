# Post-process the results
suppressPackageStartupMessages({
  library(tidyverse)
  library(rlang)
})
if(!exists('argv')) argv = commandArgs(TRUE)
out_file = argv[1] %|% 'out/dss_filtered_uncorrelated_0.25_N8.bed'
delta_limit = argv[2] |> as.numeric() %|% 0.25
N_limit = argv[3] |> as.integer() %|% 8L

get_dml = \(delta, delta_se, threshold = filt_delta) {
  # taken from DSS::callDML
  p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
  p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
  p1 + p2
}
signif_thresh = 0.995

# This is what the mixture model script does:
in_data = dir("independent_loci/dml", pattern = 'independent_X..rds', full.names = TRUE) |> map_dfr(read_rds)

metadata = read_tsv('crossing_metadata.tsv') |> 
  select(offspring = cross_name, maternal = mat_p, paternal = pat_p) |> 
  pivot_longer(-offspring, names_to = 'role', values_to = 'parent')

# Re-shape the input data to have the same format as the correlated version
out_data = in_data |>  
  filter(!is.na(delta_o), !is.na(delta_p),
         get_dml(delta_p, delta.se_p, delta_limit) >= signif_thresh) |> 
  select(-starts_with('mu')) |>
  mutate(start = pos - 1,
         diff_delta = delta_o - delta_p,
         diff_delta_se = sqrt(delta.se_o^2 + delta.se_p^2)
         ) |> 
  select(`#chrom` = chr, start, stop = pos, 
         offspring = cross, 
         min_N = N_min,
         everything()) |>
  pivot_longer(ends_with(c("_mat", "_pat")), names_to = c('.value', 'role'), 
               names_sep = '_') |> 
  rename(delta_po = delta, delta.se_po = delta.se) |> 
  mutate(role = recode(role, pat = 'paternal', mat = 'maternal')) |> 
  left_join(metadata) |> 
  relocate(1:4, parent, role, min_N) |> 
  filter(min_N >= N_limit)
  
write_tsv(out_data, out_file)

