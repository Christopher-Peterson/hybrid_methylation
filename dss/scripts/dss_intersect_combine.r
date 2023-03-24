suppressPackageStartupMessages({
  library(tidyverse)
})
cols = c('chrom', 'start', 'pos', 'delta_p', 'delta.se_p', 
         'ch2', 'st2', 'pos2', 'delta_o', 'delta.se_o')
parent_info = read_tsv('offspring_parents.tsv', col_names = c('offspring', 'dam', 'sire'))

get_dml = \(delta, delta_se, threshold = 0.5) {
  # taken from DSS::callDML
  p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
  p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
  p1 + p2
}


base_data = dir('pairwise_out/pairwise_joint', ".bed", full.names = TRUE) |>
  read_tsv(col_names = cols, col_select = c(chrom, ch2, pos, pos2, delta_p, delta.se_p, delta_o, delta.se_o),
           id = 'file') |> 
  mutate(offspring = basename(file) |> substr(1,2) ) |> 
  select(-file) |> 
  left_join(parent_info) |> 
  mutate(prob_5 = get_dml(delta_p, delta.se_p, 0.5))

write_rds(base_data, 'out/paired_loci_all.rds')

# I want details on how many samples were in each of these...
# I don't believe it as of now.
# Too many delta_0's w/ the exact same mean & sd
# Let's append the raw data there

# A task for tomorrow

write_rds(base_data |> filter(prob_5 > 0.9999), 'out/paired_loci_delta_0.5..rds')

# re-export each as cleaned up bed files
# join them w/ the methyl extraction data later

dir.create('pairwise_out/trimmed_joint', showWarnings = FALSE)
trim_bed = \(in_file) {
  local_offspring = basename(in_file) |> substr(1,2)
  
  in_data = read_tsv(in_file, col_names = cols, 
                       col_select = c(chrom, start, pos, delta_p, delta.se_p, delta_o, delta.se_o)) |> 
    mutate(offspring = local_offspring) |> 
    left_join(parent_info) |> 
    rename(`#chrom` = chrom) |> 
    mutate(prob_5 = get_dml(delta_p, delta.se_p, 0.5))
  out_file = str_replace(in_file, 'pairwise_joint', 'trimmed_joint')
  # parents = filter(parent_info, offspring == local_offspring)
  write_tsv(in_data, out_file)
  invisible()
}
dir('pairwise_out/pairwise_joint', ".bed", full.names = TRUE) |> walk(trim_bed)

# ls