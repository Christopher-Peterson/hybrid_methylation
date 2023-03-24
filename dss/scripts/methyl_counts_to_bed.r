# Take the raw methylation counts & convert to bed files
# Then create a script that will concat the bed files
suppressPackageStartupMessages({
  library(readr)
  library(rlang)
  library(stringr)
  library(dplyr)
})
if(!exists('argv')) argv = commandArgs(TRUE)
out_dir = argv[1] %|% 'count_beds/indiv'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
chrom_files = dir('chrom_data', full.names = TRUE)

# Convert count data to bed
rds_to_bed = \(in_file) {
  id = basename(in_file) |> str_remove('.rds')
  dat = read_rds(in_file) |> bind_rows() |> 
    mutate(start = pos - 1, id = id) |> 
    select(`#chrom` = chr, start, stop = pos, N, id)
  out_file = paste0(id, '.bed') |> file.path(out_dir , . = _)
  write_tsv(dat, out_file)
  invisible()
}
lapply(chrom_files, rds_to_bed)

# trim the existing bedfiles 
dir.create('pairwise_out/trimmed_joint', showWarnings = FALSE)

cols = c('chrom', 'start', 'pos', 'delta_p', 'delta.se_p', 
         'ch2', 'st2', 'pos2', 'delta_o', 'delta.se_o')
parent_info = read_tsv('offspring_parents.tsv', col_names = c('offspring', 'dam', 'sire'))

get_dml = \(delta, delta_se, threshold = 0.5) {
  # taken from DSS::callDML
  p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
  p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
  p1 + p2
}

trim_bed = \(in_file) {
  local_offspring = basename(in_file) |> substr(1,2)
  
  in_data = read_tsv(in_file, col_names = cols, 
                     col_select = c(chrom, start, pos, delta_p, 
                                    delta.se_p, delta_o, delta.se_o)) |> 
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
