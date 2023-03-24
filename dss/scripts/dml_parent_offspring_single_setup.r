suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(rlang)
  library(glue)
})
if(!exists('argv')) argv = commandArgs(TRUE)

delta = argv[1] |> as.numeric() %|% 0.5
N = argv[2] |> as.integer() %|% 5L
job_file = argv[3] %|% paste0('jobs/dss_parent_offspring_single')
max_cores = argv[4] |> as.integer() %|% 48L
data_file = argv[5] %|% as.character(glue('out/joint_with_N_delta_{delta}_N{N}.rds'))
result_dir = argv[6] %|% as.character(glue("par_off_single_out_delta_{delta}_N{N}") )


rscript = 'r-methylkit scripts/dml_parent_offspring_single.r'
chrom_list = glue("NC_0580{65+1:14}.1") |> as.character()

# First, filter the data file
if(!file.exists(data_file)) {
  full_data = 'out/joint_with_N_all.rds'
  combined_data = read_rds(out_full)
  filt_p = 0.995
  filt_delta = delta # rename temporarily 
  # Output filt_delta_N as a bed file
  get_dml = \(delta, delta_se, threshold = filt_delta) {
    # taken from DSS::callDML
    p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
    p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
    p1 + p2
  }
  
  data_file_bed  = str_replace(data_file, '.rds', '.bed')
  filtered_data = combined_data |> filter(min_N >= N, get_dml(delta_p, delta.se_p) >= filt_p)
  filtered_data |> write_tsv(data_file_bed)
  data = filtered_data |> select(-start) |> write_rds(data_file)
} else {
  data = read_rds(data_file)
}


run_data = data |>
  distinct(`#chrom`, offspring, dam, sire) |> 
  pivot_longer(dam:sire, names_to = 'type', values_to = 'parent') |> 
  mutate(chrom_num = factor(`#chrom`, levels = chrom_list) |> as.integer()) |> 
  arrange(chrom_num, parent, offspring) |> 
  mutate(job = glue('{rscript} {parent} {offspring} {chrom_num} {max_cores} {result_dir}'),
         out_dml = glue('{result_dir}/dml_test/{parent}_{offspring}_chr{chrom_num}.rds'))  

the_batch = run_data |> filter(!file.exists(out_dml))
the_batch$job |> write_lines(job_file)
