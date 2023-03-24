suppressPackageStartupMessages({
  library(tidyverse)
  library(DSS)
  library(glue)
  library(rlang)
  library(parallel)
})

# What do I want to do here?
# 1. Compare all offspring directly to their parent
# 2. Compair offspring to each other (hybrid / non)

if ( !exists('argv') ) argv = commandArgs(TRUE)

# Arguments should be:
parent_id = argv[1] %|% "A1"
chr = argv[2] |> as.integer() %|% 14L
max_cores = argv[3] |> as.integer() %|% 48L |> pmin(parallel::detectCores()/2L)
out_dir = argv[4] %|% "par_off_out"

out_files = local({
  subdir = c('bs_seq', 'dml_test')
  glue("{out_dir}/{subdir}/{parent_id}_chr{chr}.rds") |> 
    as.list() |> set_names(subdir)
})
# create directories
out_files |> unlist() |> dirname() |> walk(dir.create,recursive = TRUE, showWarnings = FALSE)
core_lst = c(60, 56, 52, 48, 44, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 19:1) |> discard(\(x) x > max_cores)



### Functions ####
# Try running DML with multicore options, reducing the number of cores required until it works

safe_dml_test = safely(\(cores, bs_dat, g1, g2) {
  mc_args = MulticoreParam(workers = cores, progressbar = FALSE)
  out_tst = DMLtest(bs_dat,  group1 = g1,  group2 = g2, 
                    smoothing = TRUE, BPPARAM=mc_args)
  out_tst
})
try_dml = \(bs_dat, g1, g2, core_opts =core_lst, no_output = FALSE) {
  if(isTRUE(no_output))  sink('/dev/null')
  for(core in core_opts) {
    out = safe_dml_test(core, bs_dat, g1, g2) 
    # sink()
    if(!is.null(out$result)) break()
    # print(glue::glue("{core} cores failed"))
    # print(out$error)
  }
  if(isTRUE(no_output)) sink()
  out$result
}
#### Run the data ####

# Get the offspring associated w/ the parent
parental_group = read_tsv("offspring_parents.tsv", 
                          col_names = c('offspring', 'g2', 'g1'), 
                          show_col_types = FALSE) |> 
  mutate(hybrid = offspring %in% paste0("X", 5:8)) |> 
  pivot_longer(c(g1, g2), names_to = 'genome', values_to = 'parent') |> 
  mutate(id = glue("{offspring}.{substr(genome, 2,2)}")) |> 
  arrange(parent) |> filter(parent == parent_id)
data_names = c(parent_id, parental_group$id)
in_data = glue('chrom_data/{data_names}.rds' )

raw_data = in_data |> map(\(x) read_rds(x)[[chr]])

bsseq_data = raw_data |>  makeBSseqData(sampleNames = data_names)

write_rds(bsseq_data, out_files$bs_seq)

dml_test = bsseq_data |> try_dml(parent_id, parental_group$id, no_output = TRUE)

write_rds(dml_test, out_files$dml_test)
# dml_call = callDML(dml_test) |> write_rds(out_files$dml_call)
# dml_call = callDMR(dml_test) |> write_rds(out_files$dmr_call)


# That's all that can be done with this one