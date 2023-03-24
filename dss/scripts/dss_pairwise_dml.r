suppressPackageStartupMessages({
  library(tidyverse)
  library(DSS)
  library(glue)
  library(rlang)
  library(parallel)
})

if ( !exists('argv') ) argv = commandArgs(TRUE)

# Arguments should be:
parent1 = argv[1] %|% "X5.1"
parent2 = argv[2] %|% "X5.2"
chr = argv[3] |> as.integer() %|% 14L
max_cores = argv[4] |> as.integer() %|% 48L |> pmin(parallel::detectCores()/2L)
out_dir = argv[5] %|% "pairwise_out"

out_files = local({
  subdir = c('bs_seq', 'dml_test', 'dml_call', 'dmr_call')
  glue("{out_dir}/{subdir}/{parent1}_{parent2}_chr{chr}.rds") |> 
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
})
try_dml = \(bs_dat, g1, g2, core_opts =core_lst) {
  for(core in core_opts) {
    out = safe_dml_test(core, bs_dat, g1, g2) 
    if(!is.null(out$result)) break()
    print(glue::glue("{core} cores failed"))
    print(out$error)
  }
  out$result
}

in_data = glue('chrom_data/{p}.rds',p = c(parent1, parent2) )
raw_data = in_data |> map(\(x) read_rds(x)[[chr]])
bsseq_data = raw_data |>  makeBSseqData(sampleNames = c(parent1, parent2))
write_rds(bsseq_data, out_files$bs_seq)
dml_test = bsseq_data |> try_dml(parent1, parent2)
write_rds(dml_test, out_files$dml_test)
dml_call = callDML(dml_test) |> write_rds(out_files$dml_call)
dml_call = callDMR(dml_test) |> write_rds(out_files$dmr_call)


# That's all that can be done with this one