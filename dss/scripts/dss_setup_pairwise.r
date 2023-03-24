suppressPackageStartupMessages({
  library(tidyverse)
  library(DSS)
  library(glue)
  library(rlang)
  library(parallel)
})

if ( !exists('argv') ) argv = commandArgs(TRUE)
jobfile = argv[1] %|% "jobs/dss_pairwise"
out_dir = argv[2] %|% "pairwise_out"
force_reset = argv[3] |> as.logical() %|% FALSE
max_cores = argv[4] |> as.integer() %|% 48L

#### Create the input files ####
read_bsseq_cov = \(file) {
  dat = read_tsv(file, col_names = c('chr', 'pos','end', 'meth_prec', 'N_meth', 'N_unmeth'), 
                 col_select = -c(end, meth_prec))
  dat |> dplyr::mutate(N = N_meth + N_unmeth) |> 
    dplyr::select(chr, pos, N, X = N_meth)
}

chrom_dir = 'chrom_data'
dir.create(chrom_dir, showWarnings = FALSE)

# Files that already exist
existing_data = dir(chrom_dir, '.rds') |> str_remove(".rds")
not_existing = ifelse(force_reset, identity, \(x) x[!x %in% existing_data])
# List input files
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
input_list = c(adult_list,offspring_list) |> not_existing() 
if(length(input_list) > 0) {
  all_data = mclapply(input_list, read_bsseq_cov, mc.cores = 7) |> set_names(names(input_list))
  # List full chromosomes
  real_chroms = all_data[[1]]$chr |> unique() |> stringr::str_subset("^NC")
  
  # Create a list, where the first level is separated by individual and the second by chrom
  data_by_chrom = all_data |>
    lapply(\(x) dplyr::filter(x, chr %in% real_chroms) |> 
             dplyr::group_by(chr) |> dplyr::group_split())
  # Save each of them, one file per bam
  iwalk(data_by_chrom, \(dat, nm) {
    file = glue("{chrom_dir}/{nm}.rds")
    write_rds(dat, file)
  })
}
#### Create the Jobfile ####
offspring_pairs =
  read_tsv("offspring_parents.tsv", 
           col_names = c('offspring', 'g2', 'g1')) |> 
  mutate(across(g2:g1, \(x) str_remove(x,"A") |> as.integer()),
         in_order = g2>g1,
         first = paste0(offspring, '.', 1 + !in_order),
         second = paste0(offspring, '.', 1 + in_order),
         ) |> 
  select(first, second) |> pmap(\(...) c(...) |> as.character())

adult_pairs = combn(names(adult_list), 2, simplify = FALSE)
# offspring_pairs = map(1:8, \(i) paste0("X", i, ".", c(1,2)))
# Need to figure out directions, though...

all_pairs = c(adult_pairs, offspring_pairs) |>
  map_dfr(\(x)  tibble(p1 = x[1], p2 = x[2], chr = seq_along(real_chroms))) |> 
  mutate(end_file = glue("{out_dir}/dml_test/{p1}_{p2}_chr{chr}.rds") )
# Remove completed data sets
if(isFALSE(force_reset)) all_pairs = all_pairs |> dplyr::filter(!file.exists(end_file))

script = 'r-methylkit scripts/dss_pairwise_dml.r'
glue_data(all_pairs, "{script} {p1} {p2} {chr} {max_cores} {out_dir}") |> 
  write_lines(jobfile)

