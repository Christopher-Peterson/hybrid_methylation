# DSS attempt
library(tidyverse)
library(readr)
library(DSS)
library(glue)


# Read data and make by-chrom  ####
adult_names = c("A1", "A2", "A3", "A4")
adult_list = dir("parents", pattern = "*cov.gz", full.names = TRUE) |> 
  set_names(adult_names)
read_bsseq_cov = \(file) {
  dat = read_tsv(file, col_names = c('chr', 'pos', 'end', 'meth_prec', 'N_meth', 'N_unmeth'), 
                 col_select = -c(end, meth_prec))
  dat |> dplyr::mutate(N = N_meth + N_unmeth) |> 
    dplyr::select(chr, pos, N, X = N_meth)
}

adult_data = map(adult_list, read_bsseq_cov)
real_chroms = adult_data[[1]]$chr |> unique() |> stringr::str_subset("^NC")

adult_by_chrom = adult_data |>
  lapply(\(x) dplyr::filter(x, chr %in% real_chroms) |> 
           dplyr::group_by(chr) |> dplyr::group_split()) |> 
  purrr::transpose()

### Convert to BSseq ####
bs_by_chrom = parallel::mclapply(adult_by_chrom, makeBSseqData, 
                   sampleNames = c('A1', 'A2', 'A3', "A4"), mc.cores = 4)
write_rds(bs_by_chrom, 'adult_bsseq_data.rds')
bs_by_chrom = read_rds('adult_bsseq_data.rds')

# Run DML tests for all 4 adults ####
mc_args = MulticoreParam(workers = 40, progressbar = FALSE)
last_iter=14L
for(i in 1:last_iter) {
  local({
    out = DMLtest(bs_by_chrom[[i]],  group1 = c("A3", "A4"), 
                         group2 = c("A1", "A2"), smoothing = TRUE,
                         BPPARAM=mc_args)
    nm = paste('chr', i, 'test.rds', sep = "_")
    write_rds(out,  nm)
  })  
}

all_dml = dir(".", "chr_[0-9]+_test.rds") |> parallel::mclapply(read_rds, mc.cores = 14)
merged_dml = all_dml |> 
  map_dfr(\(x) unclass(x) |> structure(class = 'data.frame')) |> 
  structure(class = c("data.frame", "DMLtest"))

# Next, read in the results, merge them, and use 
# call DML and Call DMR on them
called_dml = callDML(merged_dml)
write_rds(called_dml, 'adult_dml_test.rds')
called_dmr = callDMR(merged_dml)
write_rds(called_dmr, 'adult_dmr_test.rds')

called_dml = read_rds('adult_dml_test.rds')
called_dmr = read_rds('adult_dmr_test.rds')

# Try to figure out plotting #####
called_dmr_by_chrom = local({
  split = called_dmr |> group_by(chr) |> group_split()
  split |> set_names( map_chr(split, \(df) df$chr[1]))
})
pdf('dmr_test.pdf');
for(i in 2:101)  showOneDMR(called_dmr_by_chrom[[1]][i,], bs_by_chrom[[1]]);
dev.off()

### Try Pairwise comparisons ####


adult_by_chrom[[1]][c("A1", "A2")]
