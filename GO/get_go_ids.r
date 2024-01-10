# Extract all go terms
library(tidyverse)
files = dir('annotated_loci', full.names = TRUE)

bed_colnames = c("chr", "start", "end", "mu1", "mu2",
       "diff", "diff.se", "pval", "fdr", "phi1", "phi2",
       "id1", "id2" , "pair" , "generation" , 'unknown',
       'chr_2', 'gene_start', 'gene_end', 'gene_name' )

gene_names = loci = files |> 
  map(\(x) read_tsv(x, col_names = bed_colnames, col_select = gene_name) |> 
        pull(gene_name) |> unique() |> str_remove(fixed('gene-'))) |> 
  unlist() |> unique()

go_key = read_tsv('db/Amil_2.1_GO.tsv', col_names = c('gene_id', 'go_terms'))

length(gene_names)
missing_names = gene_names |> discard(\(x) x %in% go_key$gene_id)

unindexed_genes = go_key |> filter(! gene_id %in% gene_names)

go_keys_to_use = go_key |> filter(gene_id %in% gene_names) |> 
  filter(nchar(go_terms) > 0) |> 
  mutate(go_terms_list = str_split(go_terms, ',')) 
write_rds(go_keys_to_use, 'db/useful_go_keys.rds')

