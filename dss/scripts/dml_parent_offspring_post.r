suppressPackageStartupMessages({
library(tidyverse)
library(DSS)
library(glue)
})

input_data = dir('parent_offspring/dml_test', pattern = 'rds', full.names = TRUE) |> 
  map_dfr(\(x) read_rds(x) |> as_tibble() |>  mutate(file = basename(x))) |> 
  extract(file, into = c("parent", "chr"), regex = "(A.)_chr([0-9]+).rds", convert = TRUE)

as_dml = function(x) structure(x, class = c("data.frame", "DMLtest"))
to_bed = \(x) x |>   select(chr, start = pos, end = pos, type, pair, chr, diff, diff.se)


dmr_delta = 0.1
fdr_thresh = 0.05
dml_calls = input_data |> as_dml() |> 
  callDML(delta = dmr_delta) |> 
  as_tibble() |> 
  filter(fdr <= fdr_thresh)
# dml_calls only includes called dmls that meet the criteria; I want all of them
dml_indices = dml_calls |> 
  select(pos, parent, chr) |> 
  mutate(fdr_pass = paste("FDR <=", fdr_thresh))
out_dat = input_data |> left_join(dml_indices) |> 
  mutate(fdr_pass = coalesce(fdr_pass, "Above FDR"))


# dml_calls |> dplyr::slice(1) |> select(pos, chr, parent, fdr, pval) |>
  # left_join(input_data, by = c('pos', 'parent', 'chr')) |> glimpse()
# 

out_figure = out_dat |> 
  mutate(species = recode(parent, A1 = "A. sel", A2 = "A. sel", A3 = "A. mil", A4 = "A. mil") ,
         rep = recode(parent, A1 = 1, A3 = 1, A2 = 2, A4 = 2) ) |> 
  #mutate(fdr_clear = if_else(fdr <= fdr_thresh, paste('FDR <=', fdr_thresh), 'Above FDR')) |> 
  ggplot(aes(x = diff, y = -log(pval), color = fdr_pass)) +
  geom_point(alpha = .2) + 
  facet_grid(rep ~ species) + 
  scale_color_manual(values = c('black', 'red')) + 
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

ggsave(glue('figures/parent_offspring_fdr{fdr_thresh}_delta{dmr_delta}.png'), out_figure, width = 12, height = 12, dpi = 300)
write_rds(out_dat, glue('out/parent_offspring_dml_delta{dmr_delta}.rds'))


out_dat |> group_by(parent, fdr_pass) |> summarize(n = n()) |> mutate(total_N = sum(n)) |> 
  mutate(prop = n / total_N)

# Roughly 1% of loci were differentially methylated, parent vs. offspring

# Subset these datasets by which ones were significantly different in parents...



# Read them in
# DML test




# Write as bed
# then launch the bedtools scripts

# Followed by plotting...
# Volcano plots?
# 

# I need to think about this
# Look for DML's?
# 
