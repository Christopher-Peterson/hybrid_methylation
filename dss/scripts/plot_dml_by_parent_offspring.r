suppressPackageStartupMessages({
  library(tidyverse)
  library(rlang)
  library(glue)
  library(graphics)
  library(purrr)
})
if(!exists('argv')) argv = commandArgs(TRUE)

in_file = argv[1] %|% "subset_diff_meth/parents_0.5_0.05_windowBoundaries_1kb.rds"

# Get relevant info from the file name
sub_names = in_file |> basename() |> 
  str_match("parents_(0.[0-9]+)_(0.[0-9]+)_(.+).rds")


offspring_info = read_tsv('offspring_parents.tsv', col_names = c('off', 'p1', 'p2')) |> 
  mutate(pair = 1:n()) |> 
  rowwise() |> 
  mutate(type = case_when(
    all(c_across(p1:p2) %in% c("A1", "A2")) ~ "A. selago F1",
    all(c_across(p1:p2) %in% c("A3", "A4")) ~ "A. milepora F1",
    TRUE ~ "Hybrid F1"
  )) |> ungroup()

joint_data = read_rds(in_file)

offspring_info |> 
  select(p1, p2, pair) |> 
  pivot_longer(p1:p2, names_to = 'which', values_to = 'parent') |> 
  select(-which) |> chop(pair)
chroms = joint_data$chr |> unique() |> sort()

# How to subset this
# out_dat # from other script; fix

po_diff_dat = read_rds("out/parent_offspring_dml_delta0.1.rds") |> 
  select(parent, chr, pos, po_diff = diff) |> mutate(chr = chroms[chr])

combined_po_dat = joint_data |> left_join(offspring_info |> select(p1, p2, pair)) |> 
  left_join(po_diff_dat |> dplyr::rename(p1 = parent, p1_diff = po_diff) ) |> 
  left_join(po_diff_dat |> dplyr::rename(p2 = parent, p2_diff = po_diff)) |> 
  mutate(po_diff_max = pmax(abs(p1_diff), abs(p2_diff), 0, na.rm = TRUE))

plot_title = glue("Differential methylation, {sub_names[4]}")
plot_subtitle = glue("Delta = {sub_names[2]}, FDR = {sub_names[3]}")

out_file = in_file |> str_replace("subset_diff_meth/parents", 'po_deltas') |> 
                                  # "figures/diff_meth") |> 
  str_replace(".rds", ".png")

out_figure = combined_po_dat |> 
  ggplot(aes(diff_p, diff_o)) +
  facet_wrap(~type) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_abline(slope = 0, intercept = 0, linetype = 1) +
  geom_point(aes(color = po_diff_max)) + coord_fixed() + theme_classic() +
  ggtitle(plot_title, plot_subtitle) + 
  scale_color_viridis_c("Largest Parent - Offspring\ndifference in parents") + 
  # scale_color_brewer(type = 'qual', palette = "Set1") +
  # scale_color_viridis_c("Offspring\nDifference SE") +
  xlab("Parental Methylation Difference") +
  ylab("F1 Methylation Difference")


ggsave(out_file, out_figure, width = 13, height = 7, dpi = 300)
print('done')
write_rds(combined_po_dat, '04_md_out.rds')

