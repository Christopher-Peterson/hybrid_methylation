suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(glue)
  library(readr)
  library(purrr)
  library(tidyr)
  library(forcats)
})
# This should be sourced by other scripts
# data_file should be defined prior to sourcing
if(!(exists('data_file') && str_detect(data_file, '.bed$' ) ) ) {
  stop("`data_file` must be declared before sourcing this script (must be .bed format)")
}

if(!exists("N")) N = 8L
if(!exists('delta')) delta = 0.25
#### Baseline Data and figure 1 ####
pretty_parents = \(x) recode(x, 
                             A1 = 's<sub>1</sub>', A2 = 's<sub>2</sub>',
                             A3 = 'm<sub>1</sub>', A4 = 'm<sub>2</sub>')
offspring_details = tribble(~"offspring",  ~"maternal", ~"paternal",
                            "X1",      "A1",      "A2",
                            "X3",      "A3",      "A4",
                            "X5",      "A1",      "A3",
                            "X6",      "A1",      "A4",
                            "X7",      "A2",      "A3",
                            "X8",      "A2",      "A4",
                            "X2",      "A2",      "A1",
                            "X4",      "A4",      "A3") |> 
  mutate(across(maternal:paternal, pretty_parents),
         pretty_offspring = glue("{maternal} × {paternal}")
  ) |> select(offspring, pretty_offspring) |> 
  arrange(offspring) |> 
  mutate(pretty_offspring = fct_inorder(pretty_offspring))

# These data are not distinct()
# Figure out what the duplicate situation is and fix it

base_data = read_tsv(data_file, show_col_types = FALSE) |> mutate(
  
  species = case_when(
    offspring %in% c("X1", "X2") ~ "Asel",
    offspring %in% c("X3", "X4") ~ "Amil",
    offspring %in% paste0("X", 5:8) ~ "Hybrid" ) |> 
    factor(levels = c('Asel', 'Amil', 'Hybrid')),
  species_pretty = fct_recode(species,
                              "Pure <i>Acropora selago</i>" = "Asel",
                              "Pure <i>Acropora millepora</i>" = "Amil",
                              "Hybrid <i>Acropora selago</i> × <i>Acropora millepora</i>" = "Hybrid")
) |> 
  distinct() |> 
  left_join(offspring_details, by = 'offspring')

make_no_po = \(dat) {
  dat |> 
    filter(min_N >= N) |> 
    filter(role == 'maternal') |> 
    select(`#chrom`, start, stop, species, species_pretty, offspring, 
           delta_p, delta.se_p, delta_o, delta.se_o, diff_delta, diff_delta_se, min_N) |> 
    distinct()
}
no_po_point_data =  make_no_po(base_data)

### Model for consistency data ####


consist_data_all = no_po_point_data |> 
  select(chrom = `#chrom`, start, species, diff_delta, diff_delta_se, offspring, min_N, delta_p) |>
  group_by(`chrom`, start) |>
  mutate(n = n(), n_d = n_distinct(offspring)) |> 
  # ungroup() |> 
  filter(n_d >= 2, min_N >= N) |> 
  select(-min_N) |>
  mutate(diff_delta_signed = sign(delta_p) * diff_delta) |> 
  select(-n, -n_d, delta_p)

chrom_to_number = \(x) {
  if_else(str_detect(x, "NC_"), x |> as.factor() |> as.integer() |> as.character(), "Scaffolds") |> 
    fct_inorder()
}
consist_data_multi = consist_data_all |> 
  ungroup() |> 
  arrange(chrom, start) |> 
  mutate(merged_chrom = chrom_to_number(chrom)) |> 
  group_by(chrom, start) |> 
  mutate(mean = mean(diff_delta_signed)) |> 
  chop(-c(chrom, start)) |> 
  ungroup() |> 
  mutate(pos_id = 1:n()) |> 
  unchop(-c(chrom, start)) |> 
  group_by(pos_id) |>
  mutate(
    sd = sd(diff_delta_signed),
    lower  = diff_delta_signed - 1.96*diff_delta_se,
    upper = diff_delta_signed + 1.96*diff_delta_se, 
    exclude = (lower > mean) | (upper < mean)) |> 
  ungroup() |> 
  mutate(chrom_num = chrom |> as.factor() |> as.integer()) |> 
  left_join(offspring_details) |> 
  mutate(tick = mean < lower | mean > upper)

full_consist_model_data = consist_data_multi |> 
  select(chrom, start, offspring) |> 
  left_join(no_po_point_data |> rename(chrom = `#chrom`)) |> 
  select(chrom, start, offspring, delta_p, delta.se_p, delta_o, delta.se_o, diff_delta, diff_delta_se) |> 
  mutate(sgn = sign(delta_p), sdelta_p = sgn * delta_p, sdelta_o = sgn * delta_o)  |> 
  arrange(chrom, start) |> 
  group_by(chrom, start) |> 
  nest() |> 
  ungroup() |> 
  mutate(locus = 1:n()) |> 
  unnest(data) |> 
  arrange(locus, offspring)

### Data for bivariate model ####
get_recip_data = \(crosses) {
  consist_data_multi |>
    filter(offspring %in% crosses) |>
    
    left_join(no_po_point_data |> select(chrom = `#chrom`, start, offspring, delta_p, delta_o, delta.se_p, delta.se_o) ) |>
    select(chrom, start, species, offspring, delta_p, delta_o,  delta.se_p, delta.se_o) |>
    group_by(species, chrom, start) |>
    mutate(n = n()) |>
    ungroup() |> filter(n >= 2) |>
    select(-n) |>
    arrange(chrom, start) |>
    mutate(offspring = factor(offspring) |> as.integer()) |>
    pivot_wider(names_from = offspring, values_from = c( delta_p, delta_o,  delta.se_p, delta.se_o),
                names_sep = '.')
}
recip_crosses = list(Amil = c("X3", "X4"), Asel = c("X1", "X2") )
# in_data = read_tsv(data_file)
recip_consist_data = imap_dfr(recip_crosses, \(x, spp) get_recip_data(x) |> mutate(species = spp))
