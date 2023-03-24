library(tidyverse); library(ggtext)
theme_set(theme_classic())
full_dat = read_rds('dss/out/joint_with_N_all.rds')
get_dml = \(delta, delta_se, threshold = 0.5) {
  # taken from DSS::callDML
  p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
  p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
  p1 + p2
}

parent_ids = \(x) recode(x, A1 = "s<sub>1</sub>", A2 = "s<sub>2</sub>", A3 = "m<sub>1</sub>", A4 = "m<sub>2</sub>")

id_key = read_csv('metadata/dam_sire.csv') |> 
  mutate(Name = paste(parent_ids(DAM), parent_ids(SIRE), sep = " × "),
         generation = 'Offspring') |> 
  select(ID, Name, generation) |> 
  mutate(species = case_when(
    ID %in% c('X1', "X2") ~ "*A. sel*",
    ID %in% c("X3", "X4") ~ "*A. mil*",
    TRUE ~ "Hybrid"
  )) |> 
  mutate(name2 = glue::glue("{species}: {Name}"))

prob_threshold = 0.995
delta_levels = seq(0.05, 0.5, by = 0.05)

dat_by_delta = accumulate(delta_levels, \(dat, thresh) {
  dat |> filter(get_dml(delta_p, delta.se_p, threshold = thresh) >= prob_threshold) |> 
    mutate(delta_min = thresh)
}, .init = full_dat |> mutate(delta_min = 0)) |> 
  bind_rows()

dat_by_delta_n = accumulate(2:10, \(dat, n_min) {
  dat |> filter(min_N >= n_min) |> mutate(by_N = n_min)
}, .init = dat_by_delta |>  mutate(by_N = 1L)) |> bind_rows()
  
delta_n_counts = dat_by_delta_n |> group_by(offspring, by_N, delta_min) |> 
  count() |> left_join( id_key |> rename(offspring = ID) |> select(-generation)) 

# I'm considering the 0.3, n = 7 option because it seems to balance out the dominance of X7
delta_counts_figure = delta_n_counts |> 
  mutate(focal = case_when(
    delta_min == 0.5 & by_N == 5 ~ 'current',
    delta_min == 0.3 & by_N == 8 ~ 'possible',
    TRUE ~ NA_character_)) |> 
  ggplot(aes(x = by_N, y = delta_min, fill = n, color = focal)) + 
  geom_tile() + scale_fill_viridis_c('Num. Loci', trans = 'log10', breaks = c(10, 100, 1000, 10000, 100000, 1000000),
                                     labels = scales::label_number(scale_cut = scales::cut_si(''))) + 
  facet_wrap(~name2, nrow = 2) + coord_fixed(ratio = 15) + 
  theme(strip.text = element_markdown()) + 
  scale_color_manual(values = c('red', 'black'), guide = 'none', na.translate = FALSE ) + 
    
  # scale_color_manual(values = c(alpha('black', 0), 'red'), guide = 'none') + 
  scale_x_continuous(breaks = c(2,4,6,8,10)) + 
  labs(x = "Minimum Reads", y = "Minimum Parental Difference")
ggsave('figures/delta_by_N_diagnostic.png', delta_counts_figure, width = 10, height = 5, dpi = 300)


max_n = delta_n_counts |> group_by(by_N, delta_min) |> 
  filter(n == max(n)) |>  # This is the weird one 
  ungroup() |> 
  select(by_N, delta_min, n_ref = n)
delta_counts_figure_rel = delta_n_counts |> left_join(max_n) |> 
  mutate(focal = case_when(
    delta_min == 0.5 & by_N == 5 ~ 'current',
    delta_min == 0.3 & by_N == 8 ~ 'possible',
    TRUE ~ NA_character_)) |> 
  mutate(rel_n = n/n_ref) |> 
  ggplot(aes(x = by_N, y = delta_min, fill = rel_n, color = focal)) + 
  geom_tile() + scale_fill_viridis_c('Num. Loci relative \n to maximum') + 
  facet_wrap(~name2, nrow = 2) + coord_fixed(ratio = 15) + 
  theme(strip.text = element_markdown()) + 
  scale_color_manual(values = c('red', 'black'), guide = 'none', na.translate = FALSE ) + 

  # scale_color_manual(values = c(alpha('black', 0), 'red'), guide = 'none') + 
  scale_x_continuous(breaks = c(2,4,6,8,10)) + 
  labs(x = "Minimum Reads", y = "Minimum Parental Difference")
ggsave('figures/delta_by_N_diagnostic_rel.png', delta_counts_figure_rel, width = 10, height = 5, dpi = 300)

