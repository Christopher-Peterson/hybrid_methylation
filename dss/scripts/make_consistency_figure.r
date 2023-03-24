library(tidyverse)
library(withr)

data_file = "out/joint_with_N_delta_0.5_N5.rds"
in_data = read_rds(data_file)
get_dml = \(delta, delta_se, threshold = 0.5) {
  # taken from DSS::callDML
  p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
  p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
  p1 + p2
}

ch_dat = in_data |> 
  group_by(`#chrom`, start) |> 
  chop(delta_p:min_N) |> 
  ungroup() |> 
  mutate(N = map_int(prob_5, length)) |> 
  filter(N > 1) |> 
  unchop(delta_p:min_N) |> 
  # I'm including th esign of delta_p so that diff_op measures the decrease in delta
  mutate(diff_op = sign(delta_p) * (delta_o - delta_p), se_op = sqrt(delta.se_p ^2 + delta.se_o ^ 2))
  


spp_combo = function(sppa, sppb) {
  if_else(sppa == sppb, glue("{sppa} (Both)"),
          glue("{pmin(sppa, sppb)} / {pmax(sppa, sppb)}")
          )
}

make_n2_data = function(seed) {
  withr::with_seed(seed,  {
  ch_dat |> 
    filter(N == 2) |> 
    select(-c(delta_p:delta.se_o,dam:min_N)) |> 
    mutate(species = case_when(
      offspring %in% c("X1", "X2") ~ "Asel",
      offspring %in% c("X3", "X4") ~ "Amil",
      offspring %in% paste0("X", 5:8) ~ "Hybrid"
    )) |> select(-offspring) |> 
    group_by(`#chrom`, start) |> 
    mutate(rep = sample(c('a', 'b'), replace = FALSE)) |> 
    pivot_wider(values_from = c(diff_op, se_op, species),
                names_from = rep) |> 
    mutate(species = spp_combo(species_a, species_b), seed = seed) #|> 
    # select(-species_a, -species_b)
  })
}

n2_data = make_n2_data(18941)
# set.seed(1156); n2_runs = rdunif(9, b = 65536) |> 
  # map_dfr(make_n2_data)

# There's no variation caused by swapping x and y
 
ci_mult = 1 # 1.96

# library(ggforce)
  
make_comp_plot = function(data) {
  data |> 
  ggplot(aes(x = diff_op_a, y = diff_op_b, color = species)) + 
    coord_fixed() + 
    # facet_wrap(~seed) + 
    geom_vline(xintercept = 0, color = grey(.75), linetype = 3) + 
    geom_hline(yintercept = 0, color = grey(.75), linetype = 3) +
    geom_abline(slope = 1, color = grey(.75), linetype = 3) + 
    geom_linerange(aes(xmin = diff_op_a - ci_mult * se_op_a,
                       xmax = diff_op_a + ci_mult * se_op_a), alpha = .5) +
    geom_linerange(aes(ymin = diff_op_b - ci_mult * se_op_b,
                       ymax = diff_op_b + ci_mult * se_op_b), alpha = .5) +
    geom_point() +
    # geom_ellipse(aes(x0 = diff_op_a, y0 = diff_op_b,
                     # a = se_op_a, b = se_op_b, angle = 0,  alpha = .4)) + 
    scale_color_brewer(palette = "Dark2") + 
    theme_classic() + 
    # theme(plot.tit)
    # theme(strip.background = element_blank(), strip.text = element_blank()) + 
    xlab("Δ Offspring - Δ Parents, first sample") + 
    ylab("Δ Offspring - Δ Parents, second sample")  +
    labs(title = "Are Δ Parent/Offspring differences consistent within loci?",
         caption = glue("Showing {nrow(data)} differentially methylated loci that appear in two samples."),
         subtitle = "Negative values indicate reduced Δ divergence in offspring")
}
out_fig = n2_data |> # filter(species_a == species_b) |>
  make_comp_plot() 
ggsave('dss/figures/two_loci_consistency.png',out_fig, dpi = 300, width = 7, height = 7 )  



# So what should I do next?
# Jaccard? Not likely
# Look at the mu's, I think
# They can indicate gain or loss
# Filter the original DSS data by the overlap
# Do specific parent vs. offspring analysis on these
# 
# Run DSS on all of them; each offspring vs. its specific parent
# Pool results as beds

# Split the passing loci 
# This can provide gain/loss details
# 
# Try the "contingency table" style mdrad option?

