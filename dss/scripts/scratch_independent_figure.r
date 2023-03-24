# quick & dirty new plotting method
library(tidyverse)
library(ggtext)
library(glue)
library(ggpointdensity)
library(patchwork)
in_data = dir('figures/data', pattern = 'independent_X..rds', full.names = TRUE) |> map_dfr(read_rds)

get_dml = \(delta, delta_se, threshold = filt_delta) {
  # taken from DSS::callDML
  p1 <- pnorm(delta - threshold, sd = delta_se) ## Pr(delta.mu > delta)
  p2 <- pnorm(delta + threshold, sd = delta_se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
    p1 + p2
}




### Plot Elements ####

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
  ) |> select(cross = offspring, pretty_offspring)

fixed_lines = list(
  h = geom_hline(yintercept = 0, linetype = 3, color = grey(.6)),
  v = geom_vline(xintercept = 0, linetype = 3, color = grey(.6)),
  diag = geom_abline(slope = 1, intercept = 0, linetype = 3, color = grey(.6))
)
labels = list(
  point_dens = "Density of\nLoci in Area",
  delta_p = "Differential methylation between parents (<sub></sub>Δ<sub>p</sub>)",
  delta_o = "Differential methylation between offspring alleles (<sub></sub>Δ<sub>o</sub>)",
  po_delta = "Intergenerational methylation change (<sub></sub>Δ<sub>♀</sub> and Δ<sub>♂</sub>)",
  diff_delta_long = "Intergenerational differential methylation change [Sign(<sub></sub>Δ<sub>p</sub>)(<sub></sub>Δ<sub>o</sub> - Δ<sub>p</sub>)]",
  diff_delta = "Intergenerational differential methylation change",
  
  # po_delta_expr = expression(frac(plain("Intergenerational")~plain("methylation"),plain("change")~(Delta["\u2640"]~plain('and')~ Delta["\u2642"]))),
  po_delta_expr = expression(atop(plain("Intergenerational change in"), plain("methylation (")*Delta["\u2640"]~plain('and')~ Delta["\u2642"]*plain(")")) ),
  
  po_delta_wrap = list(
    maternal = "Intergenerational change in<br>methylation (<sub></sub>Δ<sub>♀</sub>)",
    paternal = "Intergenerational change in<br>methylation (<sub></sub>Δ<sub>♂</sub>)"
  ), 
  diff_delta_wrap = "Intergenerational differential<br>methylation change (<sub></sub>Δ<sub>o</sub> - Δ<sub>p</sub>)",
  delta_p_wrap = "Differential methylation<br>between parents (<sub></sub>Δ<sub>p</sub>)",
  delta_o_wrap = "Differential methylation<br>between offspring alleles (<sub></sub>Δ<sub>o</sub>)"
)


#### FIrst Figure ####

figure_1 = local({
  cross_layer_plot = \(data, regline) {
    data |> ggplot(aes(delta_p, delta_o)) +
      fixed_lines[c('h', 'diag')] + 
      facet_wrap(~offspring_herit, nrow = 2) +
      geom_pointdensity() +
      theme_classic() +
      theme(strip.background = element_blank(), 
            plot.margin = unit(c(0,0,0,0), 'mm'), 
            legend.background =  element_blank(),
            strip.text = element_markdown(size = 12), 
            legend.title = element_text(size = 10), 
            axis.title.x = element_markdown(),
            axis.title.y = element_markdown(),
            axis.title = element_blank(),
            # legend.direction = 'horizontal',
            # legend.box = 'horizontal',
            legend.position = c(0.5, 0),
            legend.just = c(0.5, 0)  ) +
      coord_fixed() + 
      geom_line(aes(y = med, group = dpar), data = regline) + 
      geom_line(aes(y = lo, group = dpar), linetype = 2, data = regline) + 
      geom_line(aes(y = hi, group = dpar), linetype = 2, data = regline) +     
      scale_color_viridis_c('Data Density',
                            guide = guide_colorbar(
                              title.hjust = 0.5,
                              title.position = 'top',
                              direction = 'horizontal',
                            )) + labs(x = NULL, y = NULL)#+ 
    #labs(x = labels$delta_p, y = labels$delta_o)
  }
  species_layer_plot = \(data) {
    ggplot(data, aes(x,y)) +  geom_blank() +
      facet_wrap(~species_pretty) + 
      annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",linewidth=1) + # Adds an underline to species
      theme(strip.background = element_blank(), 
            strip.text = element_markdown(size = 14), 
            plot.margin = unit(c(0,0,0,0), 'mm'), 
            plot.background = element_blank(), 
            panel.background = element_blank(), 
            panel.spacing = unit(c(0,0,0,0), 'mm'), 
            axis.text = element_blank(), axis.title = element_blank(),
            axis.line = element_blank(), 
            axis.ticks = element_blank())+ 
      labs(x=NULL, y = NULL)
  }
    nest_species_cross = \(sp, cross_list) {
      # browser()
      N = length(cross_list)
      wrapped_plots =   patchwork::wrap_plots(cross_list, ncol = N, nrow = 1, byrow = TRUE)
      out = sp + inset_element(wrapped_plots, 0, 0, 1, 1, on_top = FALSE) 
      out
    }
  
  # Make a plot for each cross
  cross_layer_table = point_data |> 
    arrange(species_pretty, offspring) |> 
    group_by(species, offspring) |> 
    nest() |>  ungroup() |> 
    mutate(reg_data = regline |> 
             arrange(offspring) |> 
             group_by(offspring) |> 
             group_split() ) |> 
    mutate(cross_plot = map2(data, reg_data, cross_layer_plot)) |> 
    select(-data, -reg_data)
  
  # Make a frame for each species, match it against the associated crosses, and inset the nested crosses
  species_layer_table = point_data |> arrange(species_pretty) |> 
    distinct(species, species_pretty) |> 
    mutate(x = list(c(1,-1)), y = x) |> 
    unchop(c(x,y)) |> 
    nest(data = c(x, y, species_pretty)) |> 
    mutate(species_plot = map(data, species_layer_plot)) |> 
    select(-data) |> 
    left_join(cross_layer_table |> select(-offspring) |> chop(cross_plot), by = 'species') |> 
    mutate(nested_plot = map2(species_plot, cross_plot, nest_species_cross) )
  
  plot_design = 
    'AB
     CC'
  full_fig1 = ggplot(NULL) + 
    labs(x = labels$delta_p, y = labels$delta_o) + theme_classic() + 
    theme(line = element_blank(), axis.text = element_blank(),
          axis.title.x = element_markdown(), 
          axis.title.y = element_markdown()
    ) + inset_element(
      wrap_plots(species_layer_table$nested_plot, design =plot_design ), 
      left = 0, right = 1, top = 1, bottom = 0) 
  full_fig1
})
  
  ggsave('figures/test_fig1.png',figure_1, width = 13, height = 8, dpi = 300)
  










  tmp = ggplot(iris, aes(Sepal.Width, Sepal.Length))  +geom_point() + facet_wrap(~1)
tmp2 = ggplot(iris, aes(Sepal.Width, Sepal.Length))  +geom_point() + facet_wrap(~Species)

tmp + theme(#axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.line = element_blank(), axis.ticks = element_blank()) + 
  inset_element(tmp2+labs(x = NULL, y = NULL) , 0, 0, 1, 1)



placement_data = point_data |> distinct(species_pretty, offspring_herit, offspring) |> 
  mutate(xmin = (1:n())*3, xmax = xmin + 1L, ymin = xmin, ymax = xmax) |> 
  mutate(figure = map(split_data, herit_figure))

species_scaffold = \(x)

# NEW PLAN:
# Three Stages
# First, each panel
# Then, each species group
# THEN the whole thing
# 


figure1_scaffold = placement_data |> 
  ggplot() +
  geom_blank(aes(x = xmin, y = ymin)) + 
  geom_blank(aes(x = xmax, y = ymax)) + 
  # fixed_lines[c('h', 'diag')] + 
  facet_nested_wrap(
    vars(species_pretty, offspring_herit), nrow = 2, 
    nest_line = element_line(), scales = 'free',
    strip = strip_nested(
      text_x =  list(
        element_markdown(size = 14),
        element_blank()),
      by_layer_x = TRUE)) +
  # geom_pointdensity(aes(y = delta_o)) + 
  # coord_fixed(ratio = ) + 
  # geom_line(aes(y = med, group = dpar), data = regline) + 
  # geom_line(aes(y = lo, group = dpar), linetype = 2, data = regline) + 
  # geom_line(aes(y = hi, group = dpar), linetype = 2, data = regline) +     
  # scale_color_viridis_c(labels$point_dens, trans = 'log10') + 
  theme(#legend.position = c(1,0.5),
        #legend.justification = c(1,0),
        #legend.direction = 'horizontal',
        line = element_blank(), axis.text = element_blank(),
        strip.background = element_blank(), 
        ggh4x.facet.nestline = element_line(colour = 'black') )+ 
  inset_figures +
  labs(x = labels$delta_p, y = labels$delta_o) #+ 
  
ggsave('figures/test_fig1.png',figure1_scaffold,width = 16, height = 8, dpi = 300)


# geom_text(aes(x = delta_p, y = delta_o, label = panel_id), 
  #           fontface = 'bold', data = mixture_facet_labels, 
  #           vjust = 0, inherit.aes = FALSE)


# placement_data
# split_data


mixture_model_figure = point_data |> 
  ggplot(aes(x = delta_p)) +
  fixed_lines[c('h', 'diag')] + 
  facet_nested_wrap(
    vars(species_pretty, offspring_herit), nrow = 2, 
    nest_line = element_line(),
    strip = strip_nested(
      text_x =  list(
        element_markdown(size = 14),
        element_markdown(size = 12)),
      by_layer_x = TRUE)) +
  geom_pointdensity(aes(y = delta_o)) + 
  coord_fixed() + 
  geom_line(aes(y = med, group = dpar), data = regline) + 
  geom_line(aes(y = lo, group = dpar), linetype = 2, data = regline) + 
  geom_line(aes(y = hi, group = dpar), linetype = 2, data = regline) +     
  scale_color_viridis_c(labels$point_dens, trans = 'log10') + 
  theme(legend.position = c(1,0.5),
        legend.justification = c(1,0),
        legend.direction = 'horizontal',
        strip.background = element_blank(), 
        ggh4x.facet.nestline = element_line(colour = 'black') )+ 
  labs(x = labels$delta_p, y = labels$delta_o) + 
  geom_text(aes(x = delta_p, y = delta_o, label = panel_id), 
            fontface = 'bold', data = mixture_facet_labels, 
            vjust = 0, inherit.aes = FALSE)




herit_figure(split_data[[1]])

#### Make the plot ####
mixture_facet_labels = point_data |> 
  distinct(species_pretty, offspring_herit) |> 
  arrange(species_pretty, offspring_herit) |> 
  mutate(delta_p = -1, delta_o = 1, panel_id = LETTERS[1:8])
mixture_model_figure = point_data |> 
  ggplot(aes(x = delta_p)) +
  fixed_lines[c('h', 'diag')] + 
  facet_nested_wrap(
    vars(species_pretty, offspring_herit), nrow = 2, 
    nest_line = element_line(),
    strip = strip_nested(
      text_x =  list(
        element_markdown(size = 14),
        element_markdown(size = 12)),
      by_layer_x = TRUE)) +
  geom_pointdensity(aes(y = delta_o)) + 
  coord_fixed() + 
  geom_line(aes(y = med, group = dpar), data = regline) + 
  geom_line(aes(y = lo, group = dpar), linetype = 2, data = regline) + 
  geom_line(aes(y = hi, group = dpar), linetype = 2, data = regline) +     
  scale_color_viridis_c(labels$point_dens, trans = 'log10') + 
  theme(legend.position = c(1,0.5),
        legend.justification = c(1,0),
        legend.direction = 'horizontal',
        strip.background = element_blank(), 
        ggh4x.facet.nestline = element_line(colour = 'black') )+ 
  labs(x = labels$delta_p, y = labels$delta_o) + 
  geom_text(aes(x = delta_p, y = delta_o, label = panel_id), 
            fontface = 'bold', data = mixture_facet_labels, 
            vjust = 0, inherit.aes = FALSE)
ggsave(out_files$mixture_model, mixture_model_figure, width = 15, height = 10, dpi = 300)





herit_panels = split_data |> map(herit_figure)
# The main panels
fig1_main = patchwork::wrap_plots(herit_panels, ncol = 4, nrow = 2, byrow = TRUE)
# This is a hack to make the axis titles look like they do on faceted plots instead of being repeated
fig_1_full = ggplot(NULL) + 
  labs(x = labels$delta_p, y = labels$delta_o) + theme_classic() + 
  theme(line = element_blank(), axis.text = element_blank(),
        axis.title.x = element_markdown(), 
        axis.title.y = element_markdown()
        ) + 
  inset_element(fig1_main & labs(x = NULL, y = NULL), left = 0, right = 1, top = 1, bottom = 0) 

ggsave('figures/fig1_ind.png', fig_1_full, width = 14.5, height = 8, dpi = 300)

# Figure 2 ####

chrom_to_number = \(x) {
  if_else(str_detect(x, "NC_"), x |> as.factor() |> as.integer() |> as.character(), "Scaffolds") |> 
    fct_inorder()
}

cut_data = in_data |> 
  filter(get_dml(delta_p, delta.se_p, 0.25) > 0.995, N_min >= 8) |>
  chop(delta_p:N_min) |>
  mutate(count = map_int(delta_p, length)) |> 
  filter(count >1) |> 
  unchop(delta_p:N_min) |>
  mutate(diff_delta = delta_o - delta_p,
         diff_delta_se = sqrt(delta.se_o^2 + delta.se_p^2) )

consist_data = cut_data |> 
  select(chr, pos, diff_delta, diff_delta_se, cross) |> 
  arrange(chr, pos) |>
  mutate(merged_chr = chrom_to_number(chr)) |> 
  chop(c(diff_delta, diff_delta_se, cross)) |> 
  mutate(
    pos_order = 1:n(),
    mean = map_dbl(diff_delta, mean),
    sd = map_dbl(diff_delta, sd)) |> 
  unchop(c(diff_delta, diff_delta_se, cross)) |>
  group_by(pos_order) |> 
  mutate(
    lower  = diff_delta - 1.96*diff_delta_se,
    upper = diff_delta + 1.96*diff_delta_se,
    inconsistent = any(lower > mean | upper < mean),
    ) |>
  ungroup() |> 
  left_join(offspring_details) 

multi_fig = \(dat) {
  dat |> 
  ggplot(aes(x = pos_order)) + 
  fixed_lines['h'] + 
  facet_grid(.~merged_chr, scale = 'free_x', space = 'free', switch = 'x' ) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, group = chr), fill = alpha('black', .1), color = alpha('black', 0)) + 
  geom_linerange(aes(ymin = lower, ymax = upper, color = pretty_offspring), alpha = .3,
                 position = position_dodge(width = .25)) + 
  geom_point(aes(y = diff_delta, color = pretty_offspring), position = position_dodge(width = .25)) + 
  geom_line(aes(y = mean, group = chr), color = alpha(grey(.2), .7)) + 
  theme_classic() + 
  geom_rug(data = \(x) x |> filter(inconsistent), sides = 'b', color = 'red') + 
  # geom_raster(aes(fill = inconsistent, y = 2.2, height = 0.1, width = 1)) + 
  # scale_fill_manual(guide = 'none', values = c(alpha('black', 0), 'red')) + 
  # scale_color_brewer(palette = 'Set3', guide = 'none') + 
  scale_color_viridis_d('Cross',option = 'turbo',  begin = 0, end = 1) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.direction = 'horizontal',
        legend.position = c(0,0),
        legend.background = element_blank(), 
        legend.justification = c(-0.1,-0.1),
        strip.placement = 'outside',
        legend.text = element_markdown(),
        #axis.title.x = element_blank(),
        strip.background =  element_blank(), #element_rect(fill = grey(.9)),
        # strip.text = element_blank(),
        axis.title.y = element_markdown()) +
  # scale_y_continuous(
                     # limits = range(consist_data_all$diff_delta)
  # )+
  # scale_alpha_manual(values = c(.2, .8), guide = 'none') + 
  labs(x = "Rank position in the genome",
       y = labels$diff_delta_long) ;
}
 # consist_fig_1d =  consist_data |> multi_fig()
 consist_fig_top =  consist_data |> filter(merged_chr %in% as.character(1:11)) |>  multi_fig()
 consist_fig_bot =  consist_data |> filter(!merged_chr %in% as.character(1:11)) |>  multi_fig()
 
 
 fig2_full = ggplot() + 
   theme_classic() + 
   theme(line = element_blank(),
         axis.text = element_blank(),
         axis.title.x = element_markdown(),
         axis.title.y = element_markdown()) + 
   labs(x = "Rank position in the genome",
        y = labels$diff_delta_long,) + 
   inset_element( (consist_fig_top & labs(x = NULL, y = NULL) ) / 
                  (consist_fig_bot &  labs(x = NULL, y = NULL)) , 
                  left = 0, right = 1, top = 1, bottom = 0)  
ggsave('figures/fig2_ind.png', fig2_full, dpi = 300, width = 16, height = 8)  
# There are some weird blank spots on the scaffold ends.  figure them out and fix it.