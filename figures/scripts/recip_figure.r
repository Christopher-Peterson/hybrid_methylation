# Make the reciprical model figures

# Stick w/ noarrow as best option
source('figures/scripts/setup_figures.r')

# Read the common data in ####
source('dss/scripts/standardize_data.r')


plot_recip = \(spp, add_arrows = TRUE) {
  # browser()
  recip_data = recip_consist_data |> filter(species ==spp ) 
  offspring_names = offspring_details |> filter(offspring %in% recip_crosses[[spp]]) |> pull(pretty_offspring)
  # axis_labs = glue("<sub></sub>Δ<sub>o</sub>, {offspring_names}") |> as.character()
  axis_labs = paste(labels$delta_o, offspring_names, sep = ', ')
  
  # Used for drawing rug plot at the bottom
  scale_data =  tibble(delta_p.1 = seq(-1, 1, by = .01)) |> 
    mutate(delta_o.1 = delta_p.1, delta_o.2 = -delta_p.1, group = sign(delta_p.1)) |> 
    # Cut out the middle
    filter(abs(delta_p.1) > min(abs(recip_data$delta_p.1)))
  
  arrow_layer = ifelse(isTRUE(add_arrows), list(
    geom_segment(aes(x = delta_p.1, y = delta_p.2, xend = delta_o.1, yend = delta_o.2, group = id),
                 alpha = .05, arrow = arrow(length = unit(0.03, 'npc'), type = 'closed')) 
  ), list())
  recip_data |> mutate(id = 1:n()) |> ggplot()+
    # geom_point(aes(x = delta_p.1, y = delta_p.2), color = 'red') + 
    fixed_lines[c('v', 'h', 'anti_diag')] +
    geom_line(aes(x = delta_o.1, y = delta_o.2, color = delta_p.1, group = group), data = scale_data, size = 1) + 
    arrow_layer + 
    geom_point(aes(x = delta_o.1, y = delta_o.2, fill = delta_p.1), shape = 21, color = grey(.8)) +
    # These are from the model
    # geom_abline(slope = slope_med, linetype = 1) + 
    # geom_abline(slope = slope_ci, linetype = 2, ) + 
    
    scale_color_viridis_c(guide = 'none', option = 'viridis', begin = .3, end = 1, aesthetics = c('colour', 'fill')) + 
    coord_fixed() + 
    # geom_rug(aes(x = delta_o.1, y = delta_o.2, color = delta_p.1), 
    #          data = scale_data,  length = unit(.03, 'npc')) + 
    labs(x = axis_labs[1], y = axis_labs[2]) + 
    theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())
}

recip_with_arrows = plot_recip('Asel', TRUE) + plot_recip('Amil', TRUE) 
recip_no_arrows = plot_recip('Asel', FALSE) + plot_recip('Amil', FALSE) ; recip_no_arrows

ggsave(out_files$consistency_recip, recip_with_arrows, width = 12, height = 6, dpi = 300)
ggsave(out_files$consistency_recip_noarrow, recip_no_arrows, width = 12, height = 6, dpi = 300)
