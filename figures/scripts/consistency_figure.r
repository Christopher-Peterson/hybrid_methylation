source('figures/scripts/setup_figures.r')

# Read the common data in ####
source('dss/scripts/standardize_data.r')



### Helper Functions ####
### 

# Standardize chromosome
chrom_to_number = \(x) {
  if_else(str_detect(x, "NC_"), x |> as.factor() |> as.integer() |> as.character(), "Scaffolds") |> 
    fct_inorder()
}
# Used for post-processing stan fit results
filter_and_index = \(df, var_name, index_name = 'idx') {
  var_start = paste0('^', var_name, '\\[')
  df |> filter(str_detect(variable, var_start)) |> 
    mutate( "{index_name}" := variable |> str_remove(var_start) |> str_remove('\\]') |> as.integer()) |> 
    select(-variable)
}

# the two guide args should be uncalled functions
multi_consistency_fig = \(diff_delta_data, locus_data, alpha_lims = alpha_limits, plot_name = 'add a name') {
  # Dataset for red ticks, indicating non-heritability
  diff_tick_data = locus_data |> filter(locMean_q95 < 0 | locMean_q5 > 0) 
  
  diff_delta_data |> # filter(merged_chrom %in% chroms) |> 
    ggplot(aes(x = locus)) + 
    fixed_lines['h'] + 
    facet_grid(.~merged_chrom, scale = 'free_x', space = 'free_x', drop = TRUE) +
    
    geom_ribbon(aes(ymin = locMean_q5, ymax = locMean_q95, group = chrom), 
                fill = alpha('black', .1), color = alpha('black', 0), data = locus_data) +
    
    geom_linerange(aes(ymin = diff_delta_q5, ymax = diff_delta_q95, color = pretty_offspring),
                   position = position_dodge(width = .25)) +
    geom_line(aes(y = locMean, group = chrom), color = 'black', data = locus_data) +
    geom_point(aes(y = diff_delta, color = pretty_offspring), position = position_dodge(width = .25)) +
    theme_classic() + 
    scale_color_viridis_d('Cross',option = 'turbo',  begin = 0, end = 1, guide = 'none') + 
    
    # Shows the variability at the locus based on lightness
    geom_rug(aes(x = locus, alpha = sigma^2), data = locus_data , color = 'black',  sides = 't', 
             length = unit(.08, 'npc') ) + 
    # Dummy version of the above data to force a color bar scale
    geom_tile(aes(fill = sigma^2, y = locMean), data = locus_data |> slice(1:2),
              alpha = 0, width = 0.1, height = 0.001) + 
    
    # IF red, there's a significant deviation from heritability
    geom_rug(aes(x = locus), data = diff_tick_data, color = 'red', sides = 't', 
             length = unit( 0.04, 'npc'), inherit.aes = FALSE) +
    # Color the ticks grey for variation
    # The alpha scale colors the rug; the fill scale sets the color bar (since scale_alpha doesn't like continuous legends)
    scale_alpha_continuous(limits = alpha_lims, range = c(0, 1), guide = 'none') + 
    scale_fill_gradient('Intra-locus variance', low = 'white', high = 'black', limits = alpha_lims,
                        guide = guide_colorbar(title.position = 'top', barwidth = unit(.205, 'npc'))) +
    
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          legend.direction = 'horizontal',
          legend.justification = c(0,0),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          # legend.text = element_markdown(),
          legend.background = element_blank(),
          #axis.title.x = element_blank(),
          # strip.background = element_blank(),
          # strip.text = element_blank(),
          axis.title.y.left = element_blank()) +
    scale_y_continuous(labels$diff_delta_long,expand = expansion(mult = c(0.02, 0.1)),
                       sec.axis = dup_axis(
                         # Create secondary axis for second y lab
                         name = plot_name
                       )
                       # limits = range(consist_data_all$diff_delta)
    )+ scale_x_continuous(expand = expansion(add = 1)) 
}




#### REad in the Data ####

consist_smry = read_rds(consist_smry_file)
gbm_loci = read_tsv(gbm_file) |> distinct(chrom = `#chrom`, start) |> mutate(gbm = TRUE)


full_data = full_consist_model_data |>
  mutate(index = 1:n(), merged_chrom = chrom_to_number(chrom)) |> 
  left_join(gbm_loci) |> mutate(gbm = !is.na(gbm)) |> 
  mutate(eta_hat = sgn * (delta_o - delta_p))


diff_delta_data = consist_smry  |>
  filter_and_index('diff_delta', 'index') |> 
  select(index, diff_delta = median, diff_delta_q5 = q5, diff_delta_q95 = q95) |> 
  left_join(full_data |> select(index, merged_chrom, start, offspring, locus, 
                                eta_hat, delta_o, delta_p, sgn, gbm),
            by = 'index') |> 
  left_join(offspring_details)

locus_data = left_join(
  consist_smry |> filter_and_index('sigma_sdeltao_locus', 'locus') |>
    # select(sigma = median, sigma_q5 = q5, sigma_q95 = q95, locus) ,
    select(sigma = median, sigma_q5 = q2.5, sigma_q95 = q97.5, locus) ,
  consist_smry |> filter_and_index('locus_mean', 'locus') |> 
    # select(locMean = median, locMean_q5 = q5, locMean_q95 = q95, locus) ,
    select(locMean = median, locMean_q5 = q2.5, locMean_q95 = q97.5, locus) ,
  by = 'locus') |> 
  left_join(full_data |> select(locus, start, gbm, merged_chrom, chrom), by = 'locus') |> 
  distinct()
alpha_limits = (locus_data$sigma^2) |> range()

gbm_data = lst(diff_delta_data, locus_data) |> 
  map(\(x) filter(x, gbm) |> ungroup() |> mutate(locus = factor(locus) |> as.integer())) |> 
  c(plot_name = 'Gene Body Loci')
nogene_data = lst(diff_delta_data, locus_data) |> 
  map(\(x) filter(x, !gbm) |>  ungroup() |> mutate(locus = factor(locus) |> as.integer())) |> 
  c(plot_name = "Intergenic Loci")

#### Make the plots ####

# Because ggtext has issues w/ legends, this is rendered as its own plot and frankensteined in later
color_legend = diff_delta_data |> distinct(pretty_offspring) |> 
  arrange(pretty_offspring) |> 
  mutate(col = rep((0:3)/2, 2), row = rep(2:1, each = 4)) |> 
  ggplot(aes(x = col, y = row, label = pretty_offspring)) + 
  geom_point(aes(color = pretty_offspring)) + 
  scale_color_viridis_d('Cross',option = 'turbo',  begin = 0, end = 1, guide = 'none') + 
  geom_richtext(aes(x = col + 0.05), hjust = 0, label.color = NA, fill = NA,
                 label.padding = grid::unit(rep(0, 4), 'pt')) +
  theme_void() +
  coord_fixed(xlim = c(0,2.05), ratio = c(0.125),clip = 'off') +
  annotate('text', label = 'Cross', y = 3, x = 0, hjust = 0.18)


consist_fig_gbm = do.call(multi_consistency_fig, gbm_data) +
  theme(axis.title.y.left = element_blank()) + 
  theme(legend.position = c(0.118, 0)) + # This is the alpha legend
  inset_element(color_legend, left = 0.115, right = 0.365, top = 0.48, bottom = 0.2, clip = FALSE,on_top = TRUE )
  # c(0.38,0.0))
consist_fig_nogene =  do.call(multi_consistency_fig, nogene_data) +
  theme(axis.title.y.left = element_blank()) + 
  theme(legend.position = 'none') # c(0.38,0.0))

fig2_full = ggplot() + 
  theme_classic() + 
  theme(line = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown()
        # legend.text = element_markdown(),
        # legend.position = 'bottom'
        ) + 
  labs(x = "Rank position in the genome",
       y = labels$eta ) +
       # y = 'Intergenerational Change in Differential Methylation<br>(Negative values indicate less differential methylation in offspring)') +
  # y = labels$diff_delta_long,) + 
  inset_element( (consist_fig_gbm  & labs(x = NULL)  ) / 
                   (consist_fig_nogene &  labs(x = NULL, y = NULL)) +
                   patchwork::plot_annotation(tag_levels = 'A'), 
                 left = 0, right = 1, top = 1, bottom = 0)
ggsave(out_files$consistency_1d, fig2_full, width = 16, height = 8, dpi = 1200)
    # fig2_full

# Make a figure of the summary stats

# gbm_data$

### Locus-specific model variance estimates figure ####
parm_estimates = bind_rows(
  consist_smry |> filter_and_index(c('sd_1'))|> 
    mutate(parameter = 'Among-Locus\nvariation', lab = "A"),
  consist_smry |> filter_and_index('sigma_sdeltao') |> 
    mutate(parameter = 'Within-Locus\nvariation', lab = "B"),
) |>  mutate(group = c('Intergenic Regions', 'Gene Bodies')[idx]) |> 
  select(group, parameter, median, q2.5, q97.5, lab) |> 
  mutate(lab = if_else(group == 'Gene Bodies', '', lab))

consistency_sd_plot = parm_estimates |>
  ggplot(aes(x = median, y = group)) + 
  facet_wrap(~parameter, ncol = 1, strip.position = 'right') + 
  geom_linerange(aes(xmin = q2.5, xmax = q97.5), linewidth = 1.2, color = grey(.5)) + 
  geom_point(color = 'black', pch = '|', size = 4) + 
  labs(y = '', x = 'SD(η)') + 
  geom_text(aes(label = lab), x = 0.17, vjust = -1.8, fontface = 'bold' ) + 
  theme(strip.background = element_blank())

ggsave(out_files$consistency_parameters, consistency_sd_plot, width = 5.5, height = 3, dpi = 300)
# Now, extract the locus effects
# and plot them, colored by significance

# Make a count w/ percentages
