# 

source('figures/scripts/setup_figures.r')

# Read the common data in ####
source('dss/scripts/standardize_data.r')

#### Function definitions ####
#' Format a number as a rounded percent
prec=\(x) formatC(x*100, digits = 1, format = 'f') # Format data as a precent


#' create the data structure for all of the plots
#' @param files mixed_model_files or similar, defined in setup_figures.r
#' @param point no_po_point_data or subsets thereof
#' @return a list with elements point (for the data) and regline 
read_theta_regline = \(files, point_data = no_po_point_data) {
  # browser()
  # browser()
  theta_quantiles = read_csv(files$theta, show_col_types = FALSE) |>
    # Get 95% CI
    select(offspring, med = `50%`, lo = `2.5%`, hi = `97.5%`) |> 
    mutate(across(-offspring, \(x) 1-x)) |> # Convert to theta_2, the prob. of the reg line
    # hi & low are swapped because it's 1-x
    left_join(offspring_details, by = 'offspring') |> 
    #paste0("X", pair) |> str_replace("XX", "X")) |> 
    mutate(offspring_herit = glue("{pretty_offspring} cross; Heritable Loci:<br>{prec(med)}% [{prec(hi)}%, {prec(lo)}%]")) #|> 
    # For some reason, pair is screwed up for X1; this fixes it
    # select(offspring, offspring_herit, pretty_offspring)
  main_data = point_data |> 
    left_join(theta_quantiles |>  
                select(offspring, offspring_herit, pretty_offspring), 
              by = 'offspring')
  
  # Re-format to have same form as point_data
  regline = read_csv(files$regline, show_col_types = FALSE) |> 
    filter(dpar == 'mu2') |> 
    select(offspring, delta_p, dpar, med = `50%`, lo = `2.5%`, hi = `97.5%`) |> 
    left_join(main_data |> distinct(offspring, species_pretty, offspring_herit), by = 'offspring') 
  list(
    points = main_data,
    regline = regline, 
    thetas = theta_quantiles
  )
}

# This plot is created from several nested sets of panels
.global_include_reg = TRUE
#' Create the point + regline plot for a single cross
cross_layer_plot = \(data, include_reg = .global_include_reg, .noaxis = TRUE) {
  extra_elements = list()
  if(isTRUE(include_reg)) {
    extra_elements = c(extra_elements, list(
      geom_line(aes(y = med, group = dpar), data = data$regline) ,
        geom_line(aes(y = lo, group = dpar), linetype = 2, data = data$regline), 
        geom_line(aes(y = hi, group = dpar), linetype = 2, data = data$regline) 
    ))
  }
  if(isTRUE(.noaxis)) {
    extra_elements = c(extra_elements, list(
      theme(axis.title = element_blank()),
      labs(x = NULL, y = NULL)
      ))
  } else {
    extra_elements = c(extra_elements, list(
      labs(x = labels$delta_p, y = labels$delta_o)
    ))
  }
  data$points |> ggplot(aes(delta_p, delta_o)) +
    fixed_lines[c('h', 'diag')] + 
    facet_wrap(~offspring_herit, nrow = 2) +
    # Adjust ofspring herit to remove the proportions if necessary
    # OR don't 
    geom_pointdensity() +
    theme_classic() +
    theme(strip.background = element_blank(), 
          plot.margin = unit(c(0,0,0,0), 'mm'), 
          legend.background =  element_blank(),
          strip.text = element_markdown(size = 12), 
          legend.title = element_text(size = 10), 
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          # legend.direction = 'horizontal',
          # legend.box = 'horizontal',
          legend.position = c(0.5, 0),
          legend.just = c(0.5, 0)  ) +
    coord_fixed() + 
    xlim(-1, 1) + ylim(-1, 1) + 
    extra_elements +  
    scale_color_viridis_c('Data Density',
                          guide = guide_colorbar(
                            title.hjust = 0.5,
                            title.position = 'top',
                            direction = 'horizontal',
                          )) #+ 
  #
}

#' Create a species-level frame
species_layer_plot = \(data) {
  ggplot(data$points, aes(x,y)) +  geom_blank() +
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

# Nest the crosses within species
#' @param sp species frame
#' @param cross_list list of cross_layer_plots
nest_species_cross = \(sp, cross_list) {
  # browser()
  N = length(cross_list)
  wrapped_plots =   patchwork::wrap_plots(cross_list, ncol = N, nrow = 1, byrow = TRUE)
  out = sp + inset_element(wrapped_plots, 0, 0, 1, 1, on_top = FALSE) 
  out
}

#' Make the whole plot
#' @param data_list list with names point and regline, created by read_theta_regline
make_mixture_model_figure = \(data_list) {
      # Make a plot for each cross
      # browser()
    cross_layer_table = data_list$points |> 
      arrange(species_pretty, offspring) |> 
      group_by(species, offspring, pretty_offspring) |> 
      nest() |>  ungroup() |> 
      mutate(reg_data = data_list$regline |> 
               arrange(offspring) |> 
               group_by(offspring) |> 
               group_split() ) |> 
      mutate(cross_plot = 
               map2(data, reg_data, 
                    \(x, y) list(points = x, regline = y) |> cross_layer_plot() )) |> 
      select(-data, -reg_data)
    # browser()
    # Make a frame for each species, match it against the associated crosses, and inset the nested crosses
    species_layer_table = data_list$points |> arrange(species_pretty) |> 
      distinct(species, species_pretty) |> 
      mutate(x = list(c(1,-1)), y = x) |> 
      unchop(c(x,y)) |> 
      nest(data = c(x, y, species_pretty)) |> 
      mutate(species_plot = map(data, \(x) list(points = x) |> species_layer_plot())) |> 
      select(-data) |> 
      left_join(cross_layer_table |> select(-offspring, -pretty_offspring) |>
                  chop(cross_plot), by = 'species') |> 
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
}

#### Read the data ####

# Figure this out...

gbm_data = read_tsv(gbm_file) |> select(1:2) |> distinct() |> left_join(no_po_point_data) 
nogene_data = read_tsv(nogene_file) |> select(1:2) |> distinct() |>  left_join(no_po_point_data);
# make a function that does makes no_po_point data for the othertwo

all_loci = read_theta_regline(mix_model_data_all, no_po_point_data)
gbm_loci = read_theta_regline(mix_model_data_gbm, gbm_data)
nogene_loci = read_theta_regline(mix_model_data_nogene, nogene_data)
# Ugh, fix this

# Do this for the GBM and no GBM as well

#### Make the figures ####

fig1_all_loci = make_mixture_model_figure(all_loci)
fig1_plots = list( all_loci, gbm_loci, nogene_loci ) |> 
  map(make_mixture_model_figure)
map2(out_files[c('mixture_model', 'mixture_model_gbm', 'mixture_model_nogene')], fig1_plots, 
     ggsave, width = 13, height = 8, dpi = 300)

# ggsave(out_files[['mixture_model']], fig1_plots[[1]], 
#      ggsave, width = 13, height = 8, dpi = 1200)


lumped_data = all_loci
lumped_data$points = lumped_data$points |> mutate(offspring_herit = '')

lumped_plot = cross_layer_plot(lumped_data, include_reg = FALSE, .noaxis = FALSE) 
ggsave('writing_presentations/fig_1_lumped.png', lumped_plot, width = 6, height = 6, dpi = 300)
# ggsave(out_files$mixture_model, figure_1, width = 13, height = 8, dpi = 300)
no_model_data = all_loci
no_model_data$points = no_model_data$points |> mutate(offspring_herit = str_remove(offspring_herit, ';.+$'))
.global_include_reg = FALSE
fig1_all_loci_no_model =  make_mixture_model_figure(no_model_data) ; .global_include_reg = TRUE
ggsave('writing_presentations/fig_1_no_model.png', fig1_all_loci_no_model,
       width = 13, height = 8, dpi = 300)

# all_loci$points


make_table_ci = \(dat) {
  dat |> 
    select(offspring, med, lo, hi) |>
    left_join(offspring_details) |> 
    arrange(offspring) |> 
    mutate(across(c(med, lo, hi), \(x) formatC(x, digits = 3, format = 'f'))) |> 
    mutate(.col = glue("{med} [{lo}, {hi}]")) |> 
    select(pretty_offspring, .col) 
}
make_table_rows = \(lst, nm) {
  # browser()
  reg_line = lst$regline |> 
    # To get the slope, subtract regline values at 1 from values at 0
    select(offspring, delta_p, med, lo, hi) |> 
    filter(delta_p %in% c(1, 0)) |> 
    arrange(offspring, delta_p) |> 
    group_by(offspring) |> 
    summarize(across(c(med, lo, hi), \(x) x[2] - x[1])) |> 
    make_table_ci() |> rename(regression = .col) 
  
  
  theta_details = lst$thetas |> 
    make_table_ci() |> rename(mixture_component = .col)
  sample_size = lst$points |> count(offspring, name = "N") |> 
    left_join(offspring_details) |> select(-offspring)
  # browser()
  left_join(reg_line, theta_details) |> left_join(sample_size) |> 
    mutate(Loci = nm)  |> 
    select(Cross = pretty_offspring, Loci, N, 'Heritable Proportion' = mixture_component, 'Regression Slope' = regression )
}
table_1 = bind_rows(
  make_table_rows(all_loci, 'All'),
  make_table_rows(gbm_loci, 'Coding'),
  make_table_rows(nogene_loci, 'Non-Coding')
          )
write_rds(table_1, out_files$mixture_model_table)
