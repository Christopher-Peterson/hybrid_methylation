# Make expectation graphs

source('figures/scripts/setup_figures.r')
library(tidyverse)
theme_set(theme_classic())

er_data = local({
  set.seed(23548907)
  n_per = 150
  no_herit = tibble(
    parent = c(runif(n_per, -1, -.6), runif(n_per, .6, 1)),
    offspring = rnorm(n_per*2, 0, .25) + runif(n_per*2, -.1, .1)/2,
    herit = 'Methylation is Not Heritable'
  )
  n_per = 75
  
  herit = tibble(
    parent = c(runif(n_per, -1, -.6), runif(n_per, .6, 1)),
    herit = 'Methylation is Heritable'
  ) |> mutate(
    offspring = round(parent) * .8 + rnorm(n_per * 2, 0, .1)
  )
  bind_rows(herit, no_herit) |> 
    filter(offspring |> between(-1, 1))   
})
make_er_figure = \(data, nrow = 1) {
  ggplot(data, aes(x = parent, y = offspring)) + 
    geom_abline(slope = 1, linetype = 2) + 
    geom_hline(yintercept = 0, linetype = 1) + 
    coord_fixed() + 
    facet_wrap(~herit, nrow = nrow) + 
    geom_point(shape = 21, fill = grey(.6)) + 
    labs(x = labels$delta_p, y = labels$delta_o) + 
    theme(strip.background = element_blank(),
          axis.title.x = element_markdown(size = 12), axis.title.y = element_markdown(size = 12),
          # axis.text = element_text(size = 14),
          strip.text = element_text(size = 12)
    )
}
expected_results_fig = er_data |> make_er_figure(1)

ggsave(out_files$expected_results, expected_results_fig, width = 10, height = 5.5, dpi = 300)  

expected_results_fig_v = er_data |> make_er_figure(2)
ggsave('writing_presentations/expected_v.png', expected_results_fig_v, width = 4, height = 8, dpi = 300)  
