# Define eta
library(tidyverse)
library(glue)
library(withr)
library(patchwork)
theme_set(theme_classic())
source('figures/scripts/setup_figures.r')
example_data = with_seed(12415, 
  tibble(delta_p = rnorm(100, 0, 0.1) + rep(c(0.7, -0.7), each = 50)) |> 
  mutate(delta_o = rep(c(0.65, -0.65), each = 50) * rep(c(.95,0, -.85), times = c(45, 45, 10)/5)  +
           rnorm(100, 0, 0.15) ,
          diag = delta_o - delta_p, 
          sgn = (sign(diag) * sign(delta_p)) |> as.factor(),
          eta = sign(delta_p) * diag ))
md_axes =   theme(axis.title.x = element_markdown(), axis.title.y = element_markdown())
p0 = example_data |> 
  ggplot(aes(delta_p, delta_o)) + coord_fixed(xlim = c(-1,1), ylim = c(-1, 1)) +
  fixed_lines[c('h', 'diag')] + 
  labs(x = labels$delta_p, y = labels$delta_o) + md_axes
p1 = p0 +  geom_point()  
p2 = p0 + geom_segment(aes(xend = delta_p, yend = delta_p, color = sgn), alpha = .5) + 
          scale_color_viridis_d('',option = 'turbo', begin = .2, end = .8, direction = -1, guide = 'none',
                                labels = c('Differential methylation reduced in offspring',
                                           'Differential methylation increased in offspring')) + 
  geom_point() #+ 
  # theme(legend.position = c(0, 1), legend.justification = c(0,0.55), legend.background = element_blank())
# p2
p3 = example_data |> 
  ggplot(aes(delta_p, diag)) +
  coord_fixed(xlim = c(-1,1), ylim = 1.4*c(-1, 1)) +
  fixed_lines$h + 
  labs(x = labels$delta_p, y = "<sub></sub>Δ<sub>o</sub> - <sub></sub>Δ<sub>p</sub>") +
  md_axes + 
  geom_segment(aes(xend = delta_p, yend = 0, color = sgn), alpha = .5) + 
  scale_color_viridis_d('',option = 'turbo', begin = .2, end = .8, direction = -1, guide = 'none',
                        labels = c('Differential methylation reduced in offspring',
                                   'Differential methylation increased in offspring')) + 
  geom_point() + md_axes
# p3

p4 = example_data |> 
  ggplot(aes(delta_p, eta)) +
  coord_fixed(xlim = c(-1,1), ylim = 1.4*c(-1, 1)) + 
  fixed_lines$h + 
  labs(x = labels$delta_p, y = labels$eta) +
  geom_segment(aes(xend = delta_p, yend = 0, color = sgn), alpha = .5) + 
  scale_color_viridis_d('',option = 'turbo', begin = .2, end = .8, direction = -1, 
                        labels = c('Differential methylation reduced in offspring',
                                   'Differential methylation increased in offspring')) + 
  md_axes +   geom_point() +
  theme(legend.position = c(0, 1), legend.justification = c(0,0.55), 
        legend.background = element_blank())
# p4
joint_plot = p2+p3+p4 + patchwork::plot_layout( nrow = 1, widths = c(1,1,1)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.background = element_blank()) 

ggsave(out_files$eta_diagram, joint_plot, width = 12, height = 5.5, dpi = 300)
    # tribble(
#   ~plt, ~nm, ~height,
#   p1, 'p1', 5,
#   p2,  'p2', 5,
#   p3 , 'p3', 7.5,
#   p4, , 'p4',7.5,
#   ) |> pwalk(\(plt, nm, height) ggsave(glue('writing_presentations/eta_def_{nm}.png'), plt, width = 5, dpi = 300)) 
# 
