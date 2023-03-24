# Make the Crossing design figure
# 
library(tidyverse)
source('figures/scripts/setup_figures.r')
library(ggtext)
library(ggh4x)

parents = c('s<sub>1</sub>', 's<sub>2</sub>', 'm<sub>1</sub>', 'm<sub>2</sub>') |> fct_inorder()
get_spp = function(x) x %in% parents[1:2]#, 'A. selago', 'A. millepora')

species = factor(c('<i>A. selago</i>', '<i>A. millepora</i>'))
cross_figure_data = expand_grid(
  dam = parents[c(4,3, 2,1)] |> fct_inorder(), 
  sire = parents[]) |> 
  filter(dam != sire) |> 
  mutate(d_sel = get_spp(dam), s_sel = get_spp(sire)) |> 
  mutate(F1 = case_when(
    d_sel & s_sel ~      "Pure<br><i>A. selago</i>",
    (!d_sel) & (!s_sel) ~ "Pure<br><i>A. millepora</i>",
    d_sel & (!s_sel) ~ "Hybrid",
    !d_sel & (s_sel) ~ "Inviable"
  ),
  d_spp = if_else(d_sel, species[1], species[2] ),
  s_spp = if_else(s_sel,  species[1], species[2] ),
  group = 1:n(), 
  dam_joint = paste(dam, d_spp, sep = '|') |> fct_inorder(),
  sire_joint = paste(sire, s_spp, sep  = '|') |> fct_inorder())

cross_figure = cross_figure_data |> 
  ggplot(aes(y = dam_joint, x = sire_joint, label = F1, color = F1, fill = F1, group = group)) + 
  geom_tile() +
  # facet_grid(d_spp~ s_spp, as.table = TRUE, switch = "y") + 
  scale_x_discrete("Sire", expand = c(0.01,0.01), position = 'top') +
  scale_y_discrete("Dam", expand = c(0.01,0.01)) +
  theme_classic() + 
  # geom_vline(xintercept = c(1.5, 2.5), color = 'black') +
  coord_fixed() + 
  theme(strip.background = element_blank(),
        strip.text = element_markdown(size = 15),
        panel.spacing = unit(0.03, 'cm'),
        ggh4x.axis.nesttext.x = element_markdown(size = 14),
        ggh4x.axis.nesttext.y = element_markdown(size = 14,angle = 90, hjust = 0.5),
        ggh4x.axis.nestline.x = element_line(linetype = 2, color = grey(.7)),
        ggh4x.axis.nestline.y = element_line(linetype = 2, color = grey(.7)),
        axis.text.x.top  = element_markdown(size = 11),
        axis.text.y.left = element_markdown(size = 11),
        axis.title = element_text(size = 15),
        panel.background = element_rect(fill = grey(.4)),
        strip.placement = 'outside') +
  scale_fill_manual(
    values = c("Pure<br><i>A. selago</i>" = 'lightcoral',
               "Pure<br><i>A. millepora</i>" = 'steelblue1',
               "Hybrid" = "violet",
               "Inviable" = grey(.4)),
    drop = FALSE, guide = 'none',
    
    na.value = grey(.5)
  ) + 
  geom_richtext(size = 4, fill = NA, label.colour = NA) +
  annotate('text', x = 1.5, y = 1.5, size = 8, color = 'white', label = 'Inviable\nCross') +
  annotate('tile', x = 1.5, y = 1.5, width = 2, height = 2, fill = NA, color = 'black')+
  scale_color_manual(guide = 'none',
    values = c("Pure<br><i>A. selago</i>"  ='black',
               "Pure<br><i>A. millepora</i>"  = 'black',
               "Hybrid" = "black",
               "Inviable" = alpha(grey(.4), 0))
  ) + guides(x = guide_axis_nested(delim = '|'), 
             y = guide_axis_nested(delim = '|'))
ggsave(out_files$cross_diagram,cross_figure, width = 5.5, height = 5.5, dpi = 500)

# Look into replacing facets w/ ggh4x underlines or someting? 

