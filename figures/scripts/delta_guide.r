# Delta term diagram
source('figures/scripts/setup_figures.r')

delta_diagram_data = tribble(
  ~'x', ~'y', ~'text',
  'maternal',  'parent',        'Dam',
  'paternal',  'parent',        'Sire',
  'maternal',  'offspring',     "Offspring\nDam\nAllele",
  'paternal',   'offspring',     "Offspring\nSire\nAllele",
)

arrow_text = \(between = c('parent', 'offspring', 'maternal', 'paternal')) {
  between = match.arg(between)
  # Arguments for geom_segment
  seg_lst = list(
    b = between, bend = between, s = 1.25, send = 1.75, 
    arrow = arrow(ends = 'both',length = unit(.05, 'npc')) ) 
  # Proper_names for them
  seg_nms = switch(between %in% c('parent', 'offspring') + 1,
                   mp = c('x', 'xend', 'y', 'yend', 'arrow'), # Paternal Maternal
                   po = c('y', 'yend', 'x', 'xend', 'arrow') # Parent offspring
  ) 
  
  geom_seg = do.call(geom_segment, set_names(seg_lst, seg_nms))
  sub = switch(between, parent = 'p', offspring = 'o', maternal = 'f', paternal = 'm')
  text_just = if_else(between == 'maternal', 1, 0)# switch(between, maternal = 0, )
  txt_args = list(
    fixed = between,
    mid = 1.5,
    label = glue::glue('Δ<sub>{sub}</sub>'),
    just = text_just
  ) 
  txt_nms = switch(between %in% c('parent', 'offspring') + 1,
                   mp = c('x', 'y', 'label', 'hjust'), # Paternal Maternal
                   po = c('y', 'x', 'label', 'vjust') # Parent offspring
  )
  # extra args to remove label box
  no_lab_box = list(fill = NA, label.color = NA, 
                    label.padding = grid::unit(rep(0, 4), "pt") )
  geom_txt = do.call(geom_richtext, txt_args |> set_names(txt_nms) |> c(no_lab_box))
  list(geom_seg, geom_txt)
}


delta_diagram = ggplot(delta_diagram_data) + 
  aes(x,y) +
  ggplot2::theme_void() + 
  theme(panel.background = element_rect(fill = 'white'))+
  # theme_nothing() + 
  geom_label(aes(label = text), fill = 'cadetblue1') +
  coord_fixed() + 
  arrow_text('parent') + 
  arrow_text('offspring') + 
  arrow_text('maternal') + 
  arrow_text('paternal') +
  annotate('text', x = 1.5, y = 2.3, label = "Four differential methylation\nscores per locus per cross", size = 5) 

ggsave(out_files$delta_guide, delta_diagram, dpi = 300, width = 4, height = 4)
delta_diagram