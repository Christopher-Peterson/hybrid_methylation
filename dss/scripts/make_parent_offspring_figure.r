library(tidyverse)
library(ggpointdensity)
library(patchwork)

# Make this work specifically for paternal vs. maternal
# Show gain/loss & if it differs for paternal or maternal

# how to fix the paternal/maternal confusion?
# Look at the file names
# And look at the 
# the parent/offspring tsv indicates which offspring genome each of the parents belongs to
# the 


data = 'dss/out/dml_with_po.bed'
# Note: the dam and sire columns in this dataset are inaccurately named
po_data = read_tsv(data, col_select = -c(chrom2, start2, stop2)) |> 
  mutate(diff_delta = sign(delta_p) * (delta_p - delta_o), 
         se_op = sqrt(delta.se_p ^2 + delta.se_o ^ 2)) |> 
  mutate(flip_sign = if_else(offspring %in% c("X2", "X4"), -1, 1),
         across(c(delta_o, delta_p), \(x) x * flip_sign))

maternal_paternal = tribble(
  ~'offspring',	~'maternal',	~'paternal',
  'X1',	'A1',	'A2',
  'X2',	'A2',	'A1',
  'X3',	'A3',	'A4',
  'X4',	'A4',	'A3',
  'X5',	'A1',	'A3',
  'X6',	'A1',	'A4',
  'X7',	'A2',	'A3',
  'X8',	'A2',	'A4',
) |> pivot_longer(-offspring, names_to = 'role', values_to = 'parent') 

filtered_data = maternal_paternal |> 
  left_join(po_data,  by = c('offspring', 'parent')) |> 
  filter(offspring==offspring2) |> 
  select(-dam, -sire, -prob_5, -offspring2) |> 
  mutate(parental_sign = if_else(sign(delta_p) == -1, 'Higher Paternal\nMethylation', 
                                 'Higher Maternal\nMethylation'))

make_double_delta_plot= \(prole, psign, max_color = 400) {
  filtered_data |> 
    filter(sign(delta_p) == psign, role == prole) |> 
  # filtered_data |> 
    ggplot(aes(diff_delta, po_delta)) + 
    theme_classic() +
    coord_fixed() +
    facet_wrap(~parental_sign, strip.position = 'top') + 
    # facet_grid(parental_sign~role) +
    theme(strip.background = element_blank(), strip.placement = 'outside')+
    geom_hline(yintercept = 0, color = grey(.75), linetype = 3) + 
    geom_pointdensity() + 
    theme(legend.position = 'bottom') +
    scale_color_viridis_c("Number of\nNearby Loci",trans = 'log10', limits = c(1, max_color)) + 
    scale_x_continuous('Δ Parental - Δ Offspring\nMethylation',
                       limits = range(filtered_data$diff_delta)) + 
    scale_y_continuous("Parent - Offspring\nMethylation",
                       limits = range(filtered_data$po_delta)) 
}
make_delta_plot = \(prole) {
  filtered_data |> filter(role == prole) |>
    mutate(role = recode(role, paternal = "Paternal Alleles", 
                         maternal = "Maternal Alleles")) |> 
    ggplot(aes(x = delta_p)) +
    geom_abline(slope = 0, intercept = 0, linetype = 1, color = grey(.75)) +
    geom_abline(slope = 1, intercept = 0, linetype = 1, color = grey(.75)) +
    geom_point(aes(y = delta_o, color = po_delta)) + 
    facet_wrap(~role, strip.position = 'left') +
    coord_fixed() + theme_classic() +
    theme(strip.background = element_blank(), strip.placement = 'outside')+
    scale_color_viridis_c("Parent - Offspring\nMethylation",option = 'magma',
                          limits = range(filtered_data$po_delta)) +
    xlab("Δ Parental Methylation") +
    ylab("Δ Offspring Methylation") +
    theme(legend.position = 'bottom') +
    xlim(-1, 1) + ylim(-1, 1)
}
noxt = theme(axis.title.x = element_blank())
noyt = theme(axis.title.y = element_blank())
nost = theme(strip.text = element_blank())
bigst = theme(strip.text = element_text(size = 12))

full_plot = 
  make_delta_plot('maternal') + noxt + bigst + 
  make_delta_plot('paternal') + bigst +
  make_double_delta_plot('maternal', -1) + noxt + bigst+
  make_double_delta_plot('paternal', -1)   + nost+
  make_double_delta_plot('maternal', 1) + noxt + noyt+ bigst + 
  make_double_delta_plot('paternal', 1) + noyt  + nost+
  plot_layout(guides = 'collect', nrow = 2, byrow = FALSE) &
  theme(legend.position = 'bottom'); full_plot
ggsave('dss/figures/parent_offspring_methyl.png', full_plot, dpi = 300, width = 8, height = 7)


# filtered_data |> 
#   ggplot(aes(y = delta_o, x = po_delta, color = role)) + 
#   facet_wrap(~parental_sign) + 
#   theme_classic() + 
#   geom_point()

# Make 3 versions of this data w/ different color scales & send them to misha
base_mp_plot = filtered_data |> select(offspring, `#chrom`, stop, po_delta, role, delta_p, delta_o) |> 
  pivot_wider(names_from = role, values_from = po_delta) |> 
  ggplot(aes(x = maternal, y = paternal)) +#, color = delta_p)) +# abs(delta_o))) + 
  theme_classic() + 
  coord_fixed() + 
  scale_x_continuous('Maternal - specific offspring methylation', limits = c(-1,1)) + 
  scale_y_continuous('Paternal - specific offspring methylation', limits = c(-1,1)) + 
  theme(legend.position = 'bottom')
mp_pd_plot = base_mp_plot + 
  scale_color_viridis_c("Loci Density")+
  geom_pointdensity() 
mp_o_plot = base_mp_plot +   
  scale_color_viridis_c("|Δ| Offspring", option = 'mako')+
  geom_point(aes(color = abs(delta_o))) 
mp_p_plot = base_mp_plot +   
  scale_color_viridis_c("Δ parents", option = 'turbo')+
  geom_point(aes(color = delta_p)) 

ggsave('dss/figures/parent_offspring_tripanel.png', 
       mp_p_plot + mp_o_plot + mp_pd_plot,
       dpi = 300, width = 11, height = 6)
 # 
# double_delta_plot = filtered_data |> 
#   ggplot(aes(diff_delta, po_delta)) + 
#   theme_classic() +
#   coord_fixed() +
#   facet_grid(parental_sign~role) +
#   theme(strip.background = element_blank(), strip.placement = 'outside')+
#   geom_hline(yintercept = 0, color = grey(.75), linetype = 3) + 
#   geom_pointdensity() + 
#   scale_color_viridis_c("Number of\nNearby Loci",trans = 'log10') + 
#   xlab('Δ Parental - Δ Offspring Methylation') + 
#   ylab("Parent - Offspring Methylation")
# 
# 
# delta_plot = filtered_data |> 
#   ggplot(aes(x = delta_p)) +
#   geom_abline(slope = 0, intercept = 0, linetype = 1, color = grey(.75)) +
#   geom_abline(slope = 1, intercept = 0, linetype = 1, color = grey(.75)) +
#   geom_point(aes(y = delta_o, color = po_delta)) + 
#   facet_grid("tst"~role) +
#   coord_fixed() + theme_classic() +
#   theme(strip.background = element_blank())+
#   scale_color_viridis_c("Parent - Offspring\nMethylation",option = 'magma') +
#   xlab("Δ Parental Methylation") +
#   ylab("Δ Offspring Methylation")
# cowplot::get_legend(delta_plot) + cowplot::get_legend(double_delta_plot)
# 
# cowplot::plot_grid(delta_plot + theme(legend.position = 'none'), 
#                    double_delta_plot + theme(legend.position = 'none'),
#                    align = 'vh', axis = 'lrtb', 
#                    rel_heights = c(1,2), nrow = 2, ncol = 1, greedy = FALSE)
# # 
# delta_plot + double_delta_plot +  plot_layout(guides = 'collect', ncol = 1, design = layout)
