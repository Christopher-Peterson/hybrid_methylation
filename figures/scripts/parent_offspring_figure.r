# inter-generational differences
source('figures/scripts/setup_figures.r')

# Read the common data in ####
source('dss/scripts/standardize_data.r')

# Add GBM information
gbm_loci = read_tsv(gbm_file) |> select(`#chrom`, start) |> distinct() |> mutate(gbm = TRUE)
gbm_base_data = base_data |> left_join(gbm_loci) |> 
  mutate(gbm = !is.na(gbm)) |> 
  mutate(diff_delta = diff_delta * sign(delta_p))


# Filtering based on GBM doesn't change anything
base_mp_plot = gbm_base_data |># filter(gbm) |>
  select(offspring, `#chrom`, stop, delta_po, role, delta_p, delta_o, diff_delta) |>
  pivot_wider(names_from = role, values_from = delta_po) |> 
  ggplot(aes(x = maternal, y = paternal)) +#, color = delta_p)) +# abs(delta_o))) + 
  coord_fixed() + 
  scale_x_continuous('Maternal-specific intergenerational<br>methylation change (<sub></sub>Δ<sub>f</sub>)', limits = c(-1,1)) + 
  scale_y_continuous('Paternal-specific intergenerational<br>methylation change (<sub></sub>Δ<sub>m</sub>)', limits = c(-1,1)) + 
  theme(legend.position = c(1,1),#c(.5,.95), 
        legend.justification = c(1,.4),
        legend.direction = 'horizontal', legend.background = element_blank())
mp_pd_plot = base_mp_plot + 
  scale_color_viridis_c(labels$point_dens)+
  ggtitle('A') +
  geom_pointdensity() 
mp_o_plot = base_mp_plot +   
  scale_color_viridis_c(expression(
    "|"*Delta[o]*"|"
  ), option = 'mako')+
  ggtitle('C') +
  geom_point(aes(color = abs(delta_o))) #+ theme(legend.text = element_element_textbox())
mp_p_plot = base_mp_plot +   
  scale_color_viridis_c(expression(
    Delta[p]
  ), option = 'turbo')+
  ggtitle('B') +
  geom_point(aes(color = delta_p)) #+ theme(legend.text())
mp_diff_plot = base_mp_plot +   
  scale_color_viridis_c(expression(eta), option = 'rocket', direction = -1)+
  ggtitle('D') +
  geom_point(aes(color = diff_delta))
ggsave(out_files$inter_gen_v2, 
       mp_pd_plot   + mp_p_plot+ mp_o_plot + mp_diff_plot,#+ plot_annotation(tag_levels = 'A'), 
       dpi = 300, width = 10, height = 10)

### Old version, don't use ####
#' # This is outdated
#' #' Make a single panel of the double-delta plot
#' #' @param prole parental role, maternal or paternal
#' #' @param psign parental sign, -1 or 1
#' #' @param max_color max value of the color scale, for consistency
#' make_double_delta_plot= \(prole, psign, max_color = 1000) {
#'   base_data |> 
#'     filter(sign(delta_p) == psign, role == prole) |> 
#'     mutate(parental_sign = if_else(sign(delta_p) == -1,
#'                                    'Higher Paternal\nMethylation', 'Higher Maternal\nMethylation')) |>
#'     ggplot(aes(diff_delta*sign(delta_p), delta_po)) + 
#'     coord_fixed() +
#'     facet_wrap(~parental_sign, strip.position = 'top') + 
#'     theme(strip.background = element_blank(), strip.placement = 'outside')+
#'     fixed_lines[c('h', 'v')] + 
#'     geom_pointdensity() + 
#'     theme(legend.position = 'bottom') +
#'     scale_color_viridis_c(labels$point_dens , 
#'                           trans = 'log10', limits = c(1, max_color)) + 
#'     scale_x_continuous(labels$diff_delta_wrap,
#'                        limits = base_data |> filter(sign(delta_p) == psign) |>
#'                          with(diff_delta*sign(delta_p)) |> range() ) + 
#'     scale_y_continuous(labels$delta_po_wrap[[prole]],
#'                        limits = range(base_data$delta_po)) 
#' }
#' #' @param prole parental role
#' make_delta_plot = \(prole) {
#'   base_data |> filter(role == prole) |>
#'     mutate(role = recode(role, paternal = "Paternal Alleles", 
#'                          maternal = "Maternal Alleles")) |> 
#'     ggplot(aes(x = delta_p)) +
#'     fixed_lines[c('h', 'diag')] +
#'     geom_point(aes(y = delta_o, color = delta_po)) + 
#'     facet_wrap(~role, strip.position = 'left') +
#'     coord_fixed() +
#'     theme(strip.background = element_blank(),
#'           strip.placement = 'outside',
#'           legend.text = element_blank()
#'     )+
#'     scale_color_viridis_c(labels$delta_po_expr, option = 'magma',
#'                           limits = range(base_data$delta_po)) +
#'     labs(x = labels$delta_p_wrap, y = labels$delta_o_wrap ) + 
#'     theme(legend.position = 'bottom') +
#'     xlim(-1, 1) + ylim(-1, 1)
#' }
#' # Theme elements
#' noxt = theme(axis.title.x = element_blank())
#' noyt = theme(axis.title.y = element_blank())
#' nost = theme(strip.text = element_blank())
#' bigst = theme(strip.text = element_text(size = 12))
#' intergen_plot1 = 
#'   make_delta_plot('maternal') + noxt + bigst + 
#'   make_delta_plot('paternal') + bigst +
#'   make_double_delta_plot('maternal', -1) + noxt + bigst+
#'   make_double_delta_plot('paternal', -1)   + nost+
#'   make_double_delta_plot('maternal', 1) + noxt + noyt+ bigst + 
#'   make_double_delta_plot('paternal', 1) + noyt  + nost+
#'   plot_layout(guides = 'collect', nrow = 2, byrow = FALSE, tag_level = 'new')  +
#'   plot_annotation(tag_levels = 'A') &
#'   theme(legend.position = 'bottom')
#' ggsave(out_files$inter_gen_v1, intergen_plot1, dpi = 300, width = 10, height = 8)
