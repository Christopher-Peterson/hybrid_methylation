suppressPackageStartupMessages({
  library(tidyverse)
  library(rlang)
  library(glue)
  library(graphics)
  library(purrr)
  library(ggpointdensity)
  library(patchwork)
  library(parallel)
})
if(!exists('argv')) argv = commandArgs(TRUE)

# in_file = argv[1] %|% "out/joint_full_test.rds" #  "subset_diff_meth/parents_0.3_0.05_tssBoundaries.rds"

in_files = dir("subset_diff_meth", ".rds") |> file.path('subset_diff_meth', .=_)

make_plot = \(in_file, delta = NA_character_, fdr = NA_character_, out_file = NA_character_) {
  # Get relevant info from the file name
  sub_names = in_file |> basename() |> 
    str_match("parents_(0.[0-9]+)_(0.[0-9]+)_(.+).rds")
  delta = delta %|% sub_names[2]
  fdr = fdr %|% sub_names[3]
  
  plot_title = glue("Differential methylation, {sub_names[4]}")
  plot_subtitle = glue("Δ = {delta}, FDR = {fdr}")
  
  out_file = out_file %|% in_file |> 
    str_replace("subset_diff_meth/parents", "figures/diff_meth_dens") |> 
    str_replace(".rds", ".png")
  
  joint_data = read_rds(in_file)
  
  point_figure = joint_data |> 
    ggplot(aes(diff_p, diff_o)) +
    facet_wrap(~type) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_abline(slope = 0, intercept = 0, linetype = 1) +
    geom_pointdensity() + coord_fixed() + theme_classic() +
    scale_color_viridis_c("Number of Loci") + 
    ggtitle(plot_title, plot_subtitle) + 
    # scale_color_brewer(type = 'qual', palette = "Set1") +
    # scale_color_viridis_c("Offspring\nDifference SE") +
    xlab("Parental Methylation Difference") +
    ylab("F1 Methylation Difference")
  
  disc_figure = joint_data |> 
    mutate(diff_p_disc = 
             if_else(diff_p < 0, 
                     "Negative\n(Δ < -0.5)", "Positive\n(Δ > 0.5)")) |> 
    ggplot(aes(diff_p_disc, diff_o)) +
    facet_wrap(~type) +
    coord_fixed() + 
    geom_violin(fill = grey(.7)) + theme_classic() +
    geom_abline(slope = 0, intercept = c(-1, 1), linetype = 2) +
    geom_abline(slope = 0, intercept = 0, linetype = 1) +
    xlab("Parental Methylation Difference (Binned)") +
    ylab("F1 Methylation Difference (Density)")
  out_figure = point_figure / disc_figure
  
  
  
  # Just try it a few times?
  ggsave(out_file, out_figure, width = 12, height = 10, dpi = 300)
}

mclapply(in_files, make_plot, mc.cores = 6)

make_plot2 = \(in_file, nm = "all loci", delta = 0.5, fdr = 0.05, out_file = "figures/all_split.png") {
  # Get relevant info from the file name
  
  plot_title = glue("Differential methylation, {nm}")
  plot_subtitle = glue("Δ = {delta}, FDR = {fdr}")
  
  # out_file = out_file %|% in_file |> 
  #   str_replace("subset_diff_meth/parents", "figures/diff_meth_dens") |> 
  #   str_replace(".rds", ".png")
  
  joint_data = read_rds(in_file)
  
  point_figure = joint_data |> 
    ggplot(aes(diff_p, diff_o)) +
    facet_wrap(~type + pair, nrow = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_abline(slope = 0, intercept = 0, linetype = 1) +
    geom_pointdensity() + coord_fixed() + theme_classic() +
    geom_density_2d(breaks = c(0,5, 25, 50, 100, 200, 400, 600, 800, 1000),
                    color = 'white')+
    scale_color_viridis_c("Number of Loci", trans = 'log10') + 
    ggtitle(plot_title, plot_subtitle) + 
    theme(legend.position = 'bottom') + 
    # scale_color_brewer(type = 'qual', palette = "Set1") +
    # scale_color_viridis_c("Offspring\nDifference SE") +
    xlab("Parental Methylation Difference") +
    ylab("F1 Methylation Difference")
  # 
  # disc_figure = joint_data |> 
  #   mutate(diff_p_disc = 
  #            if_else(diff_p < 0, 
  #                    "Negative\n(Δ < -0.5)", "Positive\n(Δ > 0.5)")) |> 
  #   ggplot(aes(diff_p_disc, diff_o)) +
  #   facet_wrap(~type) +
  #   coord_fixed() + 
  #   geom_violin(fill = grey(.7)) + theme_classic() +
  #   geom_abline(slope = 0, intercept = c(-1, 1), linetype = 2) +
  #   geom_abline(slope = 0, intercept = 0, linetype = 1) +
  #   xlab("Parental Methylation Difference (Binned)") +
  #   ylab("F1 Methylation Difference (Density)")
  # out_figure = point_figure / disc_figure
  
  
  
  # Just try it a few times?
  ggsave(out_file, point_figure, width = 16, height = 10, dpi = 300)
}

make_plot2(in_file)

# make_plot(in_file, '0.5', fdr = '0.05', out_file = 'figures/joint_full_0.5.png')
# make_plot = \(in_file, delta = NA_character_, fdr = NA_character_, out_file = NA_character_) {
#   # Get relevant info from the file name
#   sub_names = in_file |> basename() |> 
#     str_match("parents_(0.[0-9]+)_(0.[0-9]+)_(.+).rds")
#   delta = delta %|% sub_names[2]
#   fdr = fdr %|% sub_names[3]
#   
#   plot_title = glue("Differential methylation, {sub_names[4]}")
#   plot_subtitle = glue("Δ = {delta}, FDR = {fdr}")
#   
#   out_file = out_file %|% in_file |> 
#     str_replace("subset_diff_meth/parents", "figures/diff_meth_dens") |> 
#     str_replace(".rds", ".png")
#   
#   joint_data = read_rds(in_file)
#   
#   point_figure = joint_data |> 
#     ggplot(aes(diff_p, diff_o)) +
#     facet_wrap(~type) +
#     geom_abline(slope = 1, intercept = 0, linetype = 2) +
#     geom_abline(slope = 0, intercept = 0, linetype = 1) +
#     geom_pointdensity() + coord_fixed() + theme_classic() +
#     scale_color_viridis_c("Number of Loci") + 
#     ggtitle(plot_title, plot_subtitle) + 
#     # scale_color_brewer(type = 'qual', palette = "Set1") +
#     # scale_color_viridis_c("Offspring\nDifference SE") +
#     xlab("Parental Methylation Difference") +
#     ylab("F1 Methylation Difference")
#   
#   disc_figure = joint_data |> 
#     mutate(diff_p_disc = 
#              if_else(diff_p < 0, 
#                      "Negative\n(Δ < -0.5)", "Positive\n(Δ > 0.5)")) |> 
#     ggplot(aes(diff_p_disc, diff_o)) +
#     facet_wrap(~type) +
#     coord_fixed() + 
#     geom_violin(fill = grey(.7)) + theme_classic() +
#     geom_abline(slope = 0, intercept = c(-1, 1), linetype = 2) +
#     geom_abline(slope = 0, intercept = 0, linetype = 1) +
#     xlab("Parental Methylation Difference (Binned)") +
#     ylab("F1 Methylation Difference (Density)")
#   out_figure = point_figure / disc_figure
#   
#   
#   
#   # Just try it a few times?
#   ggsave(out_file, out_figure, width = 12, height = 10, dpi = 300)
# }

# walk(in_files, make_plot)  
  # Note: may need to swap some directions on the non-hybrids