suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggpointdensity)
  library(ggtext)
  # library(ggh4x)
  library(glue)
  # library(withr)
  
})

if(!exists('FIGURES_SETUP')) {
  # Prevent endless re-running
## Setup filepaths & such ####
  delta = 0.25
  N = 8
  # use_gbm = FALSE
  # if(isTRUE(use_gbm)) {
  fig_dir=glue('figures/delta_{delta}_N{N}')
  dir.create(fig_dir, showWarnings = FALSE)
  
  data_dir = 'dss/model_data'
  data_file = '{data_dir}/dss_filtered_uncorrelated_{delta}_N{N}.bed' |> glue()
  gbm_file = '{data_dir}/dss_filtered_uncorrelated_gbm_0.25_N8.bed' |> glue()
  nogene_file = '{data_dir}/dss_filtered_uncorrelated_nogene_0.25_N8.bed' |> glue()
  mix_model_data_all = list(
    regline = '{data_dir}/mixture_model_uncorrelated_regression_line_{delta}_N{N}.csv', 
    theta = '{data_dir}/mixture_model_uncorrelated_theta_quantiles_{delta}_N{N}.csv'
  ) |> map(glue)
  mix_model_data_gbm = list(
    regline = '{data_dir}/mixture_model_uncorrelated_gbm_regression_line_{delta}_N{N}.csv', 
    theta = '{data_dir}/mixture_model_uncorrelated_gbm_theta_quantiles_{delta}_N{N}.csv'
  ) |> map(glue)
  mix_model_data_nogene = list(
    regline = '{data_dir}/mixture_model_uncorrelated_nogene_regression_line_{delta}_N{N}.csv', 
    theta = '{data_dir}/mixture_model_uncorrelated_nogene_theta_quantiles_{delta}_N{N}.csv'
  ) |> map(glue)
  
  consist_smry_file = '{data_dir}/consist_model_smry.rds' |> glue()
  # Should we do GBM /non-gbm model fits as well?
  
  
  # Output files ####
  out_files = list(
    cross_diagram = 'cross_design.png',
    delta_guide = 'delta_guide.png',
    expected_results = 'expected_results.png',
    mixture_model = 'mixture_all.png',
    mixture_model_gbm = 'mixture_gbm.png',
    mixture_model_nogene = 'mixture_nogene.png',
    consistency_1d = 'consist.png',
    consistency_recip = 'recip_arrow.png',
    consistency_recip_noarrow = 'recip_noarrow.png',
    inter_gen_v2 = 'parent_offspring.png',
    mixture_model_table = 'table_1.rds'
  ) |> map(\(x) file.path(fig_dir, x))
  
  ## Common elements ####
  theme_set(theme_classic() + 
              theme(axis.title.x = element_markdown(),
                    axis.title.y = element_markdown()))
  fixed_lines = list(
    h = geom_hline(yintercept = 0, linetype = 3, color = grey(.6)),
    v = geom_vline(xintercept = 0, linetype = 3, color = grey(.6)),
    diag = geom_abline(slope = 1, intercept = 0, linetype = 3, color = grey(.6)),
    anti_diag = geom_abline(slope = -1, intercept = 0, linetype = 3, color = grey(.6))
  )
  # Axis labels
  # The empty <sub> tags before delta are used as spacing, as it apparently likes 
  # to bold anything directly touching the delta
  labels = list(
    point_dens = "Density of\nLoci in Area",
    delta_p = "Differential methylation between parents (<sub></sub>Δ<sub>p</sub>)",
    delta_o = "Differential methylation between offspring alleles (<sub></sub>Δ<sub>o</sub>)",
    delta_po = "Intergenerational methylation change (<sub></sub>Δ<sub>♀</sub> and Δ<sub>♂</sub>)",
    diff_delta_long = "Intergenerational differential methylation change [Sign(<sub></sub>Δ<sub>p</sub>)(<sub></sub>Δ<sub>o</sub> - Δ<sub>p</sub>)]",
    diff_delta = "Intergenerational differential methylation change",
    
    # delta_po_expr = expression(frac(plain("Intergenerational")~plain("methylation"),plain("change")~(Delta["\u2640"]~plain('and')~ Delta["\u2642"]))),
    delta_po_expr = expression(atop(plain("Intergenerational change in"), plain("methylation (")*Delta["\u2640"]~plain('and')~ Delta["\u2642"]*plain(")")) ),
    
    delta_po_wrap = list(
      maternal = "Intergenerational change in<br>methylation (<sub></sub>Δ<sub>♀</sub>)",
      paternal = "Intergenerational change in<br>methylation (<sub></sub>Δ<sub>♂</sub>)"
    ), 
    diff_delta_wrap = "Intergenerational differential<br>methylation change (<sub></sub>Δ<sub>o</sub> - Δ<sub>p</sub>)",
    delta_p_wrap = "Differential methylation<br>between parents (<sub></sub>Δ<sub>p</sub>)",
    delta_o_wrap = "Differential methylation<br>between offspring alleles (<sub></sub>Δ<sub>o</sub>)"
  )
  
  ### environmental variable to avoid re-sourcing
  FIGURES_SETUP = TRUE
}