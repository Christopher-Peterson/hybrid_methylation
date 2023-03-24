suppressPackageStartupMessages({
  library(brms)
  # library(cmdstanr)
  library(readr)
  library(dplyr)
  library(rlang)
  library(posterior)
  library(tidyr)
  library(purrr)
  library(stringr)
})
if(!exists('argv')) argv = commandArgs(TRUE)
in_file =      argv[1] %|% "out/mixture_model_uncorrelated_fit_0.25_N8.rds"
out_draws =    argv[2] %|% "out/mixture_model_uncorrelated_draws_0.25_N8.rds"
out_theta =    argv[3] %|% "out/mixture_model_uncorrelated_theta_quantiles_0.25_N8.csv"
out_regline =  argv[4] %|% "out/mixture_model_uncorrelated_regression_line_0.25_N8.csv"

#### Functions ####
inv_logit = \(x){
  p <- exp(x)/(1 + exp(x))
  p <- ifelse(is.na(p) & !is.na(x), 1, p)
  p
}
# Quantile summary
qsmry=\(x, probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) {
  q = quantile(x, probs)
  tibble(!!!q)
}
# like rename_with
rename_with_fn = \(x, .fn) {
  nms = names(x)
  fixed_nms = c(".chain", ".iteration", ".draw")
  setNames(x, if_else(nms %in% fixed_nms, nms, .fn(nms)))
}
# like mutate(across())
mutate_across_variables = \(x, .what, .fn, ...) {
  # Make sure the .fn has x in its environment
  force(.fn)
  rlang::fn_env(.fn) = list2env(x)
  # Get the column names selected in .what
  selected_cols = as_tibble(x) |> 
    select({{.what}}) |> 
    select(!any_of(c(".chain", ".iteration", ".draw"))) |> 
    names() 
  # Build expressions for each function_eval
  exprs = selected_cols |> 
    set_names() |> 
    lapply(\(z) {
      sm = rlang::sym(z)
      rlang::expr( .fn(!!sm) )
    }) #|> set_names(selected_cols)
  # Run it
  rlang::expr(mutate_variables(x, !!!exprs)) |> eval()
}

#### Draws ####
brm_fit = read_rds(in_file)
brm_draws = as_draws_df(brm_fit)
slope_draws = brm_draws |> 
  subset_draws('bsp_', regex = TRUE) |> 
  rename_with_fn(\(x) x |> str_remove(".+offspring") |> str_replace("bsp.+p", "slope")) |> 
  mutate_across_variables(!contains('slope'), \(x) x + slope) |> 
  rename_variables(X1 = slope) |> 
  rename_with_fn(\(x) paste0("slope_", x))

intercept_draws = brm_draws |> 
  subset_draws('b_mu', regex = TRUE) |> 
  rename_with_fn(\(x) x |> str_remove("b_") |> str_remove("_deltao") |> str_remove("offspring")) |> 
  mutate_across_variables(!contains('Intercept'), \(x) x + mu2_Intercept) |> 
  rename_variables(mu2_X1 = mu2_Intercept) |> 
  rename_with_fn(\(x) str_replace(x, 'mu2', 'mu2_Intercept'))

theta_draws = brm_draws |> 
  subset_draws('b_theta1', regex = TRUE) |> 
  mutate_across_variables(!contains('Intercept'), \(x) x + b_theta1_deltao_Intercept) |> 
  mutate_across_variables(everything(), \(x) inv_logit(x)) |> 
  rename_with_fn(\(x) str_remove(x, "b_theta1_deltao_")) |> 
  rename_variables(offspringX1=Intercept) |> 
  rename_with_fn(\(x) str_replace(x, "offspring", "Theta1_"))

list(slope=slope_draws,intercept=intercept_draws,theta=theta_draws) |> 
  write_rds(out_draws)

### Theta Quantiles ####
theta_quants = theta_draws |> 
  pivot_longer(starts_with("Theta1"), names_to = 'par', values_to = 'draws')|> 
  group_by(par) |> 
  summarize(quants = qsmry(draws)) |> 
  unpack(quants) |> 
  separate(par, c('param', 'offspring'), sep = '_') |>
  mutate(component = "1")
theta_quants |>  write_csv(out_theta)


# Predict slopes ####
new_dat = tibble(delta.se_p = 0.001, 
                 delta_o = 1, delta.se_o = .0001,
                 pair = 1:8, offspring = paste0("X", 1:8), pos = 1:8) |> 
  purrr::pmap_dfr(\(...) tibble(..., delta_p = seq(-1, 1, by = .05)))

get_regline = function(dpar) {
  draws = posterior_epred(brm_fit, new_dat, resp = 'deltao', dpar = dpar)
  line = as_draws_df(draws) |> 
    pivot_longer(starts_with("..."),
                 names_to = '.row', 
                 values_to = '.value')|> 
    group_by(.row) |> 
    summarize(quants = qsmry(.value)) |> 
    unpack(quants) |> 
    left_join(new_dat |> mutate(.row = 1:n(), .row = paste0("...", .row) ) ,
              by = ".row") |> select(-.row) |> 
    select(offspring, delta_p, contains("%")) |> 
    mutate(dpar = dpar)
}
herit_line = get_regline(dpar = 'mu2')
non_herit_line = get_regline(dpar = 'mu1')

bind_rows(herit_line, non_herit_line) |> write_csv(out_regline)

#### Quantiles (old) #### 
# intercept_quants = intercept_draws |> 
#   pivot_longer(everything(), names_to = 'par', values_to = 'draws')|> 
#   group_by(par) |> 
#   summarize(quants = qsmry(draws)) |> 
#   unpack(quants) |>
#   mutate(pair = str_remove(par, 'mu2_Intercept_pair_') |> recode(mu1_Intercept = NA_character_),
#          component = str_remove(par, "_.*") |> str_remove('mu'),
#          param = "Intercept"
#   ) |> select(-par)
# 
# theta_quants = theta_draws |> 
#   pivot_longer(starts_with("Theta1"), names_to = 'par', values_to = 'draws')|> 
#   group_by(par) |> 
#   summarize(quants = qsmry(draws)) |> 
#   unpack(quants) |> 
#   separate(par, c('param', 'pair'), sep = '_') |>
#   mutate(component = "1", pair = str_remove(pair, 'p'))
#   # mutate(pair = str_remove(par, "xpairp") |> recode(Intercept = "1"), 
#   #        param = "Theta", component = "1") |> select(-par)
# # 
# # slope_quants = slope_draws |> 
# #   pivot_longer(everything(), names_to = 'pair', names_prefix = 'p', values_to = 'draws')|> 
# #   group_by(pair) |> 
# #   summarize(quants = qsmry(draws)) |> 
# #   unpack(quants) |>
# #   mutate(component = '2', param = 'slope')
# theta_quants |>  write_csv(out_quantile)
# bind_rows(theta_quants, intercept_quants, slope_quants) |> 

