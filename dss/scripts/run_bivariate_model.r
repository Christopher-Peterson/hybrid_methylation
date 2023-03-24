# This models consistency of heritability

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(readr)
  library(purrr)
  library(dplyr)
  library(rlang)
  library(tidyr)
  library(posterior)
  library(stringr)
})
# Set the cmdstanr path
# Make sure to install CMDSTANR there first
Sys.getenv('SCRATCH') |> file.path('cmdstan/cmdstan-2.30.1/') |> set_cmdstan_path()

if(!exists('argv')) argv = commandArgs(TRUE)

# Input data w/ all loci
data_file = argv[1] %|%  "model_data/dss_filtered_uncorrelated_0.25_N8.bed"
# Files to be created
out_fit = argv[2]             %|% "model_data/bivariate_model_run.rds"
out_smry = argv[4]            %|% "model_data/bivariate_model_smry.rds"

seed = argv[5] |> as.integer() %|% 509341L
set.seed(seed)
stan_file = 'stan_models/bivariate_model.stan'

N = 8L
delta = 0.25

### Set up the data ####

source('scripts/standardize_data.r')




# INCORRECT FROM HERE:
# 
# 
# 
# 
# 


# full_consist_model_data # exists now

### Brms setup; ####
# We're setting up a brms model as a scaffold
# The current stan model was based off of it
# and brms can easily put the data in (mostly) the right form
model_formula = 
  bf(sdelta_o |  mi(delta.se_o) ~ mi(delta_p):sgn + (1 | locus)) +
  bf(delta_p | mi(delta.se_p) ~ 1, family = gaussian()) + 
  set_rescor(FALSE)
get_prior(model_formula, full_consist_model_data, family = gaussian())

mi_vars = c('deltap', 'sdeltao')
base_prior =
  set_prior("constant(1)", class = 'b', resp = 'sdeltao', coef = 'midelta_p:sgn') +
  set_prior('normal(0, 1.5)', class = 'Intercept', resp = mi_vars) +
  set_prior('normal(0, .75)', class = 'meanme', resp = mi_vars) +
  set_prior('student_t(7, 0, 0.25)', class = 'sdme',  lb = 0, resp = mi_vars) +
  set_prior('student_t(7, 0, 1)', class = 'sigma', resp =mi_vars, lb = 0) +
  set_prior('student_t(7, 0, 1)', class = 'sd', resp = 'sdeltao', lb = 0)

base_prior |> validate_prior(model_formula, family = gaussian(), data = full_consist_model_data)

post_proc = stanvar(scode = '
  // Diff_delta calculations
  vector[N] diff_delta = Yl_sdeltao - Intercept_sdeltao - (bsp_sdeltao[1] * Yl_sdeltap);
  vector[N] true_resid = diff_delta -  r_1_sdeltao_1[J_1_sdeltao] ;
  // Will need to do some post-processing to calculate the residual variance per locus
  real sd_diff_delta = sd(diff_delta);
  real sd_resid = sqrt(sd_diff_delta^2 - sd_1[1]^2);
', block = 'genquant')#

### Needed for the stan run ####

base_stan_data =  make_standata(model_formula, consist_model_data_gbm,
                                gaussian(), base_prior )
# Add extra data points:

extra_data = local({
  gbm_lgl = consist_model_data_gbm |> distinct(locus, is_gbm) |> 
    arrange(locus) |> pull(is_gbm)
  N_obs = nrow(consist_model_data_gbm)
  list(
    N_gene = sum(gbm_lgl),
    locus_assignment = gbm_lgl + 1L,
    which_gene = which(gbm_lgl),
    which_nogene = which(!gbm_lgl),
    sigma_concentration = 2
    # For some reason these aren't in the base data
    # Csp_sdeltao_1 = rep(1, N_obs),
    # N_deltap = N_obs,
    
  )
})


cmdstan_model = cmdstan_model(stan_file)

cmdstan_fit = cmdstan_model$sample(
  data = c(base_stan_data, extra_data),
  seed = seed,
  chains = 4,
  adapt_delta = 0.999,
  parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500
)


smry = cmdstan_fit$summary()
write_rds(cmdstan_fit, out_fit)
write_rds(smry, out_smry)
#

# 
# view_stancode = \(formula, data, family, prior, ...) {
#   tmp_file = tempfile(fileext = '.stan')
#   make_stancode(formula, data, family, prior, ...) |> write_lines(tmp_file)
#   file.edit(tmp_file)
#   invisible(tmp_file)
# }
# make_standata(model_formula, full_model_data, gaussian(), base_prior ) |> str(1)
# view_stancode(model_formula, full_model_data, gaussian(), base_prior, stanvars = post_proc)
# 
# model_run_test =
#   brm(model_formula, full_model_data, gaussian(), base_prior,
#       stanvars = post_proc,
#       chains = 4 , cores = 4, warmup = 1000, iter=2000, seed  = seed,
#       backend = 'cmdstanr' )
#
# cmdstan_data = make_standata(model_formula, data, gaussian(), base_prior )
# stan_file = 'temp/bivariate_model.stan'
# cmdstan_model = cmdstan_model(stan_file)
#
# cmdstan_fit = cmdstan_model$sample(
#   data = cmdstan_data,
#   seed = seed,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000, iter_sampling = 1000
# )
#
# cmdstan_fit$summary() -> smry
# write_rds(cmdstan_fit, 'temp/bivariate_fit.rds')
# write_rds(smry, 'temp/bivariate_smry.rds')
#
# smry |> filter(variable == 'r_corr')
# smry |> filter(str_detect(variable, 'resi')) |>
#   select(variable, median, q5, q95)
