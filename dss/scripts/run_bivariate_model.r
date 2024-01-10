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
Sys.getenv('SCRATCH') |> file.path('cmdstan2/cmdstan-2.31.0/') |> set_cmdstan_path()
  #file.path('cmdstan/cmdstan-2.30.1/') |> set_cmdstan_path()

if(!exists('argv')) argv = commandArgs(TRUE)

which_model = argv[1] %|% "deltap" |> match.arg(c('deltap', 'nop'))
# Input data w/ all loci
spp = argv[2] %|% 'Amil'
data_file = argv[3] %|%  "model_data/dss_filtered_uncorrelated_0.25_N8.bed"
# Files to be created
# out_fit = argv[3]             %|% "model_data/bivariate_model_run.rds"
# out_smry = argv[4]            %|% "model_data/bivariate_model_smry.rds"

out_fit = glue("model_data/bivariate_model_{spp}_{which_model}_run.rds")
out_smry = glue("model_data/bivariate_model_{spp}_{which_model}_smry.rds")

seed = argv[5] |> as.integer() %|% 509341L #+1
set.seed(seed)
stan_file = switch(which_model, deltap = 'stan_models/bivariate_model.stan', nop = 'stan_models/bivariate_model_nop.stan')

N = 8L
delta = 0.25

### Set up the data ####

source('scripts/standardize_data.r')
biv_data = recip_consist_data |> filter(species == spp)
# recip_consist_data


# data_file = argv[1]            %|% "temp/consistency_delta_0.25_N8.tsv"
# species = argv[2]               %|% 'Amil' |> match.arg(c('Amil', 'Asel'))
# out_file = argv[2]             %|% "temp/bivariate_model_trial.rds"
# seed = argv[3] |> as.integer() %|% 7823915L
# set.seed(seed)

# BRMS scafold setup
model_formula = # the deltap terms should be given an offset, but brms dones't like it; fixed in manual model
  bf(delta_o.1 |  mi(delta.se_o.1) ~ (mi(delta_p.1))) +
  bf(delta_o.2 |  mi(delta.se_o.2) ~ (mi(delta_p.1))) +
  set_rescor(FALSE) +
  bf( delta_p.1 | mi(delta.se_p.1) ~ 1, family = gaussian() ) #+

# get_prior(model_formula, family = gaussian(), data = data)
o_resps = c('deltao1', 'deltao2')
all_resps =  c(o_resps, 'deltap1')

base_prior = set_prior('normal(0, 1.5)', class = 'b', resp = o_resps) +
  set_prior('normal(0, 1.5)', class = 'Intercept', resp = all_resps) +
  set_prior('normal(0, .75)', class = 'meanme', resp = all_resps) +
  set_prior('student_t(7, 0, 0.25)', class = 'sdme',  lb = 0, resp = all_resps) +
  set_prior('student_t(7, 0, 1)', class = 'sigma', resp =all_resps, lb = 0)

base_prior |> validate_prior(model_formula, family = gaussian(), data = biv_data)

# get_prior(model_formula, family = gaussian(), data = recip_consist_data)


cmdstan_data = make_standata(model_formula, biv_data, gaussian(), base_prior )
cmdstan_model = cmdstan_model(stan_file)

cmdstan_fit = cmdstan_model$sample(
  data = cmdstan_data,
  seed = seed,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000
)

cmdstan_fit$summary() -> smry
write_rds(cmdstan_fit, out_fit)
write_rds(smry, out_smry)A

# smry |> filter(variable == 'r_corr')
# smry |> filter(str_detect(variable, 'resi')) |>
#   select(variable, median, q5, q95)
# The r_corr parameter isn't quite what it need
# It's looking specifically at 'residual' correlation not explained by Delta_p
#
# 
# cmdstan_fit = cmdstan_model$sample(
#   data = c(base_stan_data, extra_data),
#   seed = seed,
#   chains = 4,
#   adapt_delta = 0.999,
#   parallel_chains = 4,
#   iter_warmup = 1500, iter_sampling = 1500
# )

