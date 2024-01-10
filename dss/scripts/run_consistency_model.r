# This models consistency of heritability

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(readr)
  library(purrr)
  library(dplyr)
  library(rlang)
  library(glue)
  library(rlang)
  library(tidyr)
  library(posterior)
  library(stringr)
})
# Set the cmdstanr path
# Make sure to install CMDSTANR there first
Sys.getenv('SCRATCH') |> file.path('cmdstan/cmdstan-2.31.0/') |> set_cmdstan_path()

if(!exists('argv')) argv = commandArgs(TRUE)

# Input data w/ all loci
data_file = argv[1] %|%  "model_data/dss_filtered_uncorrelated_0.25_N8.bed"
gbm_data_file = argv[2] %|%  "model_data/gbm_loci_uncorrelated_0.25_N8.bed"
# Files to be created
out_fit = argv[3]             %|% "model_data/consist_model_autocor_run_og.rds"
out_smry = argv[4]            %|% "model_data/consist_model_autocor_smry_og.rds"

seed = argv[5] |> as.integer() %|% 18245L
posterior_cis = argv[6] %|% 'model_data/consist_model_autocor_cis_og.rds'
set.seed(seed)
stan_file = 'stan_models/consistency_varying_sigma_autocor.stan'

N = 8L
delta = 0.25

### Set up the data ####
# withr::local_dir('dss') # For local use
source('scripts/standardize_data.r')

# Get the GBM meta-data
gbm_loci = read_tsv(gbm_data_file) |> filter(role == 'maternal') |> select(chrom = `#chrom`, start) |> 
  distinct() |> mutate(gbm = TRUE)
consist_model_data_gbm = left_join(full_consist_model_data, gbm_loci, by = c('chrom', 'start')) |> 
  mutate(is_gbm = !is.na(gbm),
         sgn = sign(delta_p))
# withr::deferred_run() # Clear out temporary wd
# full_consist_model_data # exists now

### Brms setup; ####
# We're setting up a brms model as a scaffold
# The current stan model was based off of it
# and brms can easily put the data in (mostly) the right form
# reminder: sdelta_o is sgn * delta_o.  You wasted half a day trying to figure out why the math wasn't right but the model gave good eta estimates because you forgot this
model_formula = 
  bf(sdelta_o |  mi(delta.se_o) ~ mi(delta_p):sgn + (1 | locus)) +
  bf(delta_p | mi(delta.se_p) ~ 1, family = gaussian()) + 
  set_rescor(FALSE)
# get_prior(model_formula, full_consist_model_data, family = gaussian())

# 
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

### Setup the multi-GP component ####
# I'll be using a 5/2 Matern kernel with a GP to account for spatial autocorrelation
gp_data = local({
  chrom_dat = consist_model_data_gbm |> group_by(chrom) |> distinct(locus, start) |> 
    mutate(position = start / 1000) |>  select(-start) |> 
    # chop(locus:position) |> 
    mutate(N = n(), begin = min(locus)) |> 
    ungroup() 
  chopped = chop(chrom_dat, locus:position)
  list(
    N_chrom = nrow(chopped),
    positions = chrom_dat$position,
    chrom_starts = chopped$begin,
    chrom_n = chopped$N,
    
    gp_nugget = 1e-6, # small number to add to diags for numerical stability
    
    # gp_len_priors = c(0, 2) # This is bad
    gp_len_priors = c(-1, 1)
    )
# Length-scale prior:
# The arabadobsis paper has methylation autocorrelation of at least 5kbp, 
# and it looks like it stretches reasonbly beyond that
# try exponential?
# I'll be uptting a prior on log10 length scale, such that 
# This means that about 95% of the prior probability for length scale is between 100 and 100,000 bp
})

  
### Needed for the stan run ####

base_stan_data =  make_standata(model_formula, consist_model_data_gbm,
                                gaussian(), base_prior, 
                                threads = 32)
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
    data = c(base_stan_data, extra_data, gp_data),
    seed = seed,
    chains = 4,
    adapt_delta = 0.9995,
    max_treedepth = 12,
    parallel_chains = 4,
    iter_warmup = 3000, iter_sampling = 1000
  )


# smry = cmdstan_fit$summary()
write_rds(cmdstan_fit, out_fit)
# cmdstan_fit = read_rds(out_fit)
cmdstan_draws = cmdstan_fit |> as_draws()
write_rds(cmdstan_draws, 'model_data/consist_model_autocor_draws_og.rds')
# cmdstan_draws = read_rds('model_data/consist_model_autocor_draws.rds')
smry_diag = summarize_draws(cmdstan_draws, default_convergence_measures()) 
smry_quant = summarize_draws(cmdstan_draws, 'median', 
                             ~quantile2(.x, probs = c(0.025, 0.05, 0.95, 0.975)))

# Left_join by variable
smry_full = smry_quant |> left_join(smry_diag)

# smry_full$variable |> str_remove_all('\\[.+') |> unique() |> sort()


write_rds(smry_full, out_smry)
# smry_full = read_rds(out_smry)
# Some tests to run
smry_full |> filter(rhat >= 1.01)

# smry$variable |> unique() |> str_remove('\\[.+\\]') |> unique()



# Report the important posteriors
round3 = \(x) round(x, 3) |> as.character()
round3_perc = \(x) { paste0(round(x,3) *100, '%')}
round3_square = \(x) round3(x^2)
important_variables = tribble(~'name', ~'variable', ~'fn',
  'Autocorrelation length scale', 'gp_log10_length_scale', \(x) round(10^x * 1000) |> as.character(),
  'Autocorrelation variance proportion', 'gp_prop', round3_perc,
  'Within-locus variation, intergenic','sigma_sdeltao[1]', round3_square,
  'Within-locus variation, gene body','sigma_sdeltao[2]', round3_square,
  'Among-locus variation, intergenic',  'sd_1[1]', round3_square,
  'Among-locus variation, gene body',  'sd_1[2]', round3_square
)
important_variables |> left_join(smry_full, by = 'variable') |> 
  select(name, fn, median, q2.5, q97.5) |> 
  pivot_longer(c(median, q2.5, q97.5), names_to = 'var', values_to = 'val') |> 
  mutate(val = map2_chr(fn, val, \(f,v) f(v))) |>
  select(-fn) |> 
  pivot_wider(names_from = var, values_from = val) |> 
  mutate(ci = glue('{median} [95% CI: {q2.5} - {q97.5}]')) |> 
  select(name, ci) |> 
  write_tsv(posterior_cis)
# Get the counts of non-heritables
# Quick chi-squared test of significatly different locus data
# p = 0.5; perfect
# Non-herit: 89 of 231 intergenic; 110 of 231 gbm

  

# GBM vs non-GBM sigmas
smry_full |> filter(str_detect(variable, 'sigma_sdeltao\\[')) |> 
  mutate(across(where(is.numeric), \(x)round(x,3))) |> 
  glue_data('{variable} {median} [{q2.5}, {q97.5}]')
# Difference between means
smry_full |> filter(str_detect(variable, 'gene_effect'))
# average heritability
smry_full |> filter(str_detect(variable, 'Intercept_sdeltao')) |> 
  mutate(across(where(is.numeric), \(x)round(x,3))) |> 
  glue_data('{variable} {median} [{q2.5}, {q97.5}]')
smry_full |> filter(str_detect(variable, 'sd_')) |> 
  mutate(across(where(is.numeric), \(x)round(x,3))) |> 
  glue_data('{variable} {median} [{q2.5}, {q97.5}]')

# Intercept_sdeltao -0.275 [-0.34, -0.215]
# sd_1[1] 0.422 [0.373, 0.475]
# sd_1[2] 0.256 [0.205, 0.304]
# sigma_sdeltao[1] 0.059 [0.044, 0.079]
# sigma_sdeltao[2] 0.054 [0.039, 0.073]
# 