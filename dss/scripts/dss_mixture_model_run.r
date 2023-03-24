suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(readr)
  library(dplyr)
  library(rlang)
  library(posterior)
})


# This needs to be run on a full ls6 node w/ no other processes
if(!exists('argv')) argv = commandArgs(TRUE)
# Set the cmdstanr path
# Make sure to install CMDSTANR there first
Sys.getenv('SCRATCH') |> file.path('cmdstan/cmdstan-2.30.1/') |> set_cmdstan_path()

data_file = argv[1]            %|% "out/dss_filtered_uncorrelated_0.25_N8.bed" 
out_file = argv[2]             %|% "out/mixture_model_uncorrelated_fit_0.25_N8.rds"
seed = argv[3] |> as.integer() %|% 99723L
set.seed(seed)

input_data = read_tsv(data_file, col_select = c(`#chrom`, start, stop, offspring, role, 
                                                delta_p, delta.se_p, delta_o, delta.se_o)) |> 
  # Data are duplicated by role, so filter to just one
  filter(role == 'maternal') |> select(-role)

# Model likelihood
mix = mixture(gaussian, gaussian, order = FALSE)

base_prior = # c('mu1', 'mu2')) + 
  set_prior('normal(0, 1)', class = 'Intercept', resp = 'deltao', dpar = c( 'mu2')) + #,
  set_prior('normal(0, 2)', class = 'Intercept', resp = 'deltap') +
  set_prior('student_t(7, 0, 0.5)', class = 'sigma' , resp = 'deltap' ) + 
  set_prior('student_t(7, 0, 0.5)', class = 'sigma1' , resp = 'deltao' ) + 
  set_prior('normal(0, .75)', class = 'meanme', resp = c('deltap', 'deltao')) +
  set_prior('student_t(7, 0, 0.25)', class = 'sdme', resp = c('deltap', 'deltao'), lb = 0) 

### Original version
if(FALSE) {
  
  model_form  = bf(delta_o |  mi(delta.se_o) ~ 1, mu1 ~ (1|offspring), 
                   mu2 ~ mi(delta_p) * offspring , sigma2 = 'sigma1', 
                   theta1 ~  offspring, theta2 = 1,
                   family = mix) + 
    bf(delta_p | mi(delta.se_p) ~ 1, family = gaussian()) + set_rescor(FALSE)
  
  # get_prior(vtheatab_form, all_dat)
  
  
  full_prior = base_prior + 
    set_prior('normal(0, 1)', #'normal(.75, 0.25)'),
              class = 'b', resp = 'deltao', dpar = 'mu2') + 
    set_prior('normal(0,1)', class = 'Intercept', resp = 'deltao', dpar='theta1') +
    set_prior('student_t(7, 0, 1)', class = 'sd', resp = 'deltao', dpar = 'mu1') +
    set_prior('normal(0,1)', class = 'b', resp = 'deltao', dpar = 'theta1')
  # set_prior('student_t(7, 0, .25)', class = 'sd',resp = 'deltao', dpar='theta1')
  # validate_prior(full_prior, model_form, input_data, mix)
  
} else  {
  
  model_form  = bf(delta_o |  mi(delta.se_o) ~ 1, 
                   mu1 ~ 0, 
                   mu2 ~ mi(delta_p) * offspring , 
                   sigma2 = 'sigma1', 
                   theta1 ~  offspring, theta2 = 1,
                   family = mix) + 
    bf(delta_p | mi(delta.se_p) ~ 1, family = gaussian()) + set_rescor(FALSE)
  
  # get_prior(vtheatab_form, all_dat)
  
  
  full_prior = base_prior + 
    set_prior('normal(0, 1)', #'normal(.75, 0.25)'),
              class = 'b', resp = 'deltao', dpar = 'mu2') + 
    set_prior('normal(0,1)', class = 'Intercept', resp = 'deltao', dpar='theta1') +
    # set_prior('student_t(7, 0, 1)', class = 'sd', resp = 'deltao', dpar = 'mu1') +
    set_prior('normal(0,1)', class = 'b', resp = 'deltao', dpar = 'theta1')
  # set_prior('student_t(7, 0, .25)', class = 'sd',resp = 'deltao', dpar='theta1')
  validate_prior(full_prior, model_form, input_data, mix)
  
}
model_run = brm(model_form, input_data, mix, full_prior, chains = 4 , cores = 4,
                 warmup = 1500, iter=3000, seed  = seed, init = 0,
                 backend = 'cmdstanr' , threads = threading(16))

write_rds(model_run, out_file)
# model_run = read_rds(out_file)
# sink('logs/model_2_smry.txt') ; print(model_run); sink()
# vthetab_run_a

# 
# 
# vtheatab_form2 = bf(delta_o |  mi(delta.se_o) ~ 1, mu1 ~ (1|xpair), 
#                    mu2 ~ (1+ mi(delta_p) | xpair) + mi(delta_p) , sigma2 = 'sigma1', 
#                    theta1 ~  xpair, theta2 = 1,
#                    family = mix) + 
#   bf(delta_p | mi(delta.se_p) ~ 1, family = gaussian()) + set_rescor(FALSE)
# get_prior(vtheatab_form2, all_dat, mix) 
# vthetab_prior2 = base_prior + 
# #  set_prior('normal(0, 1)', class = 'b', coef = 'midelta_p', resp = 'deltao', dpar = 'mu2') +
#   set_prior('normal(0,1)', class = 'Intercept', resp = 'deltao', dpar='theta1') +
# 
#   set_prior('student_t(7, 0, 1)', class = 'sd', resp = 'deltao', 
#             coef = 'Intercept',
#             dpar = c('mu1', 'mu2')) +
#   set_prior('student_t(7, 0, 0.25)', class = 'sd',
#             coef = 'midelta_p',  resp = 'deltao', dpar = 'mu2') +
#   set_prior('normal(0,1)', class = 'b', resp = 'deltao', dpar = 'theta1')
# # set_prior('student_t(7, 0, .25)', class = 'sd',resp = 'deltao', dpar='theta1')
# validate_prior(vthetab_prior2, vtheatab_form2, all_dat, mix) #|>
#   # select(-nlpar, -lb, -ub)
#   # filter(source == 'default')
# 
# 
# vthetab_run_b = brm(vtheatab_form2, all_dat, mix, vthetab_prior2, chains = 4 , cores = 4,
#                     warmup = 1000, iter=2000, seed  = seed, init = 0,
#                     backend = 'cmdstanr' , threads = threading(16))
# file.remove("out/hybrid_mixture_model_b.rds")
