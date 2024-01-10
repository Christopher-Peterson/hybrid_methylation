// generated with brms 2.17.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_deltao1;  // number of observations
  vector[N_deltao1] Y_deltao1;  // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_deltao1] noise_deltao1;
  // information about non-missings
  int<lower=0> Nme_deltao1;
  int<lower=1> Jme_deltao1[Nme_deltao1];
  int<lower=1> Ksp_deltao1;  // number of special effects terms
  int<lower=1> N_deltao2;  // number of observations
  vector[N_deltao2] Y_deltao2;  // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_deltao2] noise_deltao2;
  // information about non-missings
  int<lower=0> Nme_deltao2;
  int<lower=1> Jme_deltao2[Nme_deltao2];
  int<lower=1> Ksp_deltao2;  // number of special effects terms
  // int<lower=1> N_deltap1;  // number of observations
  // vector[N_deltap1] Y_deltap1;  // response variable
  // data for measurement-error in the response
  // vector<lower=0>[N_deltap1] noise_deltap1;
  // information about non-missings
  // int<lower=0> Nme_deltap1;
  // int<lower=1> Jme_deltap1[Nme_deltap1];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int N_corr = 2;  // number of correlations
  vector[Ksp_deltao1] bsp_deltao1 = [1]';  // special effects coefficients
  vector[Ksp_deltao2] bsp_deltao2 = [-1]';  // special effects coefficients
}
parameters {
  vector[N_deltao1] Yl_deltao1;  // latent variable
  real Intercept_deltao1;  // temporary intercept for centered predictors
  // vector[Ksp_deltao1] bsp_deltao1;  // special effects coefficients
  real<lower=0> sigma_deltao;  // dispersion parameter
  vector[N_deltao2] Yl_deltao2;  // latent variable
  real Intercept_deltao2;  // temporary intercept for centered predictors
  // vector[Ksp_deltao2] bsp_deltao2;  // special effects coefficients
  // real<lower=0> sigma_deltao2;  // dispersion parameter
  // vector[N_deltap1] Yl_deltap1;  // latent variable
  // real Intercept_deltap1;  // temporary intercept for centered predictors
  // real<lower=0> sigma_deltap1;  // dispersion parameter
  cholesky_factor_corr[N_corr] L_corr; // Correlation
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior



  lprior += normal_lpdf(Intercept_deltao1 | 0, 1.5);
  lprior += normal_lpdf(bsp_deltao1 | 0, 1.5);
  lprior += student_t_lpdf(sigma_deltao | 7, 0, 1)
    - 1 * student_t_lccdf(0 | 7, 0, 1);
  lprior += normal_lpdf(Intercept_deltao2 | 0, 1.5);
  lprior += normal_lpdf(bsp_deltao2 | 0, 1.5);
  // lprior += student_t_lpdf(sigma_deltao2 | 7, 0, 1)
    // - 1 * student_t_lccdf(0 | 7, 0, 1);
  // lprior += normal_lpdf(Intercept_deltap1 | 0, 1.5);
  // lprior += student_t_lpdf(sigma_deltap1 | 7, 0, 1)
    // - 1 * student_t_lccdf(0 | 7, 0, 1);
  lprior += lkj_corr_cholesky_lpdf(L_corr | 1.0) ;
}
model {

  // likelihood including constants
  if (!prior_only) {
    matrix[N_corr, N_corr] L_Sigma = diag_pre_multiply(rep_vector(sigma_deltao, N_corr), L_corr);
    array[N] row_vector[N_corr] Yl_mat;
    // initialize linear predictor terms
    array[N] row_vector[N_corr] mu_deltao = rep_array([Intercept_deltao1, Intercept_deltao2], N);
    // vector[N_deltap1] mu_deltap1 = Intercept_deltap1 + rep_vector(0.0, N_deltap1);
    // vector[N_deltao1] mu_deltao1 = Intercept_deltao1 + rep_vector(0.0, N_deltao1);
    // initialize linear predictor term
    // vector[N_deltao2] mu_deltao2 = Intercept_deltao2 + rep_vector(0.0, N_deltao2);
    // initialize linear predictor term
    for (n in 1:N) {
    // for (n in 1:N_deltao1) {
      // add more terms to the linear predictor
      mu_deltao[n] += [ bsp_deltao1[1], bsp_deltao2[1] ] ;
      Yl_mat[n] = [Yl_deltao1[n], Yl_deltao2[n]];
      // mu_deltao1[n] += (bsp_deltao1[1]) * Yl_deltap1[n];
      // mu_deltao2[n] += (bsp_deltao2[1]) * Yl_deltap1[n];
    // }
    // for (n in 1:N_deltao2) {
      // add more terms to the linear predictor
    // }
    }
    target += multi_normal_cholesky_lpdf(Yl_mat | mu_deltao, L_Sigma);
    // target += normal_lpdf(Yl_deltao1 | mu_deltao1, sigma_deltao1);
    // target += normal_lpdf(Yl_deltao2 | mu_deltao2, sigma_deltao2);
    // target += normal_lpdf(Yl_deltap1 | mu_deltap1, sigma_deltap1);
  }
  // priors including constants
  target += lprior;
  target += normal_lpdf(Y_deltao1[Jme_deltao1] | Yl_deltao1[Jme_deltao1], noise_deltao1[Jme_deltao1]);
  target += normal_lpdf(Y_deltao2[Jme_deltao2] | Yl_deltao2[Jme_deltao2], noise_deltao2[Jme_deltao2]);
  // target += normal_lpdf(Y_deltap1[Jme_deltap1] | Yl_deltap1[Jme_deltap1], noise_deltap1[Jme_deltap1]);
}
generated quantities {
  // actual population-level intercept
  real b_deltao1_Intercept = Intercept_deltao1;
  // actual population-level intercept
  real b_deltao2_Intercept = Intercept_deltao2;
  // actual population-level intercept
  // real b_deltap1_Intercept = Intercept_deltap1;

  real r_corr = multiply_lower_tri_self_transpose(L_corr)[1,2]; // The correlation coefficient
  // Residual values
  // These are basically equal to diff_delta
  // The fact that works lets me make the full model
  // vector[N] resid_1 = Yl_deltao1 - bsp_deltao1[1] * Yl_deltap1;
  // vector[N] resid_2 = Yl_deltao2 - bsp_deltao2[1] * Yl_deltap1;
}

