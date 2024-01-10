// generated with brms 2.17.0
// Add hierarchical variance to sigma

functions {
  /* integer sequence of values
   * Args:
   *   start: starting integer
   *   end: ending integer
   * Returns:
   *   an integer sequence from start to end
   */
  array[] int sequence(int start, int end) {
    array[end - start + 1] int seq;
    for (n in 1 : num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_deltao_lpmf(array[] int seq_deltao, int start,
                                   int end, data vector Y_deltao,
                                   vector Yl_deltao, vector sigma1,
                                   data matrix Xc_mu2_deltao,
                                   vector b_mu2_deltao,
                                   real Intercept_mu2_deltao,
                                   vector Yl_deltap,
                                   data vector Csp_mu2_deltao_1,
                                   data vector Csp_mu2_deltao_2,
                                   data vector Csp_mu2_deltao_3,
                                   data vector Csp_mu2_deltao_4,
                                   data vector Csp_mu2_deltao_5,
                                   data vector Csp_mu2_deltao_6,
                                   data vector Csp_mu2_deltao_7,
                                   vector bsp_mu2_deltao, vector sigma2,
                                   data matrix Xc_theta1_deltao,
                                   vector b_theta1_deltao,
                                   real Intercept_theta1_deltao,
                                   array[] int cross_id) {
    real ptarget = 0;
    int N = end - start + 1;
    int N_deltao = end - start + 1;
    // initialize linear predictor term
    vector[N_deltao] mu1_deltao = rep_vector(0.0, N_deltao);
    // initialize linear predictor term
    vector[N_deltao] mu2_deltao = Intercept_mu2_deltao
                                  + Xc_mu2_deltao[start : end] * b_mu2_deltao;
    // initialize linear predictor term
    vector[N_deltao] theta1_deltao = Intercept_theta1_deltao
                                     + Xc_theta1_deltao[start : end]
                                       * b_theta1_deltao;
    vector[N_deltao] theta2_deltao = rep_vector(0.0, N_deltao);
    real log_sum_exp_theta;
    for (n in 1 : N_deltao) {
      // add more terms to the linear predictor
      int nn = n + start - 1;
      mu2_deltao[n] += bsp_mu2_deltao[1] * Yl_deltap[nn]
                       + bsp_mu2_deltao[2] * Yl_deltap[nn]
                         * Csp_mu2_deltao_1[nn]
                       + bsp_mu2_deltao[3] * Yl_deltap[nn]
                         * Csp_mu2_deltao_2[nn]
                       + bsp_mu2_deltao[4] * Yl_deltap[nn]
                         * Csp_mu2_deltao_3[nn]
                       + bsp_mu2_deltao[5] * Yl_deltap[nn]
                         * Csp_mu2_deltao_4[nn]
                       + bsp_mu2_deltao[6] * Yl_deltap[nn]
                         * Csp_mu2_deltao_5[nn]
                       + bsp_mu2_deltao[7] * Yl_deltap[nn]
                         * Csp_mu2_deltao_6[nn]
                       + bsp_mu2_deltao[8] * Yl_deltap[nn]
                         * Csp_mu2_deltao_7[nn];
    }
    for (n in 1 : N_deltao) {
      // scale theta to become a probability vector
      log_sum_exp_theta = log(exp(theta1_deltao[n]) + exp(theta2_deltao[n]));
      theta1_deltao[n] = theta1_deltao[n] - log_sum_exp_theta;
      theta2_deltao[n] = theta2_deltao[n] - log_sum_exp_theta;
    }
    // likelihood of the mixture model
    for (n in 1 : N_deltao) {
      int nn = n + start - 1;
      array[2] real ps;
      ps[1] = theta1_deltao[n]
              + normal_lpdf(Yl_deltao[nn] | mu1_deltao[n], sigma1[cross_id[n]]);
      ps[2] = theta2_deltao[n]
              + normal_lpdf(Yl_deltao[nn] | mu2_deltao[n], sigma2[cross_id[n]]);
      ptarget += log_sum_exp(ps);
    }
    return ptarget;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_deltap_lpmf(array[] int seq_deltap, int start,
                                   int end, data vector Y_deltap,
                                   vector Yl_deltap, real Intercept_deltap,
                                   real sigma_deltap) {
    real ptarget = 0;
    int N = end - start + 1;
    int N_deltap = end - start + 1;
    // initialize linear predictor term
    vector[N_deltap] mu_deltap = Intercept_deltap + rep_vector(0.0, N_deltap);
    ptarget += normal_lpdf(Yl_deltap[start : end] | mu_deltap, sigma_deltap);
    return ptarget;
  }
}
data {
  int<lower=1> N; // total number of observations
  int<lower=1> N_deltao; // number of observations
  vector[N_deltao] Y_deltao; // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_deltao] noise_deltao;
  // information about non-missings
  int<lower=0> Nme_deltao;
  array[Nme_deltao] int<lower=1> Jme_deltao;
  int<lower=1> K_mu2_deltao; // number of population-level effects
  matrix[N_deltao, K_mu2_deltao] X_mu2_deltao; // population-level design matrix
  int<lower=1> Ksp_mu2_deltao; // number of special effects terms
  // covariates of special effects terms
  vector[N_deltao] Csp_mu2_deltao_1;
  vector[N_deltao] Csp_mu2_deltao_2;
  vector[N_deltao] Csp_mu2_deltao_3;
  vector[N_deltao] Csp_mu2_deltao_4;
  vector[N_deltao] Csp_mu2_deltao_5;
  vector[N_deltao] Csp_mu2_deltao_6;
  vector[N_deltao] Csp_mu2_deltao_7;
  int<lower=1> K_theta1_deltao; // number of population-level effects
  matrix[N_deltao, K_theta1_deltao] X_theta1_deltao; // population-level design matrix
  int<lower=1> N_deltap; // number of observations
  vector[N_deltap] Y_deltap; // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_deltap] noise_deltap;
  // information about non-missings
  int<lower=0> Nme_deltap;
  array[Nme_deltap] int<lower=1> Jme_deltap;
  int grainsize; // grainsize for threading
  int prior_only; // should the likelihood be ignored?
  
  real hv_conc_prior;
  int N_cross;
  array[N_deltao] int<lower=1,upper=N_cross> cross_id; // redundant w/ X, but needed for hv prior
}
transformed data {
  int Kc_mu2_deltao = K_mu2_deltao - 1;
  matrix[N_deltao, Kc_mu2_deltao] Xc_mu2_deltao; // centered version of X_mu2_deltao without an intercept
  vector[Kc_mu2_deltao] means_X_mu2_deltao; // column means of X_mu2_deltao before centering
  int Kc_theta1_deltao = K_theta1_deltao - 1;
  matrix[N_deltao, Kc_theta1_deltao] Xc_theta1_deltao; // centered version of X_theta1_deltao without an intercept
  vector[Kc_theta1_deltao] means_X_theta1_deltao; // column means of X_theta1_deltao before centering
  array[N_deltao] int seq_deltao = sequence(1, N_deltao);
  array[N_deltap] int seq_deltap = sequence(1, N_deltap);
  
  for (i in 2 : K_mu2_deltao) {
    means_X_mu2_deltao[i - 1] = mean(X_mu2_deltao[ : , i]);
    Xc_mu2_deltao[ : , i - 1] = X_mu2_deltao[ : , i]
                                - means_X_mu2_deltao[i - 1];
  }
  for (i in 2 : K_theta1_deltao) {
    means_X_theta1_deltao[i - 1] = mean(X_theta1_deltao[ : , i]);
    Xc_theta1_deltao[ : , i - 1] = X_theta1_deltao[ : , i]
                                   - means_X_theta1_deltao[i - 1];
  }
}
parameters {
  vector[N_deltao] Yl_deltao; // latent variable
  real<lower=0> sigma1_deltao; // dispersion parameter
  vector[Kc_mu2_deltao] b_mu2_deltao; // population-level effects
  real Intercept_mu2_deltao; // temporary intercept for centered predictors
  vector[Ksp_mu2_deltao] bsp_mu2_deltao; // special effects coefficients
  real<lower=0> sigma2_deltao; // dispersion parameter
  // For HV prior:
  simplex[N_cross] sigma1_props;
  simplex[N_cross] sigma2_props;
  
  vector[Kc_theta1_deltao] b_theta1_deltao; // population-level effects
  real Intercept_theta1_deltao; // temporary intercept for centered predictors
  vector[N_deltap] Yl_deltap; // latent variable
  real Intercept_deltap; // temporary intercept for centered predictors
  real<lower=0> sigma_deltap; // dispersion parameter
}
transformed parameters {
  vector[N_cross] sigma1 = sigma1_deltao * sqrt(N_cross * sigma1_props);
  vector[N_cross] sigma2 = sigma2_deltao * sqrt(N_cross * sigma2_props);
  
  real lprior = 0; // prior contributions to the log posterior
  lprior += student_t_lpdf(sigma1_deltao | 7, 0, 2)
            - 1 * student_t_lccdf(0 | 7, 0, 2);
  lprior += normal_lpdf(b_mu2_deltao | 0, 1);
  lprior += normal_lpdf(Intercept_mu2_deltao | 0, 1);
  lprior += normal_lpdf(bsp_mu2_deltao | 0, 1);
  lprior += student_t_lpdf(sigma2_deltao | 7, 0, 0.5)
            - 1 * student_t_lccdf(0 | 7, 0, 0.5);
  lprior += normal_lpdf(b_theta1_deltao | 0, 1);
  lprior += normal_lpdf(Intercept_theta1_deltao | 0, 1);
  lprior += normal_lpdf(Intercept_deltap | 0, 2);
  lprior += student_t_lpdf(sigma_deltap | 7, 0, 2)
            - 1 * student_t_lccdf(0 | 7, 0, 2);
            
  lprior += dirichlet_lpdf(sigma1_props | rep_vector(hv_conc_prior, N_cross)) +
            dirichlet_lpdf(sigma2_props | rep_vector(hv_conc_prior, N_cross));
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += reduce_sum(partial_log_lik_deltao_lpmf, seq_deltao, grainsize,
                         Y_deltao, Yl_deltao, sigma1, Xc_mu2_deltao,
                         b_mu2_deltao, Intercept_mu2_deltao, Yl_deltap,
                         Csp_mu2_deltao_1, Csp_mu2_deltao_2,
                         Csp_mu2_deltao_3, Csp_mu2_deltao_4,
                         Csp_mu2_deltao_5, Csp_mu2_deltao_6,
                         Csp_mu2_deltao_7, bsp_mu2_deltao, sigma2,
                         Xc_theta1_deltao, b_theta1_deltao,
                         Intercept_theta1_deltao, cross_id);
    target += reduce_sum(partial_log_lik_deltap_lpmf, seq_deltap, grainsize,
                         Y_deltap, Yl_deltap, Intercept_deltap, sigma_deltap);
  }
  // priors including constants
  target += lprior;
  target += normal_lpdf(Y_deltao[Jme_deltao] | Yl_deltao[Jme_deltao], noise_deltao[Jme_deltao]);
  target += normal_lpdf(Y_deltap[Jme_deltap] | Yl_deltap[Jme_deltap], noise_deltap[Jme_deltap]);
}
generated quantities {
  // actual population-level intercept
  real b_mu2_deltao_Intercept = Intercept_mu2_deltao
                                - dot_product(means_X_mu2_deltao,
                                              b_mu2_deltao);
  // actual population-level intercept
  real b_theta1_deltao_Intercept = Intercept_theta1_deltao
                                   - dot_product(means_X_theta1_deltao,
                                                 b_theta1_deltao);
  // actual population-level intercept
  real b_deltap_Intercept = Intercept_deltap;
}


