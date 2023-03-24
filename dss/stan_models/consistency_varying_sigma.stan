// generated with brms 2.17.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_sdeltao;  // number of observations
  vector[N_sdeltao] Y_sdeltao;  // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_sdeltao] noise_sdeltao;
  // information about non-missings
  int<lower=0> Nme_sdeltao;
  int<lower=1> Jme_sdeltao[Nme_sdeltao];
  int<lower=1> Ksp_sdeltao;  // number of special effects terms
  // covariates of special effects terms
  vector[N_sdeltao] Csp_sdeltao_1;
  int<lower=1> N_deltap;  // number of observations
  vector[N_deltap] Y_deltap;  // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_deltap] noise_deltap;
  // information about non-missings
  int<lower=0> Nme_deltap;
  int<lower=1> Jme_deltap[Nme_deltap];
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels (Loci)
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1_sdeltao[N_sdeltao];  // grouping indicator per observation
  // group-level predictor values
  vector[N_sdeltao] Z_1_sdeltao_1;
  int prior_only;  // should the likelihood be ignored?
  // Manual additions
  int<lower=1,upper=N_1> N_gene; // Number of loci that are on the gene
  int<lower=1,upper=N_1> which_nogene[N_1-N_gene]; // index which loci are NOT gene body methylation
  int<lower=1,upper=N_1> which_gene[N_gene]; // index which loci are gene body methylation
  int<lower=1,upper=2> locus_assignment[N_1]; // 1 = not gene, 2 = gene
  real<lower=1> sigma_concentration;
}
transformed data {
  int N_nogene = N_1 - N_gene;
  vector[N_nogene] sigma_props_prior_nogene = rep_vector(sigma_concentration, N_nogene); // concentration paramter of 2; can change.
  vector[N_gene] sigma_props_prior_gene = rep_vector(sigma_concentration, N_gene); // concentration paramter of 2; can change.
}
parameters {
  vector[N_sdeltao] Yl_sdeltao;  // latent variable
  real Intercept_sdeltao;  // temporary intercept for centered predictors
  vector<lower=0>[2] sigma_sdeltao;  // dispersion parameter; 1 = not gene, 2 = gene
  vector[N_deltap] Yl_deltap;  // latent variable
  real Intercept_deltap;  // temporary intercept for centered predictors
  real<lower=0> sigma_deltap;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  simplex[N_nogene] sigma_sdeltao_locus_nogene_props; // proportions of total variation for each group.
  simplex[N_gene] sigma_sdeltao_locus_gene_props; // proportions of total variation for each group.
  
  real gene_effect; // change to mu for gbm loci
}
transformed parameters {
  vector[Ksp_sdeltao] bsp_sdeltao;  // special effects coefficients
  vector[N_1] r_1_sdeltao_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  vector[N_1] sigma_sdeltao_locus = sigma_sdeltao[locus_assignment]; #.* sqrt(N_1 * sigma_sdeltao_locus_props);
  bsp_sdeltao = rep_vector(1, rows(bsp_sdeltao));
  r_1_sdeltao_1 = (sd_1[1] * (z_1[1])); // mean per-locus effect
  
  // Assemble sigma_sdeltao_locus
  sigma_sdeltao_locus[which_nogene] = sigma_sdeltao_locus[which_nogene] .*
    sqrt((N_nogene) * sigma_sdeltao_locus_nogene_props);
  sigma_sdeltao_locus[which_gene] = sigma_sdeltao_locus[which_gene] .* 
    sqrt(N_gene * sigma_sdeltao_locus_gene_props);
  r_1_sdeltao_1[which_gene] += gene_effect;
  
  lprior += normal_lpdf(Intercept_sdeltao | 0, 1);
  lprior += normal_lpdf(gene_effect | 0, 0.5); // Prior N(0, 0.5) for gene effect
  lprior += student_t_lpdf(sigma_sdeltao | 7, 0, .5)
    - 1 * student_t_lccdf(0 | 7, 0, .5);
  lprior += normal_lpdf(Intercept_deltap | 0, 1);
  lprior += student_t_lpdf(sigma_deltap | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  lprior += student_t_lpdf(sd_1 | 7, 0, .5)
    - 1 * student_t_lccdf(0 | 7, 0, .5);
  lprior += dirichlet_lpdf(sigma_sdeltao_locus_gene_props | sigma_props_prior_gene);
  lprior += dirichlet_lpdf(sigma_sdeltao_locus_nogene_props | sigma_props_prior_nogene);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_sdeltao] mu_sdeltao = Intercept_sdeltao + rep_vector(0.0, N_sdeltao);
    // initialize linear predictor term
    vector[N_deltap] mu_deltap = Intercept_deltap + rep_vector(0.0, N_deltap);
    // mu_sdeltao[]
    for (n in 1:N_sdeltao) {
      // add more terms to the linear predictor
      mu_sdeltao[n] += (bsp_sdeltao[1]) * Yl_deltap[n] * Csp_sdeltao_1[n] + 
        r_1_sdeltao_1[J_1_sdeltao[n]] * Z_1_sdeltao_1[n];
    }
    target += normal_lpdf(Yl_sdeltao | mu_sdeltao, sigma_sdeltao_locus[J_1_sdeltao]);
    target += normal_lpdf(Yl_deltap | mu_deltap, sigma_deltap);
  }
  // priors including constants
  target += lprior;
  target += normal_lpdf(Y_sdeltao[Jme_sdeltao] | Yl_sdeltao[Jme_sdeltao], noise_sdeltao[Jme_sdeltao]);
  target += normal_lpdf(Y_deltap[Jme_deltap] | Yl_deltap[Jme_deltap], noise_deltap[Jme_deltap]);
  target += std_normal_lpdf(z_1[1]);
}
generated quantities {
  // actual population-level intercept
  real b_sdeltao_Intercept = Intercept_sdeltao;
  // actual population-level intercept
  real b_deltap_Intercept = Intercept_deltap;
  real diff_sigma_sdeltao = sigma_sdeltao[2] - sigma_sdeltao[1]; 
  // Diff_delta calculations
  vector[N] diff_delta = Yl_sdeltao - Intercept_sdeltao - (bsp_sdeltao[1] * Yl_deltap .* Csp_sdeltao_1);
  vector[N] true_resid = diff_delta -  r_1_sdeltao_1[J_1_sdeltao] ;
  // Will need to do some post-processing to calculate the residual variance per locus
  real sd_diff_delta = sd(diff_delta);
  real sd_resid = sqrt(sd_diff_delta^2 - sd_1[1]^2);
  
  

}

