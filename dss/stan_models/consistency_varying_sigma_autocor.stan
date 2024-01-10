// generated with brms 2.17.0
functions {
  // For a single chrom, create the gp effect
  vector assemble_gp_block(vector z_block, vector sigma_block, data vector positions, 
                     real length_scale, real gp_prop, data real nugget) {
    int N_block = rows(z_block);
    // real fake_sigma = 1;
    // Chokesly of the matern 5/2 gp, with added nugget for stabilization
    matrix[N_block, N_block] gp_cor =  gp_exponential_cov(
        // gp_matern32_cov(
        // gp_matern52_cov(
          to_array_1d(positions), gp_prop, length_scale);
          
    matrix[n_block, N_block] L_gp = cholesky_decompose(add_diag(
        gp_cor, (1 - gp_prop) + nugget ) );
    return diag_pre_multiply(sigma_block, L_gp) * z_block;
  }
  // Create the whole GP
  vector assemble_gp(vector z, vector sigmas, data vector positions, 
                     data array[] int starts, data array[] int lengths,
                     real length_scale, real gp_prop, data real nugget) {
    vector[rows(z)] out;
    for(n in 1:size(starts)) {
      int s = starts[n];
      int l = lengths[n];
      out[s:(s + l - 1)] = assemble_gp_block(
        segment(z, s, l), segment(sigmas, s,l), segment(positions, s, l),
        length_scale, gp_prop, nugget); 
      }
    return out;
    }
  // Attemptint to add reduce_sum
  /* integer sequence of values
  * Args:
    *   start: starting integer
  *   end: ending integer
  * Returns:
    *   an integer sequence from start to end
  */
    int[] sequence(int start, int end) {
      int seq[end - start + 1];
      for (n in 1:num_elements(seq)) {
        seq[n] = n + start - 1;
      }
      return seq;
    }
  // compute partial sums of the log-likelihood
  real partial_log_lik_sdeltao_lpmf(int[] seq_sdeltao, int start, int end, 
      data vector Y_sdeltao, vector Yl_sdeltao, real Intercept_sdeltao, 
      vector Yl_deltap, data vector Csp_sdeltao_1, vector bsp_sdeltao, 
      vector sigma_sdeltao_locus, data int[] J_1_sdeltao, data vector Z_1_sdeltao_1, vector r_1_sdeltao_1) {
    real ptarget = 0;
    int N = end - start + 1;
    int N_sdeltao = end - start + 1;
    // initialize linear predictor term
    vector[N_sdeltao] mu_sdeltao = Intercept_sdeltao + rep_vector(0.0, N_sdeltao);
    for (n in 1:N_sdeltao) {
      // add more terms to the linear predictor
      int nn = n + start - 1;
      mu_sdeltao[n] += (bsp_sdeltao[1]) * Yl_deltap[nn] * Csp_sdeltao_1[nn] + r_1_sdeltao_1[J_1_sdeltao[nn]] * Z_1_sdeltao_1[nn];
    }
    ptarget += normal_lpdf(Yl_sdeltao[start:end] | mu_sdeltao, sigma_sdeltao_locus[J_1_sdeltao[start:end]]);
    return ptarget;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_deltap_lpmf(int[] seq_deltap, int start, int end, data vector Y_deltap, vector Yl_deltap, real Intercept_deltap, real sigma_deltap, data int[] J_1_sdeltao, data vector Z_1_sdeltao_1, vector r_1_sdeltao_1) {
    real ptarget = 0;
    int N = end - start + 1;
    int N_deltap = end - start + 1;
    // initialize linear predictor term
    vector[N_deltap] mu_deltap = Intercept_deltap + rep_vector(0.0, N_deltap);
    ptarget += normal_lpdf(Yl_deltap[start:end] | mu_deltap, sigma_deltap);
    return ptarget;
  }
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
  int<lower=1> Ksp_sdeltao;  // number of special effects terms // smry = cmdstan_fit$summary()
  // covariates of special effects terms
  vector[N_sdeltao] Csp_sdeltao_1;
  int<lower=1> N_deltap;  // number of observations
  vector[N_deltap] Y_deltap;  // response variable
  // data for measurement-error in the response
  vector<lower=0>[N_deltap] noise_deltap;
  // information about non-missings
  int<lower=0> Nme_deltap;
  int<lower=1> Jme_deltap[Nme_deltap];
    int grainsize;  // grainsize for threading
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
  
  // GP parameters
  int N_chrom; // Number of chromosomes / linkage groups
  vector[N_1] positions; // positions in the genome, in units of 10k bp
  int chrom_starts[N_chrom]; // index of where each chrom starts in positions; no padding
  int chrom_n[N_chrom]; // number of elements from each chrom
  real gp_nugget; // small value to control noise regularization on the gp
  // length-scale priors
  real gp_len_priors[2]; // Mean, sd on log scale
}
transformed data {
  int N_nogene = N_1 - N_gene;
  vector[N_nogene] sigma_props_prior_nogene = rep_vector(sigma_concentration, N_nogene); // concentration paramter of 2; can change.
  vector[N_gene] sigma_props_prior_gene = rep_vector(sigma_concentration, N_gene); // concentration paramter of 2; can change.
  int seq_sdeltao[N_sdeltao] = sequence(1, N_sdeltao);
  int seq_deltap[N_deltap] = sequence(1, N_deltap);
}
parameters {
  vector[N_sdeltao] Yl_sdeltao;  // latent variable
  real Intercept_sdeltao;  // temporary intercept for centered predictors
  vector<lower=0>[2] sigma_sdeltao;  // dispersion parameter; 1 = not gene, 2 = gene
  vector[N_deltap] Yl_deltap;  // latent variable
  real Intercept_deltap;  // temporary intercept for centered predictors
  real<lower=0> sigma_deltap;  // dispersion parameter
  vector<lower=0>[2] sd_1;  // group-level standard deviations (1 = nogene, 2 = gene)
  vector[N_1] z_1[1];  // standardized group-level effects
  simplex[N_nogene] sigma_sdeltao_locus_nogene_props; // proportions of total variation for each group.
  simplex[N_gene] sigma_sdeltao_locus_gene_props; // proportions of total variation for each group.
  
  real gene_effect; // change to mu for gbm loci
  real gp_log10_length_scale; // length scale for GP
  real<lower=0, upper = 1> gp_prop; // Proportion of variation in random effect that is correlated; uniform prior
}
transformed parameters {
  vector[Ksp_sdeltao] bsp_sdeltao;  // special effects coefficients
  vector[N_1] r_1_sdeltao_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  vector[N_1] sigma_sdeltao_locus = sigma_sdeltao[locus_assignment]; // .* sqrt(N_1 * sigma_sdeltao_locus_props);
  bsp_sdeltao = rep_vector(1, rows(bsp_sdeltao));
  
// mean per-locus effect
  // r_1_sdeltao_1 = (sd_1[locus_assignment] .* (z_1[1])); 
  r_1_sdeltao_1 = assemble_gp(z_1[1], sd_1[locus_assignment], positions,
                              // starts, lengths, length_scale, nugget
                              chrom_starts, chrom_n, 10^gp_log10_length_scale,gp_prop, gp_nugget);
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
    target += reduce_sum(partial_log_lik_sdeltao_lpmf, seq_sdeltao, grainsize, Y_sdeltao,
            Yl_sdeltao, Intercept_sdeltao, Yl_deltap, Csp_sdeltao_1, bsp_sdeltao, 
            sigma_sdeltao_locus, 
            J_1_sdeltao, Z_1_sdeltao_1, r_1_sdeltao_1);
    target += reduce_sum(partial_log_lik_deltap_lpmf, seq_deltap, grainsize, Y_deltap, Yl_deltap, Intercept_deltap, sigma_deltap, J_1_sdeltao, Z_1_sdeltao_1, r_1_sdeltao_1);
  
    // initialize linear predictor term
    // vector[N_sdeltao] mu_sdeltao = Intercept_sdeltao + rep_vector(0.0, N_sdeltao);
    // initialize linear predictor term
    // vector[N_deltap] mu_deltap = Intercept_deltap + rep_vector(0.0, N_deltap);
    // mu_sdeltao[]
    // for (n in 1:N_sdeltao) {
      // add more terms to the linear predictor
      // mu_sdeltao[n] += (bsp_sdeltao[1]) * Yl_deltap[n] * Csp_sdeltao_1[n] + 
        // r_1_sdeltao_1[J_1_sdeltao[n]] * Z_1_sdeltao_1[n];
    // }
    // target += normal_lpdf(Yl_sdeltao | mu_sdeltao, sigma_sdeltao_locus[J_1_sdeltao]);
    // target += normal_lpdf(Yl_deltap | mu_deltap, sigma_deltap);
  }
  // priors including constants
  target += lprior;
  target += normal_lpdf(Y_sdeltao[Jme_sdeltao] | Yl_sdeltao[Jme_sdeltao], noise_sdeltao[Jme_sdeltao]);
  target += normal_lpdf(Y_deltap[Jme_deltap] | Yl_deltap[Jme_deltap], noise_deltap[Jme_deltap]);
  target += std_normal_lpdf(z_1[1]);
  // GP length scale prior
  target += normal_lpdf(gp_log10_length_scale | gp_len_priors[1], gp_len_priors[2]);
}
generated quantities {
  // actual population-level intercept
  real b_sdeltao_Intercept = Intercept_sdeltao;
  // actual population-level intercept
  real b_deltap_Intercept = Intercept_deltap;
  real diff_sigma_sdeltao = sigma_sdeltao[2] - sigma_sdeltao[1]; 
  // Diff_delta calculations
  vector[N] diff_delta = Yl_sdeltao - (bsp_sdeltao[1] * Yl_deltap .* Csp_sdeltao_1);
  vector[N] true_resid = diff_delta -  r_1_sdeltao_1[J_1_sdeltao] ;
  vector[N_1] locus_mean = Intercept_sdeltao + r_1_sdeltao_1;
  // Will need to do some post-processing to calculate the residual variance per locus
  // real sd_diff_delta = sd(diff_delta);
  // real sd_resid = sqrt(sd_diff_delta^2 - sd_1[1]^2);
  
  

}

