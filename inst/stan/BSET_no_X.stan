// BAYESIAN IMPUTATION MODEL OF POTENTIAL OUTCOMES WITHOUT COVARIATES
data {
  // Data
  int<lower=1> n;                       // Sample size
  vector[2] P_observed[n];              // Observed potential outcomes of (Y, S)
  int<lower=0,upper=1> Z[n];            // Treatment: 1 = Treatment, 0 = Control

  // Prior parameters
  vector[4] mu_0;                       // Prior mean for mu
  cov_matrix[4] Sigma_0;                // Prior covariance for mu
  vector<lower=0>[4] s;                 // Scales for Half-Normal on sigma
  real<lower=0> tau;                    // LKJ shape parameter
}

parameters {
  // Global parameters
  vector[4] z_mu;                       // Helper for non-centered mu
  cholesky_factor_corr[4] L_R;          // Cholesky factor of correlation matrix R
  vector<lower=0>[4] sigma;             // Standard deviations

  // Latent variables
  vector[2] P_unobserved[n];            // Unobserved potential outcomes of (Y, S)
}

transformed parameters {
  vector[4] mu;
  matrix[4, 4] L_Sigma;                 // Cholesky factor of Sigma

  // Non-centered parametrization of mu
  {
    matrix[4, 4] L0 = cholesky_decompose(Sigma_0);
    mu = mu_0 + L0 * z_mu;
  }

  // Cholesky factor of Sigma
  L_Sigma = diag_pre_multiply(sigma, L_R);
}

model {
  // Priors
  z_mu ~ std_normal();                  // Equivalent to mu ~ Normal(mu_0, Sigma_0)
  sigma ~ normal(0, s);                 // Half-Normal (enforced by <lower=0>)
  L_R ~ lkj_corr_cholesky(tau);

  // Likelihood
  for (i in 1:n) {
    vector[4] P_full;
    
    if (Z[i] == 1) {
      // Treatment: Observed (Y1, S1), Imputed (Y0, S0)
      P_full[1:2] = P_observed[i];
      P_full[3:4] = P_unobserved[i];
    } else {
      // Control: Imputed (Y1, S1), Observed (Y0, S0)
      P_full[1:2] = P_unobserved[i];
      P_full[3:4] = P_observed[i];
    }

    // Sampling the full vector P_full
    P_full ~ multi_normal_cholesky(mu, L_Sigma);
  }
}

generated quantities {
  // Correlation matrix
  corr_matrix[4] R = multiply_lower_tri_self_transpose(L_R);
  
  // Covariance matrix
  cov_matrix[4] Sigma = quad_form_diag(R, sigma);
  
  // Full matrix of potential outcomes
  vector[4] P[n]; 

  for (i in 1:n) {
    if (Z[i] == 1) {
     // Treatment: Observed (Y1, S1), Imputed (Y0, S0)
      P[i][1:2] = P_observed[i];
      P[i][3:4] = P_unobserved[i];
    } else {
      // Control: Imputed (Y1, S1), Observed (Y0, S0)
      P[i][1:2] = P_unobserved[i];
      P[i][3:4] = P_observed[i];
    }
  }
}
