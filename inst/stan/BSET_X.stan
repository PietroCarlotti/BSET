// BAYESIAN IMPUTATION MODEL OF POTENTIAL OUTCOMES WITH COVARIATES
data {
  // Data
  int<lower=1> n;                   // Sample size
  int<lower=1> d;                   // Number of covariates (including the intercept)
  matrix[n, d] X;                   // Covariate matrix
  vector[2] P_observed[n];          // Observed potential outcomes of (Y, S)
  int<lower=0,upper=1> Z[n];        // Treatment: 1 = Treatment, 0 = Control

  // Prior parameters
  vector[d] mu_beta;                // Prior mean for coefficients
  cov_matrix[d] Sigma_beta;         // Prior covariance for coefficients
  vector<lower=0>[4] s;             // Scales for Half-Normal on sigma
  real<lower=0> tau;                // LKJ shape parameter
}

parameters {
  // Global parameters
  matrix[4, d] z_B;                 // Helper for non-centered B
  cholesky_factor_corr[4] L_R;      // Cholesky factor of correlation matrix R
  vector<lower=0>[4] sigma;         // Standard deviations

  // Latent variables
  vector[2] P_unobserved[n];        // Unobserved potential outcomes of (Y, S)
}

transformed parameters {
  matrix[4, d] B;                   // Regression coefficients matrix
  matrix[4, 4] L_Sigma;             // Cholesky factor of Sigma

  // Non-centered parametrization of B
  // Each row beta_j ~ MultiNormal(mu_beta, Sigma_beta)
  {
    matrix[d, d] L_beta = cholesky_decompose(Sigma_beta);
    for (j in 1:4) {
      B[j] = (mu_beta + L_beta * (z_B[j]'))'; 
    }
  }

  // Cholesky factor of Sigma
  L_Sigma = diag_pre_multiply(sigma, L_R);
}

model {
  // Priors
  to_vector(z_B) ~ std_normal();    // Equivalent to B rows ~ Normal(mu_beta, Sigma_beta)
  sigma ~ normal(0, s);             // Half-Normal
  L_R ~ lkj_corr_cholesky(tau);

  // Likelihood
  for (i in 1:n) {
    vector[4] mu_i = B * (X[i]');   // Mean vector for individual i: B * Xi
    
    if (Z[i] == 1) {
      // Treatment: Observed (Y1, S1), Imputed (Y0, S0)
      target += multi_normal_cholesky_lpdf(
                  append_row(P_observed[i], P_unobserved[i]) | mu_i, L_Sigma
                );
    } else {
      // Control: Imputed (Y1, S1), Observed (Y0, S0)
      target += multi_normal_cholesky_lpdf(
                  append_row(P_unobserved[i], P_observed[i]) | mu_i, L_Sigma
                );
    }
  }
}

generated quantities {
  corr_matrix[4] R = multiply_lower_tri_self_transpose(L_R);
  cov_matrix[4] Sigma = quad_form_diag(R, sigma);
  
  // Full matrix of potential outcomes
  vector[4] P[n]; 

  for (i in 1:n) {
    if (Z[i] == 1) {
      P[i][1:2] = P_observed[i];
      P[i][3:4] = P_unobserved[i];
    } else {
      P[i][1:2] = P_unobserved[i];
      P[i][3:4] = P_observed[i];
    }
  }
}
