# Bayesian Surrogate Evaluation Test with Covariates from Carlotti and Parast (2026)

This function implements the Bayesian Surrogate Evaluation Test (BSET)
which includes covariates via a multivariate regression model for the
potential, as proposed by Carlotti and Parast (2026). The function fits
a Bayesian model using Stan to generate posterior samples for the
parameters of interest, including the regression coefficients for the
covariates. These posterior samples are then used to conduct a Bayesian
hypothesis test for evaluating the validity of the surrogate marker,
while adjusting for covariates. A frequentist test is also performed for
comparison. The Bayesian model is specified as follows: \$\$(Y_i, S_i) =
\begin{cases} P\_{1i}, & \text{if } Z_i = 1 \\ P\_{0i}, & \text{if } Z_i
= 0 \end{cases}, \quad i = 1, \ldots, n,\$\$ \$\$P_i
\overset{\text{ind.}}{\sim} \mathcal{N}\_4(B X_i, \Sigma), \quad i = 1,
\ldots, n,\$\$ \$\$B = \begin{bmatrix} \beta\_{1} & \beta\_{2} &
\beta\_{3} & \beta\_{4} \\ \end{bmatrix}^\top,\$\$ \$\$\beta\_{k} \sim
\mathcal{N}\_d (\mu\_{\beta}, \Sigma\_{\beta}), \quad k = 1, \ldots,
4,\$\$ \$\$\Sigma = \text{diag}(\sigma\_{1:4}) \\ \Omega \\
\text{diag}(\sigma\_{1:4}),\$\$ \$\$\sigma_k \sim \text{Half-Normal}(0,
s_k), \quad k = 1, \ldots, 4,\$\$ \$\$\Omega \sim \text{LKJ}(\tau).\$\$

## Usage

``` r
BSET_X(
  data,
  Y,
  S,
  Z,
  X,
  delta_true = NULL,
  theta_true = NULL,
  seed = NULL,
  n_chains = 4,
  n_iter = 2000,
  burn_in_ratio = 0.25,
  a = 1,
  b = 1,
  alpha = 0.05,
  beta = 0.2,
  V_S_zero = 0.5,
  BF_alternative = "greater",
  root_tolerance = 1e-16,
  intercept = TRUE,
  mu_beta = NULL,
  Sigma_beta = NULL,
  s = rep(1, 4),
  tau = 1,
  plot = FALSE,
  mute = TRUE,
  parallel = TRUE
)
```

## Arguments

- data:

  A data frame containing the observed data.

- Y:

  Character. Name of the outcome variable.

- S:

  Character. Name of the surrogate variable.

- Z:

  Character. Name of the treatment assignment variable.

- X:

  Character vector. Names of the covariates to include in the model.

- delta_true:

  The true value of delta, used to calculate frequentist coverage during
  simulations (optional).

- theta_true:

  The true value of theta, used to calculate Bayesian coverage during
  simulations (optional).

- seed:

  Random seed for reproducibility (optional. If not provided, a random
  seed will be generated).

- n_chains:

  Number of MCMC chains to run (default is 4).

- n_iter:

  Number of iterations per MCMC chain (default is 2000).

- burn_in_ratio:

  Proportion of iterations to discard as burn-in (default is 0.25).

- a, b:

  Prior parameters for prior Beta distribution on V_S (default is a = 1,
  b = 1).

- alpha:

  Type I error rate for the test (default is 0.05)

- beta:

  Type II error rate for the test (default is 0.2)

- V_S_zero:

  Value of V_S under the null hypothesis (default is 0.5)

- BF_alternative:

  Alternative hypothesis for the Bayes factor test ("greater", "less" or
  "two.sided").

- root_tolerance:

  Numerical tolerance for root-finding algorithms (default is 1e-16).

- intercept:

  Whether to include an intercept in the regression model (default is
  TRUE).

- mu_beta:

  Prior mean vector for regression coefficients (default is a vector of
  zeros with length equal to the number of covariates, including the
  intercept if specified).

- Sigma_beta:

  Prior covariance matrix for regression coefficients (default is a
  diagonal matrix with large variances, with dimensions equal to the
  number of covariates, including the intercept if specified).

- s:

  Prior scale parameters for the error variance (default is a vector of
  ones of length 4)

- tau:

  Prior parameter for the LKJ correlation distribution (default is 1)

- plot:

  Logical. Whether to plot the posterior distribution of theta (default
  is FALSE)

- mute:

  Logical. Whether to suppress Stan output during model fitting (default
  is TRUE)

- parallel:

  Logical. Whether to use parallel processing for MCMC sampling (default
  is TRUE)
