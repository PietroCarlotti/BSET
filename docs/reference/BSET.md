# Bayesian Surrogate Evaluation Test from Carlotti and Parast (2026)

This function implements the Bayesian Surrogate Evaluation Test (BSET)
as proposed by Carlotti and Parast (2026) . When `X` is `NULL` (the
default), the model is fit without covariates; when `X` is provided, a
multivariate regression model is used to adjust for baseline covariates.
In both cases the function fits a Bayesian model using Stan to generate
posterior samples for the parameters of interest, which are then used to
conduct a Bayesian hypothesis test for evaluating the validity of the
surrogate marker. A frequentist test is also performed for comparison.

## Usage

``` r
BSET(
  data,
  Y,
  S,
  Z,
  X = NULL,
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
  mu_0 = rep(0, 4),
  Sigma_0 = diag(4),
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

  Character vector. Names of the covariate columns to include in the
  model. When `NULL` (the default), the model is fit without covariates
  and the `mu_0` / `Sigma_0` prior parameters are used. When provided,
  the covariate-adjusted model is fit and the `intercept` / `mu_beta` /
  `Sigma_beta` prior parameters are used instead.

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

  Type I error rate for the test (default is 0.05).

- beta:

  Type II error rate for the test (default is 0.2).

- V_S_zero:

  Value of V_S under the null hypothesis (default is 0.5).

- BF_alternative:

  Alternative hypothesis for the Bayes factor test ("greater", "less" or
  "two.sided").

- root_tolerance:

  Numerical tolerance for root-finding algorithms (default is 1e-16).

- mu_0:

  Prior mean vector for the mean parameters. Used only when `X = NULL`
  (default is a vector of zeros of length 4).

- Sigma_0:

  Prior covariance matrix for the mean vector. Used only when `X = NULL`
  (default is a 4x4 identity matrix).

- intercept:

  Logical. Whether to include an intercept in the regression model. Used
  only when `X` is provided (default is `TRUE`).

- mu_beta:

  Prior mean vector for the regression coefficients. Used only when `X`
  is provided (default is a vector of zeros with length equal to the
  number of covariates, including the intercept if `intercept = TRUE`).

- Sigma_beta:

  Prior covariance matrix for the regression coefficients. Used only
  when `X` is provided (default is a diagonal matrix with variance 10,
  with dimensions equal to the number of covariates, including the
  intercept if `intercept = TRUE`).

- s:

  Prior scale parameters for the error variance (default is a vector of
  ones of length 4).

- tau:

  Prior parameter for the LKJ correlation distribution (default is 1).

- plot:

  Logical. Whether to plot the posterior distribution of theta (default
  is FALSE).

- mute:

  Logical. Whether to suppress Stan output during model fitting (default
  is TRUE).

- parallel:

  Logical. Whether to use parallel processing for MCMC sampling (default
  is TRUE).

## Value

A list containing:

- `stan_cores`: The number of CPU cores used for MCMC sampling.

- `Bayesian_test`: A list containing the results of the Bayesian test,
  including posterior samples and credible intervals.

- `frequentist_test`: A list containing the results of the frequentist
  test, including point estimates and confidence intervals.

- `theta_posterior_plot`: A `ggplot` object showing the posterior
  distribution of \\\theta\\, with vertical lines indicating the
  credible interval, the \\\eta\\ threshold, and the true values of
  \\\delta\\ and \\\theta\\ (if provided).

## Details

**Without covariates** (`X = NULL`), the Bayesian model is: \$\$(Y_i,
S_i) = \begin{cases} P\_{1i}, & \text{if } Z_i = 1 \\ P\_{0i}, &
\text{if } Z_i = 0 \end{cases}, \quad i = 1, \ldots, n,\$\$ \$\$P_i
\overset{\text{ind.}}{\sim} \mathcal{N}\_4(\mu, \Sigma), \quad i = 1,
\ldots, n,\$\$ \$\$\mu \sim \mathcal{N}\_4(\mu_0, \Sigma_0),\$\$
\$\$\Sigma = \text{diag}(\sigma\_{1:4}) \\ \Omega \\
\text{diag}(\sigma\_{1:4}),\$\$ \$\$\sigma_k \sim \text{Half-Normal}(0,
s_k), \quad k = 1, \ldots, 4,\$\$ \$\$\Omega \sim \text{LKJ}(\tau).\$\$

**With covariates** (`X` provided), the Bayesian model is: \$\$P_i
\overset{\text{ind.}}{\sim} \mathcal{N}\_4(B X_i, \Sigma), \quad i = 1,
\ldots, n,\$\$ \$\$B = \begin{bmatrix} \beta\_{1} & \beta\_{2} &
\beta\_{3} & \beta\_{4} \\ \end{bmatrix}^\top,\$\$ \$\$\beta\_{k} \sim
\mathcal{N}\_d (\mu\_{\beta}, \Sigma\_{\beta}), \quad k = 1, \ldots,
4,\$\$ with \\\Sigma\\, \\\sigma_k\\, and \\\Omega\\ as above.

This is a primary user-facing function of the package and includes
working examples below.

## References

Carlotti P, Parast L (2026). “A Bayesian Critique of Rank-Based Methods
for Surrogate Marker Evaluation.” *arXiv preprint arXiv:2603.14381*.

Parast L, Cai T, Tian L (2024). “A rank-based approach to evaluate a
surrogate marker in a small sample setting.” *Biometrics*, **80**(1),
ujad035.

## Examples

``` r
# Generate data from the perfect surrogate setting of Parast et al. (2024)
set.seed(123)
data_no_X <- DGP_no_X(
  n = 100,
  p = 0.5,
  mu_star = c(6, 6, 2.5, 2.5),
  Sigma_star = kronecker(diag(2), matrix(c(3, 3, 3, 3.1), 2, 2)),
  model = "Gaussian"
)

# Prepare the data frame
df_no_X <- data.frame(
  Y = data_no_X$P_observed[, "Y"],
  S = data_no_X$P_observed[, "S"],
  Z = data_no_X$Z
)

# Run BSET without covariates (requires Stan compilation, ~1-2 minutes)
if (FALSE) { # \dontrun{
result_no_X <- BSET(
  data = df_no_X,
  Y = "Y",
  S = "S",
  Z = "Z",
  seed = 123,
  n_chains = 2,
  n_iter = 500
)
} # }

# Generate data from the setting of Carlotti and Parast (2026) with a binary covariate
set.seed(123)
data_X <- DGP_X_binary(
  n = 100,
  p = 0.5,
  q = 0.5,
  mu_0 = c(5, 5, 0, 0),
  mu_1 = c(5, -5, 0, -10),
  Sigma_0 = kronecker(diag(2), matrix(c(1, 1, 1, 2), 2, 2)),
  Sigma_1 = kronecker(diag(2), matrix(c(1, 1, 1, 2), 2, 2))
)

# Prepare the data frame
df_X <- data.frame(
  Y = data_X$P_observed[, "Y"],
  S = data_X$P_observed[, "S"],
  Z = data_X$Z,
  X = data_X$X
)

# Run BSET with covariates (requires Stan compilation, ~1-2 minutes)
if (FALSE) { # \dontrun{
result_X <- BSET(
  data = df_X,
  Y = "Y",
  S = "S",
  Z = "Z",
  X = "X",
  seed = 123,
  n_chains = 2,
  n_iter = 500
)
} # }
```
