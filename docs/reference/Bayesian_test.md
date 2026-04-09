# Bayesian Test from Carlotti and Parast (2026)

This function performs a Bayesian test for surrogate evaluation based on
the imputation-based methodology of Carlotti and Parast (2026). It
calculates credible intervals, the \\\eta\\ threshold used, and
determines if the surrogate is valid.

## Usage

``` r
Bayesian_test(
  P_MCMC,
  alpha = 0.05,
  V_S_star,
  theta_true = NULL,
  V_Y_true = NULL
)
```

## Arguments

- P_MCMC:

  A three-dimensional array of potential outcomes obtained via MCMC
  sampling with dimensions `[n_subjects, variables, n_samples]`. The
  variables should correspond to \\(Y_1, S_1, Y_0, S_0)\\.

- alpha:

  Numeric. Significance level for the credible interval (default is
  0.05).

- V_S_star:

  Numeric. The value of \\v_S\\ that satisfies the power constraint for
  the surrogate validation test.

- theta_true:

  Numeric (optional). The true value of \\\eta\\, used to calculate
  frequentist coverage during simulations.

- V_Y_true:

  Numeric (optional). The true value of \\V_Y\\, used to calculate the
  surrogate validation threshold \\\eta\\ during simulations. If not
  provided, it will be estimated from the MCMC samples.

## Value

A list containing:

- `V_Y_MCMC`: Posterior draws for \\V_Y\\.

- `V_S_MCMC`: Posterior draws for \\V_S\\.

- `theta_MCMC`: Posterior draws for \\\theta = V_Y - V_S\\.

- `CI`: The calculated credible interval for \\\theta\\.

- `eta`: The \\\eta\\ threshold value used in the test.

- `coverage`: Logical indicating if `theta_true` falls within the `CI`
  (if `theta_true` is provided).

- `power`: Logical indicating if the upper bound of `CI` is below `eta`,
  which indicates that the test identifies the surrogate as valid.
