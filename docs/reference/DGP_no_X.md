# Data Generating Process without Baseline Covariates

This function generates potential outcomes from the simulation settings
described in Parast et al. (2024). It creates a dataset of potential
outcomes \$\$P = (Y_1, S_1, Y_0, S_0)\$\$ and observed outcomes
\$\$P\_{observed} = (Y, S)\$\$ based on a random treatment assignment
\\Z\\.

## Usage

``` r
DGP_no_X(
  n,
  p,
  mu_star = NULL,
  Sigma_star = NULL,
  model = c("Gaussian", "misspecified")
)
```

## Arguments

- n:

  Integer. Total sample size.

- p:

  Numeric. Probability of being assigned to the treatment group (Z=1).

- mu_star:

  Numeric vector. The mean vector for \\P\\. Required if
  `model = "Gaussian"`.

- Sigma_star:

  Matrix. The covariance matrix for \\P\\. Required if
  `model = "Gaussian"`.

- model:

  Character. The type of data generation: `"Gaussian"` or
  `"misspecified"`.

## Value

A list containing:

- `Z`: Treatment assignment vector.

- `n1`: Number of treated units.

- `n0`: Number of control units.

- `P`: Full matrix of potential outcomes.

- `P_observed`: Observed outcomes \\(Y, S)\\ corresponding to the
  assigned treatment \\Z\\.

- `P_unobserved`: Counterfactual outcomes under the opposite treatment.

## Details

The function supports two types of data-generating processes:

- **Gaussian model:** Potential outcomes are drawn from a multivariate
  normal distribution: \$\$P \sim \mathcal{N}\_{4}(\mu^{\*},
  \Sigma^{\*}).\$\$

- **Non-linear model:** Potential outcomes for the surrogate are
  generated from a non-Gaussian distribution, and the potential outcomes
  for the primary outcome are generated from a non-linear function of
  the surrogate plus non-Gaussian noise.
