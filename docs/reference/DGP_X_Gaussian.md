# Data Generating Process with a Gaussian Covariate

This function generates potential outcomes from a data generating
process similar to the one described in Parast et al. (2024) , but with
the addition of a Gaussian covariate X. It creates a dataset of
potential outcomes \$\$P = (Y_1, S_1, Y_0, S_0)\$\$ and observed
outcomes \$\$P\_{observed} = (Y, S)\$\$ based on a random treatment
assignment \\Z\\.

## Usage

``` r
DGP_X_Gaussian(n, p, beta, Sigma, m, s)
```

## Arguments

- n:

  Integer. Total sample size.

- p:

  Numeric. Probability of being assigned to the treatment group
  \\(Z=1)\\.

- beta:

  Numeric vector. Coefficients for the linear function of the covariate
  \\X\\.

- Sigma:

  Matrix. Covariance matrix for the potential outcomes.

- m:

  Numeric. Mean of the Gaussian covariate \\X\\.

- s:

  Numeric. Standard deviation of the Gaussian covariate \\X\\.

## Value

A list containing:

- `X`: The Gaussian covariate vector.

- `Z`: Treatment assignment vector.

- `n1`: Number of treated units.

- `n0`: Number of control units.

- `P`: Full matrix of potential outcomes.

- `P_observed`: Observed outcomes \\(Y, S)\\ corresponding to the
  assigned treatment \\Z\\.

- `P_unobserved`: Counterfactual outcomes under the opposite treatment.

This function is useful for generating synthetic data to test or explore
the method, for instance to verify the behavior of `BSET_X` under known
simulation settings.

## Details

The potential outcomes are generated from a multivariate normal
distribution with mean vector and covariance matrix that depend on the
value of \\X\\. Specifically, the mean vector is a linear function of
\\X\\: \$\$\mu(X) = x \cdot (\beta\_{Y1}, \beta\_{S1}, \beta\_{Y0},
\beta\_{S0})^T,\$\$ and the covariance matrix is constant across values
of \\X\\: \$\$\Sigma(X) = \Sigma.\$\$

## References

Carlotti P, Parast L (2026). “A Bayesian Critique of Rank-Based Methods
for Surrogate Marker Evaluation.” *arXiv preprint arXiv:2603.14381*.

## Examples

``` r
set.seed(123)
data <- DGP_X_Gaussian(
  n = 100,
  p = 0.5,
  beta = c(1, 7, 0, 6),
  Sigma = 0.5 * diag(4),
  m = 3,
  s = 1
)
```
