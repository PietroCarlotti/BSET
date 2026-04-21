# Data Generating Process with a Binary Covariate

This function generates potential outcomes from a data generating
process similar to the one described in Parast et al. (2024) , but with
the addition of a binary covariate X. It creates a dataset of potential
outcomes \$\$P = (Y_1, S_1, Y_0, S_0)\$\$ and observed outcomes
\$\$P\_{observed} = (Y, S)\$\$ based on a random treatment assignment
\\Z\\.

## Usage

``` r
DGP_X_binary(n, p, q, mu_0, mu_1, Sigma_0, Sigma_1)
```

## Arguments

- n:

  Integer. Total sample size.

- p:

  Numeric. Probability of being assigned to the treatment group
  \\(Z=1)\\.

- q:

  Numeric. Probability of the binary covariate X being 1.

- mu_0:

  Numeric vector. Mean vector for \\P\\ when \\X=0\\.

- mu_1:

  Numeric vector. Mean vector for \\P\\ when \\X=1\\.

- Sigma_0:

  Matrix. Covariance matrix for \\P\\ when \\X=0\\.

- Sigma_1:

  Matrix. Covariance matrix for \\P\\ when \\X=1\\.

## Value

A list containing:

- `X`: The binary covariate vector.

- `n_X1`: Number of units with \\X=1\\.

- `n_X0`: Number of units with \\X=0\\.

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

The potential outcomes are generated from multivariate normal
distributions with different mean vectors and covariance matrices
depending on the value of \\X\\. Specifically: \$\$P \mid X = 0 \sim
\mathcal{N}\_{4}(\mu^{0}, \Sigma^{0}), \\ P \mid X = 1 \sim
\mathcal{N}\_{4}(\mu^{1}, \Sigma^{1}).\$\$

## References

Carlotti P, Parast L (2026). “A Bayesian Critique of Rank-Based Methods
for Surrogate Marker Evaluation.” *arXiv preprint arXiv:2603.14381*.

## Examples

``` r
set.seed(123)
data <- DGP_X_binary(
  n = 100,
  p = 0.5,
  q = 0.5,
  mu_0 = c(5, 5, 0, 0),
  mu_1 = c(5, -5, 0, -10),
  Sigma_0 = kronecker(diag(2), matrix(c(1, 1, 1, 2), 2, 2)),
  Sigma_1 = kronecker(diag(2), matrix(c(1, 1, 1, 2), 2, 2))
)
```
