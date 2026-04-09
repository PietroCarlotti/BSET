# Frequentist Test from Parast et al. (2024)

This function performs a frequentist test for surrogate evaluation based
on the rank-based methodology of Parast et al. (2024). It calculates
confidence intervals, the \\\varepsilon\\ threshold used, and determines
if the surrogate is valid.

## Usage

``` r
frequentist_test(P_observed, Z, alpha = 0.05, beta = 0.2, delta_true = NULL)
```

## Arguments

- P_observed:

  Observed outcomes \\(Y, S)\\ corresponding to the assigned treatment
  \\Z\\.

- Z:

  Treatment assignment vector.

- alpha:

  Numeric. Significance level for the confidence interval (default is
  0.05).

- beta:

  Numeric. Type II error rate (default is 0.2).

- delta_true:

  Numeric (optional). The true value of \\\delta\\, used to calculate
  frequentist coverage during simulations.

## Value

A list containing:

- `CI`: The calculated confidence interval for \\\delta\\.

- `epsilon`: The \\\varepsilon\\ threshold value used in the test.

- `coverage`: Logical indicating if `delta_true` falls within the `CI`
  (if `delta_true` is provided).

- `power`: Logical indicating if the upper bound of `CI` is below
  `epsilon`, which indicates that the test identifies the surrogate as
  valid.
