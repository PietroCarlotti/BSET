# Compute \\v_S\\ from Carlotti and Parast (2026)

This function determines the value \\v_S\\ that is used to compute the
surrogate validation threshold \\\eta\\ from Carlotti and Parast (2026):
\$\$\eta = \max \\v_Y - v_S, 0\\,\$\$ where \\v_Y\\ is the hypothesized
value of the treatment effect on the primary outcome (typically set
equal to the estimate computed on the available data) and \\v_S\\ is the
value that satisfies the following power constraint: \$\$P(\text{BF}\_n
\geq \text{BF}\_{n, \alpha} \\ \| \\ V_S = v_S) = 1 - \beta,\$\$ where
\\\text{BF}\_{n, \alpha}\\ is the \\(1 - \alpha)\\ quantile of the Bayes
factor distribution under the null hypothesis \\V_S = V^0\_{S}\\, and
\\1 - \beta\\ is the desired power of the test. The computes the
distribution of the Bayes factor under the null hypothesis, derives the
critical value \\\text{BF}\_{n, \alpha}\\, and then uses a root-finding
algorithm to solve for the value of \\v_S\\ that satisfies the power
constraint.

## Usage

``` r
compute_V_S_star(
  n,
  alpha = 0.05,
  beta = 0.2,
  V_S_zero = 0.5,
  a = 1,
  b = 1,
  BF_alternative = "greater",
  root_tolerance = 1e-16
)
```

## Arguments

- n:

  Integer. Sample size.

- alpha:

  Numeric. Type I error rate (default is 0.05).

- beta:

  Numeric. Type II error rate (default is 0.2).

- V_S_zero:

  Numeric. The hypothesized value of the surrogate's treatment effect
  under the null hypothesis (default is 0.5).

- a:

  Numeric. First shape parameter alpha for the Beta prior (default is
  1).

- b:

  Numeric. Second shape parameter beta for the Beta prior (default is
  1).

- BF_alternative:

  Character. The type of alternative hypothesis: either `"two_sided"` or
  `"greater"`.

- root_tolerance:

  Numeric. Tolerance level for the root-finding algorithm (default is
  1e-16).

## Value

A list containing:

- `BF_alpha`: The critical value of the Bayes factor corresponding to
  the specified alpha level.

- `V_S_star`: The value of \\v_S\\ that satisfies the power constraint
  for the surrogate validation test.
