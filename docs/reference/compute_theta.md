# Monte Carlo Computation of the Parameter \\\theta\\ from Carlotti and Parast (2026)

This function implements a Monte Carlo approach to estimate the
parameter \\\theta\\ from Carlotti and Parast (2026). This parameter
represents the difference in treatment effects between the primary and
surrogate outcomes, both measured using the probability that the treated
outcome is larger than the control outcome.

## Usage

``` r
compute_theta(MC_data)
```

## Arguments

- MC_data:

  A list containing:

  - `P`: A matrix or data frame of potential outcomes with columns "Y1",
    "Y0", "S1", and "S0".

## Value

A list containing the true values:

- `V_Y`: The Monte Carlo estimate of \\P(Y\_{1i} \> Y\_{0i})\\ computed
  on `P`.

- `V_S`: The Monte Carlo estimate of \\P(S\_{1i} \> S\_{0i})\\ computed
  on `P`.

- `theta`: The difference `V_Y` - `V_S`.

## Details

The function processes data from a chosen data generating process,
computing the sample probabilities for both the primary outcome \\Y\\
and the surrogate \\S\\: \$\$\hat{V}\_Y = \frac{1}{n}
\sum\limits^{n}\_{i=1} I(Y\_{1i} \> Y\_{0i}),\$\$ \$\$\hat{V}\_S =
\frac{1}{n} \sum\limits^{n}\_{i=1} I(S\_{1i} \> S\_{0i}).\$\$ Then, it
calculates \$\$\hat{\theta} = \hat{V}\_Y - \hat{V}\_S.\$\$
