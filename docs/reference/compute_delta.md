# Monte Carlo Computation of the Parameter \\\delta\\ from Parast et al. (2024)

This function implements a Monte Carlo approach to estimate the
parameter \\\delta\\ from Parast et al. (2024) . This parameter
represents the difference in treatment effects between the primary and
surrogate outcomes, both measured using the Mann-Whitney statistic.

## Usage

``` r
compute_delta(MC_data)
```

## Arguments

- MC_data:

  A list containing:

  - `P_observed`: A data frame or matrix with columns "Y" and "S".

  - `Z`: Treatment assignment vector.

  - `n1`: Number of treated units.

  - `n0`: Number of control units.

## Value

A list containing:

- `U_Y`: Mann-Whitney U statistic for the primary outcome Y computed on
  `P_observed`.

- `U_S`: Mann-Whitney U statistic for the surrogate S computed on
  `P_observed`.

- `delta`: The difference `U_Y` - `U_S`.

## Details

The function processes data from a chosen data generating process,
computing the Mann-Whitney U statistic for both the primary outcome
\\Y\\ and the surrogate \\S\\: \$\$\hat{U}\_Y = \frac{1}{n_1 n_0}
\sum\limits\_{i:Z_i=1} \sum\limits\_{j:Z_j=0} I(Y_i \> Y_j),\$\$
\$\$\hat{U}\_S = \frac{1}{n_1 n_0} \sum\limits\_{i:Z_i=1}
\sum\limits\_{j:Z_j=0} I(S_i \> S_j).\$\$ Then, it calculates
\$\$\hat{\delta} = \hat{U}\_Y - \hat{U}\_S.\$\$ This function is
generally not intended to be called directly by the user and is instead
used internally within `BSET_no_X` and `BSET_X`.

## References

Parast L, Cai T, Tian L (2024). “A rank-based approach to evaluate a
surrogate marker in a small sample setting.” *Biometrics*, **80**(1),
ujad035.
