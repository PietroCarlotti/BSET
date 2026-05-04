# True estimands for Carlotti and Parast (2026)

A dataset containing the Monte Carlo estimates of the parameters of
interest for the two simulation settings considered in Carlotti and
Parast (2026): setting 1 (binary covariate) and setting 2 (Gaussian
covariate).

## Usage

``` r
estimands_Carlotti_and_Parast_2026
```

## Format

A data frame with 2 rows and 7 columns:

- setting:

  The index of the simulation setting.

- U_Y_MC:

  Monte Carlo estimate of \\U_Y = P(Y\_{1i} \> Y\_{0j})\\.

- U_S_MC:

  Monte Carlo estimate of \\U_S = P(S\_{1i} \> S\_{0j})\\.

- delta_MC:

  Monte Carlo estimate of \\\delta = U_Y - U_S\\.

- V_Y_MC:

  Monte Carlo estimate of \\V_Y = P(Y\_{1i} \> Y\_{0i})\\.

- V_S_MC:

  Monte Carlo estimate of \\V_S = P(S\_{1i} \> S\_{0i})\\.

- theta_MC:

  Monte Carlo estimate of \\\theta = V_Y - V_S\\.

## Source

Computed using the function
`compute_estimands_Carlotti_and_Parast_2026(1000000)`.
