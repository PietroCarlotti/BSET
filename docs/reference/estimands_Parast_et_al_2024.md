# True estimands for Parast et al. (2024)

A dataset containing the Monte Carlo estimates of the parameters of
interest for the simulation settings considered in Parast et al. (2024).

## Usage

``` r
estimands_Parast_et_al_2024
```

## Format

A data frame with 4 rows and 7 columns:

- setting:

  The index of the simulation setting.

- U_Y_MC:

  Numeric vector of Monte Carlo estimates.

- U_S_MC:

  Numeric vector of Monte Carlo estimates.

- delta_MC:

  Numeric vector of Monte Carlo estimates.

- V_Y_MC:

  Numeric vector of Monte Carlo estimates.

- V_S_MC:

  Numeric vector of Monte Carlo estimates.

- theta_MC:

  Numeric vector of Monte Carlo estimates.

## Source

Computed using the function
`compute_estimands_Parast_et_al_2024(1000000)`.
