# Simulation grid for Carlotti and Parast (2026)

A dataset containing the grid of simulation results for Carlotti and
Parast (2026) across two covariate settings: a binary covariate setting
(setting 1) and a Gaussian covariate setting (setting 2). Each setting
is run for 500 simulations, for a total of 1000 rows.

## Usage

``` r
Carlotti_and_Parast_2026_simulations_grid
```

## Format

A data frame with 1000 rows and 17 columns:

- setting:

  The index of the simulation setting: 1 for the binary covariate
  setting (`X_binary`) and 2 for the Gaussian covariate setting
  (`X_Gaussian`).

- simulation:

  The index of the simulation run.

- n_sample:

  The sample size used.

- n_chains:

  The number of MCMC chains used.

- n_iterations:

  The number of MCMC iterations.

- n_simulations:

  The total number of simulations per setting.

- timestamp:

  The timestamp of execution.

- seed:

  The random seed used.

- BF_alternative:

  The alternative hypothesis used for the Bayes factor.

- Bayesian_epsilon:

  The threshold for the Bayesian test.

- frequentist_epsilon:

  The threshold for the frequentist test.

- Bayesian_CI_upper:

  The upper bound of the Bayesian credible interval.

- frequentist_CI_upper:

  The upper bound of the frequentist confidence interval.

- Bayesian_coverage:

  Logical. Indicates if the true value is in the credible interval.

- frequentist_coverage:

  Logical. Indicates if the true value is in the confidence interval.

- Bayesian_power:

  Logical. Indicates if the Bayesian test rejected the null.

- frequentist_power:

  Logical. Indicates if the frequentist test rejected the null.

## Source

Generated using the code in the `simulations` folder of the BSET GitHub
repository (<https://github.com/PietroCarlotti/BSET>).
