# Run Simulations for the Settings in Parast et al. (2024)

This function executes the simulation study described in Parast et al.
(2024) across four different settings:

1.  Useless surrogate.

2.  Perfect surrogate.

3.  Imperfect surrogate.

4.  Misspecified model.

It performs the frequentist and Bayesian surrogate evaluation tests from
Parast et al. (2024) and Carlotti and Parast (2026), respectively, and
returns a structured data frame of the results. The execution is
parallelized using the `future` framework.

## Usage

``` r
Parast_et_al_2024_simulations(seed, n_simulations, parallel)
```

## Arguments

- seed:

  Numeric. A random seed for reproducibility of the simulations.

- n_simulations:

  Numeric. The number of simulations to run for each setting.

- parallel:

  Logical. Whether to run the simulations in parallel or sequentially
  (default is TRUE).

## Value

A data frame containing:

- `setting`: The index of the simulation setting (1 to 4).

- `simulation`: The individual simulation run ID.

- `Bayesian_epsilon`, `frequentist_epsilon`: Numeric. Surrogate
  validation thresholds for the Bayesian and frequentist tests,
  respectively.

- `Bayesian_coverage`, `frequentist_coverage`: Boolean. Indicates
  whether the true surrogate effect is within the credible or confidence
  interval for the Bayesian and frequentist tests, respectively.

- `Bayesian_power`, `frequentist_power`: Boolean. Indicates whether the
  upper bound of the credible or confidence interval is below the
  threshold for the Bayesian and frequentist tests, respectively.

- `n_sample`, `n_simulations`, `n_chains`, `n_iterations`: Numeric.
  Metadata about the simulation parameters.

- `timestamp`: Character. Date and time when the simulation was
  completed.

- `seed`: Numeric. The random seed used for reproducibility.

## Details

The function runs a total of 500 simulations for each of the four
settings, generating data using `DGP_no_X` and processes the results
using `BSET_no_X`. The simulation utilizes MCMC sampling via `rstan` for
the Bayesian estimation components. Note that it relies on an external
object `estimands_Parast_et_al_2024` for true parameter values, which
are computed using the function `compute_estimands_Parast_et_al_2024`.
