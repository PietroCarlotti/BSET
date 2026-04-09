# Monte Carlo Computation of the Estimands for the Simulation Study in Carlotti and Parast (2026)

This function iterates through the simulation settings defined in
Carlotti and Parast (2026) and estimates the true values of \\U_Y\\,
\\U_S\\, \\\delta\\, \\V_Y\\, \\V_S\\, and \\\theta\\ using a Monte
Carlo dataset generated according to the specified data-generating
processes. The settings are defined as follows:

- Setting 1: **X binary**

  - If \\X = 1\\: \$\$Y_1 \sim \mathcal{N}(5, 1), \quad Y_0 \sim
    \mathcal{N}(0, 1)\$\$ \$\$S_1 = Y_1 + \mathcal{N}(0, 1), \quad S_0 =
    Y_0 + \mathcal{N}(0, 1)\$\$

  - If \\X = 0\\: \$\$Y_1 \sim \mathcal{N}(5, 1), \quad Y_0 \sim
    \mathcal{N}(0, 1)\$\$ \$\$S_1 = Y_1 + \mathcal{N}(-10, 1), \quad S_0
    = Y_0 + \mathcal{N}(-10, 1)\$\$

## Usage

``` r
compute_estimands_Carlotti_and_Parast_2026(MC_samples)
```

## Arguments

- MC_samples:

  Integer. The number of Monte Carlo samples to generate per setting.

## Value

A data frame containing the Monte Carlo estimates for each setting:

- `setting`: The index of the simulation setting.

- `U_Y_MC`, `U_S_MC`, `delta_MC`: Parameters of interest from Parast et
  al. (2024).

- `V_Y_MC`, `V_S_MC`, `theta_MC`: Parameters of interest from Carlotti
  and Parast (2026).
