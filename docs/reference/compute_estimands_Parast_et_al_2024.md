# Monte Carlo Computation of the Estimands for the Simulation Study in Parast et al. (2024)

This function iterates through the four simulation settings defined in
Parast et al. (2024) and estimates the true values of \\U_Y\\, \\U_S\\,
\\\delta\\, \\V_Y\\, \\V_S\\, and \\\theta\\ using a Monte Carlo dataset
generated according to the specified data-generating processes.

## Usage

``` r
compute_estimands_Parast_et_al_2024(MC_samples)
```

## Arguments

- MC_samples:

  Integer. The number of Monte Carlo samples to generate per setting.

## Value

A data frame containing the Monte Carlo estimates for each setting:

- `setting`: The index of the simulation setting.

- `U_Y_MC`, `U_S_MC`, `delta_MC`: Parameters of interest from Parast et
  al. (2024) .

- `V_Y_MC`, `V_S_MC`, `theta_MC`: Parameters of interest from Carlotti
  and Parast (2026) .

## Details

The settings are defined as follows:

- Setting 1: **Useless surrogate (Gaussian model)** \$\$Y_1 \sim
  \mathcal{N}(5, 3), \quad Y_0 \sim \mathcal{N}(3, 3),\$\$ \$\$S_1 \sim
  \mathcal{N}(2/3, 1), \quad S_0 \sim \mathcal{N}(2/3, 1).\$\$

- Setting 2: **Perfect surrogate (Gaussian model)** \$\$Y_1 \sim
  \mathcal{N}(6, 3), \quad Y_0 \sim \mathcal{N}(5/2, 3),\$\$ \$\$S_1 =
  Y_1 + \mathcal{N}(0, 1/10), \quad S_0 = Y_0 + \mathcal{N}(0,
  1/10).\$\$

- Setting 3: **Imperfect surrogate (Gaussian model)** \$\$S_1 \sim
  \mathcal{N}(5, 3), \quad S_0 \sim \mathcal{N}(3, 3),\$\$ \$\$Y_1 =
  S_1 + \mathcal{N}(1.5, 0.6), \quad Y_0 = S_0 + \mathcal{N}(0,
  0.6).\$\$

- Setting 4: **Misspecified model (non-Gaussian model)** \$\$S_1 \sim
  \exp(\mathcal{N}(2.5, 1.5)), \quad S_0 \sim \exp(\mathcal{N}(0.5,
  1.5)),\$\$ \$\$Y_1 = 2 + 6/5 \sqrt{S_1} + 3/10 \exp(S_1 / 500) +
  \exp(\mathcal{N}(0, 0.3)),\$\$ \$\$Y_0 = 4/5 \sqrt{S_0} + 1/5 \exp(S_0
  / 50) + \exp(\mathcal{N}(0, 0.3)).\$\$

This function is generally not intended to be called directly by the
user. It is provided as a utility for computing the true parameter
values for the simulation settings described in Parast et al. (2024) .

## References

Parast L, Cai T, Tian L (2024). “A rank-based approach to evaluate a
surrogate marker in a small sample setting.” *Biometrics*, **80**(1),
ujad035.
