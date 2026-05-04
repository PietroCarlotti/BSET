# Replication Materials

## Simulations

The `simulations/` folder contains the following files:

- `Carlotti_and_Parast_2026_simulations.R` and `Parast_et_al_2024_simulations.R`: R functions implementing the simulation study for each method.
- `Carlotti_and_Parast_2026_batch.R` and `Parast_et_al_2024_batch.R`: R scripts that run a single batch of simulations for each method.
- `Carlotti_and_Parast_2026.sh` and `Parast_et_al_2024.sh`: shell scripts that submit all batches and reproduce the full simulation results reported in the paper.

To replicate the simulation study, run the shell script for the desired method from the `simulations/` directory.

## Figures

The `figures/` folder contains `generate_figures.R`, an R script that reproduces all figures in the paper.

Note that some figures require access to the DCCT dataset used in the real data application. This dataset is not distributed with the package; instructions for requesting access can be found in the Data Availability section of the paper.
