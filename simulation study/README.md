# Simulation Study

This folder contains R code, SLURM batch scripts, and result-analysis notebooks for the simulation study comparing the **EDPM** model against a **Bayesian parametric** model under two data-generating scenarios across 500 replications.

For the full project description (methodology, causal estimands, principal strata, variable definitions), see the [top-level README](../README.md).

## Scenarios

- **True DGP:** Data generated from the correctly specified parametric model
- **Misspecified DGP:** Data generated from a K=10 finite mixture (fitted model is single-component). To run the misspecified DGP with a 3-component mixture instead, manually set `K_mix_DG = 3` near the top of the data-generating section in the relevant simulation file(s).

## Contents

### Simulation Replication Scripts

| File | Description |
|------|-------------|
| `Simulation-TrueDGPEDPMFitRun500.R` | True DGP data generation + EDPM model fitting |
| `Simulation-MissDGPEDPMFitRun500.R` | Misspecified DGP data generation + EDPM model fitting |
| `Simulation-TrueDGPParametricFitRun500.R` | True DGP data generation + Parametric model fitting |
| `Simulation-MissDGPParametricFitRun500.R` | Misspecified DGP data generation + Parametric model fitting |

### SLURM Batch Scripts

| File | Description |
|------|-------------|
| `SimTrueDGPEDPMFitRun500.sbatch` | Array submission for True DGP + EDPM |
| `SimMissDGPEDPMFitRun500.sbatch` | Array submission for Misspecified DGP + EDPM |
| `SimTrueDGPParametricFitRun500.sbatch` | Array submission for True DGP + Parametric |
| `SimMissDGPParametricFitRun500.sbatch` | Array submission for Misspecified DGP + Parametric |

### Results Analysis

| File | Description |
|------|-------------|
| `SimulationTrueDGPEDPMFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for True DGP + EDPM |
| `SimulationMissDGPEDPMFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for Misspecified DGP + EDPM |
| `SimulationTrueDGPParametricFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for True DGP + Parametric |
| `SimulationMissDGPParametricFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for Misspecified DGP + Parametric |

## Running the Simulations on an HPC Cluster (SLURM)

```bash
sbatch SimTrueDGPEDPMFitRun500.sbatch
sbatch SimMissDGPEDPMFitRun500.sbatch
sbatch SimTrueDGPParametricFitRun500.sbatch
sbatch SimMissDGPParametricFitRun500.sbatch
```

Each batch script submits 500 independent array jobs. Results are collected in `EDPM_simu_results/` and `param_simu_results/` directories.

## Analyzing Results

After all simulation jobs complete, render the Rmd files to compute bias, MSE, and coverage:

```r
rmarkdown::render("SimulationTrueDGPEDPMFit-Bias,MSE,Coverage.Rmd")
rmarkdown::render("SimulationMissDGPEDPMFit-Bias,MSE,Coverage.Rmd")
rmarkdown::render("SimulationTrueDGPParametricFit-Bias,MSE,Coverage.Rmd")
rmarkdown::render("SimulationMissDGPParametricFit-Bias,MSE,Coverage.Rmd")
```

## Important Note on Simulation Sample Size

The simulation code currently generates datasets with **sample size n = 12,000** only. To run simulations at **n = 15,000**, duplicate the relevant simulation execution and results-analysis sections in each file and update `replicated_data_sample_size` to `15000`.
