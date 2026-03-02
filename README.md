# Causal mediation analysis for zero-inflated longitudinal data in the presence of treatment non-compliance and multiple mediators

This repository contains R code for the paper **Causal mediation analysis for zero-inflated longitudinal data in the presence of treatment non-compliance and multiple mediators**. The analysis investigates how treatment assignment (discount vs. non-discount promotional emails) affects customer order counts through two mediating pathways: time since the last email was opened and time since the last purchase. The methodology combines:

- **Bayesian nonparametric mixture modeling** (Enriched Dirichlet Process Mixture Model (EDPM))
- **Principal stratification** to handle noncompliance (email opened or not)
- **G-computation** for causal effect estimation via Monte Carlo integration
- **Simulation study** comparing the EDPM model against a parametric model

### Causal Effects Estimated

For each customer stratum, three principal causal effects are computed:

| Effect | Description |
|--------|-------------|
| **PIDE** | Principal Interventional Direct Effect — mediation effect through treatment receipt and mediators |
| **PJIIE** | Principal Joint Interventional Indirect Effect — effect through mediators |
| **PCE** | Principal Causal Effect — total effect (PIDE + PJIIE) |

### Principal Strata (Customer Segments)

| Stratum | Code | Behavior |
|---------|------|----------|
| **Active Customers** | AC | Always open emails regardless of assignment at t = 1|
| **Value Attentive** | VA | Open emails only when assigned (non-price) value-incentive emails at t = 1  |
| **Price Attentive** | PA | Open emails only when assigned price-incentive emails at t = 1|
| **Non-Active Customers** | NAC | Never open emails regardless of assignment at t = 1|

## Data and Variables

The analysis uses longitudinal email marketing data observed across **3 time periods**:

| Variable | Description | Type |
|----------|-------------|------|
| **Y_t** | Order count | Hurdle model (zero-inflated) |
| **M1_t** | Days since email opened | Continuous |
| **M2_t** | Days since purchase | Continuous |
| **D_t** | Email opened (treatment receipt) | Binary (probit) |
| **Z_t** | Discount email received (treatment assignment) | Binary (probit) |
| **L0** | Baseline confounders: gender (binary), account age (binary) | Baseline covariates |

## Model Specifications

### EDPM Model
- **Outer clusters :** N=10
- **Inner clusters :** M=4
- Stick-breaking priors on cluster weights
- Random intercepts for Z, D, M1, M2 across subjects

### MCMC Settings (data analysis)
- Iterations: 20,000
- Burn-in: 17,500
- Thinning: 10
- Chains: 4
- Posterior samples retained: 250 per chain (1,000 total)

## Repository Structure

### Model Fitting
| File | Description |
|------|-------------|
| `EDPLongNimbleModelFit.R` | NIMBLE implementation of the EDPM model — data preparation, prior elicitation, model definition, MCMC sampling, and posterior extraction |

### Simulation Study

The simulation study evaluates model performance under two data-generating scenarios across 500 replications:

- **True DGP:** Data generated from the correctly specified parametric model
- **Misspecified DGP:** Data generated from a K=3 finite mixture (model fits single-component)

| File | Description |
|------|-------------|
| `Simulation-TrueDGPEDPMFitRun500.R` | True DGP data generation + EDPM model fitting |
| `Simulation-MissDGPEDPMFitRun500.R` | Misspecified DGP data generation + EDPM model fitting |
| `Simulation-TrueDGPParametricFitRun500.R` | True DGP data generation + Parametric model fitting |
| `Simulation-MissDGPParametricFitRun500.R` | Misspecified DGP data generation + Parametric model fitting |

### Simulation Results Analysis

| File | Description |
|------|-------------|
| `SimulationTrueDGPEDPMFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for True DGP + EDPM |
| `SimulationMissDGPEDPMFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for Misspecified DGP + EDPM |
| `SimulationTrueDGPParametricFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for True DGP + Parametric |
| `SimulationMissDGPParametricFit-Bias,MSE,Coverage.Rmd` | Bias, MSE, and 95% CI coverage for Misspecified DGP + Parametric |

### G-Computation (Observed Data Analysis)

| File | Description |
|------|-------------|
| `Gcomp.Rmd` | G-computation code |
| `GcompPrintThetaZZstar.Rmd` | Extraction and display of counterfactual outcomes |
| `GcompFixedDtForceOpenAtTime23FixedMt0.Rmd` | G-computation with constrained treatment receipt and mediator values |
| `GcompCalibrateMediator.Rmd` | G-computation with mediator calibration |



## Important Notes on Simulation Sample Size

The simulation code currently generates datasets with **sample size n = 12,000** only. To run simulations at **n = 15,000**, duplicate the relevant simulation execution and results analysis sections in each file and update `replicated_data_sample_size` to `15000`.

## Dependencies

The following R packages are required:

- `nimble` — Bayesian nonparametric MCMC
- `coda` — MCMC diagnostics
- `data.table` — data manipulation
- `tidyverse` / `dplyr` / `tidyr` — data wrangling
- `ggplot2` — visualization
- `truncnorm` — truncated normal distributions
- `mvtnorm` — multivariate normal sampling
- `lme4` / `nlme` — mixed effects models (prior elicitation)
- `haven` — reading data files
- `bayesplot` — MCMC diagnostics plots
- `sets` — set operations

## Running the Simulations

### On an HPC Cluster (SLURM)

```bash
# Submit 500 array jobs for each simulation scenario
sbatch SimTrueDGPEDPMFitRun500.sbatch
sbatch SimMissDGPEDPMFitRun500.sbatch
sbatch SimTrueDGPParametricFitRun500.sbatch
sbatch SimMissDGPParametricFitRun500.sbatch
```

Each batch script submits 500 independent array jobs. Results are collected in `EDPM_simu_results/` and `param_simu_results/` directories.

### Analyzing Results

After all simulation jobs complete, render the Rmd files to compute bias, MSE, and coverage:

```r
rmarkdown::render("SimulationTrueDGPEDPMFit-Bias,MSE,Coverage.Rmd")
rmarkdown::render("SimulationMissDGPEDPMFit-Bias,MSE,Coverage.Rmd")
rmarkdown::render("SimulationTrueDGPParametricFit-Bias,MSE,Coverage.Rmd")
rmarkdown::render("SimulationMissDGPParametricFit-Bias,MSE,Coverage.Rmd")
```

### Observed Data Analysis

Run the G-computation Rmd files after fitting the EDPM model to observed data via `EDPLongNimbleModelFit.R`:

```r
rmarkdown::render("Gcomp.Rmd")
rmarkdown::render("GcompPrintThetaZZstar.Rmd")
rmarkdown::render("GcompFixedDtForceOpenAtTime23FixedMt0.Rmd")
rmarkdown::render("GcompCalibrateMediator.Rmd")
```
