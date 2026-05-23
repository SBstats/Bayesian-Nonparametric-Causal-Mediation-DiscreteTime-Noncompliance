# Causal mediation analysis for zero-inflated longitudinal data in the presence of treatment non-compliance and multiple mediators

This repository contains R code for the paper **Causal mediation analysis for zero-inflated longitudinal data in the presence of treatment non-compliance and multiple mediators**. The analysis investigates how treatment assignment (discount vs. non-discount promotional emails) affects customer order counts through two mediating pathways: time since the last email was opened and time since the last purchase. The methodology combines:

- **Bayesian nonparametric mixture modeling** (Enriched Dirichlet Process Mixture Model (EDPM))
- **Principal stratification** to handle noncompliance (email opened or not)
- **G-computation** for causal effect estimation via Monte Carlo integration


## Repository Structure

| Folder | Contents |
|--------|----------|
| [`real data analysis/`](real%20data%20analysis/) | EDPM model fitting (NIMBLE) and G-computation on the observed data. See its [README](real%20data%20analysis/README.md). |
| [`simulation study/`](simulation%20study/) | Simulation replication scripts, SLURM batch files, and bias/MSE/coverage analysis notebooks. See its [README](simulation%20study/README.md). |

## Causal Effects Estimated

For each customer stratum, three principal causal effects are computed:

| Effect | Description |
|--------|-------------|
| **PIDE** | Principal Interventional Direct Effect |
| **PJIIE** | Principal Joint Interventional Indirect Effect — effect through mediators |
| **PCE** | Principal Causal Effect — total effect (PIDE + PJIIE) |

## Principal Strata (Customer Segments)

| Stratum | Code | Behavior |
|---------|------|----------|
| **Active Customers** | AC | Always open emails regardless of assignment at t = 1 |
| **Value Attentive** | VA | Open emails only when assigned (non-price) value-incentive emails at t = 1 |
| **Price Attentive** | PA | Open emails only when assigned price-incentive emails at t = 1 |
| **Non-Active Customers** | NAC | Never open emails regardless of assignment at t = 1 |

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

## Model Specification

### EDPM Model
- **Outer clusters:** N = 10
- **Inner clusters:** M = 4
- Square-breaking (nested stick-breaking) priors on cluster weights
- Random intercepts for Z, D, M1, M2 across subjects

MCMC settings and run-specific details are documented in the [real data analysis README](real%20data%20analysis/README.md) and [simulation study README](simulation%20study/README.md).

## Dependencies

The following R packages are required across both subprojects:

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
