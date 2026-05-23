# Real Data Analysis

This folder contains R code for fitting the **Enriched Dirichlet Process Mixture (EDPM)** model to the observed data and for running the **G-computation** procedure used to estimate principal causal effects.


## Contents

### Model Fitting

| File | Description |
|------|-------------|
| `EDPLongNimbleModelFit.R` | NIMBLE implementation of the EDPM model — data preparation, prior elicitation, model definition, MCMC sampling, and posterior extraction |

### G-Computation (Observed Data Analysis)

| File | Description |
|------|-------------|
| `Gcomp.Rmd` | G-computation code |
| `GcompPrintThetaZZstar.Rmd` | Extraction and display of counterfactual outcomes |
| `GcompFixedDtForceOpenAtTime23FixedEmpiricalMt.Rmd` | G-computation with constrained treatment receipt and mediator values |
| `GcompCalibrateMediator.Rmd` | G-computation with mediator calibration |

## MCMC Settings (data analysis)

- Iterations: 20,000
- Burn-in: 17,500
- Thinning: 10
- Chains: 4
- Posterior samples retained: 250 per chain (1,000 total)

## Workflow

1. Fit the EDPM model to the observed data:

   ```r
   source("EDPLongNimbleModelFit.R")
   ```

2. Run the G-computation Rmd files using the saved posterior samples:

   ```r
   rmarkdown::render("Gcomp.Rmd")
   rmarkdown::render("GcompPrintThetaZZstar.Rmd")
   rmarkdown::render("GcompFixedDtForceOpenAtTime23FixedEmpiricalMt.Rmd")
   rmarkdown::render("GcompCalibrateMediator.Rmd")
   ```
