###############################################################################
## Enriched Dirichlet Process (EDP) Longitudinal Model — NIMBLE Implementation
##
## This script fits a Bayesian nonparametric mixture model (nested EDP) to
## longitudinal email-marketing data with:
##   - Outcome Y_t (order count, hurdle model)
##   - Two mediators M1_t, M2_t (days since opened/purchased,  hurdle)
##   - Treatment receipt D_t (email opened, binary probit)
##   - Treatment assignment Z_t (discount email, binary probit)
##   - Baseline confounders L0 (gender, age — both binary)
##
## Model structure:
##   Outer clusters (N=10): outcome Y parameters (beta_Y, sigma_Y, pi0_Y)
##   Inner clusters (M=4):  mediator/treatment parameters (theta_M, theta_D,
##                          theta_Z, pi0_M, sigma_M, pr_L0)
##   Stick-breaking priors on cluster weights (ksi_r, ksi_sr)
##   Random intercepts for Z, D, M1, M2 across subjects
##
## MCMC settings: Q=20000, warmup=17500, thin=10, 4 chains
###############################################################################

###############################################################################
## 1. Load libraries and data
###############################################################################
set.seed(1234)
library(sets)
library(coda)
library(haven)
library(plyr)
library(tidyverse)
library(nlme)
library(lme4)
library(dplyr)
library(data.table)
library(car)
library(ggplot2)
library(gridExtra)
library(bayesplot)

library(nimble, warn.conflicts = FALSE)

# Enable Nimble options for saving intermediates and verbose output
nimbleOptions(saveIntermediates = TRUE)  # Saves generated C++ code files
nimbleOptions(verbose = TRUE)           # Provides detailed compilation output

digCom_data <- read.csv("non-comp_dataset.csv")
digCom_df = data.frame(digCom_data)


###############################################################################
## 2. Prepare working data frame (long format)
###############################################################################

# Mode function for descriptive summaries
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

M1_t_mode = Mode(digCom_df$daysSinceOpened)
cat("The mode of M1_t is:", M1_t_mode, ".\n")

M2_t_mode = Mode(digCom_df$daysSincePurchased)
cat("The mode of M2_t is:", M2_t_mode, ".\n")

Y_t_mode = Mode(digCom_df$Order_cnt)
cat("The mode of Y_t is:", Y_t_mode, ".\n")

# Construct working data frame
digCom_workdf = data.frame(ID = digCom_df$recipient_id,
                           Recent = digCom_df$Recent,
                           GENDER = digCom_df$gender,
                           BASELINE_AGE = digCom_df$age,
                           t = digCom_df$t,
                           Z_t = digCom_df$Discount,
                           D_t = digCom_df$opened,
                           M1_t = digCom_df$daysSinceOpened,
                           M2_t = digCom_df$daysSincePurchased,
                           Y_t = digCom_df$Order_cnt)

# Binarize GENDER: M -> 1, others -> 0
digCom_workdf <- digCom_workdf %>%
  mutate(GENDER = ifelse(GENDER == "M", 1, 0))

# Binarize BASELINE_AGE: >= 4 -> 1, < 4 -> 0 (AGE means age of account/ subscription age)
digCom_workdf <- digCom_workdf %>%
  mutate(BASELINE_AGE = ifelse(BASELINE_AGE >= 4, 1, 0))

head(digCom_workdf)

# Report modes of variables
M1_t_mode_work = Mode(digCom_workdf$M1_t)
cat("The mode of M1_t is:", M1_t_mode_work, ".\n")

M2_t_mode_work = Mode(digCom_workdf$M2_t)
cat("The mode of M2_t is:", M2_t_mode_work, ".\n")

Y_t_mode_work = Mode(digCom_workdf$Y_t)
cat("The mode of Y_t is:", Y_t_mode_work, ".\n")

# Density plots for visual inspection
M1t_plot = ggplot(digCom_workdf, aes(x  = M1_t)) +
  geom_density() +
  labs(title = "M1_t") +
  theme_minimal()

M2t_plot = ggplot(digCom_workdf, aes(x  = M2_t)) +
  geom_density() +
  labs(title = "M2_t") +
  theme_minimal()

Yt_plot = ggplot(digCom_workdf, aes(x  = Y_t)) +
  geom_density() +
  labs(title = "Y_t") +
  theme_minimal()


###############################################################################
## 3. Reshape to wide format
###############################################################################
setDT(digCom_workdf)
digCom_workdf_wide = dcast(digCom_workdf, ID + GENDER+ BASELINE_AGE   ~ t,
                           value.var = c("Z_t", "D_t", "M1_t", "M2_t","Y_t"))
head(digCom_workdf_wide)

digCom_workdf_wide = na.omit(digCom_workdf_wide)


###############################################################################
## 4. Construct input matrices and compute MLE-based priors
###############################################################################

# L0 baseline confounders: Intercept + Gender + Age (all binary)
L0_intercept = rep(1,length(digCom_workdf_wide$ID))
L0_df_wide = data.frame(Intercept = L0_intercept,
                        GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE)

head(L0_df_wide)

Y_vec_wide = digCom_workdf_wide$Y_t_3
head(Y_vec_wide)


###############################################################################
## 4a. Helper functions for data-dependent base distribution parameters
###############################################################################

num_subj = length(unique(digCom_workdf_wide$ID))

# Linear regression: output ~ all covariates (for Y, M models)
fit_linear_reg_YvAllCovs = function(output_vec, all_Covs_df_wide){
  num_subjects = nrow(all_Covs_df_wide)
  lin_ref_fit = lm(output_vec~., data = all_Covs_df_wide)
  mu_beta = as.numeric(coef(lin_ref_fit))
  sig_sq_beta = as.numeric(diag(vcov(lin_ref_fit)))
  sig_sq_beta = (num_subjects/5)* sig_sq_beta
  sigma_sq_hat =  (summary(lin_ref_fit)$sigma)^2
  return_list = list(Coeff = mu_beta, Sig_sq_coeff = sig_sq_beta, Sigma_sq_hat = sigma_sq_hat)
  return(return_list)
}

# Linear regression: output ~ L0 only
fit_linear_reg = function(output_vec, L0_data_wide){
  num_subjects = nrow(L0_df_wide)
  lin_ref_fit = lm(output_vec~GENDER+BASELINE_AGE, data = L0_data_wide)
  mu_beta = as.numeric(coef(lin_ref_fit))
  sig_sq_beta = as.numeric(diag(vcov(lin_ref_fit)))
  sig_sq_beta = (num_subjects/5)* sig_sq_beta
  sigma_sq_hat =  (summary(lin_ref_fit)$sigma)^2
  return_list = list(Coeff = mu_beta, Sig_sq_coeff = sig_sq_beta, Sigma_sq_hat = sigma_sq_hat)
  return(return_list)
}

# Multivariate linear regression: cbind(outputs) ~ L0
fit_multiple_linear_reg = function(cbind_output_vecs, L0_data_wide){
  num_subjects = nrow(L0_df_wide)
  lin_ref_fit = lm(cbind_output_vecs~GENDER+BASELINE_AGE, data = L0_data_wide)
  mu_theta = coef(lin_ref_fit)
  sig_sq_theta = (num_subjects/5)* as.numeric(diag(vcov(lin_ref_fit)))
  sig_sq_theta = matrix(sig_sq_theta, nrow = ncol(L0_data_wide), ncol = ncol(cbind_output_vecs))
  sigma_hat =  cov(resid(lin_ref_fit))
  return_list = list(Coeff = mu_theta, Sig_sq_coeff = sig_sq_theta, Sigma_hat = sigma_hat)
  return(return_list)
}

# Logistic regression: output ~ all covariates (for D, Z models)
fit_logistic_reg = function(output_vec, all_Covs_df_wide){
  num_subjects = nrow(all_Covs_df_wide)
  logistic_ref_fit = glm(output_vec~., family = binomial(link = "logit"),
                         data = all_Covs_df_wide)
  mu_beta = as.numeric(coef(logistic_ref_fit))
  sig_sq_beta = as.numeric(diag(vcov(logistic_ref_fit)))
  sig_sq_beta = (num_subjects/5)* sig_sq_beta
  return_list = list(Coeff = mu_beta, Sig_sq_coeff = sig_sq_beta)
  return(return_list)
}


###############################################################################
## 4b. Fit MLE regressions for prior hyperparameters
###############################################################################

# Y ~ all covariates (L0 + Z,D,M history at all time points)
Y_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                              Z_t_1 = digCom_workdf_wide$Z_t_1, D_t_1 = digCom_workdf_wide$D_t_1,
                              M1_t_1 = digCom_workdf_wide$M1_t_1, M2_t_1 = digCom_workdf_wide$M2_t_1,
                              Z_t_2 = digCom_workdf_wide$Z_t_2, D_t_2 = digCom_workdf_wide$D_t_2,
                              M1_t_2 = digCom_workdf_wide$M1_t_2, M2_t_2 = digCom_workdf_wide$M2_t_2,
                              Z_t_3 = digCom_workdf_wide$Z_t_3, D_t_3 = digCom_workdf_wide$D_t_3,
                              M1_t_3 = digCom_workdf_wide$M1_t_3, M2_t_3 = digCom_workdf_wide$M2_t_3)

lin_reg_YvAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$Y_t_3[digCom_workdf_wide$Y_t_3 !=0],
                                              Y_reg_design_mat[digCom_workdf_wide$Y_t_3 !=0,])

# Design matrices for M regression (time-varying: L0 + Z,D history)
Mt3_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                              Z_t_1 = digCom_workdf_wide$Z_t_1, D_t_1 = digCom_workdf_wide$D_t_1,
                              Z_t_2 = digCom_workdf_wide$Z_t_2, D_t_2 = digCom_workdf_wide$D_t_2,
                              Z_t_3 = digCom_workdf_wide$Z_t_3, D_t_3 = digCom_workdf_wide$D_t_3)

Mt2_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                 Z_t_1 = digCom_workdf_wide$Z_t_1, D_t_1 = digCom_workdf_wide$D_t_1,
                                 Z_t_2 = digCom_workdf_wide$Z_t_2, D_t_2 = digCom_workdf_wide$D_t_2)

Mt1_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                 Z_t_1 = digCom_workdf_wide$Z_t_1, D_t_1 = digCom_workdf_wide$D_t_1)

# M2 ~ covariates (non-zero observations only, for hurdle continuous part)
lin_reg_M2t3vAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$M2_t_3[digCom_workdf_wide$M2_t_3 !=0],
                                              Mt3_reg_design_mat[digCom_workdf_wide$M2_t_3 !=0,])
lin_reg_M2t2vAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$M2_t_2[digCom_workdf_wide$M2_t_2 !=0],
                                              Mt2_reg_design_mat[digCom_workdf_wide$M2_t_2 !=0, ])
lin_reg_M2t1vAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$M2_t_1[digCom_workdf_wide$M2_t_1 !=0],
                                              Mt1_reg_design_mat[digCom_workdf_wide$M2_t_1 !=0, ])

# M1 ~ covariates (non-zero observations only)
lin_reg_M1t3vAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$M1_t_3[digCom_workdf_wide$M1_t_3 !=0],
                                              Mt3_reg_design_mat[digCom_workdf_wide$M1_t_3 !=0,])
lin_reg_M1t2vAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$M1_t_2[digCom_workdf_wide$M1_t_2 !=0],
                                              Mt2_reg_design_mat[digCom_workdf_wide$M1_t_2 !=0,])
lin_reg_M1t1vAllCovs = fit_linear_reg_YvAllCovs(digCom_workdf_wide$M1_t_1[digCom_workdf_wide$M1_t_1 !=0],
                                              Mt1_reg_design_mat[digCom_workdf_wide$M1_t_1 !=0,])

# Design matrices for D regression (time-varying: L0 + Z history)
Dt3_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                Z_t_1 = digCom_workdf_wide$Z_t_1,
                                Z_t_2 = digCom_workdf_wide$Z_t_2,
                                Z_t_3 = digCom_workdf_wide$Z_t_3)

Dt2_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                Z_t_1 = digCom_workdf_wide$Z_t_1,
                                Z_t_2 = digCom_workdf_wide$Z_t_2)

Dt1_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                Z_t_1 = digCom_workdf_wide$Z_t_1)

# D ~ covariates (logistic/probit)
logistic_reg_D_t3vAllCovs = fit_logistic_reg(digCom_workdf_wide$D_t_3, Dt3_reg_design_mat)
logistic_reg_D_t2vAllCovs = fit_logistic_reg(digCom_workdf_wide$D_t_2, Dt2_reg_design_mat)
logistic_reg_D_t1vAllCovs = fit_logistic_reg(digCom_workdf_wide$D_t_1, Dt1_reg_design_mat)

# Z ~ L0 (logistic/probit)
logistic_reg_Z_t3vL0 = fit_logistic_reg(digCom_workdf_wide$Z_t_3, L0_df_wide[, -1])
logistic_reg_Z_t2vL0 = fit_logistic_reg(digCom_workdf_wide$Z_t_2, L0_df_wide[, -1])
logistic_reg_Z_t1vL0 = fit_logistic_reg(digCom_workdf_wide$Z_t_1, L0_df_wide[, -1])


###############################################################################
## 4c. Extract MLE prior hyperparameters
###############################################################################

# Longitudinal data in wide format
longitudinal_data_wide_M2 = data.frame(M2_t_1 = digCom_workdf_wide$M2_t_1,
                                       M2_t_2 = digCom_workdf_wide$M2_t_2,
                                       M2_t_3 = digCom_workdf_wide$M2_t_3)
head(longitudinal_data_wide_M2)

longitudinal_data_wide_M1 = data.frame(M1_t_1 = digCom_workdf_wide$M1_t_1,
                                       M1_t_2 = digCom_workdf_wide$M1_t_2,
                                       M1_t_3 = digCom_workdf_wide$M1_t_3)
head(longitudinal_data_wide_M1)

# M2 MLE priors — separate vectors per time point (different covariate dimensions)
mle_coeff_M2_t1 = unlist(lin_reg_M2t1vAllCovs[1])
mle_coeff_M2_t2 = unlist(lin_reg_M2t2vAllCovs[1])
mle_coeff_M2_t3 = unlist(lin_reg_M2t3vAllCovs[1])

mle_sig_sq_M2_t1 = unlist(lin_reg_M2t1vAllCovs[2])
mle_sig_sq_M2_t2 = unlist(lin_reg_M2t2vAllCovs[2])
mle_sig_sq_M2_t3 = unlist(lin_reg_M2t3vAllCovs[2])

sigma_sq_hat_M2_t1 = unlist(lin_reg_M2t1vAllCovs[3])
sigma_sq_hat_M2_t2 = unlist(lin_reg_M2t2vAllCovs[3])
sigma_sq_hat_M2_t3 = unlist(lin_reg_M2t3vAllCovs[3])

# M1 MLE priors — separate vectors per time point
mle_coeff_M1_t1 = unlist(lin_reg_M1t1vAllCovs[1])
mle_coeff_M1_t2 = unlist(lin_reg_M1t2vAllCovs[1])
mle_coeff_M1_t3 = unlist(lin_reg_M1t3vAllCovs[1])

mle_sig_sq_M1_t1 = unlist(lin_reg_M1t1vAllCovs[2])
mle_sig_sq_M1_t2 = unlist(lin_reg_M1t2vAllCovs[2])
mle_sig_sq_M1_t3 = unlist(lin_reg_M1t3vAllCovs[2])

sigma_sq_hat_M1_t1 = unlist(lin_reg_M1t1vAllCovs[3])
sigma_sq_hat_M1_t2 = unlist(lin_reg_M1t2vAllCovs[3])
sigma_sq_hat_M1_t3 = unlist(lin_reg_M1t3vAllCovs[3])

# D longitudinal data and MLE priors
longitudinal_data_wide_D = data.frame(D_t_1 = digCom_workdf_wide$D_t_1,
                                      D_t_2 = digCom_workdf_wide$D_t_2,
                                      D_t_3 = digCom_workdf_wide$D_t_3)
head(longitudinal_data_wide_D)

mle_coeff_D_t1 = unlist(logistic_reg_D_t1vAllCovs[1])
mle_coeff_D_t2 = unlist(logistic_reg_D_t2vAllCovs[1])
mle_coeff_D_t3 = unlist(logistic_reg_D_t3vAllCovs[1])

mle_sig_sq_D_t1 = unlist(logistic_reg_D_t1vAllCovs[2])
mle_sig_sq_D_t2 = unlist(logistic_reg_D_t2vAllCovs[2])
mle_sig_sq_D_t3 = unlist(logistic_reg_D_t3vAllCovs[2])

# Design matrices for D fitting: prepend Intercept to MLE regression matrices
D_fit_design_mat_t1 = cbind(Intercept = L0_df_wide$Intercept, Dt1_reg_design_mat)
D_fit_design_mat_t2 = cbind(Intercept = L0_df_wide$Intercept, Dt2_reg_design_mat)
D_fit_design_mat_t3 = cbind(Intercept = L0_df_wide$Intercept, Dt3_reg_design_mat)

# Design matrices for M1/M2 fitting: prepend Intercept
M_fit_design_mat_t1 = cbind(Intercept = L0_df_wide$Intercept, Mt1_reg_design_mat)
M_fit_design_mat_t2 = cbind(Intercept = L0_df_wide$Intercept, Mt2_reg_design_mat)
M_fit_design_mat_t3 = cbind(Intercept = L0_df_wide$Intercept, Mt3_reg_design_mat)

# Z longitudinal data and MLE priors
longitudinal_data_wide_Z = data.frame(Z_t_1 = digCom_workdf_wide$Z_t_1,
                                      Z_t_2 = digCom_workdf_wide$Z_t_2,
                                      Z_t_3 = digCom_workdf_wide$Z_t_3)
head(longitudinal_data_wide_Z)

# Z coefficient priors: K * T matrix
mle_coeff_Z = data.frame(Z_t_1_coeff = unlist(logistic_reg_Z_t1vL0[1]),
                         Z_t_2_coeff = unlist(logistic_reg_Z_t2vL0[1]),
                         Z_t_3_coeff = unlist(logistic_reg_Z_t3vL0[1]))
head(mle_coeff_Z)

# Z variance priors: K * T matrix
mle_sig_sq_Z = data.frame( Z_t_1_sig_sq = unlist(logistic_reg_Z_t1vL0[2]),
                           Z_t_2_sig_sq = unlist(logistic_reg_Z_t2vL0[2]),
                           Z_t_3_sig_sq = unlist(logistic_reg_Z_t3vL0[2]))
head(mle_sig_sq_Z)

# L0 summary statistics for base distribution
num_Bin_L0 = 2  # both gender and baseline age are binary
pr_theta_Bin_L0 = as.numeric(colMeans(L0_df_wide[2:3], na.rm = TRUE))


###############################################################################
## 5. NIMBLE model specification
###############################################################################

# Custom nimbleFunction: record model sum-posterior-log-density
sumLogPostDens <- nimbleFunction(
  name = 'sumLogPostDens',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    stochNodes <- setdiff(model$getNodeNames(stochOnly = TRUE), target)
  },
  run = function() {
    model[[target]] <<- model$getLogProb(stochNodes)
  },
  methods = list( reset = function() {} )
)

# Custom zero-inflated hurdle normal distribution (density)
dzinhurdlenorm <- nimbleFunction(
  run = function(x = double(0), pi = double(0), mu = double(0), sigma = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if (sigma <= 0) stop("Sigma must be positive")
    logProb = 0
    if (x == 0) {
      logProb = log(max(1e-10, pi))  # Zero component
    } else {
      # Positive component: truncated normal (x > 0)
      truncProb = max(1e-10, 1 - pnorm(0, mean = mu, sd = sigma))
      logProb = log(max(1e-10, 1 - pi)) +
        dnorm(x, mean = mu, sd = sigma, log = TRUE) -
        log(truncProb)
    }
    if (log) return(logProb)
    return(exp(logProb))
  }
)

# Custom zero-inflated hurdle normal distribution (random generation)
rzinhurdlenorm <- nimbleFunction(
  run = function(n = integer(0), pi = double(0), mu = double(0), sigma = double(0)) {
    returnType(double(0))
    if (n != 1) stop("rzinhurdlenorm only allows n = 1")
    if (runif(1) < pi) {
      return(0)  # Zero-inflated part
    } else {
      x <- -1
      while (x <= 0) {
        x <- rnorm(1, mean = mu, sd = sigma)  # Truncated normal: resample until x > 0
      }
      return(x)
    }
  }
)

# Test custom distribution
dzinhurdlenorm(0, pi = 0.3, mu = 2, sigma = 1, log = TRUE)
rzinhurdlenorm(1, pi = 0.3, mu = 2, sigma = 1)

# Register the custom distribution in Nimble
registerDistributions(list(
  dzinhurdlenorm = list(
    BUGSdist = "dzinhurdlenorm(pi, mu, sigma)",
    types = c('value = double(0)', 'pi = double(0)', 'mu = double(0)', 'sigma = double(0)'),
    discrete = FALSE
  )
))


###############################################################################
## 5a. NIMBLE model code (nested EDP mixture)
###############################################################################

code <- nimbleCode({
  logDens ~ dnorm(0, 1)    ## placeholder for sum-posterior-log-density sampler

  # Outer-layer stick-breaking weights
  for(i in 1:(N-1)){
    breaks_r[i] ~ dbeta(1, alpha_bet)
  }
  breaks_r[N] <- 1
  ksi_r[1:N] <- stick_breaking(breaks_r[1:(N-1)])

  # Inner-layer stick-breaking weights (nested within each outer cluster)
  for(i in 1:N){
    for (j in 1:(M-1)){
      breaks_sr[i,j] ~ dbeta(1, alpha_thet[i])
    }
    breaks_sr[i,M] <- 1
    ksi_sr[i,1:M] <- stick_breaking(breaks_sr[i,1:(M-1)])
  }

  # Atoms from base distributions
  for(i in 1:N){
    # Outer cluster parameters: outcome model
    pi0_Y[i] ~ dbeta(2,2)
    sigma_sq_Yi[i] ~ dinvgamma(shape = 3, scale = 2 * sigSq_Y_prior)
    for(k in 1:num_Y_reg_coeff){
      beta_Yi[i, k]~ dnorm(mean = mu_betaY_prior[k], var = sigSq_betaY_prior[k])
    }

    for(j in 1:M){
      # Inner cluster parameters: baseline confounders
      for(k in 1:num_bin_L0){
        pr_Bin_L0i[i,j,k] ~ dbeta(2, 2)
      }

      # Mediator variance and hurdle probabilities (per time point)
      for(t in 1:max_t){
        sigma_sq_M2_t[i,j,t] ~ dinvgamma(shape = 3, scale = 2 * sigma_sq_hat_M2_t[t])
        sigma_sq_M1_t[i,j,t] ~ dinvgamma(shape = 3, scale = 2 *sigma_sq_hat_M1_t[t])
        pi0_M1_t[i,j,t]~ dbeta(2,2)
        pi0_M2_t[i,j,t]~ dbeta(2,2)
      }

      # Treatment assignment Z mean parameters (K=3 covariates at all t)
      for(k in 1:K){
        for(t in 1:max_t){
          theta_Z_t[i,j,k,t] ~ dnorm(mean = mu_theta_Z_t_prior[k,t],
                                     sd = sqrt(sigSq_theta_Z_t_prior[k,t]))
        }
      }

      # Treatment receipt D mean parameters (separate per time point)
      for(k in 1:K_D_t1){
        theta_D_t1[i,j,k] ~ dnorm(mean = mu_theta_D_t1_prior[k],
                                   sd = sqrt(sigSq_theta_D_t1_prior[k]))
      }
      for(k in 1:K_D_t2){
        theta_D_t2[i,j,k] ~ dnorm(mean = mu_theta_D_t2_prior[k],
                                   sd = sqrt(sigSq_theta_D_t2_prior[k]))
      }
      for(k in 1:K_D_t3){
        theta_D_t3[i,j,k] ~ dnorm(mean = mu_theta_D_t3_prior[k],
                                   sd = sqrt(sigSq_theta_D_t3_prior[k]))
      }

      # Mediator M2 mean parameters (separate per time point)
      for(k in 1:K_M_t1){
        theta_M2_t1[i,j,k] ~ dnorm(mean = mu_theta_M2_t1_prior[k],
                                    sd = sqrt(sigSq_theta_M2_t1_prior[k]))
      }
      for(k in 1:K_M_t2){
        theta_M2_t2[i,j,k] ~ dnorm(mean = mu_theta_M2_t2_prior[k],
                                    sd = sqrt(sigSq_theta_M2_t2_prior[k]))
      }
      for(k in 1:K_M_t3){
        theta_M2_t3[i,j,k] ~ dnorm(mean = mu_theta_M2_t3_prior[k],
                                    sd = sqrt(sigSq_theta_M2_t3_prior[k]))
      }

      # Mediator M1 mean parameters (separate per time point)
      for(k in 1:K_M_t1){
        theta_M1_t1[i,j,k] ~ dnorm(mean = mu_theta_M1_t1_prior[k],
                                    sd = sqrt(sigSq_theta_M1_t1_prior[k]))
      }
      for(k in 1:K_M_t2){
        theta_M1_t2[i,j,k] ~ dnorm(mean = mu_theta_M1_t2_prior[k],
                                    sd = sqrt(sigSq_theta_M1_t2_prior[k]))
      }
      for(k in 1:K_M_t3){
        theta_M1_t3[i,j,k] ~ dnorm(mean = mu_theta_M1_t3_prior[k],
                                    sd = sqrt(sigSq_theta_M1_t3_prior[k]))
      }

    }
  }

  # Random intercept variances
  bi_Z_sig_sq ~ dinvgamma(shape = 3, scale = 2)
  bi_D_sig_sq ~ dinvgamma(shape = 3, scale = 2)
  bi_M1_sig_sq ~ dinvgamma(shape = 3, scale = 2)
  bi_M2_sig_sq ~ dinvgamma(shape = 3, scale = 2)

  # Likelihood
  for(i in 1:n) {
    # Cluster membership indicators
    comp_num_outer[i] ~ dcat(ksi_r[1:N])
    comp_num_inner[i] ~ dcat(ksi_sr[comp_num_outer[i],1:M])

    # Subject-specific random intercepts
    bi_Z[i] ~ dnorm(0,sqrt(bi_Z_sig_sq))
    bi_D[i] ~ dnorm(0,sqrt(bi_D_sig_sq))
    bi_M1[i] ~ dnorm(0,sqrt(bi_M1_sig_sq))
    bi_M2[i] ~ dnorm(0,sqrt(bi_M2_sig_sq))

    # Binary baseline confounders
    for(k in 1:num_bin_L0){
      L0_bin[i,k] ~ dbern(pr_Bin_L0i[comp_num_outer[i], comp_num_inner[i],k])
    }

    # ========== t = 1 ==========
    # Z at t=1 (probit link)
    z_mean[i,1] <- inprod(theta_Z_t[comp_num_outer[i], comp_num_inner[i], 1:K, 1],
                          L0_mat[i, 1:K]) + bi_Z[i]
    data_wide_Z[i,1] ~ dbern(pnorm(z_mean[i,1],0,1))

    # D at t=1 (probit link)
    d_mean[i,1] <- inprod(theta_D_t1[comp_num_outer[i], comp_num_inner[i], 1:K_D_t1],
                          D_design_mat_t1[i, 1:K_D_t1]) + bi_D[i]
    data_wide_D[i,1] ~ dbern(pnorm(d_mean[i,1],0,1))

    # M1 at t=1 (hurdle model: latent indicator * truncated normal)
    Z_latent_M1[i,1] ~ dbern(1 - pi0_M1_t[comp_num_outer[i], comp_num_inner[i],1])
    m1_mean[i,1] <- inprod(theta_M1_t1[comp_num_outer[i], comp_num_inner[i], 1:K_M_t1],
                           M_design_mat_t1[i, 1:K_M_t1]) + bi_M1[i]
    m1_sd[i,1] <- sqrt(sigma_sq_M1_t[comp_num_outer[i],comp_num_inner[i],1])
    data_wide_M1_pos[i,1] ~ dnorm(m1_mean[i,1], sd = m1_sd[i,1])
    data_wide_M1[i,1] ~ dnorm(Z_latent_M1[i,1]*data_wide_M1_pos[i,1], sd = 0.25)

    # M2 at t=1 (hurdle model)
    Z_latent_M2[i,1] ~ dbern(1 - pi0_M2_t[comp_num_outer[i], comp_num_inner[i],1])
    m2_mean[i,1] <- inprod(theta_M2_t1[comp_num_outer[i], comp_num_inner[i], 1:K_M_t1],
                           M_design_mat_t1[i, 1:K_M_t1]) + bi_M2[i]
    m2_sd[i,1] <- sqrt(sigma_sq_M2_t[comp_num_outer[i],comp_num_inner[i],1])
    data_wide_M2_pos[i,1] ~ dnorm(m2_mean[i,1], sd = m2_sd[i,1])
    data_wide_M2[i,1] ~ dnorm(Z_latent_M2[i,1]*data_wide_M2_pos[i,1], sd = 0.25)

    # ========== t = 2 ==========
    z_mean[i,2] <- inprod(theta_Z_t[comp_num_outer[i], comp_num_inner[i], 1:K, 2],
                          L0_mat[i, 1:K]) + bi_Z[i]
    data_wide_Z[i,2] ~ dbern(pnorm(z_mean[i,2],0,1))

    d_mean[i,2] <- inprod(theta_D_t2[comp_num_outer[i], comp_num_inner[i], 1:K_D_t2],
                          D_design_mat_t2[i, 1:K_D_t2]) + bi_D[i]
    data_wide_D[i,2] ~ dbern(pnorm(d_mean[i,2],0,1))

    Z_latent_M1[i,2] ~ dbern(1 - pi0_M1_t[comp_num_outer[i], comp_num_inner[i],2])
    m1_mean[i,2] <- inprod(theta_M1_t2[comp_num_outer[i], comp_num_inner[i], 1:K_M_t2],
                           M_design_mat_t2[i, 1:K_M_t2]) + bi_M1[i]
    m1_sd[i,2] <- sqrt(sigma_sq_M1_t[comp_num_outer[i],comp_num_inner[i],2])
    data_wide_M1_pos[i,2] ~ dnorm(m1_mean[i,2], sd = m1_sd[i,2])
    data_wide_M1[i,2] ~ dnorm(Z_latent_M1[i,2]*data_wide_M1_pos[i,2], sd = 0.25)

    Z_latent_M2[i,2] ~ dbern(1 - pi0_M2_t[comp_num_outer[i], comp_num_inner[i],2])
    m2_mean[i,2] <- inprod(theta_M2_t2[comp_num_outer[i], comp_num_inner[i], 1:K_M_t2],
                           M_design_mat_t2[i, 1:K_M_t2]) + bi_M2[i]
    m2_sd[i,2] <- sqrt(sigma_sq_M2_t[comp_num_outer[i],comp_num_inner[i],2])
    data_wide_M2_pos[i,2] ~ dnorm(m2_mean[i,2], sd = m2_sd[i,2])
    data_wide_M2[i,2] ~ dnorm(Z_latent_M2[i,2]*data_wide_M2_pos[i,2], sd = 0.25)

    # ========== t = 3 ==========
    z_mean[i,3] <- inprod(theta_Z_t[comp_num_outer[i], comp_num_inner[i], 1:K, 3],
                          L0_mat[i, 1:K]) + bi_Z[i]
    data_wide_Z[i,3] ~ dbern(pnorm(z_mean[i,3],0,1))

    d_mean[i,3] <- inprod(theta_D_t3[comp_num_outer[i], comp_num_inner[i], 1:K_D_t3],
                          D_design_mat_t3[i, 1:K_D_t3]) + bi_D[i]
    data_wide_D[i,3] ~ dbern(pnorm(d_mean[i,3],0,1))

    Z_latent_M1[i,3] ~ dbern(1 - pi0_M1_t[comp_num_outer[i], comp_num_inner[i],3])
    m1_mean[i,3] <- inprod(theta_M1_t3[comp_num_outer[i], comp_num_inner[i], 1:K_M_t3],
                           M_design_mat_t3[i, 1:K_M_t3]) + bi_M1[i]
    m1_sd[i,3] <- sqrt(sigma_sq_M1_t[comp_num_outer[i],comp_num_inner[i],3])
    data_wide_M1_pos[i,3] ~ dnorm(m1_mean[i,3], sd = m1_sd[i,3])
    data_wide_M1[i,3] ~ dnorm(Z_latent_M1[i,3]*data_wide_M1_pos[i,3], sd = 0.25)

    Z_latent_M2[i,3] ~ dbern(1 - pi0_M2_t[comp_num_outer[i], comp_num_inner[i],3])
    m2_mean[i,3] <- inprod(theta_M2_t3[comp_num_outer[i], comp_num_inner[i], 1:K_M_t3],
                           M_design_mat_t3[i, 1:K_M_t3]) + bi_M2[i]
    m2_sd[i,3] <- sqrt(sigma_sq_M2_t[comp_num_outer[i],comp_num_inner[i],3])
    data_wide_M2_pos[i,3] ~ dnorm(m2_mean[i,3], sd = m2_sd[i,3])
    data_wide_M2[i,3] ~ dnorm(Z_latent_M2[i,3]*data_wide_M2_pos[i,3], sd = 0.25)

    # Outcome likelihood (hurdle model)
    Z_latent_Y[i] ~ dbern(1 - pi0_Y[comp_num_outer[i]])
    y_mean[i] <- inprod(beta_Yi[comp_num_outer[i],1:num_Y_reg_coeff],
                                         Y_reg_design_mat[i, 1:num_Y_reg_coeff])
    y_sd[i] <-  sqrt(sigma_sq_Yi[comp_num_outer[i]])
    Y_pos[i] ~ dnorm(y_mean[i], sd =  y_sd[i])
    Y[i] ~ dnorm(Z_latent_Y[i]*Y_pos[i], sd = 0.25)

  }
})


###############################################################################
## 6. Set up constants, data, and initial values
###############################################################################

num_outer_cluster_N = 10; num_inner_cluster_M = 4
L0_mat = as.matrix(L0_df_wide)
num_L0 =3   # one intercept and 2 L0
num_Y_reg_coeff = ncol(Y_reg_design_mat) +1 # +1 for intercept
alpha_thet_const = 0.5

constants <- list(alpha_bet = 1,
                  alpha_thet = rep(alpha_thet_const, num_outer_cluster_N),
                  N = num_outer_cluster_N,
                  M = num_inner_cluster_M,
                  n = length(unique(digCom_workdf_wide$ID)),
                  K = num_L0,
                  max_t = 3,
                  num_Y_reg_coeff = num_Y_reg_coeff,
                  num_bin_L0 =2,
                  K_D_t1 = ncol(D_fit_design_mat_t1), K_D_t2 = ncol(D_fit_design_mat_t2), K_D_t3 = ncol(D_fit_design_mat_t3),
                  K_M_t1 = ncol(M_fit_design_mat_t1), K_M_t2 = ncol(M_fit_design_mat_t2), K_M_t3 = ncol(M_fit_design_mat_t3)
)

data <- list(Y = digCom_workdf_wide$Y_t_3,
             Y_reg_design_mat = cbind(L0_df_wide$Intercept, Y_reg_design_mat),
             mu_betaY_prior = unlist(lin_reg_YvAllCovs[1]),
             sigSq_betaY_prior= unlist(lin_reg_YvAllCovs[2]),
             sigSq_Y_prior = unlist(lin_reg_YvAllCovs[3]),
             Z_latent_Y = ifelse(digCom_workdf_wide$Y_t_3 ==0, 0,1),

             L0_mat = L0_mat,
             L0_bin = L0_df_wide[2:3],

             data_wide_M2 = longitudinal_data_wide_M2,
             mu_theta_M2_t1_prior = mle_coeff_M2_t1,
             mu_theta_M2_t2_prior = mle_coeff_M2_t2,
             mu_theta_M2_t3_prior = mle_coeff_M2_t3,
             sigSq_theta_M2_t1_prior = mle_sig_sq_M2_t1,
             sigSq_theta_M2_t2_prior = mle_sig_sq_M2_t2,
             sigSq_theta_M2_t3_prior = mle_sig_sq_M2_t3,
             sigma_sq_hat_M2_t = c(sigma_sq_hat_M2_t1, sigma_sq_hat_M2_t2, sigma_sq_hat_M2_t3),
             Z_latent_M2 = ifelse(longitudinal_data_wide_M2 ==0, 0,1),
             M_design_mat_t1 = as.matrix(M_fit_design_mat_t1),
             M_design_mat_t2 = as.matrix(M_fit_design_mat_t2),
             M_design_mat_t3 = as.matrix(M_fit_design_mat_t3),

             data_wide_M1 = longitudinal_data_wide_M1,
             mu_theta_M1_t1_prior = mle_coeff_M1_t1,
             mu_theta_M1_t2_prior = mle_coeff_M1_t2,
             mu_theta_M1_t3_prior = mle_coeff_M1_t3,
             sigSq_theta_M1_t1_prior = mle_sig_sq_M1_t1,
             sigSq_theta_M1_t2_prior = mle_sig_sq_M1_t2,
             sigSq_theta_M1_t3_prior = mle_sig_sq_M1_t3,
             sigma_sq_hat_M1_t = c(sigma_sq_hat_M1_t1, sigma_sq_hat_M1_t2, sigma_sq_hat_M1_t3),
             Z_latent_M1 = ifelse(longitudinal_data_wide_M1 ==0, 0,1),

             data_wide_D = longitudinal_data_wide_D,
             mu_theta_D_t1_prior = mle_coeff_D_t1,
             mu_theta_D_t2_prior = mle_coeff_D_t2,
             mu_theta_D_t3_prior = mle_coeff_D_t3,
             sigSq_theta_D_t1_prior = mle_sig_sq_D_t1,
             sigSq_theta_D_t2_prior = mle_sig_sq_D_t2,
             sigSq_theta_D_t3_prior = mle_sig_sq_D_t3,
             D_design_mat_t1 = as.matrix(D_fit_design_mat_t1),
             D_design_mat_t2 = as.matrix(D_fit_design_mat_t2),
             D_design_mat_t3 = as.matrix(D_fit_design_mat_t3),

             data_wide_Z = longitudinal_data_wide_Z,
             mu_theta_Z_t_prior = mle_coeff_Z,
             sigSq_theta_Z_t_prior= mle_sig_sq_Z
)

inits <- list(
  # Outcome model parameters
  beta_Yi = matrix(0, nrow = constants$N, ncol = constants$num_Y_reg_coeff),
  sigma_sq_Yi = rinvgamma(constants$N, 1, 1),
  pi0_Y = rep(0.5, constants$N),
  Y_pos = rep(mean(digCom_workdf_wide$Y_t_3), constants$n),

  # Mediator M2 parameters
  bi_M2_sig_sq = rinvgamma(1, 1, 1),
  bi_M2 = rep(0,constants$n),
  theta_M2_t1 = array(0, c(constants$N, constants$M, constants$K_M_t1)),
  theta_M2_t2 = array(0, c(constants$N, constants$M, constants$K_M_t2)),
  theta_M2_t3 = array(0, c(constants$N, constants$M, constants$K_M_t3)),
  sigma_sq_M2_t = array(rinvgamma(constants$N * constants$M * constants$max_t, 1, 1),
                        c(constants$N, constants$M, constants$max_t)),
  pi0_M2_t = array(0.5, c(constants$N, constants$M, constants$max_t)),
  data_wide_M2_pos = matrix(rep(colMeans(longitudinal_data_wide_M2), constants$n),
                            nrow = constants$n, ncol = constants$max_t, byrow = TRUE),

  # Mediator M1 parameters
  bi_M1_sig_sq = rinvgamma(1, 1, 1),
  bi_M1 = rep(0,constants$n),
  theta_M1_t1 = array(0, c(constants$N, constants$M, constants$K_M_t1)),
  theta_M1_t2 = array(0, c(constants$N, constants$M, constants$K_M_t2)),
  theta_M1_t3 = array(0, c(constants$N, constants$M, constants$K_M_t3)),
  sigma_sq_M1_t = array(rinvgamma(constants$N * constants$M * constants$max_t, 1, 1),
                        c(constants$N, constants$M, constants$max_t)),
  pi0_M1_t = array(0.5, c(constants$N, constants$M, constants$max_t)),
  data_wide_M1_pos = matrix(rep(colMeans(longitudinal_data_wide_M1), constants$n),
                            nrow = constants$n, ncol = constants$max_t, byrow = TRUE),

  # Treatment receipt D parameters
  bi_D_sig_sq = rinvgamma(1, 1, 1),
  bi_D = rep(0,constants$n),
  theta_D_t1 = array(0, c(constants$N, constants$M, constants$K_D_t1)),
  theta_D_t2 = array(0, c(constants$N, constants$M, constants$K_D_t2)),
  theta_D_t3 = array(0, c(constants$N, constants$M, constants$K_D_t3)),

  # Treatment assignment Z parameters
  bi_Z_sig_sq = rinvgamma(1, 1, 1),
  bi_Z = rep(0,constants$n),
  theta_Z_t = array(0,
                    c(constants$N, constants$M, constants$K, constants$max_t)),

  # Baseline covariates
  pr_Bin_L0i = array(0.5, c(constants$N, constants$M, constants$num_bin_L0)),

  # EDP mixture parameters
  comp_num_outer = sample(1:constants$N, size = constants$n, replace = TRUE),
  comp_num_inner = sample(1:constants$M, size = constants$n, replace = TRUE),
  breaks_r = rbeta(constants$N, 1, 1),
  breaks_sr = matrix(rbeta(constants$N*constants$M, 1, 1),
                     nrow = constants$N, ncol = constants$M)
)


###############################################################################
## 7. Build and compile the model
###############################################################################

DiscEDPLong_model <- nimbleModel(code, constants, data, inits)

compile_DiscEDPLong_model <- compileNimble(DiscEDPLong_model, showCompilerOutput = TRUE)


###############################################################################
## 8. Configure MCMC samplers
###############################################################################

config_MCMC <- configureMCMC(DiscEDPLong_model, useConjugacy = TRUE,
                             enableWAIC = TRUE,
                             monitors = c("beta_Yi", "sigma_sq_Yi", "pi0_Y",
                                          "bi_M2_sig_sq",
                                          "theta_M2_t1", "theta_M2_t2", "theta_M2_t3",
                                          "sigma_sq_M2_t","pi0_M2_t",
                                          "bi_M1_sig_sq",
                                          "theta_M1_t1", "theta_M1_t2", "theta_M1_t3",
                                          "sigma_sq_M1_t","pi0_M1_t",
                                          "bi_D_sig_sq",
                                          "theta_D_t1", "theta_D_t2", "theta_D_t3",
                                          "bi_Z_sig_sq","theta_Z_t",
                                          "pr_Bin_L0i",
                                          "breaks_r", "breaks_sr",
                                          "ksi_r","ksi_sr",
                                          "logDens"))

# Replace default logDens sampler with custom sum-posterior-log-density sampler
config_MCMC$removeSamplers('logDens')
config_MCMC$addSampler(target = 'logDens', type = 'sumLogPostDens')

# Custom RW samplers for inner stick-breaking weights (breaks_sr)
config_MCMC$removeSamplers("breaks_sr")

for (k in 1:num_outer_cluster_N) {
  for (j in 1:(num_inner_cluster_M-1)) {
    # Count observations in outer cluster k, inner cluster j
    n_kj <- sum(DiscEDPLong_model$comp_num_outer == k &
                  DiscEDPLong_model$comp_num_inner == j)
    sum_n_kh <- 0
    for (h in (j+1):num_inner_cluster_M){
      sum_n_kh <- sum_n_kh + sum(DiscEDPLong_model$comp_num_outer == k &
                                   DiscEDPLong_model$comp_num_inner == h)
    }
    config_MCMC$addSampler(
      target = paste0("breaks_sr[", k, ",", j, "]"),
      type = "RW",
      control = list(
        posterior = "dbeta",
        shape1 = n_kj + 1,
        shape2 = alpha_thet_const + sum_n_kh
      )
    )
  }
}


###############################################################################
## 9. Build, compile, and run MCMC
###############################################################################

mcmc <- buildMCMC(config_MCMC)
cmcmc <- compileNimble(mcmc, project = DiscEDPLong_model, showCompilerOutput = TRUE)

num_iter = 20000;  num_burnin = 17500;  num_thin = 10;  num_chains= 4

samples_allChains <- runMCMC(cmcmc, niter = num_iter, nburnin = num_burnin,
                             nchains = num_chains,
                             setSeed = TRUE, thin = num_thin, WAIC = TRUE)

posterior_samples <- samples_allChains$samples   # list of nchains matrices

summary(samples_allChains)

samples_allChains$WAIC

saveRDS(posterior_samples, file="param_post1000_N10M4_AllData_samples.RDS")


###############################################################################
## 10. Trace plots for convergence diagnostics
###############################################################################

# Combine all chains into long-format dataframe for plotting
n_rows = (num_iter - num_burnin)/num_thin
posterior_samples_long <- bind_rows(
  lapply(names(posterior_samples), function(chain_name) {
    as.data.frame(posterior_samples[[chain_name]]) %>%
      mutate(Iteration = 1:n_rows, Chain = chain_name)
  }),
  .id = "ChainID"
)

# Helper: 4-chain color palette
chain_colors <- c("chain1" = "blue", "chain2" = "red",
                  "chain3" = "green", "chain4" = "purple")


## 10a. Sum-posterior-log-density traceplot
columns_of_interest <- paste0("logDens")
posterior_samples_long_filtered <- posterior_samples_long %>%
  select(Iteration, Chain, all_of(columns_of_interest)) %>%
  pivot_longer(cols = starts_with("logDens"), names_to = "Parameter", values_to = "Value")

traceplot_post_log_density <- ggplot(posterior_samples_long_filtered, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = chain_colors) +
  theme_minimal() +
  labs(title = "Traceplots for sum-posterior-log-density", color = "Chain",
       x = "Iteration", y = "Value") +
  theme(strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("traceplot_post_log_density.png")


## 10b. Outer cluster weights (ksi_r) traceplot
columns_of_interest <- paste0("ksi_r[", 1:num_outer_cluster_N, "]")
posterior_samples_long_filtered <- posterior_samples_long %>%
  select(Iteration, Chain, all_of(columns_of_interest)) %>%
  pivot_longer(cols = starts_with("ksi_r"), names_to = "Parameter", values_to = "Value")

traceplot_ksi_r <- ggplot(posterior_samples_long_filtered, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = chain_colors) +
  theme_minimal() +
  labs(title = "Traceplots for ksi_r", color = "Chain",
       x = "Iteration", y = "Value") +
  theme(strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("traceplot_ksi_r.png")


## 10c. Inner cluster weights (ksi_sr) traceplots — one per outer cluster
for(i in 1:num_outer_cluster_N){
  columns_of_interest <- paste0("ksi_sr[",i,", ",1:num_inner_cluster_M,"]")
  posterior_samples_long_filtered <- posterior_samples_long %>%
    select(Iteration, Chain, all_of(columns_of_interest)) %>%
    pivot_longer(cols = starts_with("ksi_sr"), names_to = "Parameter", values_to = "Value")

  traceplot_ksi_sr <- ggplot(posterior_samples_long_filtered, aes(x = Iteration, y = Value, color = Chain)) +
    geom_line() +
    facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
    scale_color_manual(values = chain_colors) +
    theme_minimal() +
    labs(title = "Traceplots for ksi_sr", color = "Chain",
         x = "Iteration", y = "Value") +
    theme(strip.text = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))

  ggsave(filename = paste0("traceplot_ksi_sr", i, ".png"), plot = traceplot_ksi_sr)
}


## 10d. Outcome hurdle probability (pi0_Y) traceplot
columns_of_interest <- paste0("pi0_Y[", 1:num_outer_cluster_N, "]")
posterior_samples_long_filtered <- posterior_samples_long %>%
  select(Iteration, Chain, all_of(columns_of_interest)) %>%
  pivot_longer(cols = starts_with("pi0_Y"), names_to = "Parameter", values_to = "Value")

traceplot_pi0_Y <- ggplot(posterior_samples_long_filtered, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = chain_colors) +
  theme_minimal() +
  labs(title = "Traceplots for pi0_Y", color = "Chain",
       x = "Iteration", y = "Value") +
  theme(strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("traceplot_pi0_Y.png")


## 10e. Outcome regression coefficients (beta_Yi) traceplots — one per outer cluster
for(i in 1:num_outer_cluster_N){
  columns_of_interest <- paste0("beta_Yi[",i,", ",1:num_Y_reg_coeff,"]")
  posterior_samples_long_filtered <- posterior_samples_long %>%
    select(Iteration, Chain, all_of(columns_of_interest)) %>%
    pivot_longer(cols = starts_with("beta_Yi"), names_to = "Parameter", values_to = "Value")

  traceplot_beta_Yi <- ggplot(posterior_samples_long_filtered, aes(x = Iteration, y = Value, color = Chain)) +
    geom_line() +
    facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
    scale_color_manual(values = chain_colors) +
    theme_minimal() +
    labs(title = "Traceplots for beta_Yi", color = "Chain",
         x = "Iteration", y = "Value") +
    theme(strip.text = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))

  ggsave(filename = paste0("traceplot_beta_Yi", i, ".png"), plot = traceplot_beta_Yi)
}


## 10f. Outcome variance (sigma_sq_Yi) traceplot
columns_of_interest <- paste0("sigma_sq_Yi[", 1:num_outer_cluster_N, "]")
posterior_samples_long_filtered <- posterior_samples_long %>%
  select(Iteration, Chain, all_of(columns_of_interest)) %>%
  pivot_longer(cols = starts_with("sigma_sq_Yi"), names_to = "Parameter", values_to = "Value")

traceplot_sigma_sq_Yi <- ggplot(posterior_samples_long_filtered, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 2) +
  scale_color_manual(values = chain_colors) +
  theme_minimal() +
  labs(title = "Traceplots for sigma_sq_Yi", color = "Chain",
       x = "Iteration", y = "Value") +
  theme(strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("traceplot_sigma_sq_Yi.png")
