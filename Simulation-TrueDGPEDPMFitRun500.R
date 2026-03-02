#-------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------Prepare for simulation----------------------------------------------#

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


txt.title = paste0("EDPM_simu_results/Results_Nimble_EDPM_Models_result.txt")
if (run_ID == 1) {
  df = data.frame(matrix(ncol = 38, nrow = 0))                 #each row contains results from each replication/dataset/run
  df_col_names = c("run_ID",
                   "PIDE_AC_post_mean", "PIDE_AC_lowerCI", "PIDE_AC_upperCI",
                   "PJIIE_AC_post_mean", "PJIIE_AC_lowerCI", "PJIIE_AC_upperCI",
                   "PCE_AC_post_mean", "PCE_AC_lowerCI", "PCE_AC_upperCI",
                   "PIDE_VA_post_mean", "PIDE_VA_lowerCI", "PIDE_VA_upperCI",
                   "PJIIE_VA_post_mean", "PJIIE_VA_lowerCI", "PJIIE_VA_upperCI",
                   "PCE_VA_post_mean", "PCE_VA_lowerCI", "PCE_VA_upperCI",
                   "PIDE_PA_post_mean", "PIDE_PA_lowerCI", "PIDE_PA_upperCI",
                   "PJIIE_PA_post_mean", "PJIIE_PA_lowerCI", "PJIIE_PA_upperCI",
                   "PCE_PA_post_mean", "PCE_PA_lowerCI", "PCE_PA_upperCI",
                   "PIDE_NAC_post_mean", "PIDE_NAC_lowerCI", "PIDE_NAC_upperCI",
                   "PJIIE_NAC_post_mean", "PJIIE_NAC_lowerCI", "PJIIE_NAC_upperCI",
                   "PCE_NAC_post_mean", "PCE_NAC_lowerCI", "PCE_NAC_upperCI",
                   "num_outer_clusters_used")
  colnames(df) = df_col_names
  write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
}


#-------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------Fit Models in the observed dataset-------------------------------------------------#


set.seed(1234)   #this will be changed later during data generation

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(truncnorm)
library(mvtnorm)

library(nimble, warn.conflicts = FALSE)

# Enable Nimble options for saving intermediates and verbose output
nimbleOptions(saveIntermediates = TRUE)  # Saves generated C++ code files
nimbleOptions(verbose = TRUE)           # Provides detailed compilation output


#############################################################################
#############Subset a smaller data for analysis##############################


digCom_data <- read.csv("non-comp_dataset.csv")
digCom_df = data.frame(digCom_data)
#############################################################################
#############################################################################


#-------------------------------------------------------------------------------------------------------------
#------------------------------------------Data in long format------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
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


# change GENDER elements from M to 1,and others to 0.
digCom_workdf <- digCom_workdf %>%
  mutate(GENDER = ifelse(GENDER == "M", 1, 0))


# change BASELINE_AGE elements from >=4 to 1,and <4 to 0.
digCom_workdf <- digCom_workdf %>%
  mutate(BASELINE_AGE = ifelse(BASELINE_AGE >= 4, 1, 0))


#-------------------------------------------------------------------------------------------------------------
#------------------------------------------Data in wide format------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
setDT(digCom_workdf)

digCom_workdf_wide <- dcast(
  digCom_workdf,
  ID + GENDER + BASELINE_AGE ~ t,
  value.var = c("Z_t", "D_t", "M1_t", "M2_t", "Y_t"),
  fun.aggregate = mean   # or unique, depending on your data
)
head(digCom_workdf_wide)

digCom_workdf_wide = digCom_workdf_wide[!duplicated(ID)]  #remove any duplicated ID rows in the wide format data

digCom_workdf_wide = na.omit(digCom_workdf_wide) #deleting the some rows when converting data from long to wide


cat("The total sample size for the analysis is:", length(unique(digCom_workdf_wide$ID)), "\n")


#-------------------------------------------------------------------------------------------------------------
#------------------------------------------Input data---------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

#Remove rows with a missing column


#L0 df order columns by intercept, binary covariates, continuous covariates
L0_intercept = rep(1,length(digCom_workdf_wide$ID))
L0_df_wide = data.frame(Intercept = L0_intercept,
                        GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE) 


Y_vec_wide = digCom_workdf_wide$Y_t_3


###############################################################################
##########Define design matrices for parametric regression#####################
###############################################################################

param_reg_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                           Z_t_1 = digCom_workdf_wide$Z_t_1, D_t_1 = digCom_workdf_wide$D_t_1, M1_t_1 = digCom_workdf_wide$M1_t_1, M2_t_1 = digCom_workdf_wide$M2_t_1,
                           Z_t_2 = digCom_workdf_wide$Z_t_2, D_t_2 = digCom_workdf_wide$D_t_2, M1_t_2 = digCom_workdf_wide$M1_t_2, M2_t_2 = digCom_workdf_wide$M2_t_2,
                           Z_t_3 = digCom_workdf_wide$Z_t_3, D_t_3 = digCom_workdf_wide$D_t_3, M1_t_3 = digCom_workdf_wide$M1_t_3, M2_t_3 = digCom_workdf_wide$M2_t_3,
                           Y_t_3 = digCom_workdf_wide$Y_t_3)


Y_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                    Z_t_3 = digCom_workdf_wide$Z_t_3, 
                                    D_t_3 = digCom_workdf_wide$D_t_3,
                                    M1_t_3 = digCom_workdf_wide$M1_t_3, 
                                    M2_t_3 = digCom_workdf_wide$M2_t_3)


M2t3_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                       M2_t_2 = digCom_workdf_wide$M2_t_2,
                                       Z_t_3 = digCom_workdf_wide$Z_t_3,
                                       D_t_3 = digCom_workdf_wide$D_t_3,
                                       M1_t_3 = digCom_workdf_wide$M1_t_3)


M1t3_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                       M1_t_2 = digCom_workdf_wide$M1_t_2,
                                       M2_t_2 = digCom_workdf_wide$M2_t_2,
                                       Z_t_3 = digCom_workdf_wide$Z_t_3,
                                       D_t_3 = digCom_workdf_wide$D_t_3)

Dt3_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                      D_t_2 = digCom_workdf_wide$D_t_2,
                                      M1_t_2 = digCom_workdf_wide$M1_t_2,
                                      M2_t_2 = digCom_workdf_wide$M2_t_2,
                                      Z_t_3 = digCom_workdf_wide$Z_t_3)


Zt3_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                      Z_t_2 = digCom_workdf_wide$Z_t_2,
                                      D_t_2 = digCom_workdf_wide$D_t_2,
                                      M1_t_2 = digCom_workdf_wide$M1_t_2,
                                      M2_t_2 = digCom_workdf_wide$M2_t_2)


M2t2_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                       M2_t_1 = digCom_workdf_wide$M2_t_1,
                                       Z_t_2 = digCom_workdf_wide$Z_t_2, 
                                       D_t_2 = digCom_workdf_wide$D_t_2,
                                       M1_t_2 = digCom_workdf_wide$M1_t_2)


M1t2_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                       M1_t_1 = digCom_workdf_wide$M1_t_1,
                                       M2_t_1 = digCom_workdf_wide$M2_t_1,
                                       Z_t_2 = digCom_workdf_wide$Z_t_2,
                                       D_t_2 = digCom_workdf_wide$D_t_2)


Dt2_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                      D_t_1 = digCom_workdf_wide$D_t_1,
                                      M1_t_1 = digCom_workdf_wide$M1_t_1, 
                                      M2_t_1 = digCom_workdf_wide$M2_t_1,
                                      Z_t_2 = digCom_workdf_wide$Z_t_2)


Zt2_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                      Z_t_1 = digCom_workdf_wide$Z_t_1, 
                                      D_t_1 = digCom_workdf_wide$D_t_1,
                                      M1_t_1 = digCom_workdf_wide$M1_t_1,
                                      M2_t_1 = digCom_workdf_wide$M2_t_1)


M2t1_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                       Z_t_1 = digCom_workdf_wide$Z_t_1,
                                       D_t_1 = digCom_workdf_wide$D_t_1, 
                                       M1_t_1 = digCom_workdf_wide$M1_t_1)


M1t1_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                       Z_t_1 = digCom_workdf_wide$Z_t_1,
                                       D_t_1 = digCom_workdf_wide$D_t_1)


Dt1_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE,
                                      Z_t_1 = digCom_workdf_wide$Z_t_1)


Zt1_param_reg_design_mat = data.frame(GENDER = digCom_workdf_wide$GENDER, BASELINE_AGE = digCom_workdf_wide$BASELINE_AGE)


longitudinal_data_wide_M2 = data.frame(M2_t_1 = digCom_workdf_wide$M2_t_1,
                                       M2_t_2 = digCom_workdf_wide$M2_t_2,
                                       M2_t_3 = digCom_workdf_wide$M2_t_3)


longitudinal_data_wide_M1 = data.frame(M1_t_1 = digCom_workdf_wide$M1_t_1, 
                                       M1_t_2 = digCom_workdf_wide$M1_t_2,
                                       M1_t_3 = digCom_workdf_wide$M1_t_3)


longitudinal_data_wide_D = data.frame(D_t_1 = digCom_workdf_wide$D_t_1,
                                      D_t_2 = digCom_workdf_wide$D_t_2,
                                      D_t_3 = digCom_workdf_wide$D_t_3)


longitudinal_data_wide_Z = data.frame(Z_t_1 = digCom_workdf_wide$Z_t_1,
                                      Z_t_2 = digCom_workdf_wide$Z_t_2,
                                      Z_t_3 = digCom_workdf_wide$Z_t_3)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
num_Bin_L0 = 2  # both gender and baseline age are binary
# sample frequency of each Binary column of L0_df_wide
pr_theta_Bin_L0 = as.numeric(colMeans(L0_df_wide[2:3], na.rm = TRUE))


# sample means of each Continuous column of L0_df_wide
# sample variances of each Continuous column of L0_df_wide


###############################################################################                                                  
##############################Nimble model#####################################
###############################################################################                                                


#Record model sum-posterior-log-density

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


# Define a custom distribution for zero-inflated hurdle normal
dzinhurdlenorm <- nimbleFunction(
  run = function(x = double(0), pi = double(0), mu = double(0), sigma = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if (sigma <= 0) stop("Sigma must be positive")
    logProb = 0
    if (x == 0) {
      logProb = log(max(1e-10, pi))  # Avoid underflow for pi
    } else {
      # Handle truncation probability
      truncProb = max(1e-10, 1 - pnorm(0, mean = mu, sd = sigma))
      logProb = log(max(1e-10, 1 - pi)) +
        dnorm(x, mean = mu, sd = sigma, log = TRUE) -
        log(truncProb)
    }
    if (log) return(logProb)
    return(exp(logProb))
  }
)


rzinhurdlenorm <- nimbleFunction(
  run = function(n = integer(0), pi = double(0), mu = double(0), sigma = double(0)) {
    returnType(double(0))
    if (n != 1) stop("rzinhurdlenorm only allows n = 1")
    if (runif(1) < pi) {
      return(0)  # Zero-inflated part
    } else {
      x <- -1  # Initialize x to an invalid value for the truncation
      while (x <= 0) {
        x <- rnorm(1, mean = mu, sd = sigma)  # Normal component
      }
      return(x)  # Return the valid sample
    }
  }
)

# Test dzinhurdlenorm
dzinhurdlenorm(0, pi = 0.3, mu = 2, sigma = 1, log = TRUE)

# Test rzinhurdlenorm
rzinhurdlenorm(1, pi = 0.3, mu = 2, sigma = 1)


# Register the Distribution in Nimble
registerDistributions(list(
  dzinhurdlenorm = list(
    BUGSdist = "dzinhurdlenorm(pi, mu, sigma)",
    types = c('value = double(0)', 'pi = double(0)', 'mu = double(0)', 'sigma = double(0)'),
    discrete = FALSE  # Set to TRUE only if the distribution produces integers
  )
))


#######################################################################################
##########################Nimble Code Function##################################################
#######################################################################################

code_DG <- nimbleCode({
  logDens ~ dnorm(0, 1)    ## this distribution does not matter: only need to compute model sum-posterior-log-density
  
  ##################################################################################
  ##############################Prior distributions#################################
  ##################################################################################
  
  
  ##############################Hurdle probabilities################################
  pi0_Y ~ dbeta(1,1)
  pi0_M1_t3~ dbeta(1,1)
  pi0_M2_t3~ dbeta(1,1)
  pi0_M1_t2~ dbeta(1,1)
  pi0_M2_t2~ dbeta(1,1)
  pi0_M1_t1~ dbeta(1,1)
  pi0_M2_t1~ dbeta(1,1)
  
  
  #####################Regression fixed effects coefficients################################
  # AT t =3
  #########Do this if you want to force coeffs for Z and D in outcome regression to be positive########
  # for (k in 1:num_Y_reg_coeff) {
  #   beta_raw[k] ~ dnorm(0, var = 10)   # unconstrained
  #   # apply positivity transformation 
  #   beta_Yi[k] <- (1 - pos_index[k]) * beta_raw[k] +
  #     pos_index[k] * abs(beta_raw[k])
  # }
  for(k in 1:num_Y_reg_coeff){beta_Yi[k]~ dnorm(mean = 0, var = 10)}  
  for(k in 1:num_M2t3_reg_coeff){theta_M2_t3[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_M1t3_reg_coeff){theta_M1_t3[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_Dt3_reg_coeff){theta_D_t3[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_Zt3_reg_coeff){theta_Z_t3[k]~ dnorm(mean = 0, var = 10)} 
  
  # AT t =2
  for(k in 1:num_M2t2_reg_coeff){theta_M2_t2[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_M1t2_reg_coeff){theta_M1_t2[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_Dt2_reg_coeff){theta_D_t2[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_Zt2_reg_coeff){theta_Z_t2[k]~ dnorm(mean = 0, var = 10)} 
  
  # AT t =1
  for(k in 1:num_M2t1_reg_coeff){theta_M2_t1[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_M1t1_reg_coeff){theta_M1_t1[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_Dt1_reg_coeff){theta_D_t1[k]~ dnorm(mean = 0, var = 10)} 
  for(k in 1:num_Zt1_reg_coeff){theta_Z_t1[k]~ dnorm(mean = 0, var = 10)} 
  
  
  #####################Regression variance coefficients################################
  sigma_sq_Yi ~ dinvgamma(shape = 3, scale = 1)
  sigma_sq_M2_t3 ~ dinvgamma(shape = 3, scale = 1)
  sigma_sq_M1_t3 ~ dinvgamma(shape = 3, scale = 1)
  sigma_sq_M2_t2 ~ dinvgamma(shape = 3, scale = 1)
  sigma_sq_M1_t2 ~ dinvgamma(shape = 3, scale = 1)
  sigma_sq_M2_t1 ~ dinvgamma(shape = 3, scale = 1)
  sigma_sq_M1_t1 ~ dinvgamma(shape = 3, scale = 1)  
  
  
  ##################################################################################
  ##############################Likelihood##########################################
  ##################################################################################
  
  
  for(i in 1:n) {
    
    
    #AT t = 1  
    
    Zt1[i] ~ dbern(pnorm(inprod(theta_Z_t1[1:num_Zt1_reg_coeff], Zt1_param_reg_design_mat[i,1:num_Zt1_reg_coeff]),0,1))
    Dt1[i] ~ dbern(pnorm(inprod(theta_D_t1[1:num_Dt1_reg_coeff], Dt1_param_reg_design_mat[i,1:num_Dt1_reg_coeff]),0,1))  
    
    Z_latent_M1t1[i] ~ dbern(1 - pi0_M1_t1)  # Indicator for positive values 
    m1t1_mean[i] <-   inprod(theta_M1_t1[1:num_M1t1_reg_coeff], 
                             M1t1_param_reg_design_mat[i, 1:num_M1t1_reg_coeff]) 
    data_wide_M1t1_pos[i] ~ dnorm(m1t1_mean[i], sd = sqrt(sigma_sq_M1_t1))
    M1t1[i] ~ dnorm(Z_latent_M1t1[i]*data_wide_M1t1_pos[i], sd = 0.25) 
    
    Z_latent_M2t1[i] ~ dbern(1 - pi0_M2_t1)  # Indicator for positive values 
    m2t1_mean[i] <-   inprod(theta_M2_t1[1:num_M2t1_reg_coeff], 
                             M2t1_param_reg_design_mat[i, 1:num_M2t1_reg_coeff]) 
    data_wide_M2t1_pos[i] ~ dnorm(m2t1_mean[i], sd = sqrt(sigma_sq_M2_t1))
    M2t1[i] ~ dnorm(Z_latent_M2t1[i]*data_wide_M2t1_pos[i], sd = 0.25) 
    
    
    #AT t = 2    
    
    Zt2[i] ~ dbern(pnorm(inprod(theta_Z_t2[1:num_Zt2_reg_coeff], Zt2_param_reg_design_mat[i,1:num_Zt2_reg_coeff]),0,1)) 
    Dt2[i] ~ dbern(pnorm(inprod(theta_D_t2[1:num_Dt2_reg_coeff], Dt2_param_reg_design_mat[i,1:num_Dt2_reg_coeff]),0,1))    
    
    Z_latent_M1t2[i] ~ dbern(1 - pi0_M1_t2)  # Indicator for positive values 
    m1t2_mean[i] <-   inprod(theta_M1_t2[1:num_M1t2_reg_coeff], 
                             M1t2_param_reg_design_mat[i, 1:num_M1t2_reg_coeff]) 
    data_wide_M1t2_pos[i] ~ dnorm(m1t2_mean[i], sd = sqrt(sigma_sq_M1_t2))
    M1t2[i] ~ dnorm(Z_latent_M1t2[i]*data_wide_M1t2_pos[i], sd = 0.25) 
    
    Z_latent_M2t2[i] ~ dbern(1 - pi0_M2_t2)  # Indicator for positive values 
    m2t2_mean[i] <-   inprod(theta_M2_t2[1:num_M2t2_reg_coeff], 
                             M2t2_param_reg_design_mat[i, 1:num_M2t2_reg_coeff]) 
    data_wide_M2t2_pos[i] ~ dnorm(m2t2_mean[i], sd = sqrt(sigma_sq_M2_t2))
    M2t2[i] ~ dnorm(Z_latent_M2t2[i]*data_wide_M2t3_pos[i], sd = 0.25) 
    
    
    #AT t = 3 
    
    Zt3[i] ~ dbern(pnorm(inprod(theta_Z_t3[1:num_Zt3_reg_coeff], Zt3_param_reg_design_mat[i,1:num_Zt3_reg_coeff]),0,1)) 
    Dt3[i] ~ dbern(pnorm(inprod(theta_D_t3[1:num_Dt3_reg_coeff], Dt3_param_reg_design_mat[i,1:num_Dt3_reg_coeff]),0,1))  
    
    Z_latent_M1t3[i] ~ dbern(1 - pi0_M1_t3)  # Indicator for positive values 
    m1t3_mean[i] <-   inprod(theta_M1_t3[1:num_M1t3_reg_coeff], 
                             M1t3_param_reg_design_mat[i, 1:num_M1t3_reg_coeff]) 
    data_wide_M1t3_pos[i] ~ dnorm(m1t3_mean[i], sd = sqrt(sigma_sq_M1_t3))
    M1t3[i] ~ dnorm(Z_latent_M1t3[i]*data_wide_M1t3_pos[i], sd = 0.25) 
    
    Z_latent_M2t3[i] ~ dbern(1 - pi0_M2_t3)  # Indicator for positive values 
    m2t3_mean[i] <-   inprod(theta_M2_t3[1:num_M2t3_reg_coeff], 
                             M2t3_param_reg_design_mat[i, 1:num_M2t3_reg_coeff]) 
    data_wide_M2t3_pos[i] ~ dnorm(m2t3_mean[i], sd = sqrt(sigma_sq_M2_t3))
    M2t3[i] ~ dnorm(Z_latent_M2t3[i]*data_wide_M2t3_pos[i], sd = 0.25) 
    
    Z_latent_Y[i] ~ dbern(1 - pi0_Y)  # Indicator for positive values
    y_mean[i] <- inprod(beta_Yi[1:num_Y_reg_coeff], Y_param_reg_design_mat[i, 1:num_Y_reg_coeff]) 
    Y_pos[i] ~ dnorm(y_mean[i], sd =  sqrt(sigma_sq_Yi))
    Y[i] ~ dnorm(Z_latent_Y[i]*Y_pos[i], sd = 0.25)
    
    
  }
})


#######################################################################################
##############Fit Parametric model and save one iter only for data generation##########
#######################################################################################

# one outer cluster and one inner cluster for fitting parametric data
L0_mat = as.matrix(L0_df_wide)
num_L0 =3   # one intercept and 2 L0


constants_DG <- list(n = length(unique(digCom_workdf_wide$ID)), # Number of individuals
                     #K = num_L0,                                # one intercept and 2 L0
                     #max_t = 3,
                     num_Y_reg_coeff = ncol(Y_param_reg_design_mat) + 1, # plus 1 to include the intercept
                     num_M2t3_reg_coeff = ncol(M2t3_param_reg_design_mat) + 1, 
                     num_M1t3_reg_coeff = ncol(M1t3_param_reg_design_mat) + 1, 
                     num_Dt3_reg_coeff = ncol(Dt3_param_reg_design_mat) + 1, 
                     num_Zt3_reg_coeff = ncol(Zt3_param_reg_design_mat) + 1, 
                     num_M2t2_reg_coeff = ncol(M2t2_param_reg_design_mat) + 1,
                     num_M1t2_reg_coeff = ncol(M1t2_param_reg_design_mat) + 1,
                     num_Dt2_reg_coeff =ncol(Dt2_param_reg_design_mat) + 1, 
                     num_Zt2_reg_coeff =ncol(Zt2_param_reg_design_mat) + 1, 
                     num_M2t1_reg_coeff =ncol(M2t1_param_reg_design_mat) + 1, 
                     num_M1t1_reg_coeff =ncol(M1t1_param_reg_design_mat) + 1, 
                     num_Dt1_reg_coeff =ncol(Dt1_param_reg_design_mat) + 1, 
                     num_Zt1_reg_coeff =ncol(Zt1_param_reg_design_mat) + 1
)


data_DG <- list(
  # regression outcomes
  #### t=3
  Y = digCom_workdf_wide$Y_t_3,
  Z_latent_Y = ifelse(digCom_workdf_wide$Y_t_3 ==0, 0,1),
  M2t3 =  digCom_workdf_wide$M2_t_3,
  M1t3 = digCom_workdf_wide$M1_t_3,
  Z_latent_M2t3 = ifelse(digCom_workdf_wide$M2_t_3 ==0, 0,1),
  Z_latent_M1t3 = ifelse(digCom_workdf_wide$M1_t_3 ==0, 0,1),
  Dt3 = digCom_workdf_wide$D_t_3,
  Zt3 = digCom_workdf_wide$Z_t_3,
  
  
  #### t=2
  M2t2 =  digCom_workdf_wide$M2_t_2,
  M1t2 = digCom_workdf_wide$M1_t_2,
  Z_latent_M2t2 = ifelse(digCom_workdf_wide$M2_t_2 ==0, 0,1),
  Z_latent_M1t2 = ifelse(digCom_workdf_wide$M1_t_2 ==0, 0,1),
  Dt2 = digCom_workdf_wide$D_t_2,
  Zt2 = digCom_workdf_wide$Z_t_2,
  
  
  #### t=1
  M2t1 =  digCom_workdf_wide$M2_t_1,
  M1t1 = digCom_workdf_wide$M1_t_1,
  Z_latent_M2t1 = ifelse(digCom_workdf_wide$M2_t_1 ==0, 0,1),
  Z_latent_M1t1 = ifelse(digCom_workdf_wide$M1_t_1 ==0, 0,1),
  Dt1 = digCom_workdf_wide$D_t_1,
  Zt1 = digCom_workdf_wide$Z_t_1,
  
  
  # data for baseline confounders
  #L0_mat = L0_mat,
  #L0_bin = L0_df_wide[2:3],
  
  
  # regression design matrices
  Y_param_reg_design_mat = cbind(L0_df_wide$Intercept, Y_param_reg_design_mat),
  M2t3_param_reg_design_mat = cbind(L0_df_wide$Intercept, M2t3_param_reg_design_mat),
  M1t3_param_reg_design_mat = cbind(L0_df_wide$Intercept, M1t3_param_reg_design_mat),
  Dt3_param_reg_design_mat = cbind(L0_df_wide$Intercept, Dt3_param_reg_design_mat),
  Zt3_param_reg_design_mat = cbind(L0_df_wide$Intercept, Zt3_param_reg_design_mat),
  M2t2_param_reg_design_mat = cbind(L0_df_wide$Intercept, M2t2_param_reg_design_mat),
  M1t2_param_reg_design_mat = cbind(L0_df_wide$Intercept, M1t2_param_reg_design_mat),
  Dt2_param_reg_design_mat = cbind(L0_df_wide$Intercept, Dt2_param_reg_design_mat),
  Zt2_param_reg_design_mat = cbind(L0_df_wide$Intercept, Zt2_param_reg_design_mat),
  M2t1_param_reg_design_mat = cbind(L0_df_wide$Intercept, M2t1_param_reg_design_mat),
  M1t1_param_reg_design_mat = cbind(L0_df_wide$Intercept, M1t1_param_reg_design_mat),
  Dt1_param_reg_design_mat = cbind(L0_df_wide$Intercept, Dt1_param_reg_design_mat),
  Zt1_param_reg_design_mat = cbind(L0_df_wide$Intercept, Zt1_param_reg_design_mat)
  
  
)


inits_DG <- list(
  # outcome model parameters
  beta_Yi = rep(0, constants_DG$num_Y_reg_coeff),
  sigma_sq_Yi = rinvgamma(1, 1, 1),
  pi0_Y = 0.5,
  Y_pos = rep(mean(digCom_workdf_wide$Y_t_3), constants_DG$n),
  #Z_latent_Y = ifelse(digCom_workdf_wide$Y_t_3 ==0, 0,1),
  
  # mediator parameters at t =3
  #bi_M2_sig_sq = rinvgamma(1, 1, 1),
  #bi_M2 = rep(0,constants_DG$n),
  theta_M2_t3 = rep(0,constants_DG$num_M2t3_reg_coeff),
  theta_M1_t3 = rep(0,constants_DG$num_M1t3_reg_coeff),
  sigma_sq_M2_t3 = rinvgamma(1, 1, 1),
  sigma_sq_M1_t3 = rinvgamma(1, 1, 1),
  pi0_M2_t3 = 0.5,
  pi0_M1_t3 = 0.5,
  data_wide_M2t3_pos = rep(mean(digCom_workdf_wide$M2_t_3), constants_DG$n),
  data_wide_M1t3_pos = rep(mean(digCom_workdf_wide$M1_t_3), constants_DG$n),
  
  
  # mediator parameters at t =2
  theta_M2_t2 = rep(0,constants_DG$num_M2t2_reg_coeff),
  theta_M1_t2 = rep(0,constants_DG$num_M1t2_reg_coeff),
  #Sigma_M1M2_t2 = matrix(rinvgamma(4, 1, 1), nrow = 2, ncol = 2, byrow = TRUE),
  sigma_sq_M2_t2 = rinvgamma(1, 1, 1),
  sigma_sq_M1_t2 = rinvgamma(1, 1, 1),
  pi0_M2_t2 = 0.5,
  pi0_M1_t2 = 0.5,
  data_wide_M2t2_pos = rep(mean(digCom_workdf_wide$M2_t_2), constants_DG$n),
  data_wide_M1t2_pos = rep(mean(digCom_workdf_wide$M1_t_2), constants_DG$n),
  
  
  # mediator parameters at t =1
  theta_M2_t1 = rep(0,constants_DG$num_M2t1_reg_coeff),
  theta_M1_t1 = rep(0,constants_DG$num_M1t1_reg_coeff),
  #Sigma_M1M2_t1 = matrix(rinvgamma(4, 1, 1), nrow = 2, ncol = 2, byrow = TRUE),
  sigma_sq_M2_t1 = rinvgamma(1, 1, 1),
  sigma_sq_M1_t1 = rinvgamma(1, 1, 1),
  pi0_M2_t1 = 0.5,
  pi0_M1_t1 = 0.5,
  data_wide_M2t1_pos = rep(mean(digCom_workdf_wide$M2_t_1), constants_DG$n),
  data_wide_M1t1_pos = rep(mean(digCom_workdf_wide$M1_t_1), constants_DG$n),
  
  
  # treatment receipt parameters 
  theta_D_t3 = rep(0,constants_DG$num_Dt3_reg_coeff),
  theta_D_t2 = rep(0,constants_DG$num_Dt2_reg_coeff),
  theta_D_t1 = rep(0,constants_DG$num_Dt1_reg_coeff),
  
  
  # treatment assignment parameters
  theta_Z_t3 = rep(0,constants_DG$num_Zt3_reg_coeff),
  theta_Z_t2 = rep(0,constants_DG$num_Zt2_reg_coeff),
  theta_Z_t1 = rep(0,constants_DG$num_Zt1_reg_coeff)
)


Param_data_generating_model <- nimbleModel(code_DG, constants_DG, data_DG, inits_DG)  # model creation


compile_Param_data_generating_model <- compileNimble(Param_data_generating_model, showCompilerOutput = TRUE)  # model compilation

# MCMC configuration
config_MCMC_DG <- configureMCMC(Param_data_generating_model, useConjugacy = TRUE,
                                enableWAIC = TRUE,
                                monitors = c("beta_Yi", "sigma_sq_Yi", "pi0_Y",  
                                             #"bi_M2_sig_sq",  
                                             "theta_M2_t3", "theta_M1_t3", "sigma_sq_M2_t3", "sigma_sq_M1_t3", "pi0_M2_t3", "pi0_M1_t3",
                                             "theta_M2_t2", "theta_M1_t2", "sigma_sq_M2_t2", "sigma_sq_M1_t2", "pi0_M2_t2", "pi0_M1_t2",
                                             "theta_M2_t1", "theta_M1_t1", "sigma_sq_M2_t1", "sigma_sq_M1_t1", "pi0_M2_t1", "pi0_M1_t1",
                                             "theta_D_t1","theta_D_t2", "theta_D_t3",
                                             "theta_Z_t1","theta_Z_t2", "theta_Z_t3",
                                             "logDens"))


config_MCMC_DG$removeSamplers('logDens')   ## remove sampler assigned to 'logDens'
config_MCMC_DG$addSampler(target = 'logDens', type = 'sumLogPostDens')   ## add our custom sampler


mcmc_DG <- buildMCMC(config_MCMC_DG)
cmcmc_DG <- compileNimble(mcmc_DG, project = Param_data_generating_model, showCompilerOutput = TRUE)

# burn-in 1000 iterations and save 1 iter for data generation
num_iter_DG = 15000;  num_burnin_DG = 10000;  num_thin_DG = 5;  num_chains_DG= 1


samples_allChains_DG <- runMCMC(cmcmc_DG, niter = num_iter_DG, nburnin = num_burnin_DG,
                                nchains = num_chains_DG, 
                                setSeed = TRUE, thin = num_thin_DG, WAIC = TRUE)

post_samples_data_generating_model_list <- samples_allChains_DG$samples   # list of nchains matrices


param_samples_data_generating_model_df <- as.data.frame(post_samples_data_generating_model_list)

# only keep the posterior means of the parameters
param_samples_data_generating_model <- as.data.frame(t(colMeans(param_samples_data_generating_model_df, 
                                                                na.rm = TRUE)))

colnames(param_samples_data_generating_model) <- colnames(param_samples_data_generating_model_df)


#######################################################################################
####################Data generation####################################################
#######################################################################################
set.seed(run_ID)   #new seed for data generation
replicated_data_sample_size = 12000
number_of_mediators = 2
zero_vec = rep(0,2)

####################Generate L0 from the empirical distribution########################

intercept_vec =  rep(1, replicated_data_sample_size)
gender_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = mean(L0_mat[,2]))
baseline_age_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = mean(L0_mat[,3]))

replicated_df = data.frame(Intercept = intercept_vec, GENDER = gender_vec, BASELINE_AGE = baseline_age_vec)


####################################Generate Z_t1#####################################
#get the col num of all cols whose name begin with theta_Z_t1                    
theta_Z_t1_col_num = grep("^theta_Z_t1\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_Z_t1 for this iter
theta_Z_t1_vec = as.numeric(param_samples_data_generating_model[theta_Z_t1_col_num]) # outcome regression mean coefficients

Z_t1_prob = pnorm(as.matrix(replicated_df) %*% theta_Z_t1_vec)

Z_t1_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = Z_t1_prob)

replicated_df$Z_t1 = Z_t1_vec


####################################Generate D_t1#####################################
#get the col num of all cols whose name begin with theta_D_t1                    
theta_D_t1_col_num = grep("^theta_D_t1\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_D_t1 for this iter
theta_D_t1_vec = as.numeric(param_samples_data_generating_model[theta_D_t1_col_num]) # outcome regression mean coefficients

D_t1_prob = pnorm(as.matrix(replicated_df) %*% theta_D_t1_vec)

D_t1_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = D_t1_prob)

replicated_df$D_t1 = D_t1_vec


####################Generate mediator 1 M1t1####################################################


#get the col num of all cols whose name begin with theta_M1_t1                      
theta_M1_t1_col_num = grep("^theta_M1_t1\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_M1_t1 for this iter
theta_M1_t1_vec = as.numeric(param_samples_data_generating_model[theta_M1_t1_col_num]) # mediator 1 regression mean coefficients

pi0_M1_t1 = as.numeric(param_samples_data_generating_model$pi0_M1_t1)

sigma_sq_M1_t1 = as.numeric(param_samples_data_generating_model$sigma_sq_M1_t1)

M1t1_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_M1_t1)
M1t1_pos_mean = (as.matrix(replicated_df) %*% theta_M1_t1_vec)  #compute dot product of design matrix and the vector of coefficients
M1t1_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M1t1_pos_mean, sd = sqrt(sigma_sq_M1_t1)) * M1t1_pos_indicator

replicated_df$M1_t1 = M1t1_vec


####################Generate mediator 2 M2t1####################################################

#get the col num of all cols whose name begin with theta_M2_t1                      
theta_M2_t1_col_num = grep("^theta_M2_t1\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_M2_t1 for this iter
theta_M2_t1_vec = as.numeric(param_samples_data_generating_model[theta_M2_t1_col_num]) # mediator 1 regression mean coefficients

pi0_M2_t1 = as.numeric(param_samples_data_generating_model$pi0_M2_t1)

sigma_sq_M2_t1 = as.numeric(param_samples_data_generating_model$sigma_sq_M2_t1)

M2t1_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_M2_t1)
M2t1_pos_mean = (as.matrix(replicated_df) %*% theta_M2_t1_vec)  #compute dot product of design matrix and the vector of coefficients
M2t1_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M2t1_pos_mean, sd = sqrt(sigma_sq_M2_t1)) * M2t1_pos_indicator

######################################################################

replicated_df$M2_t1 = M2t1_vec
######################################################################


####################################Generate Z_t2#####################################
#get the col num of all cols whose name begin with theta_Z_t2                    
theta_Z_t2_col_num = grep("^theta_Z_t2\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_Z_t2 for this iter
theta_Z_t2_vec = as.numeric(param_samples_data_generating_model[theta_Z_t2_col_num]) # outcome regression mean coefficients

Z_t2_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                BASELINE_AGE = replicated_df$BASELINE_AGE, 
                                Z_t1 = replicated_df$Z_t1,
                                D_t1 = replicated_df$D_t1,
                                M1_t1 = replicated_df$M1_t1,
                                M2_t1 = replicated_df$M2_t1)


Z_t2_prob = pnorm(as.matrix(Z_t2_design_mat_DG) %*% theta_Z_t2_vec)

Z_t2_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = Z_t2_prob)

replicated_df$Z_t2 = Z_t2_vec


####################################Generate D_t2#####################################
#get the col num of all cols whose name begin with theta_D_t2                    
theta_D_t2_col_num = grep("^theta_D_t2\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_D_t2 for this iter
theta_D_t2_vec = as.numeric(param_samples_data_generating_model[theta_D_t2_col_num]) # outcome regression mean coefficients


D_t2_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                BASELINE_AGE = replicated_df$BASELINE_AGE,
                                D_t1 = replicated_df$D_t1,
                                M1_t1 = replicated_df$M1_t1,
                                M2_t1 = replicated_df$M2_t1, 
                                Z_t2 = replicated_df$Z_t2)

D_t2_prob = pnorm(as.matrix(D_t2_design_mat_DG) %*% theta_D_t2_vec)

D_t2_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = D_t2_prob)

replicated_df$D_t2 = D_t2_vec


####################Generate mediator 1 M1t2####################################################


#get the col num of all cols whose name begin with theta_M1_t2                      
theta_M1_t2_col_num = grep("^theta_M1_t2\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_M1_t2 for this iter
theta_M1_t2_vec = as.numeric(param_samples_data_generating_model[theta_M1_t2_col_num]) # mediator 1 regression mean coefficients

pi0_M1_t2 = as.numeric(param_samples_data_generating_model$pi0_M1_t2)

sigma_sq_M1_t2 = as.numeric(param_samples_data_generating_model$sigma_sq_M1_t2)

M1t2_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_M1_t2)


M1_t2_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                 BASELINE_AGE = replicated_df$BASELINE_AGE,
                                 M1_t1 = replicated_df$M1_t1, 
                                 M2_t1 = replicated_df$M2_t1, 
                                 Z_t2 = replicated_df$Z_t2,
                                 D_t2 = replicated_df$D_t2)


M1t2_pos_mean = (as.matrix(M1_t2_design_mat_DG) %*% theta_M1_t2_vec)  #compute dot product of design matrix and the vector of coefficients
M1t2_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M1t2_pos_mean, sd = sqrt(sigma_sq_M1_t2)) * M1t2_pos_indicator


replicated_df$M1_t2 = M1t2_vec


####################Generate mediator 2 M2t2####################################################

#get the col num of all cols whose name begin with theta_M2_t2                      
theta_M2_t2_col_num = grep("^theta_M2_t2\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_M2_t2 for this iter
theta_M2_t2_vec = as.numeric(param_samples_data_generating_model[theta_M2_t2_col_num]) # mediator 2 regression mean coefficients

pi0_M2_t2 = as.numeric(param_samples_data_generating_model$pi0_M2_t2)

sigma_sq_M2_t2 = as.numeric(param_samples_data_generating_model$sigma_sq_M2_t2)

M2t2_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_M2_t2)


M2_t2_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                 BASELINE_AGE = replicated_df$BASELINE_AGE,
                                 M2_t1 = replicated_df$M2_t1, 
                                 Z_t2 = replicated_df$Z_t2,
                                 D_t2 = replicated_df$D_t2,
                                 M1_t2 = replicated_df$M1_t2)


M2t2_pos_mean = (as.matrix(M2_t2_design_mat_DG) %*% theta_M2_t2_vec)  #compute dot product of design matrix and the vector of coefficients
M2t2_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M2t2_pos_mean, sd = sqrt(sigma_sq_M2_t2)) * M2t2_pos_indicator


replicated_df$M2_t2 = M2t2_vec


####################################Generate Z_t3#####################################
#get the col num of all cols whose name begin with theta_Z_t3                    
theta_Z_t3_col_num = grep("^theta_Z_t3\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_Z_t3 for this iter
theta_Z_t3_vec = as.numeric(param_samples_data_generating_model[theta_Z_t3_col_num]) # outcome regression mean coefficients


Z_t3_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                BASELINE_AGE = replicated_df$BASELINE_AGE, 
                                Z_t2 = replicated_df$Z_t2,
                                D_t2 = replicated_df$D_t2,
                                M1_t2 = replicated_df$M1_t2,
                                M2_t2 = replicated_df$M2_t2)


Z_t3_prob = pnorm(as.matrix(Z_t3_design_mat_DG) %*% theta_Z_t3_vec)

Z_t3_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = Z_t3_prob)

replicated_df$Z_t3 = Z_t3_vec


####################################Generate D_t3#####################################
#get the col num of all cols whose name begin with theta_D_t3                    
theta_D_t3_col_num = grep("^theta_D_t3\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_D_t3 for this iter
theta_D_t3_vec = as.numeric(param_samples_data_generating_model[theta_D_t3_col_num]) # outcome regression mean coefficients


D_t3_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                BASELINE_AGE = replicated_df$BASELINE_AGE, 
                                D_t2 = replicated_df$D_t2,
                                M1_t2 = replicated_df$M1_t2,
                                M2_t2 = replicated_df$M2_t2,
                                Z_t3 = replicated_df$Z_t3)


D_t3_prob = pnorm(as.matrix(D_t3_design_mat_DG) %*% theta_D_t3_vec)

D_t3_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = D_t3_prob)

replicated_df$D_t3 = D_t3_vec


####################Generate mediator 1 M1t3####################################################


#get the col num of all cols whose name begin with theta_M1_t3                      
theta_M1_t3_col_num = grep("^theta_M1_t3\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_M1_t3 for this iter
theta_M1_t3_vec = as.numeric(param_samples_data_generating_model[theta_M1_t3_col_num]) # mediator 1 regression mean coefficients

pi0_M1_t3 = as.numeric(param_samples_data_generating_model$pi0_M1_t3)

sigma_sq_M1_t3 = as.numeric(param_samples_data_generating_model$sigma_sq_M1_t3)

M1t3_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_M1_t3)


M1_t3_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                 BASELINE_AGE = replicated_df$BASELINE_AGE, 
                                 M1_t2 = replicated_df$M1_t2,
                                 M2_t2 = replicated_df$M2_t2,
                                 Z_t3 = replicated_df$Z_t3,
                                 D_t3 = replicated_df$D_t3)


M1t3_pos_mean = (as.matrix(M1_t3_design_mat_DG) %*% theta_M1_t3_vec)  #compute dot product of design matrix and the vector of coefficients
M1t3_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M1t3_pos_mean, sd = sqrt(sigma_sq_M1_t3)) * M1t3_pos_indicator


replicated_df$M1_t3 = M1t3_vec


####################Generate mediator 2 M2t3####################################################

#get the col num of all cols whose name begin with theta_M2_t3                      
theta_M2_t3_col_num = grep("^theta_M2_t3\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with theta_M2_t3 for this iter
theta_M2_t3_vec = as.numeric(param_samples_data_generating_model[theta_M2_t3_col_num]) # mediator 2 regression mean coefficients

pi0_M2_t3 = as.numeric(param_samples_data_generating_model$pi0_M2_t3)

sigma_sq_M2_t3 = as.numeric(param_samples_data_generating_model$sigma_sq_M2_t3)

M2t3_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_M2_t3)

M2_t3_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                                 BASELINE_AGE = replicated_df$BASELINE_AGE, 
                                 M2_t2 = replicated_df$M2_t2,
                                 Z_t3 = replicated_df$Z_t3,
                                 D_t3 = replicated_df$D_t3,
                                 M1_t3 = replicated_df$M1_t3)

M2t3_pos_mean = (as.matrix(M2_t3_design_mat_DG) %*% theta_M2_t3_vec)  #compute dot product of design matrix and the vector of coefficients
M2t3_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M2t3_pos_mean, sd = sqrt(sigma_sq_M2_t3)) * M2t3_pos_indicator


replicated_df$M2_t3 = M2t3_vec


####################Generate outcome Y####################################################

#get the col num of all cols whose name begin with beta_Yi                       
beta_Yi_col_num = grep("^beta_Yi\\[\\d+\\]$", colnames(param_samples_data_generating_model))
#get the elements of all cols whose name begin with beta_Yi for this iter
beta_Yi_vec = as.numeric(param_samples_data_generating_model[beta_Yi_col_num]) # outcome regression mean coefficients

pi0_Y = as.numeric(param_samples_data_generating_model$pi0_Y)

sigma_sq_Yi = as.numeric(param_samples_data_generating_model$sigma_sq_Yi)

Y_pos_indicator = rbinom(n = replicated_data_sample_size, size = 1, prob = 1-pi0_Y)


Y_design_mat_DG = data.frame(Intercept = replicated_df$Intercept, GENDER = replicated_df$GENDER, 
                             BASELINE_AGE = replicated_df$BASELINE_AGE, 
                             Z_t3 = replicated_df$Z_t3,
                             D_t3 = replicated_df$D_t3,
                             M1_t3 = replicated_df$M1_t3,
                             M2_t3 = replicated_df$M2_t3)


Y_pos_mean = as.matrix(Y_design_mat_DG) %*% beta_Yi_vec#compute dot product of design matrix and the vector of coefficients
Y_vector = rtruncnorm(n = replicated_data_sample_size, a=0, b=Inf,
                      mean = Y_pos_mean, sd = sqrt(sigma_sq_Yi)) * Y_pos_indicator

replicated_df$Y = Y_vector


#######################################################################################
######Fit EDPM true model on generated data#########################################
#######################################################################################

L0_df_wide_fit = data.frame(Intercept = replicated_df$Intercept,
                            GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE) 

Y_vec_wide_fit = replicated_df$Y


########################################################################################
### Fit simple linear regression to find the data dependent base distribution parameters
########################################################################################

fit_linear_reg_AllCovs = function(output_vec, all_Covs_df_wide){
  num_subjects = nrow(all_Covs_df_wide)
  lin_ref_fit = lm(output_vec~., data = all_Covs_df_wide)
  # Maximum likelihood estimates (MLEs) Y ~ L0
  mu_beta = as.numeric(coef(lin_ref_fit))
  # Asymptotic variances Y ~ L0
  sig_sq_beta = as.numeric(diag(vcov(lin_ref_fit)))
  sig_sq_beta = (num_subjects/5)* sig_sq_beta
  
  # sigma_hat Y ~ L0
  sigma_sq_hat =  (summary(lin_ref_fit)$sigma)^2
  return_list = list(Coeff = mu_beta, Sig_sq_coeff = sig_sq_beta, Sigma_sq_hat = sigma_sq_hat)
  return(return_list)
}


fit_linear_reg = function(output_vec, L0_data_wide){
  num_subjects = nrow(L0_df_wide)
  lin_ref_fit = lm(output_vec~GENDER+BASELINE_AGE, data = L0_data_wide)
  # Maximum likelihood estimates (MLEs) Y ~ L0
  mu_beta = as.numeric(coef(lin_ref_fit))
  # Asymptotic variances Y ~ L0
  sig_sq_beta = as.numeric(diag(vcov(lin_ref_fit)))
  sig_sq_beta = (num_subjects/5)* sig_sq_beta
  
  # sigma_hat Y ~ L0
  sigma_sq_hat =  (summary(lin_ref_fit)$sigma)^2
  return_list = list(Coeff = mu_beta, Sig_sq_coeff = sig_sq_beta, Sigma_sq_hat = sigma_sq_hat)
  return(return_list)
}

fit_multiple_linear_reg = function(cbind_output_vecs, L0_data_wide){
  num_subjects = nrow(L0_df_wide)
  lin_ref_fit = lm(cbind_output_vecs~GENDER+BASELINE_AGE, data = L0_data_wide)
  # Maximum likelihood estimates (MLEs) 
  mu_theta = coef(lin_ref_fit)
  # Asymptotic variances Y ~ L0
  sig_sq_theta = (num_subjects/5)* as.numeric(diag(vcov(lin_ref_fit)))
  sig_sq_theta = matrix(sig_sq_theta, nrow = ncol(L0_data_wide), ncol = ncol(cbind_output_vecs))
  
  
  sigma_hat =  cov(resid(lin_ref_fit))
  return_list = list(Coeff = mu_theta, Sig_sq_coeff = sig_sq_theta, Sigma_hat = sigma_hat)
  return(return_list)
  
}


fit_logistic_reg = function(output_vec, all_Covs_df_wide){
  num_subjects = nrow(all_Covs_df_wide)
  logistic_ref_fit = glm(output_vec~., family = binomial(link = "logit"), 
                         data = all_Covs_df_wide)
  # Maximum likelihood estimates (MLEs) Y ~ L0
  mu_beta = as.numeric(coef(logistic_ref_fit))
  # Asymptotic variances Y ~ L0
  sig_sq_beta = as.numeric(diag(vcov(logistic_ref_fit)))
  sig_sq_beta = (num_subjects/5)* sig_sq_beta
  
  return_list = list(Coeff = mu_beta, Sig_sq_coeff = sig_sq_beta)
  return(return_list)
  
}


###############################################################################
######################## Do this for Y vs Z,D,L0 ###################################
Y_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                              Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1, M1_t_1 = replicated_df$M1_t1, M2_t_1 = replicated_df$M2_t1,
                              Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2, M1_t_2 = replicated_df$M1_t2, M2_t_2 = replicated_df$M2_t2,
                              Z_t_3 = replicated_df$Z_t3, D_t_3 = replicated_df$D_t3, M1_t_3 = replicated_df$M1_t3, M2_t_3 = replicated_df$M2_t3)

###############################################################################


lin_reg_YvAllCovs = fit_linear_reg_AllCovs(replicated_df$Y[replicated_df$Y !=0], 
                                           Y_reg_design_mat[replicated_df$Y !=0,])


###############################################################################
######################## Do this for M vs Z,D,L0 ###################################
Mt3_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1, 
                                Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2,
                                Z_t_3 = replicated_df$Z_t3, D_t_3 = replicated_df$D_t3)

Mt2_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1, 
                                Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2)

Mt1_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1)

###############################################################################


# Fit multiple linear regression M2 ~ All  Covs for M2
#-------------------------------------------------------------------------------------------------------------


lin_reg_M2t3vAllCovs = fit_linear_reg_AllCovs(replicated_df$M2_t3[replicated_df$M2_t3 !=0],
                                              Mt3_reg_design_mat[replicated_df$M2_t3 !=0,])
lin_reg_M2t2vAllCovs = fit_linear_reg_AllCovs(replicated_df$M2_t2[replicated_df$M2_t2 !=0],
                                              Mt2_reg_design_mat[replicated_df$M2_t2 !=0, ])
lin_reg_M2t1vAllCovs = fit_linear_reg_AllCovs(replicated_df$M2_t1[replicated_df$M2_t1 !=0],
                                              Mt1_reg_design_mat[replicated_df$M2_t1 !=0, ])


# Fit multiple linear regression M1 ~ All  Covs for M2
#-------------------------------------------------------------------------------------------------------------


lin_reg_M1t3vAllCovs = fit_linear_reg_AllCovs(replicated_df$M1_t3[replicated_df$M1_t3 !=0],
                                              Mt3_reg_design_mat[replicated_df$M1_t3 !=0,])
lin_reg_M1t2vAllCovs = fit_linear_reg_AllCovs(replicated_df$M1_t2[replicated_df$M1_t2 !=0],
                                              Mt2_reg_design_mat[replicated_df$M1_t2 !=0, ])
lin_reg_M1t1vAllCovs = fit_linear_reg_AllCovs(replicated_df$M1_t1[replicated_df$M1_t1 !=0],
                                              Mt1_reg_design_mat[replicated_df$M1_t1 !=0, ])


###############################################################################
######################## Do this for D vs Z,L0 ###################################
Dt3_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                Z_t_1 = replicated_df$Z_t1, 
                                Z_t_2 = replicated_df$Z_t2,
                                Z_t_3 = replicated_df$Z_t3)

Dt2_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                Z_t_1 = replicated_df$Z_t1, 
                                Z_t_2 = replicated_df$Z_t2)

Dt1_reg_design_mat = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                Z_t_1 = replicated_df$Z_t1)

###############################################################################


# Fit logistic regression D ~ L0
#-------------------------------------------------------------------------------------------------------------
logistic_reg_D_t3vAllCovs = fit_logistic_reg(replicated_df$D_t3, Dt3_reg_design_mat)

logistic_reg_D_t2vAllCovs = fit_logistic_reg(replicated_df$D_t2, Dt2_reg_design_mat)

logistic_reg_D_t1vAllCovs = fit_logistic_reg(replicated_df$D_t1, Dt1_reg_design_mat)


# Fit logistic regression Z ~ L0
#-------------------------------------------------------------------------------------------------------------
logistic_reg_Z_t3vL0 = fit_logistic_reg(replicated_df$Z_t3, L0_df_wide_fit[, -1])

logistic_reg_Z_t2vL0 = fit_logistic_reg(replicated_df$Z_t2, L0_df_wide_fit[, -1])

logistic_reg_Z_t1vL0 = fit_logistic_reg(replicated_df$Z_t1, L0_df_wide_fit[, -1])


longitudinal_data_wide_M2_fit = data.frame(M2_t_1 = replicated_df$M2_t1,
                                           M2_t_2 = replicated_df$M2_t2,
                                           M2_t_3 = replicated_df$M2_t3)


longitudinal_data_wide_M1_fit = data.frame(M1_t_1 = replicated_df$M1_t1, 
                                           M1_t_2 = replicated_df$M1_t2,
                                           M1_t_3 = replicated_df$M1_t3)


# M2 MLE priors -- separate vectors per time point (different lengths)
mle_coeff_M2_t1 = unlist(lin_reg_M2t1vAllCovs[1])   # length 5: Intercept, GENDER, BASELINE_AGE, Z1, D1
mle_coeff_M2_t2 = unlist(lin_reg_M2t2vAllCovs[1])   # length 7: + Z2, D2
mle_coeff_M2_t3 = unlist(lin_reg_M2t3vAllCovs[1])   # length 9: + Z3, D3

mle_sig_sq_M2_t1 = unlist(lin_reg_M2t1vAllCovs[2])  # length 5
mle_sig_sq_M2_t2 = unlist(lin_reg_M2t2vAllCovs[2])  # length 7
mle_sig_sq_M2_t3 = unlist(lin_reg_M2t3vAllCovs[2])  # length 9


sigma_sq_hat_M2_t1 = unlist(lin_reg_M2t1vAllCovs[3])  # scalar
sigma_sq_hat_M2_t2 = unlist(lin_reg_M2t2vAllCovs[3])  # scalar
sigma_sq_hat_M2_t3 = unlist(lin_reg_M2t3vAllCovs[3])  # scalar


# M1 MLE priors -- separate vectors per time point (different lengths)
mle_coeff_M1_t1 = unlist(lin_reg_M1t1vAllCovs[1])   # length 5
mle_coeff_M1_t2 = unlist(lin_reg_M1t2vAllCovs[1])   # length 7
mle_coeff_M1_t3 = unlist(lin_reg_M1t3vAllCovs[1])   # length 9

mle_sig_sq_M1_t1 = unlist(lin_reg_M1t1vAllCovs[2])  # length 5
mle_sig_sq_M1_t2 = unlist(lin_reg_M1t2vAllCovs[2])  # length 7
mle_sig_sq_M1_t3 = unlist(lin_reg_M1t3vAllCovs[2])  # length 9


sigma_sq_hat_M1_t1 = unlist(lin_reg_M1t1vAllCovs[3])  # scalar
sigma_sq_hat_M1_t2 = unlist(lin_reg_M1t2vAllCovs[3])  # scalar
sigma_sq_hat_M1_t3 = unlist(lin_reg_M1t3vAllCovs[3])  # scalar


longitudinal_data_wide_D_fit = data.frame(D_t_1 = replicated_df$D_t1,
                                          D_t_2 = replicated_df$D_t2,
                                          D_t_3 = replicated_df$D_t3)


# Design matrices for D fitting: reuse existing MLE regression matrices, prepend Intercept
D_fit_design_mat_t1 = cbind(Intercept = replicated_df$Intercept, Dt1_reg_design_mat)
D_fit_design_mat_t2 = cbind(Intercept = replicated_df$Intercept, Dt2_reg_design_mat)
D_fit_design_mat_t3 = cbind(Intercept = replicated_df$Intercept, Dt3_reg_design_mat)

# Design matrices for M1/M2 fitting: reuse existing MLE regression matrices, prepend Intercept
M_fit_design_mat_t1 = cbind(Intercept = replicated_df$Intercept, Mt1_reg_design_mat)
M_fit_design_mat_t2 = cbind(Intercept = replicated_df$Intercept, Mt2_reg_design_mat)
M_fit_design_mat_t3 = cbind(Intercept = replicated_df$Intercept, Mt3_reg_design_mat)


# D MLE priors -- separate vectors per time point (different lengths)
mle_coeff_D_t1 = unlist(logistic_reg_D_t1vAllCovs[1])   # length 4: Intercept, GENDER, BASELINE_AGE, Z1
mle_coeff_D_t2 = unlist(logistic_reg_D_t2vAllCovs[1])   # length 5: + Z2
mle_coeff_D_t3 = unlist(logistic_reg_D_t3vAllCovs[1])   # length 6: + Z3

mle_sig_sq_D_t1 = unlist(logistic_reg_D_t1vAllCovs[2])  # length 4
mle_sig_sq_D_t2 = unlist(logistic_reg_D_t2vAllCovs[2])  # length 5
mle_sig_sq_D_t3 = unlist(logistic_reg_D_t3vAllCovs[2])  # length 6


longitudinal_data_wide_Z_fit = data.frame(Z_t_1 = replicated_df$Z_t1,
                                          Z_t_2 = replicated_df$Z_t2,
                                          Z_t_3 = replicated_df$Z_t3)


# K * T matrix
mle_coeff_Z = data.frame(Z_t_1_coeff = unlist(logistic_reg_Z_t1vL0[1]),
                         Z_t_2_coeff = unlist(logistic_reg_Z_t2vL0[1]), 
                         Z_t_3_coeff = unlist(logistic_reg_Z_t3vL0[1]))


# K * T matrix
mle_sig_sq_Z = data.frame( Z_t_1_sig_sq = unlist(logistic_reg_Z_t1vL0[2]),
                           Z_t_2_sig_sq = unlist(logistic_reg_Z_t2vL0[2]),
                           Z_t_3_sig_sq = unlist(logistic_reg_Z_t3vL0[2]))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
num_Bin_L0_fit = 2  # both gender and baseline age are binary

pr_theta_Bin_L0_fit = as.numeric(colMeans(L0_df_wide_fit[2:3], na.rm = TRUE))


code_fit <- nimbleCode({
  logDens ~ dnorm(0, 1)    ## this distribution does not matter: only need to compute model sum-posterior-log-density
  # stick breaks outer layer
  # alpha_bet ~ dgamma(1, 1)
  for(i in 1:(N-1)){ #N=num outer layers
    breaks_r[i] ~ dbeta(1, alpha_bet)
  }
  breaks_r[N] <- 1  # Ensure the last stick-breaking weight closes the distribution
  # component prob outer layer
  ksi_r[1:N] <- stick_breaking(breaks_r[1:(N-1)])
  
  # stick breaks and component prob inner layer
  for(i in 1:N){
    #alpha_thet[i] ~ dgamma(1, 1)
    for (j in 1:(M-1)){
      breaks_sr[i,j] ~ dbeta(1, alpha_thet[i])
    }
    breaks_sr[i,M] <- 1  # Ensure the last stick-breaking weight closes the distribution
    ksi_sr[i,1:M] <- stick_breaking(breaks_sr[i,1:(M-1)])
  }
  
  # Atoms from the parameter base distributions
  for(i in 1:N){ #loop over outer cluster
    # outer cluster parameters
    pi0_Y[i] ~ dbeta(2,2)
    sigma_sq_Yi[i] ~ dinvgamma(shape = 3, scale = sigSq_Y_prior)
    #i^th row of the matrix beta_Yi (beta_Y is a N * K matrix)
    
    for(k in 1:num_Y_reg_coeff){beta_Yi[i, k]~ dnorm(mean = mu_betaY_prior[k],
                                                     var = sigSq_betaY_prior[k])}  
    
    for(j in 1:M){ #loop over inner cluster
      # inner cluster parameters
      ## baseline confounders paramters (each param a N*M matrix)
      ### var param of continuous L0 (d^th vector of array sigma_sq_Cont_L0i)
      #sigma_sq_Cont_L0i[i,j] ~ dinvgamma(shape = 2, scale = sigSq_contL0_prior)
      ### mean param of continuous L0 (d^th vector of array mu_Cont_L0i) 
      #mu_Cont_L0i[i,j] ~ dnorm(mean = mu_contL0_prior,
      #                         sd = sqrt(sigSq_contL0_prior))
      ### probability param of binary L0 (2 of them at the moment) 
      for(k in 1:num_bin_L0){
        pr_Bin_L0i[i,j,k] ~ dbeta(2, 2)
      }
      
      for(t in 1:max_t){ #loop over time
        ### var param for mediators (t = 1,2,...,T)
        sigma_sq_M2_t[i,j,t] ~ dinvgamma(shape = 3, scale = sigma_sq_hat_M2_t[t])
        sigma_sq_M1_t[i,j,t] ~ dinvgamma(shape = 3, scale = sigma_sq_hat_M1_t[t])
        
        ### hurdle probabilities
        pi0_M1_t[i,j,t]~ dbeta(2,2)
        pi0_M2_t[i,j,t]~ dbeta(2,2)
      }
      
      ## treatment assignment mean paramters (Z stays K=3 at all t)
      for(k in 1:K){
        for(t in 1:max_t){
          theta_Z_t[i,j,k,t] ~ dnorm(mean = mu_theta_Z_t_prior[k,t],
                                     sd = sqrt(sigSq_theta_Z_t_prior[k,t]))
        }
      }
      
      ## treatment receipt mean paramters (D: separate per time point)
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
      
      ## mediator 1 mean paramters (M1: separate per time point)
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
      
      ## mediator 2 mean paramters (M2: separate per time point)
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
      
    }
  }
  
  # likelihood
  for(i in 1:n) {
    # component numbers for this subject
    comp_num_outer[i] ~ dcat(ksi_r[1:N])
    comp_num_inner[i] ~ dcat(ksi_sr[comp_num_outer[i],1:M])
    
    
    # baseline confounders
    ## continuous baseline confounder (0 at the moment)
    ## binary baseline confounder (2 at the moment)
    for(k in 1:num_bin_L0){
      L0_bin[i,k] ~ dbern(pr_Bin_L0i[comp_num_outer[i], comp_num_inner[i],k]) 
    }
    
    
    # ========== t = 1 ==========
    # Z at t=1 (unchanged: uses L0_mat)
    z_mean[i,1] <- inprod(theta_Z_t[comp_num_outer[i], comp_num_inner[i], 1:K, 1],
                          L0_mat[i, 1:K])
    data_wide_Z[i,1] ~ dbern(pnorm(z_mean[i,1],0,1))
    
    # D at t=1 (uses D_design_mat_t1)
    d_mean[i,1] <- inprod(theta_D_t1[comp_num_outer[i], comp_num_inner[i], 1:K_D_t1],
                          D_design_mat_t1[i, 1:K_D_t1])
    data_wide_D[i,1] ~ dbern(pnorm(d_mean[i,1],0,1))
    
    # M1 at t=1 (uses M_design_mat_t1)
    Z_latent_M1[i,1] ~ dbern(1 - pi0_M1_t[comp_num_outer[i], comp_num_inner[i],1])
    m1_mean[i,1] <- inprod(theta_M1_t1[comp_num_outer[i], comp_num_inner[i], 1:K_M_t1],
                           M_design_mat_t1[i, 1:K_M_t1])
    m1_sd[i,1] <- sqrt(sigma_sq_M1_t[comp_num_outer[i],comp_num_inner[i],1])
    data_wide_M1_pos[i,1] ~ dnorm(m1_mean[i,1], sd = m1_sd[i,1])
    data_wide_M1[i,1] ~ dnorm(Z_latent_M1[i,1]*data_wide_M1_pos[i,1], sd = 0.25)
    
    # M2 at t=1 (uses M_design_mat_t1)
    Z_latent_M2[i,1] ~ dbern(1 - pi0_M2_t[comp_num_outer[i], comp_num_inner[i],1])
    m2_mean[i,1] <- inprod(theta_M2_t1[comp_num_outer[i], comp_num_inner[i], 1:K_M_t1],
                           M_design_mat_t1[i, 1:K_M_t1])
    m2_sd[i,1] <- sqrt(sigma_sq_M2_t[comp_num_outer[i],comp_num_inner[i],1])
    data_wide_M2_pos[i,1] ~ dnorm(m2_mean[i,1], sd = m2_sd[i,1])
    data_wide_M2[i,1] ~ dnorm(Z_latent_M2[i,1]*data_wide_M2_pos[i,1], sd = 0.25)
    
    
    # ========== t = 2 ==========
    z_mean[i,2] <- inprod(theta_Z_t[comp_num_outer[i], comp_num_inner[i], 1:K, 2],
                          L0_mat[i, 1:K])
    data_wide_Z[i,2] ~ dbern(pnorm(z_mean[i,2],0,1))
    
    d_mean[i,2] <- inprod(theta_D_t2[comp_num_outer[i], comp_num_inner[i], 1:K_D_t2],
                          D_design_mat_t2[i, 1:K_D_t2])
    data_wide_D[i,2] ~ dbern(pnorm(d_mean[i,2],0,1))
    
    Z_latent_M1[i,2] ~ dbern(1 - pi0_M1_t[comp_num_outer[i], comp_num_inner[i],2])
    m1_mean[i,2] <- inprod(theta_M1_t2[comp_num_outer[i], comp_num_inner[i], 1:K_M_t2],
                           M_design_mat_t2[i, 1:K_M_t2])
    m1_sd[i,2] <- sqrt(sigma_sq_M1_t[comp_num_outer[i],comp_num_inner[i],2])
    data_wide_M1_pos[i,2] ~ dnorm(m1_mean[i,2], sd = m1_sd[i,2])
    data_wide_M1[i,2] ~ dnorm(Z_latent_M1[i,2]*data_wide_M1_pos[i,2], sd = 0.25)
    
    Z_latent_M2[i,2] ~ dbern(1 - pi0_M2_t[comp_num_outer[i], comp_num_inner[i],2])
    m2_mean[i,2] <- inprod(theta_M2_t2[comp_num_outer[i], comp_num_inner[i], 1:K_M_t2],
                           M_design_mat_t2[i, 1:K_M_t2])
    m2_sd[i,2] <- sqrt(sigma_sq_M2_t[comp_num_outer[i],comp_num_inner[i],2])
    data_wide_M2_pos[i,2] ~ dnorm(m2_mean[i,2], sd = m2_sd[i,2])
    data_wide_M2[i,2] ~ dnorm(Z_latent_M2[i,2]*data_wide_M2_pos[i,2], sd = 0.25)
    
    
    # ========== t = 3 ==========
    z_mean[i,3] <- inprod(theta_Z_t[comp_num_outer[i], comp_num_inner[i], 1:K, 3],
                          L0_mat[i, 1:K])
    data_wide_Z[i,3] ~ dbern(pnorm(z_mean[i,3],0,1))
    
    d_mean[i,3] <- inprod(theta_D_t3[comp_num_outer[i], comp_num_inner[i], 1:K_D_t3],
                          D_design_mat_t3[i, 1:K_D_t3])
    data_wide_D[i,3] ~ dbern(pnorm(d_mean[i,3],0,1))
    
    Z_latent_M1[i,3] ~ dbern(1 - pi0_M1_t[comp_num_outer[i], comp_num_inner[i],3])
    m1_mean[i,3] <- inprod(theta_M1_t3[comp_num_outer[i], comp_num_inner[i], 1:K_M_t3],
                           M_design_mat_t3[i, 1:K_M_t3])
    m1_sd[i,3] <- sqrt(sigma_sq_M1_t[comp_num_outer[i],comp_num_inner[i],3])
    data_wide_M1_pos[i,3] ~ dnorm(m1_mean[i,3], sd = m1_sd[i,3])
    data_wide_M1[i,3] ~ dnorm(Z_latent_M1[i,3]*data_wide_M1_pos[i,3], sd = 0.25)
    
    Z_latent_M2[i,3] ~ dbern(1 - pi0_M2_t[comp_num_outer[i], comp_num_inner[i],3])
    m2_mean[i,3] <- inprod(theta_M2_t3[comp_num_outer[i], comp_num_inner[i], 1:K_M_t3],
                           M_design_mat_t3[i, 1:K_M_t3])
    m2_sd[i,3] <- sqrt(sigma_sq_M2_t[comp_num_outer[i],comp_num_inner[i],3])
    data_wide_M2_pos[i,3] ~ dnorm(m2_mean[i,3], sd = m2_sd[i,3])
    data_wide_M2[i,3] ~ dnorm(Z_latent_M2[i,3]*data_wide_M2_pos[i,3], sd = 0.25)
    
    
    #outcome likelihood
    # Hurdle model for the outcome
    Z_latent_Y[i] ~ dbern(1 - pi0_Y[comp_num_outer[i]])  # Indicator for positive values
    y_mean[i] <- inprod(beta_Yi[comp_num_outer[i],1:num_Y_reg_coeff], 
                        Y_reg_design_mat[i, 1:num_Y_reg_coeff]) 
    y_sd[i] <-  sqrt(sigma_sq_Yi[comp_num_outer[i]])
    # follows a normal distribution between 0 and Inf (inclusive of 0 and Inf)
    #Y[i] ~ T(dnorm(y_mean[i], sd =  y_sd[i]),0,)
    
    
    Y_pos[i] ~ dnorm(y_mean[i], sd =  y_sd[i])
    
    Y[i] ~ dnorm(Z_latent_Y[i]*Y_pos[i], sd = 0.25)
    
    #Y[i] ~ dconstraint(Y[i] == Z_latent_Y[i]* Y_pos[i])
    
  }
})


num_outer_cluster_N = 10; num_inner_cluster_M = 4
num_Y_reg_coeff = ncol(Y_reg_design_mat) +1 #ncol(Y_reg_design_mat) does not include intercept at this point
alpha_thet_const = 0.5

L0_mat_fit = as.matrix(L0_df_wide_fit)
num_L0_fit =3   # one intercept and 2 L0


constants_fit <- list(alpha_bet = 1,
                      alpha_thet = rep(alpha_thet_const, num_outer_cluster_N),
                      N = num_outer_cluster_N,
                      M = num_inner_cluster_M, 
                      n = replicated_data_sample_size, # Number of individuals
                      K = num_L0_fit,                                # one intercept and 2 L0
                      max_t = 3,
                      num_Y_reg_coeff = num_Y_reg_coeff, #ncol(Y_reg_design_mat) does not include intercept at this point
                      num_bin_L0 =2,           # number of  binary baseline confounders
                      #num_cont_L0 =1               # num of  cont baseline confounders
                      K_D_t1 = ncol(D_fit_design_mat_t1), K_D_t2 = ncol(D_fit_design_mat_t2), K_D_t3 = ncol(D_fit_design_mat_t3),
                      K_M_t1 = ncol(M_fit_design_mat_t1), K_M_t2 = ncol(M_fit_design_mat_t2), K_M_t3 = ncol(M_fit_design_mat_t3)
)


data_fit <- list(Y = replicated_df$Y,
                 Y_reg_design_mat = cbind(replicated_df$Intercept, Y_reg_design_mat),
                 mu_betaY_prior = unlist(lin_reg_YvAllCovs[1]), # vector of means for betaY base distributions
                 sigSq_betaY_prior= unlist(lin_reg_YvAllCovs[2]),               # vector of variances for betaY base distributions
                 sigSq_Y_prior = unlist(lin_reg_YvAllCovs[3]),
                 Z_latent_Y = ifelse(replicated_df$Y ==0, 0,1),
                 
                 # --------------------theta level data---------------------------
                 #rr = rep(1:num_outer_cluster_N, each = num_inner_cluster_M),
                 #ss = rep(1:num_inner_cluster_M, num_outer_cluster_N),
                 # data for baseline confounders
                 
                 L0_mat = L0_mat_fit,
                 L0_bin = L0_df_wide_fit[2:3],
                 
                 #pr_binL0_prior = pr_theta_Bin_L0,
                 #mu_contL0_prior = mu_theta_L0,   
                 #sigSq_contL0_prior = sig_sq_theta_L0,
                 
                 data_wide_M2 = longitudinal_data_wide_M2_fit,   # long data in col order: M2_1,M2_2, ..., M2_T
                 mu_theta_M2_t1_prior = mle_coeff_M2_t1,
                 mu_theta_M2_t2_prior = mle_coeff_M2_t2,
                 mu_theta_M2_t3_prior = mle_coeff_M2_t3,
                 sigSq_theta_M2_t1_prior = mle_sig_sq_M2_t1,
                 sigSq_theta_M2_t2_prior = mle_sig_sq_M2_t2,
                 sigSq_theta_M2_t3_prior = mle_sig_sq_M2_t3,
                 sigma_sq_hat_M2_t = c(sigma_sq_hat_M2_t1, sigma_sq_hat_M2_t2, sigma_sq_hat_M2_t3),
                 Z_latent_M2 = ifelse(longitudinal_data_wide_M2_fit ==0, 0,1),
                 M_design_mat_t1 = M_fit_design_mat_t1,
                 M_design_mat_t2 = M_fit_design_mat_t2,
                 M_design_mat_t3 = M_fit_design_mat_t3,
                 
                 data_wide_M1 = longitudinal_data_wide_M1_fit, # long data in col order: M1_1,M1_2, ..., M1_T
                 mu_theta_M1_t1_prior = mle_coeff_M1_t1,
                 mu_theta_M1_t2_prior = mle_coeff_M1_t2,
                 mu_theta_M1_t3_prior = mle_coeff_M1_t3,
                 sigSq_theta_M1_t1_prior = mle_sig_sq_M1_t1,
                 sigSq_theta_M1_t2_prior = mle_sig_sq_M1_t2,
                 sigSq_theta_M1_t3_prior = mle_sig_sq_M1_t3,
                 sigma_sq_hat_M1_t = c(sigma_sq_hat_M1_t1, sigma_sq_hat_M1_t2, sigma_sq_hat_M1_t3),
                 Z_latent_M1 = ifelse(longitudinal_data_wide_M1_fit ==0, 0,1),
                 
                 data_wide_D = longitudinal_data_wide_D_fit,  # long data in col order: D_1, ...,D_T
                 mu_theta_D_t1_prior = mle_coeff_D_t1,
                 mu_theta_D_t2_prior = mle_coeff_D_t2,
                 mu_theta_D_t3_prior = mle_coeff_D_t3,
                 sigSq_theta_D_t1_prior = mle_sig_sq_D_t1,
                 sigSq_theta_D_t2_prior = mle_sig_sq_D_t2,
                 sigSq_theta_D_t3_prior = mle_sig_sq_D_t3,
                 D_design_mat_t1 = D_fit_design_mat_t1,
                 D_design_mat_t2 = D_fit_design_mat_t2,
                 D_design_mat_t3 = D_fit_design_mat_t3,
                 
                 data_wide_Z = longitudinal_data_wide_Z_fit, # long data in col order: Z_1, ...,Z_T
                 mu_theta_Z_t_prior = mle_coeff_Z,       # num_L0 * T matrix of means for theta_Z_t base distn
                 sigSq_theta_Z_t_prior= mle_sig_sq_Z     # num_L0 * T matrix of vars for theta_Z_t base
                 
)


inits_fit <-list(
  # outcome model parameters
  beta_Yi = matrix(0,
                   nrow = constants_fit$N, ncol = constants_fit$num_Y_reg_coeff),
  sigma_sq_Yi = rinvgamma(constants_fit$N, 1, 1),
  pi0_Y = rep(0.5, constants_fit$N),
  Y_pos = rep(mean(replicated_df$Y), constants_fit$n),
  #Z_latent_Y = ifelse(digCom_workdf_wide$Y_t_3 ==0, 0,1),
  
  # mediator parameters M2
  #bi_M2_sig_sq = rinvgamma(1, 1, 1),
  #bi_M2 = rep(0,constants_fit$n),
  theta_M2_t1 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_M_t1)),
  theta_M2_t2 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_M_t2)),
  theta_M2_t3 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_M_t3)),
  sigma_sq_M2_t = array(rinvgamma(constants_fit$N * constants_fit$M * constants_fit$max_t, 1, 1),
                        c(constants_fit$N, constants_fit$M, constants_fit$max_t)),
  pi0_M2_t = array(0.5, c(constants_fit$N, constants_fit$M, constants_fit$max_t)),
  data_wide_M2_pos = matrix(rep(colMeans(longitudinal_data_wide_M2_fit), constants_fit$n),
                            nrow = constants_fit$n, ncol = constants_fit$max_t, byrow = TRUE),
  #Z_latent_M2 = ifelse(longitudinal_data_wide_M2 ==0, 0,1),
  
  # mediator parameters M1
  #bi_M1_sig_sq = rinvgamma(1, 1, 1),
  #bi_M1 = rep(0,constants_fit$n),
  theta_M1_t1 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_M_t1)),
  theta_M1_t2 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_M_t2)),
  theta_M1_t3 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_M_t3)),
  sigma_sq_M1_t = array(rinvgamma(constants_fit$N * constants_fit$M * constants_fit$max_t, 1, 1),
                        c(constants_fit$N, constants_fit$M, constants_fit$max_t)),
  pi0_M1_t = array(0.5, c(constants_fit$N, constants_fit$M, constants_fit$max_t)),
  data_wide_M1_pos = matrix(rep(colMeans(longitudinal_data_wide_M1_fit), constants_fit$n),
                            nrow = constants_fit$n, ncol = constants_fit$max_t, byrow = TRUE),
  #Z_latent_M1 = ifelse(longitudinal_data_wide_M1 ==0, 0,1),
  
  # treatment receipt parameters
  #bi_D_sig_sq = rinvgamma(1, 1, 1),
  #bi_D = rep(0,constants_fit$n),
  theta_D_t1 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_D_t1)),
  theta_D_t2 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_D_t2)),
  theta_D_t3 = array(0, c(constants_fit$N, constants_fit$M, constants_fit$K_D_t3)),
  
  
  # treatment assignment parameters
  #bi_Z_sig_sq = rinvgamma(1, 1, 1),
  #bi_Z = rep(0,constants_fit$n),
  theta_Z_t = array(0,
                    c(constants_fit$N, constants_fit$M, constants_fit$K, constants_fit$max_t)),
  
  
  pr_Bin_L0i = array(0.5, c(constants_fit$N, constants_fit$M, constants_fit$num_bin_L0)),
  # EDPM parameters
  comp_num_outer = sample(1:constants_fit$N, size = constants_fit$n, 
                          replace = TRUE),
  comp_num_inner = sample(1:constants_fit$M, size = constants_fit$n, 
                          replace = TRUE),
  breaks_r = rbeta(constants_fit$N, 1, 1),
  breaks_sr = matrix(rbeta(constants_fit$N*constants_fit$M, 1, 1),
                     nrow = constants_fit$N, ncol = constants_fit$M)
  #alpha_thet = rep(1, constants$N)
)


EDP_model_fit <- nimbleModel(code_fit, constants_fit, data_fit, inits_fit)  # model creation


compile_EDP_model_fit <- compileNimble(EDP_model_fit, showCompilerOutput = TRUE)  # model compilation

# MCMC configuration
config_MCMC_fit <- configureMCMC(EDP_model_fit, useConjugacy = TRUE,
                                 enableWAIC = TRUE,
                                 monitors = c("beta_Yi", "sigma_sq_Yi", "pi0_Y",
                                              #"bi_M2_sig_sq",
                                              "theta_M2_t1", "theta_M2_t2", "theta_M2_t3",
                                              "sigma_sq_M2_t","pi0_M2_t",
                                              #"bi_M1_sig_sq",
                                              "theta_M1_t1", "theta_M1_t2", "theta_M1_t3",
                                              "sigma_sq_M1_t","pi0_M1_t",
                                              "theta_D_t1", "theta_D_t2", "theta_D_t3",
                                              "theta_Z_t",
                                              "pr_Bin_L0i",
                                              "breaks_r", "breaks_sr",
                                              "ksi_r","ksi_sr",
                                              "logDens"))


config_MCMC_fit$removeSamplers('logDens')   ## remove sampler assigned to 'logDens'
config_MCMC_fit$addSampler(target = 'logDens', type = 'sumLogPostDens')   ## add our custom sampler


# List of conjugate parameters: these parameters will be sampled from their conjugate posterior distributions
# Remove default samplers for the target parameters
config_MCMC_fit$removeSamplers("breaks_sr")


for (k in 1:num_outer_cluster_N) {
  for (j in 1:(num_inner_cluster_M-1)) {
    # number of obs in the k^th outer cluster and j^th inner cluster
    n_kj <- sum(EDP_model_fit$comp_num_outer == k &
                  EDP_model_fit$comp_num_inner == j)
    sum_n_kh <- 0
    for (h in (j+1):num_inner_cluster_M){
      sum_n_kh <- sum_n_kh + sum(EDP_model_fit$comp_num_outer == k &
                                   EDP_model_fit$comp_num_inner == h)
    }
    # Use n_kj in your custom sampler or computations
    config_MCMC_fit$addSampler(
      target = paste0("breaks_sr[", k, ",", j, "]"),
      type = "RW",
      control = list(
        posterior = "dbeta",
        shape1 = n_kj + 1,
        shape2 = alpha_thet_const + sum_n_kh
        #shape2 = EDP_model_fit$alpha_thet[k] + sum_n_kh  # Add other terms as needed
      )
    )
  }
}


mcmc_fit <- buildMCMC(config_MCMC_fit)
cmcmc_fit <- compileNimble(mcmc_fit, project = EDP_model_fit, showCompilerOutput = TRUE)

# burn-in 10000 iterations, thin by 5 to reduce autocorrelation
num_iter_fit = 15000;  num_burnin_fit = 10000;  num_thin_fit = 5;  num_chains_fit= 1


samples_allChains_fit <- runMCMC(cmcmc_fit, niter = num_iter_fit, nburnin = num_burnin_fit,
                                 nchains = num_chains_fit, 
                                 setSeed = TRUE, thin = num_thin_fit, WAIC = TRUE)

post_samples_model_fit <- samples_allChains_fit$samples   # list of nchains matrices


post_samples_model_fit <- as.data.frame(post_samples_model_fit)


# =====================================================================
# COMPUTE POSTERIOR MEAN OF NUMBER OF OUTER CLUSTERS USED
# =====================================================================
# Extract ksi_r columns (outer cluster probabilities)
ksi_r_col_names = grep("^ksi_r\\[\\d+\\]$", colnames(post_samples_model_fit), value = TRUE)
ksi_r_samples = post_samples_model_fit[, ksi_r_col_names]

# For each posterior sample (row), count clusters with weight >= 0.1
num_outer_used_per_sample = apply(ksi_r_samples, 1, function(row) sum(row >= 0.1))

# Compute posterior mean
mean_num_outer_used = mean(num_outer_used_per_sample)
cat("Posterior mean of outer clusters used (ksi_r >= 0.1):", mean_num_outer_used, "\n")
# =====================================================================


#######################################################################################
###############G-computation###########################################################
#######################################################################################


# =====================================================================
# PRE-COMPUTE COLUMN INDICES FOR FASTER LOOKUPS
# =====================================================================
# These are computed ONCE and reused across all MCInteg iterations

precompute_column_indices = function(BNPModels_df, N, M, num_Bin_L0) {
  # Pre-compute ksi_r column indices
  ksi_r_col_idx = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  
  
  # Pre-compute ksi_sr column indices as a matrix (M x N)
  ksi_sr_col_idx = matrix(NA, nrow = M, ncol = N)
  for(i in 1:N) {
    for(j in 1:M) {
      col_name = paste0("ksi_sr[", i, ", ", j, "]")
      ksi_sr_col_idx[j, i] = which(colnames(BNPModels_df) == col_name)
    }
  }
  
  # Pre-compute pr_Bin_L0i column indices as an array (N x M x num_Bin_L0)
  pr_Bin_L0_col_idx = array(NA, dim = c(N, M, num_Bin_L0))
  for(i in 1:N) {
    for(j in 1:M) {
      for(k in 1:num_Bin_L0) {
        col_name = paste0("pr_Bin_L0i[", i, ", ", j, ", ", k, "]")
        idx = which(colnames(BNPModels_df) == col_name)
        if(length(idx) > 0) pr_Bin_L0_col_idx[i, j, k] = idx
      }
    }
  }
  
  return(list(
    ksi_r_col_idx = ksi_r_col_idx,
    ksi_sr_col_idx = ksi_sr_col_idx,
    pr_Bin_L0_col_idx = pr_Bin_L0_col_idx
  ))
}

# =====================================================================
# END PRE-COMPUTE SECTION
# =====================================================================


# returns matrix of dimensions num_MC_samples by 2
compute_EDPM_cluster_component_number = function(
    BNPModels_df,
    num_MC_samples,# number of MC samples to be generated
    N,
    M,
    it,   # an integer from 1:Q (number of post iterations)
    precomp_idx = NULL  # pre-computed column indices (optional for backwards compatibility)
){
  # Use pre-computed indices if available, otherwise compute (slower)
  if(!is.null(precomp_idx)) {
    ksi_r_vec = as.numeric(BNPModels_df[it, precomp_idx$ksi_r_col_idx])
  } else {
    ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
    ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num])
  }
  
  outer_comp_vec = sample(x = 1:N, size = num_MC_samples,
                          replace = TRUE, prob = ksi_r_vec)
  
  # Use pre-computed indices for ksi_sr lookup
  if(!is.null(precomp_idx)) {
    # Fast path: use numeric indices directly
    ksi_sr_mat_iter = sapply(outer_comp_vec, function(outer_comp) {
      as.numeric(BNPModels_df[it, precomp_idx$ksi_sr_col_idx[, outer_comp]])
    })
  } else {
    # Original slow path
    ksi_sr_col_name = sapply(outer_comp_vec, function(outer_comp) {
      paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")
    })
    ksi_sr_mat_iter = apply(ksi_sr_col_name, c(1,2),
                            function(col_name){BNPModels_df[[col_name]][it]})
  }
  
  # vector of length num_MC_samples
  inner_comp_vec = apply(ksi_sr_mat_iter, 2,
                         function(ksi_sr){sample(x = c(1:M),size =1,
                                                 replace = TRUE,
                                                 prob = ksi_sr)})
  
  return_df = data.frame(outer_cluster = outer_comp_vec, inner_cluster = inner_comp_vec)
  return(return_df)
}


sample_L0_bin_single = function(MCdensity= NULL,
                                BNPModels_df,
                                EDPM_cluster_component_number_df,   # df of dim num_MC_samples times 2
                                num_MC_samples,                     # number of MC samples to be generated
                                N,
                                M,
                                it                                  # an integer from 1:Q (number of post iterations)
){
  
  
  #get the col num of all cols whose name begin with ksi_r                       
  ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  #get the elements of all cols whose name begin with ksi_r for this iter
  ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num]) # outer comp prob vec of length N
  
  # matrix of dimension M by num_MC_samples: a column for each MC sample
  ksi_sr_col_name = sapply(1:N, function(outer_comp) {
    paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")  # Column names for ksi_sr[outer_comp,1] to ksi_sr[outer_comp,M]
  })
  
  # matrix of dimension M by N: a column for each outer cluster
  ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                     function(col_name){BNPModels_df[[col_name]][it]}) 
  
  
  pr_Bin_L0_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("pr_Bin_L0i[", i, ", ", 1:M, "]"))
  )
  
  pr_Bin_L0_colnames = pr_Bin_L0_colnames_ordered[pr_Bin_L0_colnames_ordered %in%
                                                    colnames(BNPModels_df)]
  #vector of size N*M
  pr_Bin_L0_vec_all_clusters = as.numeric(sapply(pr_Bin_L0_colnames, 
                                                 function(col_name) BNPModels_df[[col_name]][it]))
  
  
  # Treat the nested finite mixture model as a single-level mixture model with N*M components and compute component probs
  total_component_prob = c()  #vector of size N*M
  for(i in 1:N){
    total_component_prob = c(total_component_prob, (ksi_r_vec[i] * as.numeric(ksi_sr_mat[,i])))
  }
  
  
  if (is.null(MCdensity)) {
    # Extract coefficients for the specific iteration: vector of length num_MCsamples
    param_coeff_pr_Bin_L0_vec = as.numeric(apply(EDPM_cluster_component_number_df,1,
                                                 function(x){BNPModels_df[[paste0("pr_Bin_L0i[", x[1], ", ", x[2], "]")]][it]}))
  }else{
    ### Sample component numbers (one for each MC sample) based on the total number of clusters (N*M)
    ## Compute weights/new mixture probabilities
    # Compute weights/new mixture probabilities
    weights_num = sweep(MCdensity, 2, total_component_prob, `*`)
    #cat("The number of 0s in  weights_num  vector is:", sum(weights_num == 0) , ".\n")
    weights_num = as.data.frame(weights_num)
    weights_denom=  rowSums(weights_num)
    
    zero_rows <- is.na(weights_denom) | (weights_denom == 0)
    
    if (any(zero_rows)) {
      ncomp <- ncol(weights_num)                    # N*M
      # set each problematic row to 1 for all comps (then will normalize below)
      weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
      # recompute denom
      weights_denom <- rowSums(weights_num)
    }
    
    weights_df = weights_num / weights_denom
    component_vec =  apply(weights_df, 1, function(probabilities) {
      sample(c(1:(N*M)), size = 1, prob = probabilities)
    })
    #print(component_vec)
    param_coeff_pr_Bin_L0_vec = pr_Bin_L0_vec_all_clusters[component_vec]
    #print(param_coeff_pr_Bin_L0_vec)
  }
  
  
  # Sample the binary outcome from a Bernoulli distribution
  random_samples =  rbinom(num_MC_samples, 1, param_coeff_pr_Bin_L0_vec)
  
  
  #density_dataframe: matrix if dim  num_MC_samples by (N*M) 
  density_dataframe = as.data.frame(sapply(pr_Bin_L0_vec_all_clusters, function(p) {
    dbinom(random_samples, size = 1, prob = p)
  }))
  
  
  return_list = list(
    samples =random_samples,
    density_df = density_dataframe
  )
  
  #return_df = data.frame(samples =random_samples, density =  density_values)
  
  return(return_list)
  
}


sample_L0_bin_multiple = function(MCdensity,
                                  BNPModels_df,
                                  EDPM_cluster_component_number_df,   # df of dim num_MC_samples times 2
                                  num_MC_samples,                     # number of MC samples to be generated
                                  N,
                                  M,
                                  L0_bin_number,
                                  it,                                 # an integer from 1:Q (number of post iterations)
                                  precomp_idx = NULL  # pre-computed column indices
){
  # Use pre-computed indices if available
  if(!is.null(precomp_idx)) {
    ksi_r_vec = as.numeric(BNPModels_df[it, precomp_idx$ksi_r_col_idx])
    # matrix of dimension M by N: use pre-computed indices
    ksi_sr_mat = matrix(as.numeric(BNPModels_df[it, as.vector(precomp_idx$ksi_sr_col_idx)]),
                        nrow = M, ncol = N)
    # vector of size N*M for pr_Bin_L0
    pr_Bin_L0_idx_vec = as.vector(t(precomp_idx$pr_Bin_L0_col_idx[,,L0_bin_number]))
    pr_Bin_L0_idx_vec = pr_Bin_L0_idx_vec[!is.na(pr_Bin_L0_idx_vec)]
    pr_Bin_L0_vec_all_clusters = as.numeric(BNPModels_df[it, pr_Bin_L0_idx_vec])
  } else {
    # Original slow path
    ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
    ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num])
    ksi_sr_col_name = sapply(1:N, function(outer_comp) {
      paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")
    })
    ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                       function(col_name){BNPModels_df[[col_name]][it]})
    pr_Bin_L0_colnames_ordered = unlist(
      lapply(1:N, function(i) paste0("pr_Bin_L0i[", i, ", ", 1:M,", ", L0_bin_number, "]"))
    )
    pr_Bin_L0_colnames = pr_Bin_L0_colnames_ordered[pr_Bin_L0_colnames_ordered %in%
                                                      colnames(BNPModels_df)]
    pr_Bin_L0_vec_all_clusters = as.numeric(sapply(pr_Bin_L0_colnames,
                                                   function(col_name) BNPModels_df[[col_name]][it]))
  }
  
  
  # Treat the nested finite mixture model as a single-level mixture model with N*M components and compute component probs
  total_component_prob = c()  #vector of size N*M
  for(i in 1:N){
    total_component_prob = c(total_component_prob, (ksi_r_vec[i] * as.numeric(ksi_sr_mat[,i])))
  }
  
  
  if (L0_bin_number == 1) {
    # Extract coefficients for the specific iteration: vector of length num_MCsamples
    param_coeff_pr_Bin_L0_vec = as.numeric(apply(EDPM_cluster_component_number_df,1,
                                                 function(x){
                                                   BNPModels_df[[paste0("pr_Bin_L0i[", x[1], ", ", x[2],", ", 1, "]")]][it]}))
  }else{
    ### Sample component numbers (one for each MC sample) based on the total number of clusters (N*M)
    ## Compute weights/new mixture probabilities
    weights_num = sweep(MCdensity, 2, total_component_prob, `*`)
    weights_num = as.data.frame(weights_num)
    weights_denom=  rowSums(weights_num)
    
    zero_rows <- is.na(weights_denom) | (weights_denom == 0)
    
    if (any(zero_rows)) {
      ncomp <- ncol(weights_num)                    # N*M
      # set each problematic row to 1 for all comps (then will normalize below)
      weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
      # recompute denom
      weights_denom <- rowSums(weights_num)
    }
    
    weights_df = weights_num / weights_denom
    
    
    component_vec =  apply(weights_df, 1, function(probabilities) {
      sample(c(1:(N*M)), size = 1, prob = probabilities)
    })
    param_coeff_pr_Bin_L0_vec = pr_Bin_L0_vec_all_clusters[component_vec]
  }
  
  
  # Sample the binary outcome from a Bernoulli distribution
  random_samples =  rbinom(num_MC_samples, 1, param_coeff_pr_Bin_L0_vec)
  
  
  #density_dataframe: matrix if dim  num_MC_samples by (N*M) 
  density_dataframe = as.data.frame(sapply(pr_Bin_L0_vec_all_clusters, function(p) {
    dbinom(random_samples, size = 1, prob = p)
  }))
  
  
  return_list = list(
    samples =random_samples,
    density_df = density_dataframe
  )
  
  
  return(return_list)
  
}


sample_L0_cont_single = function(MCdensity= NULL,
                                 BNPModels_df,
                                 EDPM_cluster_component_number_df,   # df of dim num_MC_samples times 2
                                 num_MC_samples,                     # number of MC samples to be generated
                                 N,
                                 M,
                                 number_of_L0_bin,
                                 it                                  # an integer from 1:Q (number of post iterations)
){
  
  
  #get the col num of all cols whose name begin with ksi_r                       
  ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  #get the elements of all cols whose name begin with ksi_r for this iter
  ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num]) # outer comp prob vec of length N
  
  # matrix of dimension M by num_MC_samples: a column for each MC sample
  ksi_sr_col_name = sapply(1:N, function(outer_comp) {
    paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")  # Column names for ksi_sr[outer_comp,1] to ksi_sr[outer_comp,M]
  })
  
  # matrix of dimension M by N: a column for each outer cluster
  ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                     function(col_name){BNPModels_df[[col_name]][it]}) 
  
  
  mu_Cont_L0_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("mu_Cont_L0i[", i, ", ", 1:M, "]"))
  )
  
  mu_Cont_L0_colnames = mu_Cont_L0_colnames_ordered[mu_Cont_L0_colnames_ordered %in%
                                                      colnames(BNPModels_df)]
  #vector of size N*M
  mu_Cont_L0_vec_all_clusters = as.numeric(sapply(mu_Cont_L0_colnames, 
                                                  function(col_name) BNPModels_df[[col_name]][it]))
  
  sigma_sq_Cont_L0_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("sigma_sq_Cont_L0i[", i, ", ", 1:M, "]"))
  )
  
  sigma_sq_Cont_L0_colnames =  sigma_sq_Cont_L0_colnames_ordered[
    sigma_sq_Cont_L0_colnames_ordered %in%
      colnames(BNPModels_df)]
  #vector of size N*M
  sigma_sq_Cont_L0_vec_all_clusters = as.numeric(sapply(sigma_sq_Cont_L0_colnames, 
                                                        function(col_name) BNPModels_df[[col_name]][it]))
  
  
  # treat the nested finite mixture model as a single-level mixture model with N*M  components and compute component probs
  total_component_prob = c()  #vector of size N*M
  for(i in 1:N){
    total_component_prob = c(total_component_prob, (ksi_r_vec[i] *
                                                      as.numeric(ksi_sr_mat[,i])))
  }
  
  
  if (number_of_L0_bin ==0) {
    # Extract coefficients for the specific iteration: vectors of length num_MCsamples
    param_coeff_mean_vec = as.numeric(apply(EDPM_cluster_component_number_df,1,
                                            function(x){BNPModels_df[[paste0("mu_Cont_L0i[", x[1], ", ", x[2], "]")]][it]}))
    
    param_coeff_var_vec = as.numeric(apply(EDPM_cluster_component_number_df,1,
                                           function(x){BNPModels_df[[paste0("sigma_sq_Cont_L0i[", x[1], ", ", x[2], "]")]][it]}))
    
  }else{
    ### Sample component numbers (one for each MC sample) based on the total number of clusters (N*M)
    ## Compute weights/new mixture probabilities
    # weights_num = as.data.frame(unlist(lapply(1:ncol(MCdensity),
    #                             function(i) MCdensity[,i] * total_component_prob[i])))
    # Compute weights/new mixture probabilities
    weights_num = sweep(MCdensity, 2, total_component_prob, `*`)
    #cat("The number of 0s in  weights_num  vector is:", sum(weights_num == 0) , ".\n")
    weights_num = as.data.frame(weights_num)
    weights_denom=  rowSums(weights_num)
    
    zero_rows <- is.na(weights_denom) | (weights_denom == 0)
    
    if (any(zero_rows)) {
      ncomp <- ncol(weights_num)                    # N*M
      # set each problematic row to 1 for all comps (then will normalize below)
      weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
      # recompute denom
      weights_denom <- rowSums(weights_num)
    }
    
    weights_df = weights_num / weights_denom
    # Handle rows where weights_denom was zero
    #weights[zero_rows, ] = 1 / ncol(weights)  # Assign equal weights_dfacross columns
    component_vec =  apply(weights_df, 1, function(probabilities) {
      sample(c(1:(N*M)), size = 1, prob = probabilities)
    })
    param_coeff_mean_vec = mu_Cont_L0_vec_all_clusters[component_vec]
    #print(param_coeff_mean_vec)
    param_coeff_var_vec = sigma_sq_Cont_L0_vec_all_clusters[component_vec]
    #print(param_coeff_var_vec)
  }
  
  
  # Sample the continuous outcome from a Normal distribution
  random_samples = rtruncnorm(n = num_MC_samples, a = 0, b = Inf, 
                              mean = param_coeff_mean_vec, sd = sqrt(param_coeff_var_vec))
  
  
  density_dataframe=  as.data.frame(sapply(1:(N*M), function(i) {
    dtruncnorm(x = random_samples,a = 0, b = Inf,
               mean = mu_Cont_L0_vec_all_clusters[i], 
               sd = sqrt(sigma_sq_Cont_L0_vec_all_clusters[i]))
  }))
  
  
  return_list = list(
    samples =random_samples,
    density_df = density_dataframe 
  )
  return(return_list)
  
}


sample_L0_cont_multiple = function(MCdensity= NULL,
                                   BNPModels_df,
                                   EDPM_cluster_component_number_df,   # df of dim num_MC_samples times 2
                                   num_MC_samples,                     # number of MC samples to be generated
                                   N,
                                   M,
                                   L0_cont_number,
                                   number_of_L0_bin,
                                   it                                  # an integer from 1:Q (number of post iterations)
){
  
  
  #get the col num of all cols whose name begin with ksi_r                       
  ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  #get the elements of all cols whose name begin with ksi_r for this iter
  ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num]) # outer comp prob vec of length N
  
  # matrix of dimension M by num_MC_samples: a column for each MC sample
  ksi_sr_col_name = sapply(1:N, function(outer_comp) {
    paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")  # Column names for ksi_sr[outer_comp,1] to ksi_sr[outer_comp,M]
  })
  
  # matrix of dimension M by N: a column for each outer cluster
  ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                     function(col_name){BNPModels_df[[col_name]][it]}) 
  
  
  mu_Cont_L0_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("mu_Cont_L0i[", i, ", ", 1:M,", ", L0_cont_number, "]"))
  )
  
  mu_Cont_L0_colnames = mu_Cont_L0_colnames_ordered[mu_Cont_L0_colnames_ordered %in%
                                                      colnames(BNPModels_df)]
  #vector of size N*M
  mu_Cont_L0_vec_all_clusters = as.numeric(sapply(mu_Cont_L0_colnames, 
                                                  function(col_name) BNPModels_df[[col_name]][it]))
  
  sigma_sq_Cont_L0_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("sigma_sq_Cont_L0i[", i, ", ", 1:M,", ", L0_cont_number, "]"))
  )
  
  sigma_sq_Cont_L0_colnames =  sigma_sq_Cont_L0_colnames_ordered[
    sigma_sq_Cont_L0_colnames_ordered %in%
      colnames(BNPModels_df)]
  #vector of size N*M
  sigma_sq_Cont_L0_vec_all_clusters = as.numeric(sapply(sigma_sq_Cont_L0_colnames, 
                                                        function(col_name) BNPModels_df[[col_name]][it]))
  
  
  # treat the nested finite mixture model as a single-level mixture model with N*M  components and compute component probs
  total_component_prob = c()  #vector of size N*M
  for(i in 1:N){
    total_component_prob = c(total_component_prob, (ksi_r_vec[i] *
                                                      as.numeric(ksi_sr_mat[,i])))
  }
  
  
  if (number_of_L0_bin ==0) {
    # Extract coefficients for the specific iteration: vectors of length num_MCsamples
    param_coeff_mean_vec = as.numeric(apply(EDPM_cluster_component_number_df,1,
                                            function(x){
                                              BNPModels_df[[paste0("mu_Cont_L0i[", x[1], ", ", x[2],", ", L0_cont_number, "]")]][it]}))
    
    param_coeff_var_vec = as.numeric(apply(EDPM_cluster_component_number_df,1,
                                           function(x){BNPModels_df[[paste0("sigma_sq_Cont_L0i[", x[1], ", ", x[2],", ", L0_cont_number, "]")]][it]}))
    
  }else{
    ### Sample component numbers (one for each MC sample) based on the total number of clusters (N*M)
    ## Compute weights/new mixture probabilities
    
    # Compute weights/new mixture probabilities
    weights_num = sweep(MCdensity, 2, total_component_prob, `*`)
    #cat("The number of 0s in  weights_num  vector is:", sum(weights_num == 0) , ".\n")
    weights_num = as.data.frame(weights_num)
    weights_denom=  rowSums(weights_num)
    
    zero_rows <- is.na(weights_denom) | (weights_denom == 0)
    
    if (any(zero_rows)) {
      ncomp <- ncol(weights_num)                    # N*M
      # set each problematic row to 1 for all comps (then will normalize below)
      weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
      # recompute denom
      weights_denom <- rowSums(weights_num)
    }
    
    weights_df = weights_num / weights_denom
    
    # Handle rows where weights_denom was zero
    #weights[zero_rows, ] = 1 / ncol(weights)  # Assign equal weights_dfacross columns
    component_vec =  apply(weights_df, 1, function(probabilities) {
      sample(c(1:(N*M)), size = 1, prob = probabilities)
    })
    param_coeff_mean_vec = mu_Cont_L0_vec_all_clusters[component_vec]
    #print(param_coeff_mean_vec)
    param_coeff_var_vec = sigma_sq_Cont_L0_vec_all_clusters[component_vec]
    #print(param_coeff_var_vec)
  }
  
  
  # Sample the continuous outcome from a Normal distribution
  random_samples = rtruncnorm(n = num_MC_samples, a = 0, b = Inf, 
                              mean = param_coeff_mean_vec, sd = sqrt(param_coeff_var_vec))
  
  
  density_dataframe=  as.data.frame(sapply(1:(N*M), function(i) {
    dtruncnorm(x = random_samples,a = 0, b = Inf,
               mean = mu_Cont_L0_vec_all_clusters[i], 
               sd = sqrt(sigma_sq_Cont_L0_vec_all_clusters[i]))
  }))
  
  
  return_list = list(
    samples =random_samples,
    density_df = density_dataframe 
  )
  return(return_list)
  
}


get_Zt_density_df = function(x,                                  #vector at which we want to compute the density 
                             MCdata,
                             BNPModels_df,
                             num_MC_samples,                     # number of MC samples to be generated
                             N,
                             M,
                             num_L0,                             # number of L0 INCLUDING intercept
                             t,                                  # t = {1, ..., T}
                             it                                  # an integer from 1:Q (number of post iterations)
){
  
  #extract coefficients for the specific iteration ( param_coeff_theta_Zt_df is a num_components * num_L0 dim. matrix)
  param_coeff_theta_Zt_df = matrix(nrow = N*M, ncol = num_L0)
  for(k in 1:num_L0){
    theta_Zt_colnames_ordered = unlist(
      lapply(1:N, function(i) paste0("theta_Z_t[", i, ", ", 1:M,", ", k,", ", t, "]"))
    )
    
    theta_Zt_colnames =  theta_Zt_colnames_ordered[theta_Zt_colnames_ordered %in%
                                                     colnames(BNPModels_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_theta_Zt_df[,k] = as.numeric(sapply(theta_Zt_colnames, 
                                                    function(col_name) BNPModels_df[[col_name]][it]))
  }
  
  
  # sample random intercept for all MCsamples
  #param_coeff_bi_Zt_sig_sq = BNPModels_df$bi_Z_sig_sq[it]
  #bi_Zt = rnorm(num_MC_samples, mean = 0, sd = sqrt(param_coeff_bi_Zt_sig_sq))
  
  linpred_df = (as.matrix(MCdata[,1:num_L0]) %*% t(as.matrix(param_coeff_theta_Zt_df))) #+ bi_Zt # this is num_MCsample * MN dim df
  
  #density_df_prob = 1 / (1 + exp(-linpred_df))      #do this for logit regression
  density_df_prob = pnorm(linpred_df)                #do this for probit regression
  
  
  density_dataframe = as.data.frame(apply(density_df_prob,2, function(p) {
    dbinom(x, size = 1, prob = p)              # compute num_MCsample many density values at each col of density_df_prob
  }))
  
  
  return(density_dataframe)
  
  
}


sample_Dt = function(MCdata,                             # MC data upto this point in the temporal order
                     D1_vec,                             # vector of fixed D1 values
                     MCdensity,                          # MC density upto this point in the temporal order
                     BNPModels_df,
                     num_MC_samples,                     # number of MC samples to be generated
                     N,
                     M,
                     num_D_covs,                         # number of D covariates INCLUDING intercept (varies by t)
                     t,                                  # t = {1, ..., T}
                     it                                  # an integer from 1:Q (number of post iterations)
){
  
  #extract coefficients for the specific iteration ( param_coeff_theta_Dt_df is a num_components * num_D_covs dim. matrix)
  param_coeff_theta_Dt_df = matrix(nrow = N*M, ncol = num_D_covs)
  for(k in 1:num_D_covs){
    theta_Dt_colnames_ordered = unlist(
      lapply(1:N, function(i) paste0("theta_D_t", t, "[", i, ", ", 1:M,", ", k, "]"))
    )
    
    theta_Dt_colnames =  theta_Dt_colnames_ordered[theta_Dt_colnames_ordered %in%
                                                     colnames(BNPModels_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_theta_Dt_df[,k] = as.numeric(sapply(theta_Dt_colnames, 
                                                    function(col_name) BNPModels_df[[col_name]][it]))
  }
  
  
  #get the col num of all cols whose name begin with ksi_r                       
  ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  #get the elements of all cols whose name begin with ksi_r for this iter
  ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num]) # outer comp prob vec of length N
  
  # matrix of dimension M by num_MC_samples: a column for each MC sample
  ksi_sr_col_name = sapply(1:N, function(outer_comp) {
    paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")  # Column names for ksi_sr[outer_comp,1] to ksi_sr[outer_comp,M]
  })
  
  # matrix of dimension M by N: a column for each outer cluster, sum(col) =1
  ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                     function(col_name){BNPModels_df[[col_name]][it]}) 
  #ksi_sr_df = as.data.frame(t(ksi_sr_mat)) # nrow = N, ncol = M, sum(row) =1
  
  
  # sample random intercept for all MCsamples
  #param_coeff_bi_Dt_sig_sq = BNPModels_df$bi_D_sig_sq[it]
  #bi_Dt = rnorm(num_MC_samples, mean = 0, sd = sqrt(param_coeff_bi_Dt_sig_sq))
  
  # treat the nested finite mixture model as a single-level mixture model with N*M components and compute component probs
  total_component_prob = c()  #vector of size N*M
  for(i in 1:N){
    total_component_prob = c(total_component_prob, (ksi_r_vec[i] *
                                                      as.numeric(ksi_sr_mat[,i])))
  }
  
  
  # Compute weights/new mixture probabilities
  weights_num = sweep(MCdensity, 2, total_component_prob, `*`)
  weights_num = as.data.frame(weights_num)
  weights_denom=  rowSums(weights_num)
  
  zero_rows <- is.na(weights_denom) | (weights_denom == 0)
  
  if (any(zero_rows)) {
    ncomp <- ncol(weights_num)                    # N*M
    # set each problematic row to 1 for all comps (then will normalize below)
    weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
    # recompute denom
    weights_denom <- rowSums(weights_num)
  }
  
  weights_df = weights_num / weights_denom
  
  
  component_vec =  apply(weights_df, 1, function(probabilities) {
    sample(c(1:(N*M)), size = 1, prob = probabilities)
  })
  
  
  # Extract sampling probability params vector (size num_D_covs) based on the component number drawn by using the component probabilities
  theta_Dt_sample_df = param_coeff_theta_Dt_df[component_vec,]  # this is a df with nrows = num_MC_samples, ncol = num_D_covs at iter it and time t
  # Compute sampling probability
  
  # do this for logit regression
  #theta_Dt_prob = 1 / (1 + exp(-(rowSums(MCdata[,1:num_D_covs] *theta_Dt_sample_df) + bi_Dt)))
  
  #message("reached Dt before sample generation \n")
  # do this for probit regression
  theta_Dt_prob = pnorm(rowSums(MCdata[,1:num_D_covs] *theta_Dt_sample_df)) # + bi_Dt)
  
  
  # Sample the binary outcome from a Bernoulli distribution
  random_samples =  rbinom(num_MC_samples, 1, theta_Dt_prob)
  
  linpred_df = (as.matrix(MCdata[,1:num_D_covs]) %*% t(as.matrix(param_coeff_theta_Dt_df))) #+ bi_Dt # this is num_MCsample * MN dim df
  
  
  density_df_prob = pnorm(linpred_df)  # do this for probit regression
  
  if(t==1){
    density_dataframe = as.data.frame(apply(density_df_prob,2, function(p) {
      dbinom(D1_vec, size = 1, prob = p)              # compute num_MCsample many density values at each col of density_df_prob
    }))
    return_list = list(
      density_df = density_dataframe
    )
    return(return_list)
    
  }else{
    density_dataframe = as.data.frame(apply(density_df_prob,2, function(p) {
      dbinom(random_samples, size = 1, prob = p)              # compute num_MCsample many density values at each col of density_df_prob
    }))
    
    return_list = list(
      samples =random_samples,
      density_df = density_dataframe
    )
    return(return_list)
  }
  
  
}


rnorm_hurdle = function(n,mean,sd, hurdle_prob){
  hurdle_comp = rbinom(length(hurdle_prob),1,hurdle_prob)
  return_sample = (1-hurdle_comp)* rtruncnorm(n = length(hurdle_comp), a = 1, b=10, mean =mean, sd = sd )
  return(as.numeric(return_sample))
}

dnorm_hurdle = function(x,mean,sd, hurdle_prob){
  density = ifelse(x ==0, hurdle_prob,
                   ((1-hurdle_prob)*dtruncnorm(x = x,a = 1, b=10, mean =mean, sd = sd)))
  return(as.numeric(density))
}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


sample_hurdle_data = function(MCdata,     # MC data upto this point in the temporal order
                              MCdensity,           # MC density upto this point in the temporal order
                              BNPModels_df,
                              Mediator_num,        # mediator location in temporal order 1,...,J
                              num_MC_samples,      # number of MC samples to be generated
                              N,
                              M,
                              num_M_covs,           # number of M covariates INCLUDING intercept (varies by t)
                              t,                   # t = {1, ..., T}
                              it                   # an integer from 1:Q (number of post iterations)
){
  
  
  #extract coefficients for the specific iteration ( param_coeff_mean_df is a num_components * num_M_covs dim. matrix)
  param_coeff_mean_df = matrix(nrow = N*M, ncol = num_M_covs)
  for(k in 1:num_M_covs){
    theta_Mt_colnames_ordered = unlist(
      lapply(1:N, function(i) paste0("theta_M",Mediator_num,"_t", t, "[",
                                     i, ", ", 1:M,", ", k, "]"))
    )
    
    theta_Mt_colnames =  theta_Mt_colnames_ordered[theta_Mt_colnames_ordered %in%
                                                     colnames(BNPModels_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_mean_df[,k] = as.numeric(sapply(theta_Mt_colnames, 
                                                function(col_name) BNPModels_df[[col_name]][it]))
  }
  
  
  sigSq_Mt_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("sigma_sq_M",Mediator_num,"_t[", 
                                   i, ", ", 1:M,", ", t, "]"))
  )
  
  sigSq_Mt_colnames =  sigSq_Mt_colnames_ordered[sigSq_Mt_colnames_ordered %in%
                                                   colnames(BNPModels_df)]
  # vector of  length NM
  param_coeff_var_vec = as.numeric(sapply(sigSq_Mt_colnames, 
                                          function(col_name) BNPModels_df[[col_name]][it]))
  
  
  pi0_Mt_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("pi0_M",Mediator_num,"_t[", 
                                   i, ", ", 1:M,", ", t, "]"))
  )
  
  pi0_Mt_colnames =  pi0_Mt_colnames_ordered[pi0_Mt_colnames_ordered %in%
                                               colnames(BNPModels_df)]
  # vector of  length NM
  param_coeff_hurdle_prob_vec = as.numeric(sapply(pi0_Mt_colnames, 
                                                  function(col_name) BNPModels_df[[col_name]][it]))
  
  
  #get the col num of all cols whose name begin with ksi_r                       
  ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  #get the elements of all cols whose name begin with ksi_r for this iter
  ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num]) # outer comp prob vec of length N
  
  # matrix of dimension M by num_MC_samples: a column for each MC sample
  ksi_sr_col_name = sapply(1:N, function(outer_comp) {
    paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")  # Column names for ksi_sr[outer_comp,1] to ksi_sr[outer_comp,M]
  })
  
  # matrix of dimension M by N: a column for each outer cluster, sum(col) =1
  ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                     function(col_name){BNPModels_df[[col_name]][it]}) 
  #ksi_sr_df = as.data.frame(t(ksi_sr_mat)) # nrow = N, ncol = M, sum(row) =1
  
  
  # sample random intercept for all MCsamples 
  #bi_Mt_col_name = paste0("bi_M",Mediator_num,"_sig_sq")
  #param_coeff_bi_Mt_sig_sq = BNPModels_df[[bi_Mt_col_name]][it]
  #bi_vec = rnorm(num_MC_samples, mean = 0, sd = sqrt(param_coeff_bi_Mt_sig_sq))
  
  
  # treat the nested finite mixture model as a single-level mixture model with N*M components and compute component probs
  total_component_prob = c()  #vector of size N*M
  for(i in 1:N){
    total_component_prob = c(total_component_prob, (ksi_r_vec[i] *
                                                      as.numeric(ksi_sr_mat[,i])))
  }
  
  
  # Compute weights/new mixture probabilities
  weights_num = sweep(MCdensity, 2, total_component_prob, `*`)
  #cat("The number of 0s in  weights_num  vector is:", sum(weights_num == 0) , ".\n")
  weights_num = as.data.frame(weights_num)
  weights_denom=  rowSums(weights_num)
  
  zero_rows <- is.na(weights_denom) | (weights_denom == 0)
  
  if (any(zero_rows)) {
    ncomp <- ncol(weights_num)                    # N*M
    # set each problematic row to 1 for all comps (then will normalize below)
    weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
    # recompute denom
    weights_denom <- rowSums(weights_num)
  }
  
  weights_df = weights_num / weights_denom
  
  
  component_vec =  apply(weights_df, 1, function(probabilities) {
    sample(c(1:(N*M)), size = 1, prob = probabilities)
  })
  
  
  # Extract sampling mean params vector (size num_M_covs) based on the component number drawn by using the component probabilities
  mean_sample_df = param_coeff_mean_df[component_vec,]  # this is a df with nrows = num_MC_samples, ncol = num_M_covs at iter it and time t
  # Compute sampling mean
  mean_sample_vec = rowSums(MCdata[,1:num_M_covs] *mean_sample_df) #+ bi_vec
  
  # Extract sampling var param based on the component number drawn by using the component probabilities
  var_sample_vec = param_coeff_var_vec[component_vec]
  
  # Extract sampling hurdle prob param based on the component number drawn by using the component probabilities
  
  hurdle_prob_sample_vec = param_coeff_hurdle_prob_vec[component_vec]
  
  
  random_samples =  rnorm_hurdle(num_MC_samples, mean_sample_vec , sqrt(var_sample_vec),hurdle_prob_sample_vec)
  
  
  density_df_mean = as.matrix(MCdata[,1:num_M_covs]) %*% t(as.matrix(param_coeff_mean_df)) #+ bi_vec # this is num_MCsample * MN dim df
  
  
  density_dataframe=  as.data.frame(sapply(1:(N*M), function(i) {
    dnorm_hurdle(random_samples, mean = density_df_mean[,i], 
                 sd = sqrt(param_coeff_var_vec[i]),
                 hurdle_prob =param_coeff_hurdle_prob_vec[i])
  }))
  
  return_list = list(
    samples =random_samples,
    density_df = density_dataframe
  )
  
  return(return_list)
}


compute_Y_mean = function(MCdata,                             # MC data upto this point in the temporal order
                          MCdensity,                      # MC density upto this point in the temporal order
                          BNPModels_df,
                          num_MC_samples,                     # number of MC samples to be generated
                          N,
                          M,
                          num_covariates,               # number of local regression INCLUDING intercept
                          it                                  # an integer from 1:Q (number of post iterations)
){
  
  #extract coefficients for the specific iteration ( param_coeff_mean_df is a num_components * num_covariates dim. matrix)
  beta_Y_colnames_ordered = unlist(lapply(1:N, function(i) paste0("beta_Yi[", 
                                                                  i, ", ", 1:num_covariates, "]")))
  beta_Y_colnames =  beta_Y_colnames_ordered[beta_Y_colnames_ordered %in%
                                               colnames(BNPModels_df)]
  #df with nrow = N, ncol = num_covariates
  param_coeff_mean_df = matrix(sapply(beta_Y_colnames, 
                                      function(col_name) BNPModels_df[[col_name]][it]),
                               nrow = N, ncol = num_covariates, byrow = TRUE)
  
  
  sigSq_Y_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("sigma_sq_Yi[", 
                                   i, "]"))
  )
  
  sigSq_Y_colnames =  sigSq_Y_colnames_ordered[sigSq_Y_colnames_ordered %in%
                                                 colnames(BNPModels_df)]
  # vector of  length N
  param_coeff_var_vec = as.numeric(sapply(sigSq_Y_colnames, 
                                          function(col_name) BNPModels_df[[col_name]][it]))
  
  
  pi0_Y_colnames_ordered = unlist(
    lapply(1:N, function(i) paste0("pi0_Y[", i,  "]"))
  )
  
  pi0_Y_colnames =  pi0_Y_colnames_ordered[pi0_Y_colnames_ordered %in%
                                             colnames(BNPModels_df)]
  # vector of  length N
  param_coeff_hurdle_prob_vec = as.numeric(sapply(pi0_Y_colnames, 
                                                  function(col_name) BNPModels_df[[col_name]][it]))
  
  
  #get the col num of all cols whose name begin with ksi_r                       
  ksi_r_col_num = grep("^ksi_r\\[\\d+\\]$", colnames(BNPModels_df))
  #get the elements of all cols whose name begin with ksi_r for this iter
  ksi_r_vec = as.numeric(BNPModels_df[it, ksi_r_col_num]) # outer comp prob vec of length N
  
  # matrix of dimension M by num_MC_samples: a column for each MC sample
  ksi_sr_col_name = sapply(1:N, function(outer_comp) {
    paste0("ksi_sr[", outer_comp, ", ", 1:M, "]")  # Column names for ksi_sr[outer_comp,1] to ksi_sr[outer_comp,M]
  })
  
  # matrix of dimension M by N: a column for each outer cluster, sum(col) =1
  ksi_sr_mat = apply(ksi_sr_col_name, c(1,2),
                     function(col_name){BNPModels_df[[col_name]][it]})
  
  ksi_sr_df = as.data.frame(t(ksi_sr_mat)) # nrow = N, ncol = M, sum(row) =1
  
  
  # Compute weights/new mixture probabilities
  weights_num = as.data.frame(matrix(nrow = num_MC_samples, ncol = N))
  for (i in 1:N) {
    start_col = (i - 1) * M + 1           #{1, M + 1, 2M + 1, 3M + 1, ... }
    end_col = i * M                       #{M, 2M, 3M, 4M, ...}
    #MCdensity_times_ksi_sr[,start_col:end_col] = MCdensity[, start_col:end_col] * ksi_sr_df[i,]
    # Sum the rows for each block of M columns
    weights_num[, i] = rowSums(MCdensity[, start_col:end_col]* 
                                 as.numeric(ksi_sr_df[i,])) * 
      ksi_r_vec[i]
  }
  
  
  weights_denom=  rowSums(weights_num)
  
  zero_rows <- is.na(weights_denom) | (weights_denom == 0)
  
  if (any(zero_rows)) {
    ncomp <- ncol(weights_num)                    # N*M
    # set each problematic row to 1 for all comps (then will normalize below)
    weights_num[zero_rows, ] <- matrix(1, nrow = sum(zero_rows), ncol = ncomp, byrow = TRUE)
    # recompute denom
    weights_denom <- rowSums(weights_num)
  }
  weights_df = weights_num/weights_denom
  # Handle rows where weights_denom was zero
  #weights[zero_rows, ] = 1 / ncol(weights)  # Assign equal weights_dfacross columns
  
  
  #sample the outer component number for each num_MC_samples subjects
  outer_comp_num = apply(weights_df, 1, function(w) sample(1:ncol(weights_df), size = 1, prob = w))
  # vectors of size num_MC_samples
  param_coeff_hurdle_prob_MC_samples = param_coeff_hurdle_prob_vec[outer_comp_num]
  # binary variable with 0 indicating Y =0 and 1 indicating Y~Normal()
  # vector of size num_MC_samples
  
  
  # this is num_MCsample * N dim df of untruncated normal mean
  density_df_mean =as.matrix(MCdata[,1:num_covariates]) %*% t(as.matrix(param_coeff_mean_df)) 
  
  
  # compute mean for truncated Normal
  z = sweep(-density_df_mean, 2, sqrt(param_coeff_var_vec), `/`)
  phi_z = dnorm(z)  # PDF of the standard normal at z
  Phi_z = pnorm(z)  # CDF of the standard normal at z
  trunc_norm_param_denom = pmax(1 - Phi_z, 1e-10)
  density_df_mean_trunc_norm = density_df_mean + sweep(phi_z / trunc_norm_param_denom, 2, 
                                                       sqrt(param_coeff_var_vec), FUN = "*")
  
  
  # E[Y|...] = sum_j(w_j * (1-pi0_j) * E[trunc_norm_j])
  density_df_mean_Y_hurdle = sweep(density_df_mean_trunc_norm, 2, 1-param_coeff_hurdle_prob_vec, FUN = "*")
  mean_Y = rowSums(density_df_mean_Y_hurdle * weights_df)

  return(mean_Y)
}


####################################################################
######.............Function for MC integration .............########
MCInteg <- function(var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), D=trt received, M1= mediator 1,  M2= mediator 2, Y=outcome. 
                    BNPModels, #fitted BNP model parameter samples list
                    fixed.regime.trt, # A vector specifying the Exposure  regime Z
                    fixed.regime.control, # A vector specifying the Exposure  regime Z_star
                    U1,       # could be one of: "c" for VALUE ATTENTIVE, "d" for PRICE ATTENTIVE, "n" for NON-ACTIVE CUSTOMERS and "a" for ACTIVE CUSTOMERS
                    J=10000, # Size of pseudo data. Default is set to 2,000.
                    Ndraws=1000, # Number of posterior draws. Default is set to 200.
                    num_outer_clusters,
                    num_inner_clusters,
                    num_Bin_L0,
                    num_Cont_L0,
                    ...
){
  
  n_Y <- length(which(var.type=="Y")) # Number of outcome variables, should be equal to 1
  n_L0 <- length(which(var.type=="L0")) # Number of baseline confounder variables INCLUDING the intercept, >1 allowed
  n_Fi <- length(which(var.type=="Fi")) # Number of fixed variables, >1 allowed
  if(!is.null(fixed.regime.trt)) {
    n_Reg_trt <- max(ifelse(!is.null(fixed.regime.trt), length(fixed.regime.trt), 0))
  } else {
    n_Reg_trt <- 1
  }
  
  if(!is.null(fixed.regime.control)) {
    n_Reg_control <- max(ifelse(!is.null(fixed.regime.control), length(fixed.regime.control), 0))
  } else {
    n_Reg_control <- 1
  }
  
  if (n_Y != 1) stop("The number of outcome variables is not equal to 1 ")
  #if (ncol(data) != length(var.type)) stop("The number of columns in data is not equal to the length of var.type ")
  
  if(!is.null(fixed.regime.trt)){
    if (length(fixed.regime.trt) != n_Fi) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
  }
  
  if(!is.null(fixed.regime.control)){
    if (length(fixed.regime.control) != n_Fi) stop("Warning: The number of fixed variables is not equal to the length of the fixed regime(s)")
  }
  
  if (!(U1 %in% c("c", "d", "n", "a"))) {
    stop("U1 must be one of 'c', 'd', 'n', or 'a'")
  }
  
  
  ########################################################
  ## Create matrices to store the output.
  
  
  Y_zzstar_mat <- matrix(nrow =  1 , ncol = Ndraws)
  MC_error_vec = c()
  
  # =====================================================================
  # PRE-COMPUTE COLUMN INDICES ONCE (before the loop)
  # =====================================================================
  precomp_idx = precompute_column_indices(BNPModels, num_outer_clusters,
                                          num_inner_clusters, num_Bin_L0)
  
  
  #L0data <- subset_data[id==1,c(1:n_L0)]
  for(it in 1:Ndraws) {
    #m <- 2                               # m indexes rows of the outcome matrix P_zzstar
    zt <- 1                               # zt indexes FIXED longitudinal trt regime in G-comp
    dt <- 1
    mt <- 1
    intercept = rep(1, J)
    # initialize a matrix to store density values used for computing cluster component probs
    density_prod_df = matrix(1, nrow = J, ncol = num_outer_clusters* num_inner_clusters)
    # sample data for the first binary baseline confounder (Add a loop for L0 data if the number of L0 data increase)
    EDPM_cluster_component_number = compute_EDPM_cluster_component_number(
      BNPModels_df = BNPModels,
      num_MC_samples = J,
      N = num_outer_clusters,
      M = num_inner_clusters,
      it = it,
      precomp_idx = precomp_idx
    )
    
    if(num_Bin_L0 ==1){
      L0_bin_sampled_vec_iter = sample_L0_bin_single(MCdensity= NULL,
                                                     BNPModels_df= BNPModels,
                                                     EDPM_cluster_component_number_df = EDPM_cluster_component_number,
                                                     num_MC_samples = J,
                                                     N = num_outer_clusters,
                                                     M = num_inner_clusters,
                                                     it = it)
      L0_data = cbind(intercept, L0_bin_sampled_vec_iter$samples)
      density_prod_df = density_prod_df * L0_bin_sampled_vec_iter$density_df
      
    } else if(num_Bin_L0 >1){
      L0_bin_sampled_mat_iter = NULL # initialize a matrix to store samples
      for (k in 1: num_Bin_L0){
        L0_bin_sampled_vec_iter_k = sample_L0_bin_multiple(
          MCdensity= density_prod_df,
          BNPModels_df= BNPModels,
          EDPM_cluster_component_number_df = EDPM_cluster_component_number,
          num_MC_samples = J,
          N = num_outer_clusters,
          M = num_inner_clusters,
          L0_bin_number = k,
          it = it,
          precomp_idx = precomp_idx)
        L0_bin_sampled_mat_iter = cbind(L0_bin_sampled_mat_iter,
                                        L0_bin_sampled_vec_iter_k$samples)
        density_prod_df = density_prod_df * L0_bin_sampled_vec_iter_k$density_df
      }
      
      L0_data = cbind(intercept, L0_bin_sampled_mat_iter)
    }
    
    
    # sample data for the continuous baseline confounder
    if(num_Cont_L0 ==1){
      L0_cont_sampled_vec_iter = sample_L0_cont_single(MCdensity= density_prod_df,
                                                       BNPModels_df= BNPModels,
                                                       EDPM_cluster_component_number_df = EDPM_cluster_component_number,
                                                       num_MC_samples = J,
                                                       N = num_outer_clusters,
                                                       M = num_inner_clusters,
                                                       number_of_L0_bin = num_Bin_L0,
                                                       it = it)
      L0_data = cbind(L0_data, L0_cont_sampled_vec_iter$samples)
      density_prod_df = density_prod_df * L0_cont_sampled_vec_iter$density_df
      
    } else if(num_Cont_L0 >1){
      L0_cont_sampled_mat_iter = NULL # initialize a matrix to store samples
      for (k in 1: num_Cont_L0){
        L0_cont_sampled_vec_iter_k = sample_L0_cont_multiple(
          MCdensity= density_prod_df,
          BNPModels_df= BNPModels,
          EDPM_cluster_component_number_df = EDPM_cluster_component_number,
          num_MC_samples = J,
          N = num_outer_clusters,
          M = num_inner_clusters,
          L0_cont_number = k,
          number_of_L0_bin = num_Bin_L0,
          it = it)
        L0_cont_sampled_mat_iter = cbind(L0_cont_sampled_mat_iter,
                                         L0_cont_sampled_vec_iter_k$samples)
        density_prod_df = density_prod_df * L0_cont_sampled_vec_iter_k$density_df
      }
      
      L0_data = cbind(L0_data, L0_cont_sampled_mat_iter)
    }
    
    
    # design matrix for Y regression (initialize with L0_data)
    Y_design_mat = L0_data
    # Zt_df = as.data.frame(matrix(rep(fixed.regime.trt, each = J), nrow = J, ncol = length(fixed.regime.trt)))  
    # Y_design_mat = cbind(L0_data, Zt_df)   
    
    
    Z1_trt = rep(fixed.regime.trt[1], J)
    Z1_control = rep(fixed.regime.control[1], J)
    if (U1 == "c") {
      D1_trt = Z1_trt
      D1_control = Z1_control
    } else if (U1 == "d") {
      D1_trt = 1- Z1_trt
      D1_control = 1-Z1_control
    } else if (U1 == "n") {
      D1_trt = rep(0,J)
      D1_control = rep(0,J)
    } else if (U1 == "a") {
      D1_trt = rep(1,J)
      D1_control = rep(1,J)
    }
    
    
    # Derive K vectors from var.type: count preceding Fi/Dt entries for each Dt/M position
    dt_positions = which(var.type == "Dt")
    m2t_positions = which(var.type == "M2t")
    K_D_vec = sapply(dt_positions, function(pos) n_L0 + sum(var.type[1:pos] == "Fi"))
    K_M_vec = sapply(m2t_positions, function(pos) n_L0 + sum(var.type[1:pos] == "Fi") + sum(var.type[1:pos] == "Dt"))
    
    # History tracking matrices (built up as we iterate through time)
    Z_history_trt = NULL      # will accumulate Z columns for trt regime
    Z_history_control = NULL  # will accumulate Z columns for control regime
    D_history_trt = NULL      # will accumulate D columns for trt regime
    D_history_control = NULL  # will accumulate D columns for control regime
    
    for (j in (n_L0+1):length(var.type)) {
      
      if(var.type[j] == "Fi") {
        Zt_density_df_trt = get_Zt_density_df(x = rep(fixed.regime.trt[zt], J),
                                              MCdata = L0_data,
                                              BNPModels_df = BNPModels,
                                              num_MC_samples = J,
                                              N = num_outer_clusters,
                                              M = num_inner_clusters,
                                              num_L0= n_L0,
                                              t =zt,
                                              it= it )
        
        
        Zt_density_df_control = get_Zt_density_df(x = rep(fixed.regime.control[zt], J),
                                                  MCdata =L0_data,
                                                  BNPModels_df= BNPModels,
                                                  num_MC_samples = J,
                                                  N = num_outer_clusters,
                                                  M = num_inner_clusters,
                                                  num_L0= n_L0,
                                                  t = zt,
                                                  it = it)
        
        density_prod_df_trt = density_prod_df * Zt_density_df_trt
        density_prod_df_control = density_prod_df * Zt_density_df_control
        Y_design_mat = cbind(Y_design_mat, rep(fixed.regime.trt[zt], J))
        
        # Track Z history for building D and M covariate sets
        Z_history_trt = cbind(Z_history_trt, rep(fixed.regime.trt[zt], J))
        Z_history_control = cbind(Z_history_control, rep(fixed.regime.control[zt], J))
        zt = zt + 1
      }
      else if(var.type[j] == "Dt") {
        # Build D covariate data: L0 + Z history up to time dt
        D_MCdata_trt = cbind(L0_data, Z_history_trt[, 1:dt, drop=FALSE])
        D_MCdata_control = cbind(L0_data, Z_history_control[, 1:dt, drop=FALSE])
        
        Dt_sampled_vec_iter_trt = sample_Dt(MCdata = D_MCdata_trt,
                                            D1_vec = D1_trt,
                                            MCdensity = density_prod_df_trt,
                                            BNPModels_df= BNPModels,
                                            num_MC_samples = J,
                                            N=num_outer_clusters,
                                            M=num_inner_clusters,
                                            num_D_covs = K_D_vec[dt],
                                            t = dt,
                                            it =it)
        
        
        Dt_sampled_vec_iter_control = sample_Dt(MCdata = D_MCdata_control,
                                                D1_vec = D1_control,
                                                MCdensity= density_prod_df_control,
                                                BNPModels_df= BNPModels,
                                                num_MC_samples = J,
                                                N=num_outer_clusters,
                                                M=num_inner_clusters,
                                                num_D_covs = K_D_vec[dt],
                                                t = dt,
                                                it =it)
        
        
        density_prod_df_trt = density_prod_df_trt * Dt_sampled_vec_iter_trt$density_df
        density_prod_df_control = density_prod_df_control * Dt_sampled_vec_iter_control$density_df
        if(dt == 1){
          Y_design_mat = cbind(Y_design_mat, D1_trt)
          # Track D history
          D_history_trt = cbind(D_history_trt, D1_trt)
          D_history_control = cbind(D_history_control, D1_control)
        }else{
          Y_design_mat = cbind(Y_design_mat, Dt_sampled_vec_iter_trt$samples)
          # Track D history
          D_history_trt = cbind(D_history_trt, Dt_sampled_vec_iter_trt$samples)
          D_history_control = cbind(D_history_control, Dt_sampled_vec_iter_control$samples)
        }
        dt = dt +1
      }
      else if(var.type[j] == "M1t") {
        # Build M covariate data: L0 + interleaved (Z1,D1,Z2,D2,...) up to time mt
        M_history_control = NULL
        for(tt in 1:mt){
          M_history_control = cbind(M_history_control,
                                    Z_history_control[, tt, drop=FALSE],
                                    D_history_control[, tt, drop=FALSE])
        }
        M_MCdata_control = cbind(L0_data, M_history_control)
        
        M1t_sampled_vec_iter_control = sample_hurdle_data(MCdata= M_MCdata_control,
                                                          MCdensity = density_prod_df_control,
                                                          BNPModels_df= BNPModels,
                                                          Mediator_num = 1,
                                                          num_MC_samples = J,
                                                          N=num_outer_clusters,
                                                          M=num_inner_clusters,
                                                          num_M_covs = K_M_vec[mt],
                                                          t = mt,
                                                          it =it)
        
        density_prod_df_trt = density_prod_df_trt * M1t_sampled_vec_iter_control$density_df    # Mt density arises from the control regime
        density_prod_df_control = density_prod_df_control * M1t_sampled_vec_iter_control$density_df
        Y_design_mat = cbind(Y_design_mat, M1t_sampled_vec_iter_control$samples)
        
        
      }
      else if(var.type[j] == "M2t") {
        # Build M covariate data: L0 + interleaved (Z1,D1,Z2,D2,...) up to time mt
        M_history_control_M2 = NULL
        for(tt in 1:mt){
          M_history_control_M2 = cbind(M_history_control_M2,
                                       Z_history_control[, tt, drop=FALSE],
                                       D_history_control[, tt, drop=FALSE])
        }
        M_MCdata_control_M2 = cbind(L0_data, M_history_control_M2)
        
        M2t_sampled_vec_iter_control=sample_hurdle_data(MCdata= M_MCdata_control_M2,
                                                        MCdensity = density_prod_df_control,
                                                        BNPModels_df= BNPModels,
                                                        Mediator_num = 2,
                                                        num_MC_samples = J,
                                                        N=num_outer_clusters,
                                                        M=num_inner_clusters,
                                                        num_M_covs = K_M_vec[mt],
                                                        t = mt,
                                                        it =it)
        
        density_prod_df_trt = density_prod_df_trt * M2t_sampled_vec_iter_control$density_df    # Mt density arises from the control regime
        density_prod_df_control = density_prod_df_control * M2t_sampled_vec_iter_control$density_df
        Y_design_mat = cbind(Y_design_mat, M2t_sampled_vec_iter_control$samples)
        mt = mt + 1
      }
      else if(var.type[j] == "Y") {
        Y_mean_vec = compute_Y_mean(MCdata= Y_design_mat, 
                                    MCdensity= density_prod_df_trt, 
                                    BNPModels_df= BNPModels,
                                    num_MC_samples= nrow(Y_design_mat),
                                    N =num_outer_clusters,
                                    M=num_inner_clusters,
                                    num_covariates= ncol(Y_design_mat),
                                    it = it )
        
        Y_zzstar_mat[1,it] = mean(Y_mean_vec)
        MC_error_vec[it] = sd(Y_mean_vec) / sqrt(length(Y_mean_vec))
        if(it == 500){
          mc_error = sd(Y_mean_vec) / sqrt(length(Y_mean_vec))
          cat("The MC Error at iteration ", it, "is: ", mc_error, ".\n")
          cat("MC estimate at iteration ", it, "is: ", Y_zzstar_mat[1,it], ".\n")
        }
      }
    }
  }
  MC_error_to_post_SD = MC_error_vec/ sd(as.numeric(Y_zzstar_mat[1,]))
  return(Y_zzstar_mat)
}


compute_PIDE = function(theta_10, theta_00){
  return(theta_10 - theta_00)
}

compute_PJIIE = function(theta_11, theta_10){
  return(theta_11 - theta_10)
}


compute_PCE = function(PIDE, PJIIE){
  return(PIDE+PJIIE)
}

compute_total_effects = function(theta_11, theta_00){
  return(theta_11 - theta_00)
}


calculate_95bayesian_credible_interval  = function(data_row){
  posterior_samples = as.numeric(data_row)
  posterior_mean =  as.numeric(format(round(mean(posterior_samples) , 4), nsmall = 4))
  #posterior_median =  as.numeric(format(round(median(posterior_samples) , 4), nsmall = 4))
  quantiles = quantile(posterior_samples, c(0.025, 0.975))
  interval = as.numeric(format(round(quantiles, 4), nsmall = 4))
  return_df = data.frame(posterior_mean = posterior_mean, 
                         lower_95CI = interval[1],   upper_95CI = interval[2])
  return(return_df)
}


plot_ConfInt <-function(conf_int_mat) {
  conf_int_df = data.frame(conf_int_mat)
  conf_int_df$CI <- 1:nrow(conf_int_df)
  ci_long <- gather(conf_int_df, key = "Bound", value = "Value", -CI)
  ggplot(ci_long, aes(x = CI, ymin = Value, ymax = Value, color = Bound )) +
    geom_errorbar(width = 0.2)+
    labs(x = "Visits", y = "Confidence/Credible Intervals")+
    scale_color_manual(values = c("lower" = "red", "upper" = "blue")) +
    theme_minimal()
  
}


#######################################################################################
###############Compute causal effects for all compliance classes#######################
#######################################################################################


input_varType1 = c("L0","L0","L0","Fi", "Dt", "M1t","M2t", "Fi", "Dt", "M1t","M2t","Fi", "Dt", "M1t","M2t", "Y")


#########################################ACTIVE CUSTOMERS#################################################


active_customer_theta_zzstar = MCInteg(var.type = input_varType1,
                                    BNPModels = post_samples_model_fit,
                                    fixed.regime.trt = rep(1,3),
                                    fixed.regime.control = rep(0,3),
                                    U1 = "a",            
                                    J=10000, 
                                    Ndraws= 1000, 
                                    num_outer_clusters =10,
                                    num_inner_clusters=4,
                                    num_Bin_L0 = 2, 
                                    num_Cont_L0 = 0)


active_customer_theta_zz = MCInteg(var.type = input_varType1,
                                BNPModels = post_samples_model_fit,
                                fixed.regime.trt = rep(1,3),
                                fixed.regime.control = rep(1,3),
                                U1 = "a",            
                                J=10000, 
                                Ndraws= 1000, 
                                num_outer_clusters =10,
                                num_inner_clusters=4,
                                num_Bin_L0 = 2, 
                                num_Cont_L0 = 0)


active_customer_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                        BNPModels = post_samples_model_fit,
                                        fixed.regime.trt = rep(0,3),
                                        fixed.regime.control = rep(0,3),
                                        U1 = "a",            
                                        J=10000, 
                                        Ndraws= 1000, 
                                        num_outer_clusters =10,
                                        num_inner_clusters=4,
                                        num_Bin_L0 = 2, 
                                        num_Cont_L0 = 0)


active_customer_PIDE = compute_PIDE(active_customer_theta_zzstar, active_customer_theta_zstarzstar)
active_customer_PIDE_CI = calculate_95bayesian_credible_interval(active_customer_PIDE) 

active_customer_PJIIE = compute_PJIIE(active_customer_theta_zz, active_customer_theta_zzstar)
active_customer_PJIIE_CI = calculate_95bayesian_credible_interval(active_customer_PJIIE) 

active_customer_PCE = compute_total_effects(active_customer_theta_zz,active_customer_theta_zstarzstar)
active_customer_PCE_CI = calculate_95bayesian_credible_interval(active_customer_PCE) 


#########################################Value Attentive#################################################


value_attentive_theta_zzstar = MCInteg(var.type = input_varType1,
                                BNPModels = post_samples_model_fit,
                                fixed.regime.trt = rep(1,3),
                                fixed.regime.control = rep(0,3),
                                U1 = "c",            
                                J=10000, 
                                Ndraws= 1000, 
                                num_outer_clusters =10,
                                num_inner_clusters=4,
                                num_Bin_L0 = 2, 
                                num_Cont_L0 = 0)


value_attentive_theta_zz = MCInteg(var.type = input_varType1,
                            BNPModels = post_samples_model_fit,
                            fixed.regime.trt = rep(1,3),
                            fixed.regime.control = rep(1,3),
                            U1 = "c",            
                            J=10000, 
                            Ndraws= 1000, 
                            num_outer_clusters =10,
                            num_inner_clusters=4,
                            num_Bin_L0 = 2, 
                            num_Cont_L0 = 0)


value_attentive_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                    BNPModels = post_samples_model_fit,
                                    fixed.regime.trt = rep(0,3),
                                    fixed.regime.control = rep(0,3),
                                    U1 = "c",            
                                    J=10000, 
                                    Ndraws= 1000, 
                                    num_outer_clusters =10,
                                    num_inner_clusters=4,
                                    num_Bin_L0 = 2, 
                                    num_Cont_L0 = 0)


value_attentive_PIDE = compute_PIDE(value_attentive_theta_zzstar, value_attentive_theta_zstarzstar)
value_attentive_PIDE_CI = calculate_95bayesian_credible_interval(value_attentive_PIDE) 

value_attentive_PJIIE = compute_PJIIE(value_attentive_theta_zz, value_attentive_theta_zzstar)
value_attentive_PJIIE_CI = calculate_95bayesian_credible_interval(value_attentive_PJIIE) 

value_attentive_PCE = compute_total_effects(value_attentive_theta_zz,value_attentive_theta_zstarzstar)
value_attentive_PCE_CI = calculate_95bayesian_credible_interval(value_attentive_PCE) 


#########################################Price Attentive#################################################


price_attentive_theta_zzstar = MCInteg(var.type = input_varType1,
                              BNPModels = post_samples_model_fit,
                              fixed.regime.trt = rep(1,3),
                              fixed.regime.control = rep(0,3),
                              U1 = "d",            
                              J=10000, 
                              Ndraws= 1000, 
                              num_outer_clusters =10,
                              num_inner_clusters=4,
                              num_Bin_L0 = 2, 
                              num_Cont_L0 = 0)


price_attentive_theta_zz = MCInteg(var.type = input_varType1,
                          BNPModels = post_samples_model_fit,
                          fixed.regime.trt = rep(1,3),
                          fixed.regime.control = rep(1,3),
                          U1 = "d",            
                          J=10000, 
                          Ndraws= 1000, 
                          num_outer_clusters =10,
                          num_inner_clusters=4,
                          num_Bin_L0 = 2, 
                          num_Cont_L0 = 0)


price_attentive_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                  BNPModels = post_samples_model_fit,
                                  fixed.regime.trt = rep(0,3),
                                  fixed.regime.control = rep(0,3),
                                  U1 = "d",            
                                  J=10000, 
                                  Ndraws= 1000, 
                                  num_outer_clusters =10,
                                  num_inner_clusters=4,
                                  num_Bin_L0 = 2, 
                                  num_Cont_L0 = 0)


price_attentive_PIDE = compute_PIDE(price_attentive_theta_zzstar, price_attentive_theta_zstarzstar)
price_attentive_PIDE_CI = calculate_95bayesian_credible_interval(price_attentive_PIDE) 

price_attentive_PJIIE = compute_PJIIE(price_attentive_theta_zz, price_attentive_theta_zzstar)
price_attentive_PJIIE_CI = calculate_95bayesian_credible_interval(price_attentive_PJIIE) 

price_attentive_PCE = compute_total_effects(price_attentive_theta_zz,price_attentive_theta_zstarzstar)
price_attentive_PCE_CI = calculate_95bayesian_credible_interval(price_attentive_PCE) 


#########################################Non-Active Customers#################################################


non_active_customer_theta_zzstar = MCInteg(var.type = input_varType1,
                                   BNPModels = post_samples_model_fit,
                                   fixed.regime.trt = rep(1,3),
                                   fixed.regime.control = rep(0,3),
                                   U1 = "n",            
                                   J=10000, 
                                   Ndraws= 1000, 
                                   num_outer_clusters =10,
                                   num_inner_clusters=4,
                                   num_Bin_L0 = 2, 
                                   num_Cont_L0 = 0)


non_active_customer_theta_zz = MCInteg(var.type = input_varType1,
                               BNPModels = post_samples_model_fit,
                               fixed.regime.trt = rep(1,3),
                               fixed.regime.control = rep(1,3),
                               U1 = "n",            
                               J=10000, 
                               Ndraws= 1000, 
                               num_outer_clusters =10,
                               num_inner_clusters=4,
                               num_Bin_L0 = 2, 
                               num_Cont_L0 = 0)


non_active_customer_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                       BNPModels = post_samples_model_fit,
                                       fixed.regime.trt = rep(0,3),
                                       fixed.regime.control = rep(0,3),
                                       U1 = "n",            
                                       J=10000, 
                                       Ndraws= 1000, 
                                       num_outer_clusters =10,
                                       num_inner_clusters=4,
                                       num_Bin_L0 = 2, 
                                       num_Cont_L0 = 0)


non_active_customer_PIDE = compute_PIDE(non_active_customer_theta_zzstar, non_active_customer_theta_zstarzstar)
non_active_customer_PIDE_CI = calculate_95bayesian_credible_interval(non_active_customer_PIDE) 

non_active_customer_PJIIE = compute_PJIIE(non_active_customer_theta_zz, non_active_customer_theta_zzstar)
non_active_customer_PJIIE_CI = calculate_95bayesian_credible_interval(non_active_customer_PJIIE) 

non_active_customer_PCE = compute_total_effects(non_active_customer_theta_zz,non_active_customer_theta_zstarzstar)
non_active_customer_PCE_CI = calculate_95bayesian_credible_interval(non_active_customer_PCE) 


#######################################################################################
#####################Save results######################################################
#######################################################################################


#-------------------------------------------------------------------------------------------------------------------#
#--------------------------------WRITE RESULTS ON A DF--------------------------------------------------------------#

allinfo <- data.frame(
  runID = run_ID,
  
  PIDE_AC_post_mean  = as.numeric(active_customer_PIDE_CI[1]),
  PIDE_AC_lowerCI    = as.numeric(active_customer_PIDE_CI[2]),
  PIDE_AC_upperCI    = as.numeric(active_customer_PIDE_CI[3]),
  
  PJIIE_AC_post_mean = as.numeric(active_customer_PJIIE_CI[1]),
  PJIIE_AC_lowerCI   = as.numeric(active_customer_PJIIE_CI[2]),
  PJIIE_AC_upperCI   = as.numeric(active_customer_PJIIE_CI[3]),
  
  PCE_AC_post_mean   = as.numeric(active_customer_PCE_CI[1]),
  PCE_AC_lowerCI     = as.numeric(active_customer_PCE_CI[2]),
  PCE_AC_upperCI     = as.numeric(active_customer_PCE_CI[3]),
  
  PIDE_VA_post_mean   = as.numeric(value_attentive_PIDE_CI[1]),
  PIDE_VA_lowerCI     = as.numeric(value_attentive_PIDE_CI[2]),
  PIDE_VA_upperCI     = as.numeric(value_attentive_PIDE_CI[3]),
  
  PJIIE_VA_post_mean  = as.numeric(value_attentive_PJIIE_CI[1]),
  PJIIE_VA_lowerCI    = as.numeric(value_attentive_PJIIE_CI[2]),
  PJIIE_VA_upperCI    = as.numeric(value_attentive_PJIIE_CI[3]),
  
  PCE_VA_post_mean    = as.numeric(value_attentive_PCE_CI[1]),
  PCE_VA_lowerCI      = as.numeric(value_attentive_PCE_CI[2]),
  PCE_VA_upperCI      = as.numeric(value_attentive_PCE_CI[3]),
  
  PIDE_PA_post_mean   = as.numeric(price_attentive_PIDE_CI[1]),
  PIDE_PA_lowerCI     = as.numeric(price_attentive_PIDE_CI[2]),
  PIDE_PA_upperCI     = as.numeric(price_attentive_PIDE_CI[3]),
  
  PJIIE_PA_post_mean  = as.numeric(price_attentive_PJIIE_CI[1]),
  PJIIE_PA_lowerCI    = as.numeric(price_attentive_PJIIE_CI[2]),
  PJIIE_PA_upperCI    = as.numeric(price_attentive_PJIIE_CI[3]),
  
  PCE_PA_post_mean    = as.numeric(price_attentive_PCE_CI[1]),
  PCE_PA_lowerCI      = as.numeric(price_attentive_PCE_CI[2]),
  PCE_PA_upperCI      = as.numeric(price_attentive_PCE_CI[3]),
  
  PIDE_NAC_post_mean  = as.numeric(non_active_customer_PIDE_CI[1]),
  PIDE_NAC_lowerCI    = as.numeric(non_active_customer_PIDE_CI[2]),
  PIDE_NAC_upperCI    = as.numeric(non_active_customer_PIDE_CI[3]),
  
  PJIIE_NAC_post_mean = as.numeric(non_active_customer_PJIIE_CI[1]),
  PJIIE_NAC_lowerCI   = as.numeric(non_active_customer_PJIIE_CI[2]),
  PJIIE_NAC_upperCI   = as.numeric(non_active_customer_PJIIE_CI[3]),
  
  PCE_NAC_post_mean   = as.numeric(non_active_customer_PCE_CI[1]),
  PCE_NAC_lowerCI     = as.numeric(non_active_customer_PCE_CI[2]),
  PCE_NAC_upperCI     = as.numeric(non_active_customer_PCE_CI[3]),
  
  num_outer_clusters_used = as.numeric(mean_num_outer_used)
)


write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
gc()

