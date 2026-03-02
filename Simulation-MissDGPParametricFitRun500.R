#-------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------Prepare for simulation----------------------------------------------#

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


txt.title = paste0("param_simu_results/Results_Nimble_Param_Models_result.txt")
if (run_ID == 1) {
  df = data.frame(matrix(ncol = 37, nrow = 0))                 #each row contains results from each replication/dataset/run
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
                   "PCE_NAC_post_mean", "PCE_NAC_lowerCI", "PCE_NAC_upperCI")
  colnames(df) = df_col_names
  write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
}


#-------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------Fit Models in the observed dataset-------------------------------------------------#


set.seed(1234)   #this will be changed later during data generation

# Set library path for compute node packages
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# Debug: Print library paths and R version (remove after debugging)
cat("R version:", R.version.string, "\n")
cat("Library paths:\n")
print(.libPaths())
cat("Nimble installed:", "nimble" %in% rownames(installed.packages()), "\n")

# Load packages (install on compute node first using install_packages.sbatch)
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


# scale the mediators


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
digCom_workdf$GENDER <- ifelse(digCom_workdf$GENDER == "M", 1, 0)


# change BASELINE_AGE elements from >=4 to 1,and <4 to 0.
digCom_workdf$BASELINE_AGE <- ifelse(digCom_workdf$BASELINE_AGE >= 4, 1, 0)


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


num_Bin_L0 = 2  # both gender and baseline age are binary

# sample frequency of each Binary column of L0_df_wide
pr_theta_Bin_L0 = as.numeric(colMeans(L0_df_wide[2:3], na.rm = TRUE))

# sample means of each Continuous column of L0_df_wide

# sample variances of each Continuous column of L0_df_wide


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


code_DG <- nimbleCode({
  logDens ~ dnorm(0, 1)    ## this distribution does not matter: only need to compute model sum-posterior-log-density
  
  ##################################################################################
  ##############################Prior distributions#################################
  ##################################################################################
  
  
  ##############################Mixture component probabilities######################
  pi_mix[1:K] ~ ddirch(alpha[1:K])
  
  
  ##############################Hurdle probabilities################################
  for(c in 1:K){
    pi0_Y[c] ~ dbeta(1,1)
    pi0_M1_t3[c] ~ dbeta(1,1)
    pi0_M2_t3[c] ~ dbeta(1,1)
    pi0_M1_t2[c] ~ dbeta(1,1)
    pi0_M2_t2[c] ~ dbeta(1,1)
    pi0_M1_t1[c] ~ dbeta(1,1)
    pi0_M2_t1[c] ~ dbeta(1,1)
  }
  
  
  #####################Regression fixed effects coefficients################################
  for(c in 1:K){
    # AT t =3
    for(k in 1:num_Y_reg_coeff){beta_Yi[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_M2t3_reg_coeff){theta_M2_t3[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_M1t3_reg_coeff){theta_M1_t3[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_Dt3_reg_coeff){theta_D_t3[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_Zt3_reg_coeff){theta_Z_t3[c,k] ~ dnorm(mean = 0, var = 10)}
    
    # AT t =2
    for(k in 1:num_M2t2_reg_coeff){theta_M2_t2[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_M1t2_reg_coeff){theta_M1_t2[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_Dt2_reg_coeff){theta_D_t2[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_Zt2_reg_coeff){theta_Z_t2[c,k] ~ dnorm(mean = 0, var = 10)}
    
    # AT t =1
    for(k in 1:num_M2t1_reg_coeff){theta_M2_t1[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_M1t1_reg_coeff){theta_M1_t1[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_Dt1_reg_coeff){theta_D_t1[c,k] ~ dnorm(mean = 0, var = 10)}
    for(k in 1:num_Zt1_reg_coeff){theta_Z_t1[c,k] ~ dnorm(mean = 0, var = 10)}
  }
  
  
  #####################Regression variance coefficients################################
  for(c in 1:K){
    sigma_sq_Yi[c] ~ dinvgamma(shape = 3, scale = 1)
    sigma_sq_M2_t3[c] ~ dinvgamma(shape = 3, scale = 1)
    sigma_sq_M1_t3[c] ~ dinvgamma(shape = 3, scale = 1)
    sigma_sq_M2_t2[c] ~ dinvgamma(shape = 3, scale = 1)
    sigma_sq_M1_t2[c] ~ dinvgamma(shape = 3, scale = 1)
    sigma_sq_M2_t1[c] ~ dinvgamma(shape = 3, scale = 1)
    sigma_sq_M1_t1[c] ~ dinvgamma(shape = 3, scale = 1)
  }
  
  
  ##################################################################################
  ##############################Likelihood##########################################
  ##################################################################################
  
  
  for(i in 1:n) {
    
    
    # Mixture component assignment
    s[i] ~ dcat(pi_mix[1:K])
    
    
    #AT t = 1
    
    Zt1[i] ~ dbern(pnorm(inprod(theta_Z_t1[s[i], 1:num_Zt1_reg_coeff], Zt1_param_reg_design_mat[i,1:num_Zt1_reg_coeff]),0,1))
    Dt1[i] ~ dbern(pnorm(inprod(theta_D_t1[s[i], 1:num_Dt1_reg_coeff], Dt1_param_reg_design_mat[i,1:num_Dt1_reg_coeff]),0,1))
    
    Z_latent_M1t1[i] ~ dbern(1 - pi0_M1_t1[s[i]])  # Indicator for positive values
    m1t1_mean[i] <-   inprod(theta_M1_t1[s[i], 1:num_M1t1_reg_coeff],
                             M1t1_param_reg_design_mat[i, 1:num_M1t1_reg_coeff])
    data_wide_M1t1_pos[i] ~ dnorm(m1t1_mean[i], sd = sqrt(sigma_sq_M1_t1[s[i]]))
    M1t1[i] ~ dnorm(Z_latent_M1t1[i]*data_wide_M1t1_pos[i], sd = 0.25)
    
    Z_latent_M2t1[i] ~ dbern(1 - pi0_M2_t1[s[i]])  # Indicator for positive values
    m2t1_mean[i] <-   inprod(theta_M2_t1[s[i], 1:num_M2t1_reg_coeff],
                             M2t1_param_reg_design_mat[i, 1:num_M2t1_reg_coeff])
    data_wide_M2t1_pos[i] ~ dnorm(m2t1_mean[i], sd = sqrt(sigma_sq_M2_t1[s[i]]))
    M2t1[i] ~ dnorm(Z_latent_M2t1[i]*data_wide_M2t1_pos[i], sd = 0.25)
    
    
    #AT t = 2
    
    Zt2[i] ~ dbern(pnorm(inprod(theta_Z_t2[s[i], 1:num_Zt2_reg_coeff], Zt2_param_reg_design_mat[i,1:num_Zt2_reg_coeff]),0,1))
    Dt2[i] ~ dbern(pnorm(inprod(theta_D_t2[s[i], 1:num_Dt2_reg_coeff], Dt2_param_reg_design_mat[i,1:num_Dt2_reg_coeff]),0,1))
    
    Z_latent_M1t2[i] ~ dbern(1 - pi0_M1_t2[s[i]])  # Indicator for positive values
    m1t2_mean[i] <-   inprod(theta_M1_t2[s[i], 1:num_M1t2_reg_coeff],
                             M1t2_param_reg_design_mat[i, 1:num_M1t2_reg_coeff])
    data_wide_M1t2_pos[i] ~ dnorm(m1t2_mean[i], sd = sqrt(sigma_sq_M1_t2[s[i]]))
    M1t2[i] ~ dnorm(Z_latent_M1t2[i]*data_wide_M1t2_pos[i], sd = 0.25)
    
    Z_latent_M2t2[i] ~ dbern(1 - pi0_M2_t2[s[i]])  # Indicator for positive values
    m2t2_mean[i] <-   inprod(theta_M2_t2[s[i], 1:num_M2t2_reg_coeff],
                             M2t2_param_reg_design_mat[i, 1:num_M2t2_reg_coeff])
    data_wide_M2t2_pos[i] ~ dnorm(m2t2_mean[i], sd = sqrt(sigma_sq_M2_t2[s[i]]))
    M2t2[i] ~ dnorm(Z_latent_M2t2[i]*data_wide_M2t3_pos[i], sd = 0.25)
    
    
    #AT t = 3
    
    Zt3[i] ~ dbern(pnorm(inprod(theta_Z_t3[s[i], 1:num_Zt3_reg_coeff], Zt3_param_reg_design_mat[i,1:num_Zt3_reg_coeff]),0,1))
    Dt3[i] ~ dbern(pnorm(inprod(theta_D_t3[s[i], 1:num_Dt3_reg_coeff], Dt3_param_reg_design_mat[i,1:num_Dt3_reg_coeff]),0,1))
    
    Z_latent_M1t3[i] ~ dbern(1 - pi0_M1_t3[s[i]])  # Indicator for positive values
    m1t3_mean[i] <-   inprod(theta_M1_t3[s[i], 1:num_M1t3_reg_coeff],
                             M1t3_param_reg_design_mat[i, 1:num_M1t3_reg_coeff])
    data_wide_M1t3_pos[i] ~ dnorm(m1t3_mean[i], sd = sqrt(sigma_sq_M1_t3[s[i]]))
    M1t3[i] ~ dnorm(Z_latent_M1t3[i]*data_wide_M1t3_pos[i], sd = 0.25)
    
    Z_latent_M2t3[i] ~ dbern(1 - pi0_M2_t3[s[i]])  # Indicator for positive values
    m2t3_mean[i] <-   inprod(theta_M2_t3[s[i], 1:num_M2t3_reg_coeff],
                             M2t3_param_reg_design_mat[i, 1:num_M2t3_reg_coeff])
    data_wide_M2t3_pos[i] ~ dnorm(m2t3_mean[i], sd = sqrt(sigma_sq_M2_t3[s[i]]))
    M2t3[i] ~ dnorm(Z_latent_M2t3[i]*data_wide_M2t3_pos[i], sd = 0.25)
    
    Z_latent_Y[i] ~ dbern(1 - pi0_Y[s[i]])  # Indicator for positive values
    y_mean[i] <- inprod(beta_Yi[s[i], 1:num_Y_reg_coeff], Y_param_reg_design_mat[i, 1:num_Y_reg_coeff])
    Y_pos[i] ~ dnorm(y_mean[i], sd =  sqrt(sigma_sq_Yi[s[i]]))
    Y[i] ~ dnorm(Z_latent_Y[i]*Y_pos[i], sd = 0.25)
    
    
  }
})


# one outer cluster and one inner cluster for fitting parametric data
L0_mat = as.matrix(L0_df_wide)
num_L0 =3   # one intercept and 2 L0

constants_DG <- list(n = length(unique(digCom_workdf_wide$ID)), # Number of individuals
                     K = 3,                                      # Number of mixture components
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
  # Dirichlet concentration parameter
  alpha = rep(1, 3),
  
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
  # mixture parameters
  pi_mix = rep(1/3, 3),
  s = sample(1:3, constants_DG$n, replace = TRUE),
  
  # outcome model parameters (K x num_coeff matrices)
  beta_Yi = matrix(0, nrow = 3, ncol = constants_DG$num_Y_reg_coeff),
  sigma_sq_Yi = rinvgamma(3, 1, 1),
  pi0_Y = rep(0.5, 3),
  Y_pos = rep(mean(digCom_workdf_wide$Y_t_3), constants_DG$n),
  
  # mediator parameters at t =3
  theta_M2_t3 = matrix(0, nrow = 3, ncol = constants_DG$num_M2t3_reg_coeff),
  theta_M1_t3 = matrix(0, nrow = 3, ncol = constants_DG$num_M1t3_reg_coeff),
  sigma_sq_M2_t3 = rinvgamma(3, 1, 1),
  sigma_sq_M1_t3 = rinvgamma(3, 1, 1),
  pi0_M2_t3 = rep(0.5, 3),
  pi0_M1_t3 = rep(0.5, 3),
  data_wide_M2t3_pos = rep(mean(digCom_workdf_wide$M2_t_3), constants_DG$n),
  data_wide_M1t3_pos = rep(mean(digCom_workdf_wide$M1_t_3), constants_DG$n),
  
  # mediator parameters at t =2
  theta_M2_t2 = matrix(0, nrow = 3, ncol = constants_DG$num_M2t2_reg_coeff),
  theta_M1_t2 = matrix(0, nrow = 3, ncol = constants_DG$num_M1t2_reg_coeff),
  sigma_sq_M2_t2 = rinvgamma(3, 1, 1),
  sigma_sq_M1_t2 = rinvgamma(3, 1, 1),
  pi0_M2_t2 = rep(0.5, 3),
  pi0_M1_t2 = rep(0.5, 3),
  data_wide_M2t2_pos = rep(mean(digCom_workdf_wide$M2_t_2), constants_DG$n),
  data_wide_M1t2_pos = rep(mean(digCom_workdf_wide$M1_t_2), constants_DG$n),
  
  # mediator parameters at t =1
  theta_M2_t1 = matrix(0, nrow = 3, ncol = constants_DG$num_M2t1_reg_coeff),
  theta_M1_t1 = matrix(0, nrow = 3, ncol = constants_DG$num_M1t1_reg_coeff),
  sigma_sq_M2_t1 = rinvgamma(3, 1, 1),
  sigma_sq_M1_t1 = rinvgamma(3, 1, 1),
  pi0_M2_t1 = rep(0.5, 3),
  pi0_M1_t1 = rep(0.5, 3),
  data_wide_M2t1_pos = rep(mean(digCom_workdf_wide$M2_t_1), constants_DG$n),
  data_wide_M1t1_pos = rep(mean(digCom_workdf_wide$M1_t_1), constants_DG$n),
  
  # treatment receipt parameters (K x num_coeff matrices)
  theta_D_t3 = matrix(0, nrow = 3, ncol = constants_DG$num_Dt3_reg_coeff),
  theta_D_t2 = matrix(0, nrow = 3, ncol = constants_DG$num_Dt2_reg_coeff),
  theta_D_t1 = matrix(0, nrow = 3, ncol = constants_DG$num_Dt1_reg_coeff),
  
  # treatment assignment parameters (K x num_coeff matrices)
  theta_Z_t3 = matrix(0, nrow = 3, ncol = constants_DG$num_Zt3_reg_coeff),
  theta_Z_t2 = matrix(0, nrow = 3, ncol = constants_DG$num_Zt2_reg_coeff),
  theta_Z_t1 = matrix(0, nrow = 3, ncol = constants_DG$num_Zt1_reg_coeff)
)


Param_data_generating_model <- nimbleModel(code_DG, constants_DG, data_DG, inits_DG)  # model creation


compile_Param_data_generating_model <- compileNimble(Param_data_generating_model, showCompilerOutput = TRUE)  # model compilation

# MCMC configuration
config_MCMC_DG <- configureMCMC(Param_data_generating_model, useConjugacy = TRUE,
                                enableWAIC = TRUE,
                                monitors = c("pi_mix", "s",
                                             "beta_Yi", "sigma_sq_Yi", "pi0_Y",
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


# Override specific coefficients for all K=3 mixture components
# Column names are now in format "param[c, k]" where c=component, k=covariate index


set.seed(run_ID)   #new seed for data generation
replicated_data_sample_size = 12000
number_of_mediators = 2
zero_vec = rep(0,2)

####################Generate L0 from the empirical distribution########################

intercept_vec =  rep(1, replicated_data_sample_size)
gender_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = mean(L0_mat[,2]))
baseline_age_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = mean(L0_mat[,3]))

replicated_df = data.frame(Intercept = intercept_vec, GENDER = gender_vec, BASELINE_AGE = baseline_age_vec)


####################Generate mixture component assignments s[i]######################
# Extract posterior mean of pi_mix for each component
K_mix = 3
pi_mix_cols = paste0("pi_mix[", 1:K_mix, "]")
pi_mix_vec = as.numeric(param_samples_data_generating_model[, pi_mix_cols])
# Assign each individual to a mixture component
s_vec = sample(1:K_mix, size = replicated_data_sample_size, replace = TRUE, prob = pi_mix_vec)
replicated_df$s = s_vec


####################Helper function to extract component-specific coefficients########
# For a parameter like "theta_Z_t1[c, k]", extract coefficients for component c
get_component_coeffs = function(param_name, component, num_coeffs, param_df) {
  col_names = paste0(param_name, "[", component, ", ", 1:num_coeffs, "]")
  return(as.numeric(param_df[, col_names]))
}


####################################Generate Z_t1#####################################
num_Zt1_coeffs = constants_DG$num_Zt1_reg_coeff
Z_t1_prob = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_Z_t1_cc = get_component_coeffs("theta_Z_t1", cc, num_Zt1_coeffs, param_samples_data_generating_model)
  
  Z_t1_prob[idx] = pnorm(as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE")]) %*% theta_Z_t1_cc)
}
Z_t1_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = Z_t1_prob)
replicated_df$Z_t1 = Z_t1_vec


####################################Generate D_t1#####################################
num_Dt1_coeffs = constants_DG$num_Dt1_reg_coeff
D_t1_prob = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_D_t1_cc = get_component_coeffs("theta_D_t1", cc, num_Dt1_coeffs, param_samples_data_generating_model)
  D_t1_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","Z_t1")])
  D_t1_prob[idx] = pnorm(D_t1_design %*% theta_D_t1_cc)
}
D_t1_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = D_t1_prob)
replicated_df$D_t1 = D_t1_vec


####################Generate mediator 1 M1t1####################################################
num_M1t1_coeffs = constants_DG$num_M1t1_reg_coeff
M1t1_vec = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  n_idx = length(idx)
  theta_M1_t1_cc = get_component_coeffs("theta_M1_t1", cc, num_M1t1_coeffs, param_samples_data_generating_model)
  pi0_M1_t1_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_M1_t1[", cc, "]")])
  sigma_sq_M1_t1_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_M1_t1[", cc, "]")])
  M1t1_pos_indicator = rbinom(n = n_idx, size = 1, prob = 1 - pi0_M1_t1_cc)
  M1t1_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","Z_t1","D_t1")])
  M1t1_pos_mean = M1t1_design %*% theta_M1_t1_cc
  M1t1_vec[idx] = rtruncnorm(n = n_idx, a=0, b=10, mean = M1t1_pos_mean, sd = sqrt(sigma_sq_M1_t1_cc)) * M1t1_pos_indicator
}
replicated_df$M1_t1 = M1t1_vec


####################Generate mediator 2 M2t1####################################################
num_M2t1_coeffs = constants_DG$num_M2t1_reg_coeff
M2t1_vec = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  n_idx = length(idx)
  theta_M2_t1_cc = get_component_coeffs("theta_M2_t1", cc, num_M2t1_coeffs, param_samples_data_generating_model)
  pi0_M2_t1_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_M2_t1[", cc, "]")])
  sigma_sq_M2_t1_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_M2_t1[", cc, "]")])
  M2t1_pos_indicator = rbinom(n = n_idx, size = 1, prob = 1 - pi0_M2_t1_cc)
  M2t1_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","Z_t1","D_t1","M1_t1")])
  M2t1_pos_mean = M2t1_design %*% theta_M2_t1_cc
  M2t1_vec[idx] = rtruncnorm(n = n_idx, a=0, b=10, mean = M2t1_pos_mean, sd = sqrt(sigma_sq_M2_t1_cc)) * M2t1_pos_indicator
}
replicated_df$M2_t1 = M2t1_vec


####################################Generate Z_t2#####################################
num_Zt2_coeffs = constants_DG$num_Zt2_reg_coeff
Z_t2_prob = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_Z_t2_cc = get_component_coeffs("theta_Z_t2", cc, num_Zt2_coeffs, param_samples_data_generating_model)
  Z_t2_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","Z_t1","D_t1","M1_t1","M2_t1")])
  Z_t2_prob[idx] = pnorm(Z_t2_design %*% theta_Z_t2_cc)
}
Z_t2_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = Z_t2_prob)
replicated_df$Z_t2 = Z_t2_vec


####################################Generate D_t2#####################################
num_Dt2_coeffs = constants_DG$num_Dt2_reg_coeff
D_t2_prob = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_D_t2_cc = get_component_coeffs("theta_D_t2", cc, num_Dt2_coeffs, param_samples_data_generating_model)
  D_t2_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","D_t1","M1_t1","M2_t1","Z_t2")])
  D_t2_prob[idx] = pnorm(D_t2_design %*% theta_D_t2_cc)
}
D_t2_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = D_t2_prob)
replicated_df$D_t2 = D_t2_vec


####################Generate mediator 1 M1t2####################################################
num_M1t2_coeffs = constants_DG$num_M1t2_reg_coeff
M1t2_vec = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  n_idx = length(idx)
  theta_M1_t2_cc = get_component_coeffs("theta_M1_t2", cc, num_M1t2_coeffs, param_samples_data_generating_model)
  pi0_M1_t2_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_M1_t2[", cc, "]")])
  sigma_sq_M1_t2_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_M1_t2[", cc, "]")])
  M1t2_pos_indicator = rbinom(n = n_idx, size = 1, prob = 1 - pi0_M1_t2_cc)
  M1t2_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","M1_t1","M2_t1","Z_t2","D_t2")])
  M1t2_pos_mean = M1t2_design %*% theta_M1_t2_cc
  M1t2_vec[idx] = rtruncnorm(n = n_idx, a=0, b=10, mean = M1t2_pos_mean, sd = sqrt(sigma_sq_M1_t2_cc)) * M1t2_pos_indicator
}
replicated_df$M1_t2 = M1t2_vec


####################Generate mediator 2 M2t2####################################################
num_M2t2_coeffs = constants_DG$num_M2t2_reg_coeff
M2t2_vec = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  n_idx = length(idx)
  theta_M2_t2_cc = get_component_coeffs("theta_M2_t2", cc, num_M2t2_coeffs, param_samples_data_generating_model)
  pi0_M2_t2_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_M2_t2[", cc, "]")])
  sigma_sq_M2_t2_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_M2_t2[", cc, "]")])
  M2t2_pos_indicator = rbinom(n = n_idx, size = 1, prob = 1 - pi0_M2_t2_cc)
  M2t2_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","M2_t1","Z_t2","D_t2","M1_t2")])
  M2t2_pos_mean = M2t2_design %*% theta_M2_t2_cc
  M2t2_vec[idx] = rtruncnorm(n = n_idx, a=0, b=10, mean = M2t2_pos_mean, sd = sqrt(sigma_sq_M2_t2_cc)) * M2t2_pos_indicator
}
replicated_df$M2_t2 = M2t2_vec


####################################Generate Z_t3#####################################
num_Zt3_coeffs = constants_DG$num_Zt3_reg_coeff
Z_t3_prob = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_Z_t3_cc = get_component_coeffs("theta_Z_t3", cc, num_Zt3_coeffs, param_samples_data_generating_model)
  Z_t3_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","Z_t2","D_t2","M1_t2","M2_t2")])
  Z_t3_prob[idx] = pnorm(Z_t3_design %*% theta_Z_t3_cc)
}
Z_t3_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = Z_t3_prob)
replicated_df$Z_t3 = Z_t3_vec
####################################Generate D_t3#####################################
num_Dt3_coeffs = constants_DG$num_Dt3_reg_coeff
D_t3_prob = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_D_t3_cc = get_component_coeffs("theta_D_t3", cc, num_Dt3_coeffs, param_samples_data_generating_model)
  D_t3_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","D_t2","M1_t2","M2_t2","Z_t3")])
  D_t3_prob[idx] = pnorm(D_t3_design %*% theta_D_t3_cc)
}
D_t3_vec = rbinom(n = replicated_data_sample_size, size = 1, prob = D_t3_prob)
replicated_df$D_t3 = D_t3_vec

####################Generate mediator 1 M1t3####################################################
num_M1t3_coeffs = constants_DG$num_M1t3_reg_coeff
M1t3_pos_mean = rep(NA, replicated_data_sample_size)
M1t3_pos_indicator = rep(NA, replicated_data_sample_size)
M1t3_sigma_sq = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_M1_t3_cc = get_component_coeffs("theta_M1_t3", cc, num_M1t3_coeffs, param_samples_data_generating_model)
  pi0_M1_t3_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_M1_t3[", cc, "]")])
  sigma_sq_M1_t3_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_M1_t3[", cc, "]")])
  M1_t3_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","M1_t2","M2_t2","Z_t3","D_t3")])
  M1t3_pos_mean[idx] = M1_t3_design %*% theta_M1_t3_cc
  M1t3_pos_indicator[idx] = rbinom(n = length(idx), size = 1, prob = 1 - pi0_M1_t3_cc)
  M1t3_sigma_sq[idx] = sigma_sq_M1_t3_cc
}
M1t3_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M1t3_pos_mean, sd = sqrt(M1t3_sigma_sq)) * M1t3_pos_indicator
replicated_df$M1_t3 = M1t3_vec

####################Generate mediator 2 M2t3####################################################
num_M2t3_coeffs = constants_DG$num_M2t3_reg_coeff
M2t3_pos_mean = rep(NA, replicated_data_sample_size)
M2t3_pos_indicator = rep(NA, replicated_data_sample_size)
M2t3_sigma_sq = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  theta_M2_t3_cc = get_component_coeffs("theta_M2_t3", cc, num_M2t3_coeffs, param_samples_data_generating_model)
  pi0_M2_t3_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_M2_t3[", cc, "]")])
  sigma_sq_M2_t3_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_M2_t3[", cc, "]")])
  M2_t3_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","M2_t2","Z_t3","D_t3","M1_t3")])
  M2t3_pos_mean[idx] = M2_t3_design %*% theta_M2_t3_cc
  M2t3_pos_indicator[idx] = rbinom(n = length(idx), size = 1, prob = 1 - pi0_M2_t3_cc)
  M2t3_sigma_sq[idx] = sigma_sq_M2_t3_cc
}
M2t3_vec = rtruncnorm(n = replicated_data_sample_size, a=0, b=10,
                      mean = M2t3_pos_mean, sd = sqrt(M2t3_sigma_sq)) * M2t3_pos_indicator
replicated_df$M2_t3 = M2t3_vec

####################Generate outcome Y####################################################
num_Y_coeffs = constants_DG$num_Y_reg_coeff
Y_pos_mean = rep(NA, replicated_data_sample_size)
Y_pos_indicator = rep(NA, replicated_data_sample_size)
Y_sigma_sq = rep(NA, replicated_data_sample_size)
for(cc in 1:K_mix){
  idx = which(s_vec == cc)
  beta_Yi_cc = get_component_coeffs("beta_Yi", cc, num_Y_coeffs, param_samples_data_generating_model)
  pi0_Y_cc = as.numeric(param_samples_data_generating_model[, paste0("pi0_Y[", cc, "]")])
  sigma_sq_Yi_cc = as.numeric(param_samples_data_generating_model[, paste0("sigma_sq_Yi[", cc, "]")])
  Y_design = as.matrix(replicated_df[idx, c("Intercept","GENDER","BASELINE_AGE","Z_t3","D_t3","M1_t3","M2_t3")])
  Y_pos_mean[idx] = Y_design %*% beta_Yi_cc
  Y_pos_indicator[idx] = rbinom(n = length(idx), size = 1, prob = 1 - pi0_Y_cc)
  Y_sigma_sq[idx] = sigma_sq_Yi_cc
}
Y_vector = rtruncnorm(n = replicated_data_sample_size, a=0, b=Inf,
                      mean = Y_pos_mean, sd = sqrt(Y_sigma_sq)) * Y_pos_indicator
replicated_df$Y = Y_vector


L0_df_wide_fit = data.frame(Intercept = replicated_df$Intercept,
                            GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE) 

Y_vec_wide_fit = replicated_df$Y


param_reg_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                               Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1, M1_t_1 = replicated_df$M1_t1, M2_t_1 = replicated_df$M2_t1,
                               Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2, M1_t_2 = replicated_df$M1_t2, M2_t_2 = replicated_df$M2_t2,
                               Z_t_3 = replicated_df$Z_t3, D_t_3 = replicated_df$D_t3, M1_t_3 = replicated_df$M1_t3, M2_t_3 = replicated_df$M2_t3,
                               Y_t_3 = replicated_df$Y)


Y_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                        Z_t_3 = replicated_df$Z_t3, D_t_3 = replicated_df$D_t3, M1_t_3 = replicated_df$M1_t3, M2_t_3 = replicated_df$M2_t3)


M2t3_param_reg_design_mat_fit  = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                            M2_t_2 = replicated_df$M2_t2,
                                            Z_t_3 = replicated_df$Z_t3, D_t_3 = replicated_df$D_t3, M1_t_3 = replicated_df$M1_t3)


M1t3_param_reg_design_mat_fit  = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                            M1_t_2 = replicated_df$M1_t2, M2_t_2 = replicated_df$M2_t2,
                                            Z_t_3 = replicated_df$Z_t3, D_t_3 = replicated_df$D_t3)


Dt3_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                          D_t_2 = replicated_df$D_t2, M1_t_2 = replicated_df$M1_t2, M2_t_2 = replicated_df$M2_t2,
                                          Z_t_3 = replicated_df$Z_t3)


Zt3_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                          Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2, M1_t_2 = replicated_df$M1_t2, M2_t_2 = replicated_df$M2_t2)


M2t2_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                           M2_t_1 = replicated_df$M2_t1,
                                           Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2, M1_t_2 = replicated_df$M1_t2)


M1t2_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                           M1_t_1 = replicated_df$M1_t1, M2_t_1 = replicated_df$M2_t1,
                                           Z_t_2 = replicated_df$Z_t2, D_t_2 = replicated_df$D_t2)


Dt2_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                          D_t_1 = replicated_df$D_t1, M1_t_1 = replicated_df$M1_t1, M2_t_1 = replicated_df$M2_t1,
                                          Z_t_2 = replicated_df$Z_t2)


Zt2_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                          Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1, M1_t_1 = replicated_df$M1_t1, M2_t_1 = replicated_df$M2_t1)


M2t1_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                           Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1, M1_t_1 = replicated_df$M1_t1)


M1t1_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                           Z_t_1 = replicated_df$Z_t1, D_t_1 = replicated_df$D_t1)


Dt1_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE,
                                          Z_t_1 = replicated_df$Z_t1)


Zt1_param_reg_design_mat_fit = data.frame(GENDER = replicated_df$GENDER, BASELINE_AGE = replicated_df$BASELINE_AGE)


num_Bin_L0_fit = 2  # both gender and baseline age are binary

pr_theta_Bin_L0_fit = as.numeric(colMeans(L0_df_wide_fit[2:3], na.rm = TRUE))


L0_mat_fit = as.matrix(L0_df_wide_fit)
num_L0_fit =3   # one intercept and 2 L0


code_fit <- nimbleCode({
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
######################################################################################
#######################################################################################


constants_fit <- list(n = replicated_data_sample_size, # Number of individuals
                      num_Y_reg_coeff = ncol(Y_param_reg_design_mat_fit) + 1, # plus 1 to include the intercept
                      num_M2t3_reg_coeff = ncol(M2t3_param_reg_design_mat_fit) + 1,
                      num_M1t3_reg_coeff = ncol(M1t3_param_reg_design_mat_fit) + 1,
                      num_Dt3_reg_coeff = ncol(Dt3_param_reg_design_mat_fit) + 1, 
                      num_Zt3_reg_coeff = ncol(Zt3_param_reg_design_mat_fit) + 1, 
                      num_M2t2_reg_coeff = ncol(M2t2_param_reg_design_mat_fit) + 1, 
                      num_M1t2_reg_coeff = ncol(M1t2_param_reg_design_mat_fit) + 1, 
                      num_Dt2_reg_coeff =ncol(Dt2_param_reg_design_mat_fit) + 1, 
                      num_Zt2_reg_coeff =ncol(Zt2_param_reg_design_mat_fit) + 1, 
                      num_M2t1_reg_coeff =ncol(M2t1_param_reg_design_mat_fit) + 1,
                      num_M1t1_reg_coeff =ncol(M1t1_param_reg_design_mat_fit) + 1,
                      num_Dt1_reg_coeff =ncol(Dt1_param_reg_design_mat_fit) + 1, 
                      num_Zt1_reg_coeff =ncol(Zt1_param_reg_design_mat_fit) + 1
)


data_fit <- list(
  # regression outcomes
  #### t=3
  Y = param_reg_mat_fit$Y_t_3,
  Z_latent_Y = ifelse(param_reg_mat_fit$Y_t_3 ==0, 0,1),
  M2t3 =  param_reg_mat_fit$M2_t_3,
  Z_latent_M2t3 = ifelse(param_reg_mat_fit$M2_t_3 ==0, 0,1),
  M1t3 = param_reg_mat_fit$M1_t_3,
  Z_latent_M1t3 = ifelse(param_reg_mat_fit$M1_t_3 ==0, 0,1),
  Dt3 = param_reg_mat_fit$D_t_3,
  Zt3 = param_reg_mat_fit$Z_t_3,
  
  
  #### t=2
  M2t2 =  param_reg_mat_fit$M2_t_2,
  Z_latent_M2t2 = ifelse(param_reg_mat_fit$M2_t_2 ==0, 0,1),
  M1t2 = param_reg_mat_fit$M1_t_2,
  Z_latent_M1t2 = ifelse(param_reg_mat_fit$M1_t_2 ==0, 0,1),
  Dt2 = param_reg_mat_fit$D_t_2,
  Zt2 = param_reg_mat_fit$Z_t_2,
  
  
  #### t=1
  M2t1 =  param_reg_mat_fit$M2_t_1,
  Z_latent_M2t1 = ifelse(param_reg_mat_fit$M2_t_1 ==0, 0,1),
  M1t1 = param_reg_mat_fit$M1_t_1,
  Z_latent_M1t1 = ifelse(param_reg_mat_fit$M1_t_1 ==0, 0,1),
  Dt1 = param_reg_mat_fit$D_t_1,
  Zt1 = param_reg_mat_fit$Z_t_1,
  
  
  # data for baseline confounders
  
  
  # regression design matrices
  Y_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Y_param_reg_design_mat_fit),
  M2t3_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, M2t3_param_reg_design_mat_fit),
  M1t3_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, M1t3_param_reg_design_mat_fit),
  Dt3_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Dt3_param_reg_design_mat_fit),
  Zt3_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Zt3_param_reg_design_mat_fit),
  M2t2_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, M2t2_param_reg_design_mat_fit),
  M1t2_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, M1t2_param_reg_design_mat_fit),
  Dt2_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Dt2_param_reg_design_mat_fit),
  Zt2_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Zt2_param_reg_design_mat_fit),
  M2t1_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, M2t1_param_reg_design_mat_fit),
  M1t1_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, M1t1_param_reg_design_mat_fit),
  Dt1_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Dt1_param_reg_design_mat_fit),
  Zt1_param_reg_design_mat = cbind(L0_df_wide_fit$Intercept, Zt1_param_reg_design_mat_fit)
  
  
)


inits_fit <- list(
  # outcome model parameters
  beta_Yi = rep(0, constants_fit$num_Y_reg_coeff),
  sigma_sq_Yi = rinvgamma(1, 1, 1),
  pi0_Y = 0.5,
  Y_pos = rep(mean(param_reg_mat_fit$Y_t_3), constants_fit$n),
  
  # mediator parameters at t =3
  theta_M2_t3 = rep(0,constants_fit$num_M2t3_reg_coeff),
  theta_M1_t3 = rep(0,constants_fit$num_M1t3_reg_coeff),
  sigma_sq_M2_t3 = rinvgamma(1, 1, 1),
  sigma_sq_M1_t3 = rinvgamma(1, 1, 1),
  pi0_M2_t3 = 0.5,
  pi0_M1_t3 = 0.5,
  data_wide_M2t3_pos = rep(mean(param_reg_mat_fit$M2_t_3), constants_fit$n),
  data_wide_M1t3_pos = rep(mean(param_reg_mat_fit$M1_t_3), constants_fit$n),
  
  
  # mediator parameters at t =2
  theta_M2_t2 = rep(0,constants_fit$num_M2t2_reg_coeff),
  theta_M1_t2 = rep(0,constants_fit$num_M1t2_reg_coeff),
  sigma_sq_M2_t2 = rinvgamma(1, 1, 1),
  sigma_sq_M1_t2 = rinvgamma(1, 1, 1),
  pi0_M2_t2 = 0.5,
  pi0_M1_t2 = 0.5,
  data_wide_M2t2_pos = rep(mean(param_reg_mat_fit$M2_t_2), constants_fit$n),
  data_wide_M1t2_pos = rep(mean(param_reg_mat_fit$M1_t_2), constants_fit$n),
  
  
  # mediator parameters at t =1
  theta_M2_t1 = rep(0,constants_fit$num_M2t1_reg_coeff),
  theta_M1_t1 = rep(0,constants_fit$num_M1t1_reg_coeff),
  sigma_sq_M2_t1 = rinvgamma(1, 1, 1),
  sigma_sq_M1_t1 = rinvgamma(1, 1, 1),
  pi0_M2_t1 = 0.5,
  pi0_M1_t1 = 0.5,
  data_wide_M2t1_pos = rep(mean(param_reg_mat_fit$M2_t_1), constants_fit$n),
  data_wide_M1t1_pos = rep(mean(param_reg_mat_fit$M1_t_1), constants_fit$n),
  
  
  # treatment receipt parameters 
  theta_D_t3 = rep(0,constants_fit$num_Dt3_reg_coeff),
  theta_D_t2 = rep(0,constants_fit$num_Dt2_reg_coeff),
  theta_D_t1 = rep(0,constants_fit$num_Dt1_reg_coeff),
  
  
  # treatment assignment parameters
  theta_Z_t3 = rep(0,constants_fit$num_Zt3_reg_coeff),
  theta_Z_t2 = rep(0,constants_fit$num_Zt2_reg_coeff),
  theta_Z_t1 = rep(0,constants_fit$num_Zt1_reg_coeff)
)


Param_model_fit <- nimbleModel(code_fit, constants_fit, data_fit, inits_fit)  # model creation


compile_Param_model_fit <- compileNimble(Param_model_fit, showCompilerOutput = TRUE)  # model compilation

# MCMC configuration
config_MCMC_fit <- configureMCMC(Param_model_fit, useConjugacy = TRUE,
                                 enableWAIC = TRUE,
                                 monitors = c("beta_Yi", "sigma_sq_Yi", "pi0_Y",  
                                              "theta_M2_t3", "theta_M1_t3", "sigma_sq_M2_t3", "sigma_sq_M1_t3", "pi0_M2_t3", "pi0_M1_t3",
                                              "theta_M2_t2", "theta_M1_t2", "sigma_sq_M2_t2", "sigma_sq_M1_t2", "pi0_M2_t2", "pi0_M1_t2",
                                              "theta_M2_t1", "theta_M1_t1", "sigma_sq_M2_t1", "sigma_sq_M1_t1", "pi0_M2_t1", "pi0_M1_t1",
                                              "theta_D_t1","theta_D_t2", "theta_D_t3",
                                              "theta_Z_t1","theta_Z_t2", "theta_Z_t3",
                                              "logDens"))


config_MCMC_fit$removeSamplers('logDens')   ## remove sampler assigned to 'logDens'
config_MCMC_fit$addSampler(target = 'logDens', type = 'sumLogPostDens')   ## add our custom sampler


mcmc_fit <- buildMCMC(config_MCMC_fit)
cmcmc_fit <- compileNimble(mcmc_fit, project = Param_model_fit, showCompilerOutput = TRUE)

# burn-in 1000 iterations and save 1 iter for data generation
num_iter_fit = 15000;  num_burnin_fit = 10000;  num_thin_fit = 5;  num_chains_fit= 1


samples_allChains_fit <- runMCMC(cmcmc_fit, niter = num_iter_fit, nburnin = num_burnin_fit,
                                 nchains = num_chains_fit, 
                                 setSeed = TRUE, thin = num_thin_fit, WAIC = TRUE)

post_samples_model_fit <- samples_allChains_fit$samples   # list of nchains matrices


post_samples_model_fit <- as.data.frame(post_samples_model_fit)


sample_L0_bin_multiple = function(L0_obs_bin_data_no_intercept, # dataframe of all binary baseline confounders
                                  num_MC_samples
){
  
  probs = colMeans(L0_obs_bin_data_no_intercept)  
  
  return_df = as.data.frame(
    matrix(rbinom(num_MC_samples * length(probs), size = 1, prob = rep(probs, each = num_MC_samples)), 
           nrow = num_MC_samples, ncol = length(probs))
  )
  
  
  return(return_df)
  
}


sample_Zt = function(MCdata,                             # MC data upto this point in the temporal order
                     Param_coeffs_df,
                     num_MC_samples,                     # number of MC samples to be generated
                     #num_L0,                             # number of L0 INCLUDING intercept
                     t,                                  # t = {1, ..., T}
                     it                                  # an integer from 1:Q (number of post iterations)
){
  
  # get parameter samples for this iteration at time = t
  param_coeff_theta_Zt_vec = c()
  for(k in 1:ncol(MCdata)){
    theta_Zt_colnames_ordered = paste0("theta_Z_t",t,"[",k,"]")
    
    theta_Zt_colnames =  theta_Zt_colnames_ordered[theta_Zt_colnames_ordered %in%
                                                     colnames(Param_coeffs_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_theta_Zt_vec[k] = as.numeric(sapply(theta_Zt_colnames, 
                                                    function(col_name) Param_coeffs_df[[col_name]][it]))
  }
  
  
  Z_t_prob = pnorm(as.matrix(MCdata) %*% param_coeff_theta_Zt_vec)
  
  Z_t_vec = rbinom(n = num_MC_samples, size = 1, prob = Z_t_prob)
  
  return(Z_t_vec)
  
}


sample_Dt = function(MCdata,                             # MC data upto this point in the temporal order
                     Param_coeffs_df,
                     num_MC_samples,                     # number of MC samples to be generated
                     #num_L0,                             # number of L0 INCLUDING intercept
                     t,                                  # t = {1, ..., T}
                     it                                  # an integer from 1:Q (number of post iterations)
){
  
  # get parameter samples for this iteration at time = t
  param_coeff_theta_Dt_vec = c()
  for(k in 1:ncol(MCdata)){
    theta_Dt_colnames_ordered = paste0("theta_D_t",t,"[",k,"]")
    
    theta_Dt_colnames =  theta_Dt_colnames_ordered[theta_Dt_colnames_ordered %in%
                                                     colnames(Param_coeffs_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_theta_Dt_vec[k] = as.numeric(sapply(theta_Dt_colnames, 
                                                    function(col_name) Param_coeffs_df[[col_name]][it]))
  }
  
  
  D_t_prob = pnorm(as.matrix(MCdata) %*% param_coeff_theta_Dt_vec)
  
  D_t_vec = rbinom(n = num_MC_samples, size = 1, prob = D_t_prob)
  
  return(D_t_vec)
  
}


sample_allMediators_Mt = function(MCdata,                # MC data upto this point in the temporal order
                                  Param_coeffs_df,
                                  num_MC_samples,                     # number of MC samples to be generated
                                  number_of_mediators,                # number of mediators
                                  t,                                  # t = {1, ..., T}
                                  it                                  # an integer from 1:Q (number of post iterations)
){
  
  param_coeff_theta_Mt_mat = matrix(NA, nrow = ncol(MCdata), ncol = number_of_mediators) # each column contains regression mean coeffs for a mediator
  param_coeff_sigSq_Mt_vec = c()
  param_coeff_pi0_Mt_vec = c()
  for(j in 1:number_of_mediators){
    theta_Mt_colnames_ordered = unlist(
      lapply(1:ncol(MCdata), function(k) paste0("theta_M",j,"_t",t,"[",k,"]"))
    )
    
    theta_Mt_colnames =  theta_Mt_colnames_ordered[theta_Mt_colnames_ordered %in%
                                                     colnames(Param_coeffs_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_theta_Mt_mat[,j] = as.numeric(sapply(theta_Mt_colnames, 
                                                     function(col_name) Param_coeffs_df[[col_name]][it]))
    
    
    sigSq_Mt_colnames_ordered = paste0("sigma_sq_M",j,"_t",t)
    
    sigSq_Mt_colnames =  sigSq_Mt_colnames_ordered[sigSq_Mt_colnames_ordered %in%
                                                     colnames(Param_coeffs_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_sigSq_Mt_vec[j] = as.numeric(sapply(sigSq_Mt_colnames, 
                                                    function(col_name) Param_coeffs_df[[col_name]][it]))
    
    
    pi0_Mt_colnames_ordered = paste0("pi0_M",j,"_t",t)
    
    pi0_Mt_colnames =  pi0_Mt_colnames_ordered[pi0_Mt_colnames_ordered %in%
                                                 colnames(Param_coeffs_df)]
    
    param_coeff_pi0_Mt_vec[j] = as.numeric(sapply(pi0_Mt_colnames, 
                                                  function(col_name) Param_coeffs_df[[col_name]][it]))
    
    
  }
  
  
  ##################Get coeffs for the random effects###############################
  Sigma_b_M1M2_t_col_name = sapply(1:number_of_mediators, function(mediator_num) {
    paste0("Sigma_b_M1M2_t",t,"[", mediator_num, ", ", 1:number_of_mediators, "]")  # Column names for Sigma_b_M1M2_t1[1,1] to Sigma_b_M1M2_t1[2,2]
  })
  
  # matrix of dimension J by J for J mediators
  Sigma_b_M1M2_t_mat = apply(Sigma_b_M1M2_t_col_name, c(1,2),
                             function(col_name){Param_coeffs_df[[col_name]][it]}) 
  
  #b_M1M2_t_df is a df of dimensions num_MC_samples by length(zero_vec)
  b_M1M2_t_mat = as.matrix(rmvnorm(n = num_MC_samples, mean = zero_vec, sigma = Sigma_b_M1M2_t_mat))   #random intercepts
  
  
  # Generate both M1 and M2
  
  #M1M2_pos_indicator_df is a df of dimensions num_MC_samples by length(zero_vec)
  M1M2_pos_indicator_df =  as.data.frame(
    matrix(rbinom(num_MC_samples * length(param_coeff_pi0_Mt_vec), 
                  size = 1, 
                  prob = rep(1-param_coeff_pi0_Mt_vec, each = num_MC_samples)),
           nrow = num_MC_samples, ncol = length(param_coeff_pi0_Mt_vec))
  )
  
  M1M2_pos_mean_df = (as.matrix(MCdata) %*% param_coeff_theta_Mt_mat) +  b_M1M2_t_mat # Result: [num_MC_samples × num_of_mediators]
  
  truncnorm_samples = matrix(NA, nrow = num_MC_samples, ncol = number_of_mediators)
  for (j in 1:number_of_mediators) {
    truncnorm_samples[, j] = rtruncnorm(
      n = num_MC_samples,
      a = 0,                
      b = Inf,                 
      mean = M1M2_pos_mean_df[, j],     
      sd = sqrt(param_coeff_sigSq_Mt_vec[j])           
    ) * M1M2_pos_indicator_df[,j]
  }
  
  
  return_df = as.data.frame(truncnorm_samples)
  return(return_df)
}


sample_Mediator_Mt = function(MCdata,                # MC data upto this point in the temporal order
                              Param_coeffs_df,
                              num_MC_samples,                     # number of MC samples to be generated
                              mediator_number,                   # 1 to J (num of mediators)
                              t,                                  # t = {1, ..., T}
                              it                                  # an integer from 1:Q (number of post iterations)
){
  
  
  param_coeff_theta_Mt_vec = c()
  for(k in 1:ncol(MCdata)){
    theta_Mt_colnames_ordered = paste0("theta_M",mediator_number,"_t",t,"[",k,"]")
    
    
    theta_Mt_colnames =  theta_Mt_colnames_ordered[theta_Mt_colnames_ordered %in%
                                                     colnames(Param_coeffs_df)]
    #df with nrow = N.M, ncol = num_L0
    param_coeff_theta_Mt_vec[k] = as.numeric(sapply(theta_Mt_colnames, 
                                                    function(col_name) Param_coeffs_df[[col_name]][it]))
  }
  
  
  sigSq_Mt_colnames_ordered = paste0("sigma_sq_M",mediator_number,"_t",t)
  
  sigSq_Mt_colnames =  sigSq_Mt_colnames_ordered[sigSq_Mt_colnames_ordered %in%
                                                   colnames(Param_coeffs_df)]
  #df with nrow = N.M, ncol = num_L0
  param_coeff_sigSq_Mt = as.numeric(sapply(sigSq_Mt_colnames, 
                                           function(col_name) Param_coeffs_df[[col_name]][it]))
  
  
  pi0_Mt_colnames_ordered = paste0("pi0_M",mediator_number,"_t",t)
  
  pi0_Mt_colnames =  pi0_Mt_colnames_ordered[pi0_Mt_colnames_ordered %in%
                                               colnames(Param_coeffs_df)]
  
  param_coeff_pi0_Mt = as.numeric(sapply(pi0_Mt_colnames, 
                                         function(col_name) Param_coeffs_df[[col_name]][it]))
  
  
  # Generate both M1 and M2
  
  #Mt_pos_indicator_df is a df of dimensions num_MC_samples by length(zero_vec)
  Mt_pos_indicator =  rbinom(num_MC_samples, 
                             size = 1, 
                             prob = rep(1-param_coeff_pi0_Mt))
  
  Mt_pos_mean_df = (as.matrix(MCdata) %*% param_coeff_theta_Mt_vec)  # Result: [num_MC_samples × 1]
  
  
  truncnorm_samples = rtruncnorm(
    n = num_MC_samples,
    a = 0,                
    b = 10,                 
    mean = Mt_pos_mean_df,     
    sd = sqrt(param_coeff_sigSq_Mt)           
  ) * Mt_pos_indicator
  
  
  return(as.numeric(truncnorm_samples))
}


compute_Y_mean = function(MCdata,                # MC data upto this point in the temporal order
                          Param_coeffs_df,
                          num_MC_samples,                     # number of MC samples to be generated
                          it                                  # an integer from 1:Q (number of post iterations)
){
  
  param_coeff_beta_Yi_vec = c() 
  for(k in 1:ncol(MCdata)){
    beta_Yi_colnames_ordered = paste0("beta_Yi[",k,"]")
    
    beta_Yi_colnames =  beta_Yi_colnames_ordered[beta_Yi_colnames_ordered %in%
                                                   colnames(Param_coeffs_df)]
    param_coeff_beta_Yi_vec[k] = as.numeric(sapply(beta_Yi_colnames, 
                                                   function(col_name) Param_coeffs_df[[col_name]][it]))
  }
  
  
  sigSq_Yi_colnames_ordered = paste0("sigma_sq_Yi")
  sigSq_Yi_colnames =  sigSq_Yi_colnames_ordered[sigSq_Yi_colnames_ordered %in%
                                                   colnames(Param_coeffs_df)]
  param_coeff_sigSq_Yi = as.numeric(sapply(sigSq_Yi_colnames, 
                                           function(col_name) Param_coeffs_df[[col_name]][it]))
  
  
  pi0_Y_colnames_ordered = paste0("pi0_Y")
  pi0_Y_colnames =  pi0_Y_colnames_ordered[pi0_Y_colnames_ordered %in%
                                             colnames(Param_coeffs_df)]
  param_coeff_pi0_Y = as.numeric(sapply(pi0_Y_colnames, 
                                        function(col_name) Param_coeffs_df[[col_name]][it]))
  
  
  # Generate Yi
  #Y_pos_indicator_vec is a vec of length num_MC_samples 
  Y_pos_indicator_vec =  rbinom(n = num_MC_samples,
                                size = 1,
                                prob = 1-param_coeff_pi0_Y)
  Y_pos_mean_vec = (as.matrix(MCdata) %*% param_coeff_beta_Yi_vec)
  
  
  #compute mean of truncated Normal
  z = -Y_pos_mean_vec/sqrt(param_coeff_sigSq_Yi)
  phi_z = dnorm(z)  # PDF of the standard normal at z
  Phi_z = pnorm(z)  # CDF of the standard normal at z
  trunc_norm_param_denom = pmax(1 - Phi_z, 1e-10)
  trunc_norm_mean_vec = Y_pos_mean_vec + (phi_z / trunc_norm_param_denom) * sqrt(param_coeff_sigSq_Yi)
  
  
  # E[Y|...] = (1-pi0) * E[trunc_norm]
  mean_Y_hurdle = trunc_norm_mean_vec * (1 - param_coeff_pi0_Y)

  return(mean_Y_hurdle)
}


MCInteg <- function(var.type, # Vector of variable specifications for data.  Fi=fixed (e.g. the exposure), D=trt received, M1= mediator 1,  M2= mediator 2, Y=outcome. 
                    Param_coeffs, #fitted BNP model parameter samples list
                    fixed.regime.trt, # A vector specifying the Exposure  regime Z
                    fixed.regime.control, # A vector specifying the Exposure  regime Z_star
                    U1,       # could be one of: "c" for VALUE ATTENTIVE, "d" for PRICE ATTENTIVE, "n" for NON-ACTIVE CUSTOMERS and "a" for ACTIVE CUSTOMERS
                    J=10000, # Size of pseudo data. Default is set to 2,000.
                    Ndraws=1000, # Number of posterior draws. Default is set to 200.
                    L0_obs_data_no_intercept, #matrix of all baseline confounders excluding the intercept
                    num_Bin_L0,
                    num_Cont_L0,
                    max_t,
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
  
  
  for(it in 1:Ndraws) {
    zt <- 1                               # zt indexes FIXED longitudinal trt regime in G-comp
    dt <- 1
    mt <- 1
    intercept = rep(1, J)
    
    
    ################### We only have binary L0 in replicated dataframe###########################################   
    L0_bin_sampled_mat_iter = NULL # initialize a matrix to store samples
    L0_bin_sampled_mat_iter = sample_L0_bin_multiple(L0_obs_bin_data_no_intercept = L0_obs_data_no_intercept[, 1:num_Bin_L0], 
                                                     num_MC_samples= J)
    L0_data = cbind(intercept, L0_bin_sampled_mat_iter)
    
    
    Zt_mat_MCdata_trt = matrix(NA, nrow = J, ncol = max_t)
    Dt_mat_MCdata_trt = matrix(NA, nrow = J, ncol = max_t)
    M1t_mat_MCdata_trt = matrix(NA, nrow = J, ncol = max_t)
    M2t_mat_MCdata_trt = matrix(NA, nrow = J, ncol = max_t)
    
    
    Zt_mat_MCdata_control = matrix(NA, nrow = J, ncol = max_t)
    Dt_mat_MCdata_control = matrix(NA, nrow = J, ncol = max_t)
    M1t_mat_MCdata_control = matrix(NA, nrow = J, ncol = max_t)
    M2t_mat_MCdata_control = matrix(NA, nrow = J, ncol = max_t)
    
    
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
    
    
    for (j in (n_L0+1):length(var.type)) {                                         
      
      if(var.type[j] == "Fi") {
        Zt_mat_MCdata_trt[,zt] = rep(fixed.regime.trt[zt], J)
        Zt_mat_MCdata_control[,zt] = rep(fixed.regime.control[zt], J)
        zt = zt + 1
      }
      else if(var.type[j] == "Dt") {
        
        if(dt == 1){
          Dt_mat_MCdata_trt[,dt] = D1_trt
          Dt_mat_MCdata_control[,dt] = D1_control
        }else{
          Dt_sampled_vec_iter_trt =sample_Dt(MCdata = cbind(L0_data, Dt_mat_MCdata_trt[,dt-1],
                                                            M1t_mat_MCdata_trt[,dt-1], M2t_mat_MCdata_trt[,dt-1],
                                                            Zt_mat_MCdata_trt[,dt]),                             
                                             Param_coeffs_df = Param_coeffs,
                                             num_MC_samples = J,                     
                                             t = dt,                                  
                                             it = it)
          Dt_sampled_vec_iter_control =sample_Dt(MCdata = cbind(L0_data, Dt_mat_MCdata_control[,dt-1],
                                                                M1t_mat_MCdata_control[,dt-1], M2t_mat_MCdata_control[,dt-1],
                                                                Zt_mat_MCdata_control[,dt]),                            
                                                 Param_coeffs_df = Param_coeffs,
                                                 num_MC_samples = J,                     
                                                 t = dt,                                  
                                                 it = it)
          
          
          Dt_mat_MCdata_trt[,dt] =  Dt_sampled_vec_iter_trt
          Dt_mat_MCdata_control[,dt] =  Dt_sampled_vec_iter_control
        }
        dt = dt +1
      }
      else if(var.type[j] == "M1t") { 
        if(mt == 1){
          M1t_sampled_vec_iter_trt = sample_Mediator_Mt(MCdata = cbind(L0_data, Zt_mat_MCdata_trt[,1], Dt_mat_MCdata_trt[,1]),
                                                        Param_coeffs_df = Param_coeffs,
                                                        num_MC_samples = J,                     
                                                        mediator_number =  1,                
                                                        t = mt,                                  
                                                        it = it)
          
          
          M1t_sampled_vec_iter_control =sample_Mediator_Mt(MCdata = cbind(L0_data, Zt_mat_MCdata_control[,1], Dt_mat_MCdata_control[,1]),
                                                           Param_coeffs_df = Param_coeffs,
                                                           num_MC_samples = J,                     
                                                           mediator_number =  1,                
                                                           t = mt,                                  
                                                           it = it)
          
          
          M1t_mat_MCdata_trt[,1] =  M1t_sampled_vec_iter_trt
          M1t_mat_MCdata_control[,1] =  M1t_sampled_vec_iter_control
          
        }else{
          M1t_sampled_vec_iter_trt = sample_Mediator_Mt(MCdata = cbind(L0_data, M1t_mat_MCdata_trt[,mt-1], M2t_mat_MCdata_trt[,mt-1],
                                                                       Zt_mat_MCdata_trt[,mt], Dt_mat_MCdata_trt[,mt]),
                                                        Param_coeffs_df = Param_coeffs,
                                                        num_MC_samples = J,                     
                                                        mediator_number =  1,                
                                                        t = mt,                                  
                                                        it = it)
          
          
          M1t_sampled_vec_iter_control =sample_Mediator_Mt(MCdata = cbind(L0_data, M1t_mat_MCdata_control[,mt-1], M2t_mat_MCdata_control[,mt-1],
                                                                          Zt_mat_MCdata_control[,mt], Dt_mat_MCdata_control[,mt]),
                                                           Param_coeffs_df = Param_coeffs,
                                                           num_MC_samples = J,                     
                                                           mediator_number =  1,                
                                                           t = mt,                                  
                                                           it = it)
          
          
          M1t_mat_MCdata_trt[,mt] =  M1t_sampled_vec_iter_trt
          M1t_mat_MCdata_control[,mt] =  M1t_sampled_vec_iter_control  
        }
        
        
      }
      else if(var.type[j] == "M2t") { 
        if(mt == 1){
          M2t_sampled_vec_iter_trt = sample_Mediator_Mt(MCdata = cbind(L0_data, Zt_mat_MCdata_trt[,1], Dt_mat_MCdata_trt[,1],
                                                                       M1t_mat_MCdata_trt[,1]),
                                                        Param_coeffs_df = Param_coeffs,
                                                        num_MC_samples = J,                     
                                                        mediator_number =  2,                
                                                        t = mt,                                  
                                                        it = it)
          
          
          M2t_sampled_vec_iter_control =sample_Mediator_Mt(MCdata = cbind(L0_data, Zt_mat_MCdata_control[,1], Dt_mat_MCdata_control[,1],
                                                                          M1t_mat_MCdata_control[,1]),
                                                           Param_coeffs_df = Param_coeffs,
                                                           num_MC_samples = J,                     
                                                           mediator_number =  2,                
                                                           t = mt,                                  
                                                           it = it)
          
          M2t_mat_MCdata_trt[,1] =  M2t_sampled_vec_iter_trt
          M2t_mat_MCdata_control[,1] =  M2t_sampled_vec_iter_control
          
          
        }else{
          M2t_sampled_vec_iter_trt = sample_Mediator_Mt(MCdata = cbind(L0_data, M2t_mat_MCdata_trt[,mt-1],
                                                                       Zt_mat_MCdata_trt[,mt], Dt_mat_MCdata_trt[,mt],
                                                                       M1t_mat_MCdata_trt[,mt]),
                                                        Param_coeffs_df = Param_coeffs,
                                                        num_MC_samples = J,                     
                                                        mediator_number =  2,                
                                                        t = mt,                                  
                                                        it = it)
          
          
          M2t_sampled_vec_iter_control =sample_Mediator_Mt(MCdata = cbind(L0_data,  M2t_mat_MCdata_control[,mt-1],
                                                                          Zt_mat_MCdata_control[,mt], Dt_mat_MCdata_control[,mt],
                                                                          M1t_mat_MCdata_control[,mt]),
                                                           Param_coeffs_df = Param_coeffs,
                                                           num_MC_samples = J,                     
                                                           mediator_number =  2,                
                                                           t = mt,                                  
                                                           it = it)
          
          M2t_mat_MCdata_trt[,mt] =  M2t_sampled_vec_iter_trt
          M2t_mat_MCdata_control[,mt] =  M2t_sampled_vec_iter_control 
          
        }
        
        
        mt = mt + 1
        
      }
      
      else if(var.type[j] == "Y") {
        Y_mean_vec = compute_Y_mean(MCdata = cbind(L0_data,  
                                                   Zt_mat_MCdata_trt[,ncol(Zt_mat_MCdata_trt)],
                                                   Dt_mat_MCdata_trt[,ncol(Dt_mat_MCdata_trt)],
                                                   M1t_mat_MCdata_control[,ncol(M1t_mat_MCdata_control)],
                                                   M2t_mat_MCdata_control[,ncol(M2t_mat_MCdata_control)]),                
                                    Param_coeffs_df = Param_coeffs,
                                    num_MC_samples = J,                     
                                    it = it)
        
        
        Y_zzstar_mat[1,it] = mean(Y_mean_vec)
        MC_error_vec[it] = sd(Y_mean_vec) / sqrt(length(Y_mean_vec))
        if(it == 500){
          mc_error = sd(Y_mean_vec) / sqrt(length(Y_mean_vec))
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
  quantiles = quantile(posterior_samples, c(0.025, 0.975))
  interval = as.numeric(format(round(quantiles, 4), nsmall = 4))
  return_df = data.frame(posterior_mean = posterior_mean, 
                         lower_95CI = interval[1],   upper_95CI = interval[2])
  return(return_df)
}


plot_ConfInt <-function(conf_int_mat) {
  conf_int_df = data.frame(conf_int_mat)
  conf_int_df$CI <- 1:nrow(conf_int_df)
  ci_long <- reshape(conf_int_df,
                     direction = "long",
                     varying = setdiff(names(conf_int_df), "CI"),
                     v.names = "Value",
                     timevar = "Bound",
                     times = setdiff(names(conf_int_df), "CI"))
  ggplot(ci_long, aes(x = CI, ymin = Value, ymax = Value, color = Bound )) +
    geom_errorbar(width = 0.2)+
    labs(x = "Visits", y = "Confidence/Credible Intervals")+
    scale_color_manual(values = c("lower" = "red", "upper" = "blue")) +
    theme_minimal()
  
}


input_varType1 = c("L0","L0","L0","Fi", "Dt", "M1t","M2t", "Fi", "Dt", "M1t","M2t","Fi", "Dt", "M1t","M2t", "Y")


#########################################ACTIVE CUSTOMERS#################################################


active_customer_theta_zzstar = MCInteg(var.type = input_varType1,
                                    Param_coeffs = post_samples_model_fit,
                                    fixed.regime.trt = rep(1,3),
                                    fixed.regime.control = rep(0,3),
                                    U1 = "a",            
                                    J=10000, 
                                    Ndraws= 1000, 
                                    L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                    num_Bin_L0 = 2, 
                                    num_Cont_L0 = 0,
                                    max_t = 3)


active_customer_theta_zz = MCInteg(var.type = input_varType1,
                                Param_coeffs = post_samples_model_fit,
                                fixed.regime.trt = rep(1,3),
                                fixed.regime.control = rep(1,3),
                                U1 = "a",
                                J=10000, 
                                Ndraws= 1000, 
                                L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                num_Bin_L0 = 2, 
                                num_Cont_L0 = 0,
                                max_t = 3)


active_customer_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                        Param_coeffs = post_samples_model_fit,
                                        fixed.regime.trt = rep(0,3),
                                        fixed.regime.control = rep(0,3),
                                        U1 = "a",
                                        J=10000, 
                                        Ndraws= 1000, 
                                        L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                        num_Bin_L0 = 2, 
                                        num_Cont_L0 = 0,
                                        max_t = 3)


active_customer_PIDE = compute_PIDE(active_customer_theta_zzstar, active_customer_theta_zstarzstar)
active_customer_PIDE_CI = calculate_95bayesian_credible_interval(active_customer_PIDE) 

active_customer_PJIIE = compute_PJIIE(active_customer_theta_zz, active_customer_theta_zzstar)
active_customer_PJIIE_CI = calculate_95bayesian_credible_interval(active_customer_PJIIE) 

active_customer_PCE = compute_PCE(active_customer_PIDE,active_customer_PJIIE)
active_customer_PCE_CI = calculate_95bayesian_credible_interval(active_customer_PCE) 


#########################################Value Attentive#################################################


value_attentive_theta_zzstar = MCInteg(var.type = input_varType1,
                                Param_coeffs = post_samples_model_fit,
                                fixed.regime.trt = rep(1,3),
                                fixed.regime.control = rep(0,3),
                                U1 = "c",            
                                J=10000, 
                                Ndraws= 1000, 
                                L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                num_Bin_L0 = 2, 
                                num_Cont_L0 = 0,
                                max_t = 3)


value_attentive_theta_zz = MCInteg(var.type = input_varType1,
                            Param_coeffs = post_samples_model_fit,
                            fixed.regime.trt = rep(1,3),
                            fixed.regime.control = rep(1,3),
                            U1 = "c",
                            J=10000, 
                            Ndraws= 1000, 
                            L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                            num_Bin_L0 = 2, 
                            num_Cont_L0 = 0,
                            max_t = 3)


value_attentive_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                    Param_coeffs = post_samples_model_fit,
                                    fixed.regime.trt = rep(0,3),
                                    fixed.regime.control = rep(0,3),
                                    U1 = "c",
                                    J=10000, 
                                    Ndraws= 1000, 
                                    L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                    num_Bin_L0 = 2, 
                                    num_Cont_L0 = 0,
                                    max_t = 3)


value_attentive_PIDE = compute_PIDE(value_attentive_theta_zzstar, value_attentive_theta_zstarzstar)
value_attentive_PIDE_CI = calculate_95bayesian_credible_interval(value_attentive_PIDE) 

value_attentive_PJIIE = compute_PJIIE(value_attentive_theta_zz, value_attentive_theta_zzstar)
value_attentive_PJIIE_CI = calculate_95bayesian_credible_interval(value_attentive_PJIIE) 

value_attentive_PCE = compute_PCE(value_attentive_PIDE,value_attentive_PJIIE)
value_attentive_PCE_CI = calculate_95bayesian_credible_interval(value_attentive_PCE) 


#########################################Price Attentive#################################################


price_attentive_theta_zzstar = MCInteg(var.type = input_varType1,
                              Param_coeffs = post_samples_model_fit,
                              fixed.regime.trt = rep(1,3),
                              fixed.regime.control = rep(0,3),
                              U1 = "d",            
                              J=10000, 
                              Ndraws= 1000, 
                              L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                              num_Bin_L0 = 2, 
                              num_Cont_L0 = 0,
                              max_t = 3)


price_attentive_theta_zz = MCInteg(var.type = input_varType1,
                          Param_coeffs = post_samples_model_fit,
                          fixed.regime.trt = rep(1,3),
                          fixed.regime.control = rep(1,3),
                          U1 = "d",
                          J=10000, 
                          Ndraws= 1000, 
                          L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                          num_Bin_L0 = 2, 
                          num_Cont_L0 = 0,
                          max_t = 3)


price_attentive_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                  Param_coeffs = post_samples_model_fit,
                                  fixed.regime.trt = rep(0,3),
                                  fixed.regime.control = rep(0,3),
                                  U1 = "d",
                                  J=10000, 
                                  Ndraws= 1000, 
                                  L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                  num_Bin_L0 = 2, 
                                  num_Cont_L0 = 0,
                                  max_t = 3)


price_attentive_PIDE = compute_PIDE(price_attentive_theta_zzstar, price_attentive_theta_zstarzstar)
price_attentive_PIDE_CI = calculate_95bayesian_credible_interval(price_attentive_PIDE) 

price_attentive_PJIIE = compute_PJIIE(price_attentive_theta_zz, price_attentive_theta_zzstar)
price_attentive_PJIIE_CI = calculate_95bayesian_credible_interval(price_attentive_PJIIE) 

price_attentive_PCE = compute_PCE(price_attentive_PIDE,price_attentive_PJIIE)
price_attentive_PCE_CI = calculate_95bayesian_credible_interval(price_attentive_PCE) 


#########################################Non-Active Customers#################################################


non_active_customer_theta_zzstar = MCInteg(var.type = input_varType1,
                                   Param_coeffs = post_samples_model_fit,
                                   fixed.regime.trt = rep(1,3),
                                   fixed.regime.control = rep(0,3),
                                   U1 = "n",            
                                   J=10000, 
                                   Ndraws= 1000, 
                                   L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                   num_Bin_L0 = 2, 
                                   num_Cont_L0 = 0,
                                   max_t = 3)


non_active_customer_theta_zz = MCInteg(var.type = input_varType1,
                               Param_coeffs = post_samples_model_fit,
                               fixed.regime.trt = rep(1,3),
                               fixed.regime.control = rep(1,3),
                               U1 = "n",
                               J=10000, 
                               Ndraws= 1000, 
                               L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                               num_Bin_L0 = 2, 
                               num_Cont_L0 = 0,
                               max_t = 3)


non_active_customer_theta_zstarzstar = MCInteg(var.type = input_varType1,
                                       Param_coeffs = post_samples_model_fit,
                                       fixed.regime.trt = rep(0,3),
                                       fixed.regime.control = rep(0,3),
                                       U1 = "n",
                                       J=10000, 
                                       Ndraws= 1000, 
                                       L0_obs_data_no_intercept = L0_df_wide_fit[2:3],
                                       num_Bin_L0 = 2, 
                                       num_Cont_L0 = 0,
                                       max_t = 3)


non_active_customer_PIDE = compute_PIDE(non_active_customer_theta_zzstar, non_active_customer_theta_zstarzstar)
non_active_customer_PIDE_CI = calculate_95bayesian_credible_interval(non_active_customer_PIDE) 

non_active_customer_PJIIE = compute_PJIIE(non_active_customer_theta_zz, non_active_customer_theta_zzstar)
non_active_customer_PJIIE_CI = calculate_95bayesian_credible_interval(non_active_customer_PJIIE) 

non_active_customer_PCE = compute_PCE(non_active_customer_PIDE,non_active_customer_PJIIE)
non_active_customer_PCE_CI = calculate_95bayesian_credible_interval(non_active_customer_PCE) 


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
  PCE_NAC_upperCI     = as.numeric(non_active_customer_PCE_CI[3])
)


write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
gc()


#-------------------------------------------------------------------------------------------------------------------#
