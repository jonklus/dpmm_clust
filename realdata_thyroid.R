## DPMM Real Data Example

# load necessary functions
source("./Multivariate_DPMM_unknownvar_DEV.R")
# source("./Multivariate_DPMM_unknownvar_DEE.R")
source("./Multivariate_DPMM_unknownvar_UVV.R")
source("./posterior_helper_fxns.R")
source("./post_processing_inf.R")

# load R libraries
library(ggplot2)
library(gridExtra)
library(plotly)
library(label.switching)
library(LaplacesDemon)
library(parallel)
library(stringr)
library(mclust)

raw_data = mclust::thyroid
exposure = scale(log(raw_data[,3:5])) # data.frame(scale(raw_data[,2:6])) # center and scale exposure data
# exposure$Diagnosis = factor(raw_data$Diagnosis)
data = data.frame(exposure)
data$Diagnosis = raw_data$Diagnosis

# put exposure data into format that model will accept
y = lapply(X = 1:nrow(exposure), 
           FUN = function(x){
             matrix(exposure[x,], nrow = length(exposure[x,]))
           })

SLURM_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

model = c("conjDEV", "conjDEE", "conjUVV")[1]
SM = c(TRUE,FALSE)[SLURM_ID] # do split merge
cat("\n", model, SM, "\n")



# fit model

if(model == "conjDEV"){
  
  mod1 = MVN_CRP_sampler_DEV(
    S = 12000, seed = 516, y = y,
    alpha = 1, r = 10, g = 1, h = 10,
    sigma_hyperprior = FALSE, fix_r = FALSE,
    mu0 = matrix(rep(0, nrow(y[[1]])), ncol = 1),
    a = 1, b = 10,
    k_init = 1, diag_weights = FALSE,
    verbose = TRUE, split_merge = SM)
  
} else if(model == "conjUVV"){
  
  lambda_mat = diag(x = 5, nrow = nrow(y[[1]]))
  offdiag_mag = 0.2*5
  lambda_mat[1,2] = offdiag_mag
  lambda_mat[1,3] = -offdiag_mag
  lambda_mat[2,1] = offdiag_mag
  lambda_mat[2,3] = -offdiag_mag
  lambda_mat[3,1] = -offdiag_mag
  lambda_mat[3,2] = -offdiag_mag
  
  mod1 = MVN_CRP_sampler_UVV(
    S = 12000, seed = 516, y = y,
    alpha = 1, r = 10, g = 1, h = 10, nu = 4, 
    fix_r = FALSE,
    mu0 = matrix(rep(0, nrow(y[[1]])), ncol = 1),
    lambda0 = lambda_mat,
    k_init = 1, diag_weights = FALSE,
    verbose = TRUE, split_merge = SM)
  
}



saveRDS(object = mod1, file = paste0("./RealData/thyroid_modfit_", model, 
                                     "_", ifelse(SM == TRUE, "withSM", "noSM"),".rds"))
# 
# 
# 
# saveRDS(object = mod2, file = "../MCMC_Runs/thyroid_UVV_modfit.rds")