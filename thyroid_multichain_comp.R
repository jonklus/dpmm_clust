## thyroid data model fit

# load necessary functions
source("./Multivariate_DPMM_unknownvar_DEV.R")
# source("./Multivariate_DPMM_unknownvar_DEE.R")
source("./Multivariate_DPMM_unknownvar_UVV.R")
source("./posterior_helper_fxns.R")
source("./post_processing_inf.R")

# load R libraries
library(dplyr)
# library(ggplot2)
# library(gridExtra)
# library(ggpubr)
# library(ggforce) # draw circles
# library(stringr)
library(mclust)

raw_data = mclust::thyroid
exposure = scale(log(raw_data[,3:5])) # data.frame(scale(raw_data[,2:6])) # center and scale exposure data
# exposure$Diagnosis = factor(raw_data$Diagnosis)
data = data.frame(exposure)
y = lapply(X = 1:nrow(exposure), 
           FUN = function(x){
             matrix(exposure[x,], nrow = length(exposure[x,]))
           })

# sample 4 seeds for multichain run
seeds = sample(x = 1:10^4, size = 4, replace = FALSE)

# sbatch-fed settings
i = as.numeric(Sys.getenv("i")) # mod type
mod = c("DEV_noSM", "DEV_withSM", "UVV_noSM", "UVV_withSM")[i]
SLURM_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # seed & start

dir_name = paste0("./Thyroid")
# stringr::str_replace_all(string = Sys.Date(),
#                          pattern = "-",
#                          replacement = "_"))

if(dir.exists(dir_name) == FALSE){
  dir.create(dir_name)
} 


# fit model
seed = seeds[i]

if(mod == "DEV_noSM"){
  
  mod1 = MVN_CRP_sampler_DEV(
    S = 12000, seed = seed, y = y,
    alpha = 1, r = 10, g = 1, h = 10,
    sigma_hyperprior = FALSE, fix_r = FALSE,
    mu0 = matrix(rep(0, nrow(y[[1]])), ncol = 1),
    a = 1, b = 10,
    k_init = 1, diag_weights = FALSE,
    verbose = TRUE, split_merge = FALSE)
  
  mod1_sum = dpmm_summary(output = mod1,
                          print_phi_sum = TRUE,
                          print_k_sum = TRUE,
                          make_traceplot = FALSE,
                          burn_in = 2000, t_hold = 250,
                          num_dims = 3,
                          calc_perf = FALSE,
                          equal_var = FALSE)
  
  # rm(mod1) # free up memory
  
  mod1_sum$settings
  mod1_sum$mean_summary
  mod1_sum$var_summary
  mod1_sum$splitmerge_accept
  mod1_sum$k_relfreqtab
  mod1_sum$fit_runtime
  
} else if(mod == "DEV_withSM"){
  
  mod1 = MVN_CRP_sampler_DEV(
    S = 12000, seed = seed, y = y,
    alpha = 1, r = 10, g = 1, h = 10,
    sigma_hyperprior = FALSE, fix_r = FALSE,
    mu0 = matrix(rep(0, nrow(y[[1]])), ncol = 1),
    a = 1, b = 10,
    k_init = 1, diag_weights = FALSE,
    verbose = TRUE, split_merge = TRUE)
  
  mod1_sum = dpmm_summary(output = mod1,
                          print_phi_sum = TRUE,
                          print_k_sum = TRUE,
                          make_traceplot = FALSE,
                          burn_in = 2000, t_hold = 250,
                          num_dims = 3,
                          calc_perf = FALSE,
                          equal_var = FALSE)
  
  
  # rm(mod1) # free up memory
  
  mod1_sum$settings
  mod1_sum$mean_summary
  mod1_sum$var_summary
  mod1_sum$splitmerge_accept
  mod1_sum$k_relfreqtab
  mod1_sum$fit_runtime
  
} else if(mod == "UVV_noSM"){
  
  lambda_mat = diag(x = 5, nrow = nrow(y[[1]]))
  offdiag_mag = 0.2*5
  lambda_mat[1,2] = offdiag_mag
  lambda_mat[1,3] = -offdiag_mag
  lambda_mat[2,1] = offdiag_mag
  lambda_mat[2,3] = -offdiag_mag
  lambda_mat[3,1] = -offdiag_mag
  lambda_mat[3,2] = -offdiag_mag
  
  mod1 = MVN_CRP_sampler_UVV(
    S = 12000, seed = seed, y = y,
    alpha = 1, r = 10, g = 1, h = 10, nu = 4, 
    fix_r = FALSE,
    mu0 = matrix(rep(0, nrow(y[[1]])), ncol = 1),
    lambda0 = lambda_mat,
    k_init = 1, diag_weights = FALSE,
    verbose = TRUE, split_merge = FALSE)
  
  mod1_sum = dpmm_summary(output = mod1,
                          print_phi_sum = TRUE,
                          print_k_sum = TRUE,
                          make_traceplot = FALSE,
                          burn_in = 2000, t_hold = 250,
                          num_dims = 3,
                          calc_perf = FALSE,
                          equal_var = FALSE)
  
} else if(mod == "UVV_withSM"){
  
  lambda_mat = diag(x = 5, nrow = nrow(y[[1]]))
  offdiag_mag = 0.2*5
  lambda_mat[1,2] = offdiag_mag
  lambda_mat[1,3] = -offdiag_mag
  lambda_mat[2,1] = offdiag_mag
  lambda_mat[2,3] = -offdiag_mag
  lambda_mat[3,1] = -offdiag_mag
  lambda_mat[3,2] = -offdiag_mag
  
  mod1 = MVN_CRP_sampler_UVV(
    S = 12000, seed = seed, y = y,
    alpha = 1, r = 10, g = 1, h = 10, nu = 4, 
    fix_r = FALSE,
    mu0 = matrix(rep(0, nrow(y[[1]])), ncol = 1),
    lambda0 = lambda_mat,
    k_init = 1, diag_weights = FALSE,
    verbose = TRUE, split_merge = TRUE)
  
  mod1_sum = dpmm_summary(output = mod1,
                          print_phi_sum = TRUE,
                          print_k_sum = TRUE,
                          make_traceplot = FALSE,
                          burn_in = 2000, t_hold = 250,
                          num_dims = 3,
                          calc_perf = FALSE,
                          equal_var = FALSE)
  
  
}

saveRDS(object = mod1, file = paste0("./MODFIT_", "_", mod, "_", SLURM_ID))
saveRDS(object = mod1_sum, file = paste0("./MODSUM_", "_", mod, "_", SLURM_ID))