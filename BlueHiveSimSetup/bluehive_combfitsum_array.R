##############################################################
############## BLUEHIVE SIM COMBINED SCRIPT ##################
######################  DPM MODELS ###########################
##############################################################

## Author: Jonathan Klus
## Date: 24 February 2024
## Description: Call to functions to fit and then summarize simulation study results
## run on the BlueHive server. Does not save full model fit information in the
## interest of conserving hard disk space.

###################### load packages and set seed ##############################
source("./posterior_helper_fxns.R")
source("./post_processing_inf.R")

# DEFINE INPUTS --- USER DEFINED IN R SCRIPT
model = c("conjDEV", "conjDEE", "conjUVV", "nonconjDEV", "nonconjUVV")[1]
scenario = c("3close", "3wellsep", "5grp3d")[3]
SM = TRUE # do split merge
cat("\n", model, scenario, SM, "\n")

# sbatch-fed settings
n_array = c(30,100,300)
i = as.numeric(Sys.getenv("i"))
n = n_array[i]
cat("\n n=", n, "\n")

dir_name = paste0("./SummaryProposal/MODSUM_", model, "_", scenario, "_n", n, 
                  ifelse(SM == TRUE, "_withSM", "_noSM"),"_sim_results_",
                  stringr::str_replace_all(string = Sys.Date(),
                                           pattern = "-",
                                           replacement = "_"))

if(dir.exists(dir_name) == FALSE){
  dir.create(dir_name)
} 

# extract seed
seeds = readRDS("./BHsimseeds.rds") # file with 100 random seeds

SLURM_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(SLURM_ID)

# simulate data
set.seed(seeds[SLURM_ID])

if(scenario == "3close"){
  
  w = c(.35, .25, .4)
  means = list(
    #c(10, 10),
    c(0, -5),
    c(0, 8),
    c(10,0)
  )
  
  var = diag(5, length(means[[1]])) # variances diagonal, and equal
  
  # set hyperparameters
  a = 1; b = 50
  d = 1; f = 1
  g = 1; h = 50
  
} else if(scenario == "3wellsep"){
  
  w = c(0.4, 0.3, 0.3)
  
  means = list(
    c(-20, 20),
    c(20, -20),
    c(0, 0)
  )
  
  var = diag(10, length(means[[1]])) # variances diagonal, and equal
  
  # set hyperparameters
  a = 1; b = 50
  d = 1; f = 1
  g = 1; h = 50
  
} else if(scenario == "5grp3d"){
  
  w = rep(0.2, 5)
  
  means = list(
    c(4, 4, -2),
    c(2, 3, 2),
    c(3, 2, 5),
    c(8, 8, 0),
    c(9, 7, 3)
  )
  
  var = diag(0.25, length(means[[1]]))
  
  # set hyperparameters
  a = 1; b = 1
  d = 1; f = 1
  g = 1; h = 10
  
}


assign = sample(x = 1:length(means), size = n, replace = TRUE, prob = w)
y = lapply(X = assign,
           FUN = function(x){
             t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
           })


# fit model
p = nrow(y[[1]])


  
if(model == "conjDEE"){
  
  ################################# DEE #########################################
  source("./Multivariate_DPMM_unknownvar_DEE.R")
  
  try(expr = {
    
    output = MVN_CRP_sampler_DEE(
      S = 12000, seed = seeds[SLURM_ID], y = y,
      alpha = 1, a = a, b = b, r = 10, d = d, f = f, g = g, h = h,
      sigma_hyperprior = FALSE, fix_r = FALSE,
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      k_init = 1, diag_weights = FALSE,
      verbose = FALSE, split_merge = SM, sm_iter = 5)
    
  },
  
  outFile = stdout()
  
  )
  
} else if(model == "conjDEV"){
  
  ################################# DEV #########################################
  source("./Multivariate_DPMM_unknownvar_DEV.R")
  
  try(expr = {
    
    output = MVN_CRP_sampler_DEV(
      S = 12000, seed = seeds[SLURM_ID], y = y,
      alpha = 1, a = a, b = b, r = 10, d = d, f = f, g = g, h = h,
      sigma_hyperprior = FALSE, fix_r = FALSE,
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      k_init = 1, diag_weights = FALSE,
      verbose = FALSE, split_merge = SM, sm_iter = 5)
    
  },
  
  outFile = stdout()
  
  )
  
} else if(model == "conjUVV"){
  
  ################################# UVV #########################################
  source("./Multivariate_DPMM_unknownvar_UVV.R")
  
  # PRIOR SETTINGS
  ## FOR CLOSE TOGETHER: 
  ## FOR WELL SEPARATED: lambda0 = diag(15), nu = 2 
  
  try(expr = {
    
    output = MVN_CRP_sampler_UVV(
      S = 12000, seed = seeds[SLURM_ID], y = y,
      alpha = 1, r = 10, g = g, h = h, nu = 2, fix_r = FALSE,
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      lambda0 = diag(x = 15, nrow = p),
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      k_init = 1, diag_weights = FALSE,
      verbose = FALSE, split_merge = SM, sm_iter = 5)
  },
  
  outFile = stdout()
  
  )
  
} else if(model == "nonconjDEV"){
  
  ############################ nonconj DEV #####################################
  source("./Multivariate_DPMM_nonconj_DEV.R")
  
  try(expr = {
    
    output = MVN_CRP_nonconj_DEV(
      S = 12000, seed = seeds[SLURM_ID], y = y, alpha = 1, 
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      sigma0 = 25,
      a = a, b = b, sigma_hyperprior = FALSE,
      k_init = 1, diag_weights = FALSE,
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      verbose = FALSE, split_merge = SM)
  },
  
  
  outFile = stdout()
  
  )
  
} else if(model == "nonconjUVV"){
  
  ############################ nonconj UVV #####################################
  source("./Multivariate_DPMM_nonconj_UVV.R")
  
  try(expr = {
    
    output = MVN_CRP_nonconj_UVV(
      S = 12000, seed = seeds[SLURM_ID], y = y, alpha = 1, 
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      Sigma0 = diag(25, p), Lambda0 = diag(10,p), nu = 2,
      k_init = 1, diag_weights = FALSE,
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      verbose = FALSE, split_merge = SM)
  },
  
  
  outFile = stdout()
  
  )
  
}



########################### summarize results ##################################

# for DEE model - pooled variance

if(model == "conjDEE"){
  
  try(expr = {
    
    mod_sum = dpmm_summary(output = output, 
                           print_phi_sum = TRUE,
                           print_k_sum = TRUE, 
                           make_traceplot = FALSE,
                           burn_in = 2000, t_hold = 250, 
                           num_dims = p, 
                           calc_perf = TRUE, 
                           mu_true = output$truth$mu_true, 
                           var_true = output$truth$var_true, 
                           assign_true = output$truth$assign_true, 
                           equal_var = TRUE)
    
  },
  
  outFile = stdout()
  
  )
  
} else{
  
  try(expr = {
    
    mod_sum = dpmm_summary(output = output,
                           print_phi_sum = TRUE,
                           print_k_sum = TRUE,
                           make_traceplot = FALSE,
                           burn_in = 2000, t_hold = 250,
                           num_dims = p,
                           calc_perf = TRUE,
                           mu_true = output$truth$mu_true,
                           var_true = lapply(X = 1:length(output$truth$mu_true),
                                             FUN = function(x){output$truth$var_true}),
                           assign_true = output$truth$assign_true,
                           equal_var = FALSE)
    
  },
  
  outFile = stdout()
  
  )
  
}


# save summary
saveRDS(object = mod_sum, 
        file = paste0(dir_name,"/sum_", SLURM_ID,".rds"))
