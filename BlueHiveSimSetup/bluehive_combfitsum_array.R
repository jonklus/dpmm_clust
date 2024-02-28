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
# source("./Multivariate_DPMM_unknownvar_DEV.R")
source("./Multivariate_DPMM_unknownvar_DEE.R")
# source("./Multivariate_DPMM_unknownvar_UVV.R")
source("./posterior_helper_fxns.R")
source("./post_processing_inf.R")

# DEFINE INPUTS --- USER DEFINED IN R SCRIPT
# model = "DEV" # "DEE" "UVV"
# scenario = "close"  # "wellsep"
# SM = TRUE # do split merge

# sbatch-fed settings
n_array = c(30,100,300)
i = as.numeric(Sys.getenv("i"))
n = n_array[i]
cat("\n n=", n, "\n")

SLURM_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

dir_name = paste0("./Summary/MODSUM_conjDEE_close_n", n, "_withSM_sim_results_",
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

w = c(.35, .25, .4)
means = list(
  #c(10, 10),
  c(0, -5),
  c(0, 8),
  c(10,0)
)

var = diag(5, length(means[[1]])) # variances diagonal, and equal

assign = sample(x = 1:length(means), size = n, replace = TRUE, prob = w)
y = lapply(X = assign,
           FUN = function(x){
             t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
           })


# fit model

################################# DEE #########################################

try(expr = {

  output = MVN_CRP_sampler_DEE(
    S = 12000, seed = seeds[SLURM_ID], y = y,
    alpha = 1, r = 10, g = 1, h = 50,
    sigma_hyperprior = FALSE, fix_r = FALSE,
    mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
    a = 1, b = 50,
    truth = list(mu_true = means, var_true = var, assign_true = assign),
    k_init = 1, diag_weights = FALSE,
    verbose = FALSE, split_merge = TRUE, sm_iter = 5)

},

outFile = stdout()

)


################################# UVV #########################################

# try(expr = {
#   
#   output = MVN_CRP_sampler_UVV(
#     S = 12000, seed = seeds[SLURM_ID], y = y,
#     alpha = 1, r = 10, g = 1, h = 50, nu = 2, fix_r = FALSE,
#     mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
#     lambda0 = diag(x = 15, nrow = 2),
#     truth = list(mu_true = means, var_true = var, assign_true = assign),
#     k_init = 1, diag_weights = FALSE,
#     verbose = FALSE, split_merge = TRUE, sm_iter = 5)
# },
# 
# outFile = stdout()
# 
# )

################################# DEV #########################################

# try(expr = {
#   
#   output = MVN_CRP_sampler_DEV(
#     S = 12000, seed = seeds[SLURM_ID], y = y,
#     alpha = 1, r = 10, g = 1, h = 50,
#     sigma_hyperprior = FALSE, fix_r = FALSE,
#     mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
#     a = 1, b = 50,
#     truth = list(mu_true = means, var_true = var, assign_true = assign),
#     k_init = 1, diag_weights = FALSE,
#     verbose = FALSE, split_merge = TRUE, sm_iter = 5)
#   
# },
# 
# outFile = stdout()
# 
# )
# 
# rm(output) # get rid of output in memory to make room for next stage

########################### summarize results ##################################

# for DEE model - pooled variance

try(expr = {
  
  mod_sum = dpmm_summary(output = output, 
                         print_phi_sum = TRUE,
                         print_k_sum = TRUE, 
                         make_traceplot = FALSE,
                         burn_in = 2000, t_hold = 250, 
                         num_dims = 2, 
                         calc_perf = TRUE, 
                         mu_true = output$truth$mu_true, 
                         var_true = output$truth$var_true, 
                         assign_true = output$truth$assign_true, 
                         equal_var = TRUE)
  
},

outFile = stdout()

)

# for all other models

# try(expr = {
#   
#   mod_sum = dpmm_summary(output = output, 
#                 print_phi_sum = TRUE,
#                 print_k_sum = TRUE, 
#                 make_traceplot = FALSE,
#                 burn_in = 2000, t_hold = 250, 
#                 num_dims = 2, 
#                 calc_perf = TRUE, 
#                 mu_true = output$truth$mu_true, 
#                 var_true = lapply(X = 1:length(output$truth$mu_true), 
#                                   FUN = function(x){output$truth$var_true}), 
#                 assign_true = output$truth$assign_true, 
#                 equal_var = FALSE)
# 
# },
# 
# outFile = stdout()
# 
# )


# save summary
saveRDS(object = mod_sum, 
        file = paste0(dir_name,"/sum_",  SLURM_ID,".rds"))
