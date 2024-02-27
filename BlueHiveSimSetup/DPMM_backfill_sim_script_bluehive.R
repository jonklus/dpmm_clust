################################################################################
################## SIMULATION STUDY SCRIPT #####################################
################## FOR BLUEHIVE SIMS       #####################################
################## CLOSE TOGETHER          #####################################
################################################################################
###### Jonathan Klus
###### 16 February 2024
################################################################################

source("./Multivariate_DPMM_unknownvar_DEV.R")
# source("./Multivariate_DPMM_unknownvar_DEE.R")
# source("./Multivariate_DPMM_unknownvar_UVV.R")
source("./posterior_helper_fxns.R")

# sim settings
# n_array = c(30,100,300)
# i = as.numeric(Sys.getenv("i"))
# n = n_array[i]
n=300
cat("\n n=", n, "\n")

# store results
dir_name = paste0("./conjDEV_close_n", n, "_withSM_sim_results_",
                  "2024_02_20")
                  # stringr::str_replace_all(string = Sys.Date(),
                  #                          pattern = "-",
                  #                          replacement = "_"))
if(dir.exists(dir_name) == FALSE){
  dir.create(dir_name)
}

# dir2_name = paste0("./conjUVV_close_n", n, "_withSM_sim_results_",
#                    stringr::str_replace_all(string = Sys.Date(), 
#                                             pattern = "-", 
#                                             replacement = "_"))
# if(dir.exists(dir2_name) == FALSE){
#   dir.create(dir2_name)
# } 


# extract seed
seeds = readRDS("./BHsimseeds.rds") # file with 100 random seeds
sim_num = c(1,23,37,38,65)
i = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
SLURM_ID = sim_num[i]
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

# try(expr = {
#   
#   mod1 = MVN_CRP_sampler_DEE(
#     S = 12000, seed = seeds[SLURM_ID], y = y,
#     alpha = 1, r = 10, g = 1, h = 50,
#     sigma_hyperprior = FALSE, fix_r = FALSE,
#     mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
#     a = 1, b = 50,
#     truth = list(mu_true = means, var_true = var, assign_true = assign),
#     k_init = 1, diag_weights = FALSE,
#     verbose = FALSE, split_merge = FALSE)
#   
#   filename = paste0(dir_name, "/output_", SLURM_ID, ".rds")
#   saveRDS(object = mod1, file = filename)
#   
# },
# 
# outFile = stdout()
# 
# )


################################# UVV #########################################

# try(expr = {
#   
#   mod1 = MVN_CRP_sampler_UVV(
#     S = 12000, seed = seeds[SLURM_ID], y = y,
#     alpha = 1, r = 10, g = 1, h = 50, nu = 2, fix_r = FALSE,
#     mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
#     lambda0 = diag(x = 15, nrow = 2),
#     truth = list(mu_true = means, var_true = var, assign_true = assign),
#     k_init = 1, diag_weights = FALSE,
#     verbose = FALSE, split_merge = TRUE, sm_iter = 5
#   )
#   
#   filename = paste0(dir_name, "/output_", SLURM_ID, ".rds")
#   saveRDS(object = mod1, file = filename)
#   
# },
# 
# outFile = stdout()
# 
# )

################################# DEV #########################################

### no SM

# try(expr = {
#   
#   mod1 = MVN_CRP_sampler_DEV(
#     S = 12000, seed = seeds[SLURM_ID], y = y,
#     alpha = 1, r = 10, g = 1, h = 50,
#     sigma_hyperprior = FALSE, fix_r = FALSE,
#     mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
#     a = 1, b = 50,
#     truth = list(mu_true = means, var_true = var, assign_true = assign),
#     k_init = 1, diag_weights = FALSE,
#     verbose = FALSE, split_merge = FALSE)
#   
#   filename = paste0(dir_name, "/output_", SLURM_ID, ".rds")
#   saveRDS(object = mod1, file = filename)
#   
# },
# 
# outFile = stdout()
# 
# )
# 
# rm(mod1) # get rid of mod1 in memory to make room for next stage

### with SM

try(expr = {

  mod1 = MVN_CRP_sampler_DEV(
    S = 12000, seed = seeds[SLURM_ID], y = y,
    alpha = 1, r = 10, g = 1, h = 50,
    sigma_hyperprior = FALSE, fix_r = FALSE,
    mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = 2))),0), ncol = 1),
    a = 1, b = 50,
    truth = list(mu_true = means, var_true = var, assign_true = assign),
    k_init = 1, diag_weights = FALSE,
    verbose = FALSE, split_merge = TRUE, sm_iter = 5)

  filename = paste0(dir_name, "/output_", SLURM_ID, ".rds")
  saveRDS(object = mod1, file = filename)

},

outFile = stdout()

)