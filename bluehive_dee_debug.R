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
model = c("conjDEV", "conjDEE", "conjUVV", "nonconjDEV", "nonconjUVV")[2]
scenario = c("3close", "3wellsep", "5grp3d")[1]
SM = TRUE # do split merge
cat("\n", model, scenario, SM, "\n")

# sbatch-fed settings
n_array = c(30,100) #,300)
i = 1 # as.numeric(Sys.getenv("i"))
n = n_array[i]
cat("\n n=", n, "\n")

dir_name = paste0("../MCMC_Runs/SummaryProposal/MODSUM_", model, "_", scenario, "_n", n, 
                  ifelse(SM == TRUE, "_withSM", "_noSM"),"_sim_results_",
                  stringr::str_replace_all(string = Sys.Date(),
                                           pattern = "-",
                                           replacement = "_"))

if(dir.exists(dir_name) == FALSE){
  dir.create(dir_name)
} 

# extract seed
seeds = readRDS("./BHsimseeds.rds") # file with 100 random seeds

# other settings
off_diag = FALSE # unless reset to TRUE below

SLURM_ID = 1 # as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
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
  a = 1; b = 25
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
  sigma0 = 200
  lambda0 = diag(x = 25, nrow = length(means[[1]])) # needed for UVV
  nu = length(means[[1]])  # must have nu >= p...Otherwise rinvwish will fail!
  
  
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
  sigma0 = 100
  lambda0 = diag(x = 1, nrow = length(means[[1]])) # needed for UVV
  nu = length(means[[1]]) # must have nu >= p!! Otherwise rinvwish will fail!
  
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
      alpha = 1, r = 10, g = g, h = h, nu = nu, fix_r = FALSE,
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      lambda0 = lambda0,
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      k_init = 1, diag_weights = FALSE,
      verbose = FALSE, split_merge = SM, sm_iter = 5)
  },
  
  outFile = stdout()
  
  )
  
  off_diag = TRUE
  
} else if(model == "nonconjDEV"){
  
  ############################ nonconj DEV #####################################
  source("./Multivariate_DPMM_nonconj_DEV.R")
  
  try(expr = {
    
    output = MVN_CRP_nonconj_DEV(
      S = 12000, seed = seeds[SLURM_ID], y = y, alpha = 1, 
      mu0 = matrix(round((colMeans(matrix(unlist(y), ncol = p))),0), ncol = 1),
      sigma0 = sigma0,
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
      Sigma0 = diag(sigma0, p), Lambda0 = lambda0, nu = nu,
      k_init = 1, diag_weights = FALSE,
      truth = list(mu_true = means, var_true = var, assign_true = assign),
      verbose = FALSE, split_merge = SM)
  },
  
  
  outFile = stdout()
  
  )
  
  off_diag = TRUE
  
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
                           equal_var = FALSE, off_diag = off_diag)
    
  },
  
  outFile = stdout()
  
  )
  
}


# save summary
saveRDS(object = mod_sum, 
        file = paste0(dir_name,"/sum_", SLURM_ID,".rds"))




# ### doing post sum fxn step by step to debug
# 
# ######################## DEFINE NEW FUNCTIONS ##################################
# 
# 
# print_phi_sum = TRUE
# print_k_sum = TRUE
# make_traceplot = FALSE
# burn_in = 2000
# t_hold = 250
# num_dims = p
# calc_perf = TRUE
# mu_true = output$truth$mu_true
# var_true = output$truth$var_true
# assign_true = output$truth$assign_true
# equal_var = TRUE
# off_diag = FALSE
# 
#   # output is the list of results from the DPMM simulation study function
#   # dataset is a numeric argument to summarize a specific result in the output, the desired index
#   # sum_all is a logical argument to provide a summary of all data sets --- NOT CURRENTLY IMPLEMENTED
#   # print_k_sum is a logical, if TRUE print summary of no groups found 
#   # print_phi_sum is a logical, if TRUE print summary of estimated model parameters
#   # burn_in is # of burn in iterations to discard
#   # t_hold is the threshold # of iterations for a given k in order to report results
#   # num_dims is the dimensionality of the problem (i.e. a bivariate normal is dim 2)
#   # calc_perf is a logical - whether to calculate performance metrics (currently KLD and ARI)
#   # mu_true and var_true are list arguments with the true values of the model parameters
#   # from a simulation study used to calculate the KL divergence 
#   # equal_var is a logical argument for whether the equal variance assumption was made
#   # in the model. The function will then expect a scalar variance instead of a var-covar matrix 
#   
#   start = Sys.time()
#   
#   # show basic summary
#   k_freqtab = table(output$k)
#   k_relfreqtab = round((table(output$k)/sum(table(output$k)))*100,1)
#   
#   if(print_k_sum == TRUE){
#     cat("\n Raw summary: \n")
#     
#     cat("\n Frequency of MCMC iterations finding k groups:")
#     print(k_freqtab)
#     
#     cat("\n Percentage of MCMC iterations finding k groups:")
#     print(k_relfreqtab)
#     
#     cat("\n *Note that above frequency summaries of MCMC iterations were made before burn-in or thresholds were applied. 
#           All inference on phi will be made after accounting for burn-in and thresholding. \n")
#   }
#   
#   # filter by number of iterations for each k and address label switching
#   prob_list_by_k = get_probs_by_k(probs = output$group_probs, 
#                                   n_groups = output$k, 
#                                   burn_in = burn_in, 
#                                   iter_threshold = t_hold)
#   
#   group_assign_list_by_k = get_assign_by_k(assign = output$group_assign, 
#                                            n_groups = output$k, 
#                                            burn_in = burn_in, 
#                                            iter_threshold = t_hold)
#   
#   # correct label switching 
#   stephens_result = get_stephens_result(group_assign_list_by_k = group_assign_list_by_k, 
#                                         prob_list_by_k = prob_list_by_k$prob_list)
#   
#   # summary table by k after burn-in
#   if(print_k_sum == TRUE){
#     cat("\n Summary aftern burn-in and thresholding: \n")
#     
#     num_groups = unlist(sapply(X = 1:length(group_assign_list_by_k),
#                         FUN = function(x){
#                           k = x # list entry
#                           sapply(X = 1:nrow(group_assign_list_by_k[[x]]), 
#                                  FUN = function(x){
#                                    length(unique(group_assign_list_by_k[[k]][x,]))
#                                  })
#                         }))
#     
#     
#     adj_k_freqtab = table(num_groups)
#     adj_k_relfreqtab = round((table(num_groups)/sum(table(num_groups)))*100,1)
#     
#     cat("\n Frequency of MCMC iterations finding k groups after burn-in:")
#     print(adj_k_freqtab)
#     
#     cat("\n Percentage of MCMC iterations finding k groups after burn-in:")
#     print(adj_k_relfreqtab)
#     
#   }
#   
#   # summarize means & variances
#   mean_list_by_k_stephens = list_params_by_k(draws = output$means, 
#                                              k_vec = output$k,
#                                              # burn_in = burn_in, 
#                                              # iter_threshold = thold,
#                                              iter_list = prob_list_by_k$iter_list,
#                                              relabel = TRUE,
#                                              permutation = stephens_result, 
#                                              param_type = "Mean")
#   
#   if(off_diag == FALSE){
#     
#     var_list_by_k_stephens = list_params_by_k(draws = output$vars, 
#                                               iter_list = prob_list_by_k$iter_list,
#                                               k_vec = output$k,
#                                               relabel = TRUE, equal_var = equal_var,
#                                               permutation = stephens_result,
#                                               param_type = "Var")
# 
#     
#     
#   } else{
#     
#     var_list_by_k_stephens = list_params_by_k(draws = output$vars, 
#                                               iter_list = prob_list_by_k$iter_list,
#                                               k_vec = output$k,
#                                               relabel = TRUE, equal_var = equal_var,
#                                               permutation = stephens_result,
#                                               param_type = "Covar")
#   }
#   
#   
#   # compute KL divergence
#   group_assign_list_by_k_corr = correct_group_assign(
#     group_assign_list_by_k = group_assign_list_by_k, 
#     stephens_result = stephens_result)
#   
#   if(calc_perf == TRUE){
#     
#     # KL divergence
#     kl_res = calc_KL_diverg(y = output$data,
#                             mu_est = mean_list_by_k_stephens,
#                             Sigma_est = var_list_by_k_stephens,
#                             group_assign = group_assign_list_by_k_corr,
#                             true_assign = assign_true,
#                             mu_true = mu_true,
#                             Sigma_true = var_true,
#                             equal_var_assump = equal_var,
#                             off_diag = off_diag)
#     
#     # py is truth, px is estimate
#     kl_div = kl_res$sum.KLD.py.px # how far is estimate px from truth py
#     
#     # Adjusted RAND index
#     ARI = lapply(X = 1:length(group_assign_list_by_k), 
#                  FUN = function(x){
#                    k = x 
#                    sapply(X = 1:nrow(group_assign_list_by_k[[k]]), 
#                           FUN = function(x){
#                             mclust::adjustedRandIndex(
#                               x = assign_true, 
#                               y = group_assign_list_by_k[[k]][x,])
#                           })
#                  })
#     
#     
#   }
#   
#   mean_summary = vector(mode = "list", length = length(mean_list_by_k_stephens))
#   var_summary = vector(mode = "list", length = length(mean_list_by_k_stephens))
#   # total_iter = sapply(X = 1:length(mean_list_by_k_stephens), 
#   #                     FUN = function(x){
#   #                       
#   #                     })
#   for(k in 1:length(mean_list_by_k_stephens)){
#     
#     # make mean summary table
#     mean_summary[[k]] = make_postsum(mcmc_df = mean_list_by_k_stephens[[k]], digits = 2)
#     
#     # make variance summary table
#     var_summary[[k]] = make_postsum(mcmc_df = var_list_by_k_stephens[[k]], digits = 2)
#     
#     k_i = max(as.numeric(stringr::str_extract_all(
#       string = stringr::str_extract_all(string = row.names(mean_summary[[k]]), 
#                                         pattern = "_[:digit:]_"), pattern = "[:digit:]")))
#     if(print_phi_sum == TRUE){
#       # give summary of counts after thresholding
#       cat("\n k =", k_i, " n_iter =", nrow(mean_list_by_k_stephens[[k]]), "after burn-in and thresholding\n")
#       
#       cat("\n Mean Summary: \n")
#       print(mean_summary[[k]])
#       cat("\n (Co)variance Summary: \n")
#       print(var_summary[[k]])
#     }
#     
#     if(make_traceplot == TRUE){
#       for(dim_i in 1:num_dims){
#         make_traceplot(param_list_by_k = mean_list_by_k_stephens, 
#                        k = k_i, 
#                        component_no = dim_i,
#                        p= num_dims,
#                        # title_note = cat("K=", k_i),
#                        param_type = "Mean")
#       }
#     }
#   }
