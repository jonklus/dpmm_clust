################################################################################
############ POST PROCESSING AND FORMAL INFERENCE FUNCTIONS ####################
##################               VERSION 1            ##########################
################################################################################
## AUTHOR: JONATHAN KLUS
## DATE: 10 JANUARY 2024
## DESCRIPTION: CONVENIENCE FUNCTIONS FOR FINDING POSTERIOR SUMMARIES FOR INFERENCE
## AND DIAGNOSTICS FROM MCMC DRAWS IN A SIMULATION STUDY WITH MULTIPLE DATA SETS. 
## CAN PROVIDE OVERALL SUMMARY OR BY DATA SET. 

####################### LOAD ANY REQUIRED PACKAGES #############################

# homemade functions
source("./posterior_helper_fxns.R")

# load R libraries
library(ggplot2)
library(gridExtra)
library(plotly)
library(label.switching)
library(LaplacesDemon)
library(parallel)
library(stringr)


######################## DEFINE NEW FUNCTIONS ##################################

dpmm_summary <- function(output, dataset_ind = 1, print_phi_sum = FALSE,
                         print_k_sum = TRUE, make_traceplot = TRUE,
                         burn_in = 1000, t_hold = 0, num_dims = 2, equal_var = FALSE
                         ){
  # output is the list of results from the DPMM simulation study function
  # dataset is a numeric argument to summarize a specific result in the output, the desired index
  # sum_all is a logical argument to provide a summary of all data sets --- NOT CURRENTLY IMPLEMENTED
  # print_k_sum is a logical, if TRUE print summary of no groups found 
  # print_phi_sum is a logical, if TRUE print summary of estimated model parameters
  # burn_in is # of burn in iterations to discard
  # t_hold is the threshold # of iterations for a given k in order to report results
  # num_dims is the dimensionality of the problem (i.e. a bivariate normal is dim 2)
  
  # show basic summary
  if(print_k_sum == TRUE){
    cat("\n Frequency of MCMC iterations finding K groups:")
    print(table(output[[dataset_ind]]$k))
    
    cat("\n Percentage of MCMC iterations finding K groups:")
    print(round((table(output[[dataset_ind]]$k)/sum(table(output[[dataset_ind]]$k)))*100,1))
    
    cat("\n *Note that above frequency summaries of MCMC iterations were made before burn-in or thresholds were applied. 
  All inference on phi will be made after accounting for burn-in and thresholding. \n")
  }
  
  # filter by number of iterations for each k and address label switching
    prob_list_by_k = get_probs_by_k(probs = output[[dataset_ind]]$group_probs, 
                                    n_groups = output[[dataset_ind]]$k, 
                                    burn_in = burn_in, 
                                    iter_threshold = t_hold)
    group_assign_list_by_k = get_assign_by_k(assign = output[[dataset_ind]]$group_assign, 
                                             n_groups = output[[dataset_ind]]$k, 
                                             burn_in = burn_in, 
                                             iter_threshold = t_hold)
    
    # correct label switching 
    stephens_result = get_stephens_result(group_assign_list_by_k = group_assign_list_by_k, 
                                          prob_list_by_k = prob_list_by_k$prob_list)
    
    # summarize means & variances
    mean_list_by_k_stephens = list_params_by_k(draws = output[[dataset_ind]]$means, 
                                               # k_vec = output[[1]]$k,
                                               # burn_in = burn_in, 
                                               # iter_threshold = thold,
                                               iter_list = prob_list_by_k$iter_list,
                                               relabel = TRUE,
                                               permutation = stephens_result, 
                                               param_type = "Mean")
    
    var_list_by_k_stephens = list_params_by_k(draws = output[[dataset_ind]]$vars, 
                                              iter_list = prob_list_by_k$iter_list,
                                              relabel = TRUE, equal_var = equal_var,
                                              permutation = stephens_result,
                                              param_type = "Var")
    
    # compute KL divergence
    group_assign_list_by_k_corr = correct_group_assign(
      group_assign_list_by_k = group_assign_list_by_k, 
      stephens_result = stephens_result)
    
    kl_res = calc_KL_diverg(y = output[[dataset_ind]]$data,
                            mu_est = mean_list_by_k_stephens,
                            Sigma_est = var_list_by_k_stephens,
                            group_assign = group_assign_list_by_k_corr,
                            true_assign = yreps[[dataset_ind]]$assign,
                            mu_true = mu_true,
                            Sigma_true = var_true,
                            equal_var_assump = TRUE)
    # py is truth, px is estimate
    kl_div = kl_res$sum.KLD.py.px # how far is estimate px from truth py
    
    
    mean_summary = vector(mode = "list", length = length(mean_list_by_k_stephens))
    var_summary = vector(mode = "list", length = length(mean_list_by_k_stephens))
    for(k in 1:length(mean_list_by_k_stephens)){
      
      # make mean summary table
      mean_summary[[k]] = make_postsum(mcmc_df = mean_list_by_k_stephens[[k]], digits = 2)
      
      # make variance summary table
      var_summary[[k]] = make_postsum(mcmc_df = var_list_by_k_stephens[[k]], digits = 2)
      
      k_i = ncol(var_list_by_k_stephens[[k]])
      if(print_phi_sum == TRUE){
        # give summary of counts after thresholding
        cat("\n K =", k_i, " n_k =", nrow(mean_list_by_k_stephens[[k]]), "after burn-in and thresholding\n")
        print(mean_summary[[k]])
        print(var_summary[[k]])
      }
      
      if(make_traceplot == TRUE){
        for(dim_i in 1:num_dims){
          make_traceplot(param_list_by_k = mean_list_by_k_stephens, 
                         k = k_i, 
                         component_no = dim_i,
                         p= num_dims,
                         # title_note = cat("K=", k_i),
                         param_type = "Mean")
        }
      }
    }
    
    # KL divergence for entire model --- across all k
    if(print_k_sum == TRUE){
      cat("\n KL Divergence: KL(p_est||p_true)=", round(kl_div, 3), "\n")
    }
    
  # return summary of all results
  return(list(
    mean_list_by_k_stephens = mean_list_by_k_stephens,
    var_list_by_k_stephens = var_list_by_k_stephens,
    mean_summary = mean_summary,
    var_summary = var_summary,
    kl_div = round(kl_div, 4)
  ))
  
}

# 
# # test function
# output = readRDS("../MCMC_Runs/conjDEVsamp_minisimstudy_close_withSM_2024_01_16.rds")
# # test = dpmm_summary(output = output, dataset = 1, 
# #                       print_result = TRUE, make_traceplot = TRUE,
# #                       burn_in = 1000, t_hold = 100, num_dims = 2)
# # 
# output2 = readRDS("../MCMC_Runs/conjDEVsamp_minisimstudy_close_noSM_2024_01_10.rds")
# 
# dataset_ind = 1 
# print_result = TRUE
# burn_in = 1000
# t_hold = 100
# num_dims = 2
# 
# prob_list_by_k = get_probs_by_k(probs = output[[dataset_ind]]$group_probs, 
#                                 n_groups = output[[dataset_ind]]$k, 
#                                 burn_in = burn_in, 
#                                 iter_threshold = t_hold)
# group_assign_list_by_k = get_assign_by_k(assign = output[[dataset_ind]]$group_assign, 
#                                          n_groups = output[[dataset_ind]]$k, 
#                                          burn_in = burn_in, 
#                                          iter_threshold = t_hold)
# 
# # correct label switching 
# stephens_result = get_stephens_result(group_assign_list_by_k = group_assign_list_by_k, 
#                                       prob_list_by_k = prob_list_by_k$prob_list)
# 
# # summarize means & variances
# mean_list_by_k_stephens = list_params_by_k(draws = output[[dataset_ind]]$means, 
#                                            # k_vec = output[[1]]$k,
#                                            # burn_in = burn_in, 
#                                            # iter_threshold = thold,
#                                            iter_list = prob_list_by_k$iter_list,
#                                            relabel = TRUE,
#                                            permutation = stephens_result, 
#                                            param_type = "Mean")
# 
# var_list_by_k_stephens = list_params_by_k(draws = output[[dataset_ind]]$vars, 
#                                           iter_list = prob_list_by_k$iter_list,
#                                           relabel = TRUE,
#                                           permutation = stephens_result,
#                                           param_type = "Var")
# 
# calc_KL_diverg(y = output[[1]]$data, 
#                mu_est = mean_list_by_k_stephens, 
#                Sigma_est = var_list_by_k_stephens, 
#                group_assign = stephens_result, 
#                true_assign = yreps[[1]]$assign, 
#                mu_true = mu_true, 
#                Sigma_true = var_true, 
#                equal_var_assump = TRUE)