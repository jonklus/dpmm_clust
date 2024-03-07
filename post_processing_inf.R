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
library(label.switching)
library(LaplacesDemon)
library(parallel)
library(stringr)
library(dplyr)
library(mclust)


######################## DEFINE NEW FUNCTIONS ##################################

dpmm_summary <- function(output, print_phi_sum = FALSE,
                         print_k_sum = TRUE, make_traceplot = TRUE,
                         burn_in = 1000, t_hold = 0, num_dims = 2, 
                         calc_perf = FALSE, mu_true = NULL, var_true = NULL, 
                         assign_true = NULL, equal_var = FALSE
                         ){
  # output is the list of results from the DPMM simulation study function
  # dataset is a numeric argument to summarize a specific result in the output, the desired index
  # sum_all is a logical argument to provide a summary of all data sets --- NOT CURRENTLY IMPLEMENTED
  # print_k_sum is a logical, if TRUE print summary of no groups found 
  # print_phi_sum is a logical, if TRUE print summary of estimated model parameters
  # burn_in is # of burn in iterations to discard
  # t_hold is the threshold # of iterations for a given k in order to report results
  # num_dims is the dimensionality of the problem (i.e. a bivariate normal is dim 2)
  # calc_perf is a logical - whether to calculate performance metrics (currently KLD and ARI)
  # mu_true and var_true are list arguments with the true values of the model parameters
  # from a simulation study used to calculate the KL divergence 
  # equal_var is a logical argument for whether the equal variance assumption was made
  # in the model. The function will then expect a scalar variance instead of a var-covar matrix 
  
  start = Sys.time()
  
  # show basic summary
  k_freqtab = table(output$k)
  k_relfreqtab = round((table(output$k)/sum(table(output$k)))*100,1)
  
  if(print_k_sum == TRUE){
    cat("\n Frequency of MCMC iterations finding K groups:")
    print(k_freqtab)
    
    cat("\n Percentage of MCMC iterations finding K groups:")
    print(k_relfreqtab)
    
    cat("\n *Note that above frequency summaries of MCMC iterations were made before burn-in or thresholds were applied. 
          All inference on phi will be made after accounting for burn-in and thresholding. \n")
  }
  
  # filter by number of iterations for each k and address label switching
    prob_list_by_k = get_probs_by_k(probs = output$group_probs, 
                                    n_groups = output$k, 
                                    burn_in = burn_in, 
                                    iter_threshold = t_hold)
    group_assign_list_by_k = get_assign_by_k(assign = output$group_assign, 
                                             n_groups = output$k, 
                                             burn_in = burn_in, 
                                             iter_threshold = t_hold)
    
    # correct label switching 
    stephens_result = get_stephens_result(group_assign_list_by_k = group_assign_list_by_k, 
                                          prob_list_by_k = prob_list_by_k$prob_list)
    
    # summarize means & variances
    mean_list_by_k_stephens = list_params_by_k(draws = output$means, 
                                               k_vec = output$k,
                                               # burn_in = burn_in, 
                                               # iter_threshold = thold,
                                               iter_list = prob_list_by_k$iter_list,
                                               relabel = TRUE,
                                               permutation = stephens_result, 
                                               param_type = "Mean")
    
    var_list_by_k_stephens = list_params_by_k(draws = output$vars, 
                                              iter_list = prob_list_by_k$iter_list,
                                              k_vec = output$k,
                                              relabel = TRUE, equal_var = equal_var,
                                              permutation = stephens_result,
                                              param_type = "Var")
    
    # compute KL divergence
    group_assign_list_by_k_corr = correct_group_assign(
      group_assign_list_by_k = group_assign_list_by_k, 
      stephens_result = stephens_result)
    
    if(calc_perf == TRUE){
      
      # KL divergence
      kl_res = calc_KL_diverg(y = output$data,
                              mu_est = mean_list_by_k_stephens,
                              Sigma_est = var_list_by_k_stephens,
                              group_assign = group_assign_list_by_k_corr,
                              true_assign = assign_true,
                              mu_true = mu_true,
                              Sigma_true = var_true,
                              equal_var_assump = equal_var)
      # py is truth, px is estimate
      kl_div = kl_res$sum.KLD.py.px # how far is estimate px from truth py
      
      # Adjusted RAND index
      ARI = lapply(X = 1:length(group_assign_list_by_k), 
                   FUN = function(x){
                     k = x 
                     sapply(X = 1:nrow(group_assign_list_by_k[[k]]), 
                            FUN = function(x){
                              mclust::adjustedRandIndex(
                                x = assign_true, 
                                y = group_assign_list_by_k[[k]][x,])
                            })
                   })

      
    }
    
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
    if(print_k_sum == TRUE & calc_perf == TRUE){
      cat("\n KL Divergence: KL(p_est||p_true) = ", round(kl_div, 3), "\n")
      cat("\n Adj RAND Index = ", round(mean(unlist(ARI)), 3), "\n")
    }
    
    # MH acceptance probs for split merge or non conj MH alg
    # need to check if SM results exists bc not all samplers will do this...
    # also check if return is NA
    if(is.null(output$sm_results) == FALSE){
      # check that sm_results slot in list output exists
      if(nrow(output$sm_results) > 1){
        # check that more than 1 row exists in sm_results...default is a single
        # row of NA when split_merge = FALSE in sampler...don't want to summarize this
        sm_df = data.frame(output$sm_results) %>%
          dplyr::filter(is.na(s) == FALSE) %>%
          dplyr::mutate(
            move_type = factor(move_type),
            s = as.integer(s),
            sm_iter = as.integer(sm_iter),
            accept = as.integer(accept),
            prob = as.numeric(prob)
          ) %>%
          dplyr::group_by(move_type) %>%
          dplyr::summarise(Accept_Prob = round(mean(accept),3), 
                           Count = dplyr::n())
        
        if(print_k_sum == TRUE){
          cat("\n Split/Merge MH Steps: \n")
          print(sm_df)
        }
        
      } else{
        sm_df = NULL
      }
    } else{
      sm_df = NULL
    }

    end = Sys.time()
    cat("\n Summary function runtime is", difftime(end, start, units = "m"), "mins \n")
    
  # return summary of all results
    if(calc_perf == TRUE){
      
      return(list(
        settings = output$settings,
        mean_list_by_k_stephens = mean_list_by_k_stephens,
        var_list_by_k_stephens = var_list_by_k_stephens,
        mean_summary = mean_summary,
        var_summary = var_summary,
        splitmerge_accept = sm_df, 
        kl_div = round(kl_div, 4),
        mean_ARI = mean(unlist(ARI)),
        fit_runtime = output$runtime,
        truth = output$truth,
        avg_pairwise_mat = output$pairwise_mats$avg_adj,
        k_freqtab = k_freqtab,
        k_relfreqtab = k_relfreqtab
      ))
      
    } else{
      
      return(list(
        settings = output$settings,
        mean_list_by_k_stephens = mean_list_by_k_stephens,
        var_list_by_k_stephens = var_list_by_k_stephens,
        mean_summary = mean_summary,
        var_summary = var_summary,
        splitmerge_accept = sm_df,
        fit_runtime = output$runtime,
        avg_pairwise_mat = output$pairwise_mats$avg_adj,
        k_freqtab = k_freqtab,
        k_relfreqtab = k_relfreqtab
      ))
      
    }

  
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
# prob_list_by_k = get_probs_by_k(probs = output$group_probs, 
#                                 n_groups = output$k, 
#                                 burn_in = burn_in, 
#                                 iter_threshold = t_hold)
# group_assign_list_by_k = get_assign_by_k(assign = output$group_assign, 
#                                          n_groups = output$k, 
#                                          burn_in = burn_in, 
#                                          iter_threshold = t_hold)
# 
# # correct label switching 
# stephens_result = get_stephens_result(group_assign_list_by_k = group_assign_list_by_k, 
#                                       prob_list_by_k = prob_list_by_k$prob_list)
# 
# # summarize means & variances
# mean_list_by_k_stephens = list_params_by_k(draws = output$means, 
#                                            # k_vec = output[[1]]$k,
#                                            # burn_in = burn_in, 
#                                            # iter_threshold = thold,
#                                            iter_list = prob_list_by_k$iter_list,
#                                            relabel = TRUE,
#                                            permutation = stephens_result, 
#                                            param_type = "Mean")
# 
# var_list_by_k_stephens = list_params_by_k(draws = output$vars, 
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