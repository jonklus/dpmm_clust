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

simstudy_summary <- function(output, dataset = NULL, sum_all = TRUE, print_result = TRUE, burn_in = 1000){
  # output is the list of results from the DPMM simulation study function
  # dataset is a numeric argument to summarize a specific result in the output, the desired index
  # sum_all is a logical argument to provide a summary of all data sets
  
  # filter by number of iterations for each k and address label switching
  
  
  # summarize means
  
  
  # summarize variances
  
  
  # compute KL divergence
  
  
}

output = readRDS(file = "../MCMC_Runs/conjDEVsamp_minisimstudy_noSM_2024_01_09.rds")
burn_in = 1000
thold = 1000
prob_list_by_k = get_probs_by_k(probs = output[[1]]$group_probs, 
                                n_groups = output[[1]]$k, 
                                burn_in = burn_in, iter_threshold = thold)
group_assign_list_by_k = get_assign_by_k(assign = output[[1]]$group_assign, 
                                         n_groups = output[[1]]$k, 
                                         burn_in = burn_in, iter_threshold = thold)
## deal with label switching for each k

## in this case need to skip first result bc only one group was found -- can't 
## have label switching when there is only 1 group so the function fails
# stephens_result = lapply(X = 1:length(group_assign_list_by_k), 
#                          FUN = function(x){
#                            label.switching(method = "STEPHENS", 
#                                            z = group_assign_list_by_k[[x]],
#                                            p = prob_list_by_k$prob_list[[x]])
#                          }) 

stephens_result = get_stephens_result(group_assign_list_by_k = group_assign_list_by_k, 
                                      prob_list_by_k = prob_list_by_k$prob_list)

## now reorder means to deal with label switching and redo posterior inference and 
## traceplots
mean_list_by_k_stephens = list_params_by_k(draws = output[[1]]$means, 
                                           # k_vec = output[[1]]$k,
                                           # burn_in = burn_in, 
                                           # iter_threshold = thold,
                                           iter_list = prob_list_by_k$iter_list,
                                           relabel = TRUE,
                                           permutation = stephens_result, 
                                           param_type = "Mean")

var_list_by_k_stephens = list_params_by_k(draws = output[[1]]$vars, 
                                          iter_list = prob_list_by_k$iter_list,
                                          relabel = TRUE,
                                          permutation = stephens_result,
                                          param_type = "Var")

make_postsum(mcmc_df = mean_list_by_k_stephens[[1]], digits = 2)
make_postsum(mcmc_df = mean_list_by_k_stephens[[2]], digits = 2)
make_postsum(mcmc_df = mean_list_by_k_stephens[[3]], digits = 2)

make_postsum(mcmc_df = var_list_by_k_stephens[[1]], digits = 2)
make_postsum(mcmc_df = var_list_by_k_stephens[[2]], digits = 2)
make_postsum(mcmc_df = var_list_by_k_stephens[[3]], digits = 2)

make_traceplot(param_list_by_k = mean_list_by_k_stephens, 
               k_index = 3, 
               component_no = 1,
               title_note = "K=3",
               param_type = "Mean")
make_traceplot(param_list_by_k = mean_list_by_k_stephens, 
               k_index = 3, 
               component_no = 2,
               title_note = "K=3",
               param_type = "Mean")

make_traceplot(param_list_by_k = var_list_by_k_stephens, 
               k_index = 3, 
               component_no = 1,
               title_note = "K=3",
               param_type = "Var")
make_traceplot(param_list_by_k = var_list_by_k_stephens, 
               k_index = 3, 
               component_no = 2,
               title_note = "K=3",
               param_type = "Var")