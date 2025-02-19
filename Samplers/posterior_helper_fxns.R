################################################################################
################## POSTERIOR SUMMARY HELPER FUNCTIONS ##########################
##################               VERSION 1            ##########################
################################################################################
## AUTHOR: JONATHAN KLUS
## DATE: 4 APRIL 2023
## DESCRIPTION: CONVENIENCE FUNCTIONS FOR FINDING POSTERIOR SUMMARIES FOR INFERENCE
## AND DIAGNOSTICS FROM MCMC DRAWS. DOES NOT REQUIRE INPUTS TO BE MCMC OBJECTS
## R, THOUGH THIS FUNCTIONALITY (OPTIONAL) MAY BE ADDED AT A LATER DATE.

####################### LOAD ANY REQUIRED PACKAGES #############################
library(ggplot2)
library(tidyr)
library(dplyr)
library(gtools) # mixedsort function

source("./Samplers/labelswitching_with_tol_sourceCode.R") # modified stephens
# method from label.switching with tol argument

######################### COMPUTATIONAL TOOLS ##################################
sum_log <- function(n_tot, n_1, n_2){
  # function to compute the factorial component of the second term of the Metropolis-
  # Hastings acceptance probability of the conjugate Split-Merge algorithm
  # written in terms of the MERGE step...put a negative in front of result
  # for equivalent split term
  
  # check for 1s...if either equals 1 set =2 so since 
  # factorial(0) = 1 ---> exp(log(1)) = exp(0) = 1
  if(n_1 == 1){
    n_1 = 2
  } 
  
  if(n_2 == 1){
    n_2 = 2
  }
  
  term1 = sum(log(1:(n_tot-1)))
  term2 = sum(log(1:(n_1-1)))
  term3 = sum(log(1:(n_2-1)))
  
  return(term1-(term2+term3))
  # can exponentiate to get back to original form
}

####################### FORMAT MCMC OUTPUT, SORT BY K ##########################


# this is problematic because is doesn't drop non-component columns when
# pivoted -- there needs to be some functionality to do this 
# seem to have solved this above -- we'll see...

get_probs_by_k <- function(probs, n_groups, burn_in = 50, iter_threshold = 0){
  ## Function takes a list of length no. MCMC iterations and outputs a list of 
  ## data frames where each list element is the contains the draws for iterations
  ## That found a specified number of groups k. columns in output are each parameter.
  
  # probs in a list of length # of MCMC iterations where each list element
  # contains a matrix of MCMC output with dimension (no. obs)*(no. groups) = n*k
  # n_groups is an array of length no. of MCMC iterations that describes the no. 
  # of groups found at the end of the iteration (i.e. the final k for the ith iteration)
  # burn-in is the number of initial MCMC iterations to discard
  # iter_threshold is the minimum number of iterations where a particular k was visited
  # need to cause the function to return results
  
  
  # for debugging purposes
  original_index = 1:length(probs)
  
  if(any(table(n_groups == 1))){
    singleton_labs = as.numeric(names(which(table(n_groups) == 1)))
    singleton_iters = which(n_groups %in% singleton_labs)
  } else{
    singleton_iters = NULL
  }
  
  drop_iters = unique(c(1:burn_in, singleton_iters))
  
  # drop burn-in AND any singleton iterations before proceeding
  prob_list = probs[-drop_iters]
  final_k = n_groups[-drop_iters]
  original_index = original_index[-drop_iters]
  
  # check again for singletons because some may occur after burn in dropped (i.e.
  # dropping burn in may create singletons)
  if(any(table(final_k) == 1)){
    singleton_labs = as.numeric(names(which(table(final_k) == 1)))
    singleton_iters = which(final_k %in% singleton_labs)
    final_k = final_k[-singleton_iters]
    prob_list = prob_list[-singleton_iters]
    original_index = original_index[-singleton_iters]
  } 
  
  
  # print(table(final_k))
  # print(length(prob_list))
  
  # filter out any k with less than the threshold number of iterations
  k_count = as.numeric(table(final_k))
  if(any(k_count < iter_threshold)){
    if(sum(k_count < iter_threshold) == length(k_count)){
      # if all k have n_iter < iter_threshold
      stop("Error: iter_threshold has not been met by any component k. Choose a 
           lower threshold or set iter_threshold = 0 for results.")
    } else{
      
      threshold_labs = as.numeric(names(which(table(final_k) < iter_threshold)))
      threshold_iters = which(final_k %in% threshold_labs)
      final_k = final_k[-threshold_iters]
      prob_list = prob_list[-threshold_iters]
      original_index = original_index[-threshold_iters]
    }
  }
  
  unique_k = sort(unique(final_k))
  #print(unique_k)
  # need to create separate MCMC object for each # of groups k 
  prob_list_by_k = vector(mode = "list", length = length(unique_k))
  iter_list_by_k = vector(mode = "list", length = length(unique_k))
  
  #print(unique_k)
  
  # iterate through unique k
  for(i in 1:length(unique_k)){
    
    k = unique_k[i]

    k_index = which(final_k == unique_k[i]) # indices of all iters with k probs
    # print(k_index)
    
    # preallocate list of arrays for results
    prob_list_by_k[[i]] = array(data = NA, 
                                dim = c(length(k_index), # no. iterations with k groups
                                        nrow(prob_list[[1]]), # no. observations 
                                        k)  # no. groups in iteration
                                )
    # cat("\n dim prob_list", dim(prob_list_by_k[[i]]))
    
    # also save iterations for each element to allow for easier referencing later in label switching!
    iter_list_by_k[[i]] = original_index[k_index]
    
    #print(dim(prob_list_by_k[[i]]))
    # iterate through all iterations that found the same unique k
    for(j in 1:length(k_index)){
      
      # need to output a n_iter*n_obs*k array (n_iter*n_obs*n_groups)
      # print(dim(prob_list[[k_index[j]]]))
      # print(k_index[j])
      probs_kij = prob_list[[k_index[j]]]
      
      # check if prob_list has col for new group prob. if so, rectify situation
      # by dropping "new" column
      if(ncol(prob_list[[k_index[j]]]) == (k+1)){
        
        probs_kij = probs_kij[,-(k+1)] # remove last column
        
      } else if(ncol(prob_list[[k_index[j]]]) > (k+1)){
        # print(k_index[j])
        # print(original_index[k_index[j]])
        # print(head(prob_list[[k_index[j]]]))
        stop("Error: number of group membership probs provided > (k+1)")
        
      } # else continue
      
      # cat("\n dim(probs_kij) #2", dim(probs_kij), "\n")
      
      prob_list_by_k[[i]][j,,] = probs_kij
      
    }
    
    # kgroup_probs[[i]] = array(data = unlist(ex1$probs[k_index]), 
    #                           dim = c(nrow(ex1$probs[k_index][[1]]),  
    #                                   # dimension of problem - i.e. = 2 if bivariate normal
    #                                   ncol(ex1$probs[k_index][[1]]), # number of groups (i.e. k)
    #                                   length(ex1$probs[k_index]))  # number of draws in index
    #           )
  }
  
  return(list(iter_list = iter_list_by_k, prob_list = prob_list_by_k))
}

get_assign_by_k <- function(assign, n_groups, burn_in = 50, iter_threshold = 0){
  ## Function takes a list of length no. MCMC iterations and outputs a list of 
  ## data frames where each list element contains the draws for iterations
  ## that found a specified number of groups k. columns in output are each  observation.
  
  # assign in a matrix of length # of MCMC iterations where each row is an
  # iteration and each column is an observation
  # n_groups is an array of length no. of MCMC iterations that describes the no. 
  # of groups found at the end of the iteration (i.e. the final k for the ith iteration)
  # burn-in is the number of initial MCMC iterations to discard
  
  k_vec_bi = n_groups # [(burn_in+1):nrow(assign)]
  if(any(table(k_vec_bi == 1))){
    singleton_labs = as.numeric(names(which(table(k_vec_bi) == 1)))
    singleton_iters = which(k_vec_bi %in% singleton_labs)
    # print(singleton_labs)
    # print(singleton_iters)
  } else{
    singleton_iters = NULL
  }
  
  drop_iters = unique(c(1:burn_in, singleton_iters))
  
  # drop burn-in AND any singleton iterations before proceeding
  assign_mat = assign[-drop_iters,]
  final_k = n_groups[-drop_iters]
  
  # check again for singletons because some may occur after burn in dropped (i.e.
  # dropping burn in may create singletons)
  if(any(table(final_k) == 1)){
    singleton_labs = as.numeric(names(which(table(final_k) == 1)))
    singleton_iters = which(final_k %in% singleton_labs)
    final_k = final_k[-singleton_iters]
    assign_mat = assign_mat[-singleton_iters,]
  } 
  
  unique_k = sort(unique(final_k))
  threshold_index = c()
  for(x in 1:length(unique_k)){
    k_index = which(final_k == unique_k[x])
    if(length(k_index) >= iter_threshold){
      threshold_index = c(threshold_index, x)
    }
  }
  
  if(is.null(threshold_index) == TRUE){
    # if all groups are below threshold
    stop("Error: iter_threshold has not been met by any component k. Choose a 
          lower threshold or set iter_threshold = 0 for results.")
  }
 
  assign_list_by_k = vector(mode = "list", length = length(threshold_index))
  for(index in 1:length(threshold_index)){
    
    x = threshold_index[index]
    # grab indices of all iters with k probs
    k_index = which(final_k == unique_k[x])

    testlabs = assign_mat[k_index,]
    # relabel so 1,...,k instead of skipping some numbers
    # label.switching package gets made if the max label is greater than k
    #print(dim(testlabs))
    
    # causing problems if n_iter=1 for a particular k
    # need rows = no. iterations, cols = no. obs
    newlab_mat = matrix(data = 0, nrow = length(k_index), ncol = ncol(assign_mat))
    if(length(k_index) == 1){
      # if there is only a single iteration with k groups
      
      curr_labs = as.numeric(names(table(testlabs)))
      #new_labs = 1:length(curr_labs)
      for(j in 1:length(curr_labs)){
        # relabel in order from 1:k
        newlab_mat[1,] = newlab_mat[1,] + (testlabs == curr_labs[j])*j
      }
      
    } else{
      
      for(i in 1:length(k_index)){
        curr_labs = as.numeric(names(table(testlabs[i,])))
        #new_labs = 1:length(curr_labs)
        for(j in 1:length(curr_labs)){
          # relabel in order from 1:k
          newlab_mat[i,] = newlab_mat[i,] + (testlabs[i,] == curr_labs[j])*j
        }
      }
      
    }
    
    
    assign_list_by_k[[index]] = newlab_mat
  }
  
  return(assign_list_by_k)
}

############################# LABEL SWITCHING ##################################

get_stephens_result <- function(group_assign_list_by_k, prob_list_by_k){
  # function to call modified label.switching function and generate permutation matrices
  # to fix label switching problem in MCMC output 
  
  num_params = sapply(X = 1:length(group_assign_list_by_k),
                      FUN = function(x){
                        length(unique(group_assign_list_by_k[[x]][1,]))
                      })
  
  if(1 %in% num_params){
    start_apply = 2 # stephens function will fail if you input clustering w/ 1 group 
  } else{
    start_apply = 1
  }
  
  # cat("\n num_params:", num_params, "\n")
  # cat("\n start_apply:", start_apply, "\n")
  # cat("\n length(group_assign_list_by_k):", length(group_assign_list_by_k), "\n")
  # cat("\n length(prob_list_by_k):", length(prob_list_by_k), "\n")
  
  stephens_result = lapply(X = start_apply:length(group_assign_list_by_k), 
                           FUN = function(x){
                             # use capture.output to suppress print and cat messages
                             # from label.switching function
                             # from https://stat.ethz.ch/pipermail/r-help/2008-January/151471.html
                             # note modified with lp_tol to address rare edge cases
                             log <- capture.output({
                               res <- label.switching(method = "STEPHENS", 
                                                      z = group_assign_list_by_k[[x]],
                                                      p = prob_list_by_k[[x]])
                             })
                             return(res)
                           }) 
  
  return(stephens_result)
  
}

relabel_groups <- function(curr_label_mat, permutations, means){
  
  # function to relabel group assignments and reorder means so that posterior
  # inference may be performed
  
  # curr_labels in a n_iter*n_obs matrix. it is a single list entry of the result from 
  # the get_assign_by_k function  
  # permutations is a result of the label.switching function, an n_iter*k matrix
  # that describes the relabeling scheme
  

  new_label_mat = matrix(data = 0, nrow = nrow(curr_label_mat), ncol = ncol(curr_label_mat))
  for(i in 1:nrow(curr_label_mat)){
    curr_labs = as.numeric(names(table(curr_label_mat[i,])))
    new_labs = permutations[i,]
    for(j in 1:length(curr_labs)){
      # relabel in order from 1:k
      index = which(curr_label_mat[i,] == curr_labs[j])
      new_label_mat[i,index] = new_labs[j]
    }
  }
  
  return(new_label_mat)
}


# list_means_by_k <- function(means, burn_in = 50, relabel = FALSE, permutation = NULL){
#   
#   ### DEPRECATED FUNCTION ### 
#   
#   ## Function takes a list of length no. MCMC iterations and outputs a list of 
#   ## data frames where each list element is the contains the draws for iterations
#   ## That found a specified number of groups k. columns in output are each parameter.
#   
#   # means in a list of length # of MCMC iterations where each list element
#   # contains a matrix of MCMC output with dimension (no. obs)*(no. groups) = n*k
#   # burn-in is the number of initial MCMC iterations to discard
#   # relabel is logical, if TRUE then use permutation to relabel means to solve
#   # the label switching problems
#   # permutation is a list sorted by k (# groups) resulting from the label.switching
#   # package
#   
#   mean_list = means[(burn_in+1):length(means)]
#   num_means = sapply(X = 1:length(mean_list),
#                      FUN = function(x){
#                        ncol(mean_list[[x]])
#                      })
#   unique_k = sort(unique(num_means))
#   #print(unique_k)
#   npar = nrow(mean_list[[1]]) # number of parameters per group, does not change with k
#   #print(npar)
#   # need to create separate object for each # of groups k 
#   mean_list_by_k = vector(mode = "list", length = length(unique_k))
#   
#   # if label switching solution is given, also reorder means to address this
#   if(relabel == FALSE){
#     
#     for(i in 1:length(unique_k)){
#       k_index = which(num_means == unique_k[i]) # indices of all iters with k means
#       mean_list_by_k[[i]] = data.frame(matrix(data = unlist(mean_list[k_index]), 
#                                               ncol = npar*unique_k[i],
#                                               byrow = TRUE))
#       # name columns i.e. mu23 is mean for group 2, 3rd component 
#       # (i.e. from a length 3 mean vector)
#       col_header_names = unlist(lapply(X = 1:npar, 
#                                        FUN = function(x){
#                                          paste0("mu", 1:unique_k[i], x)
#                                        }))
#       names(mean_list_by_k[[i]]) = sort(col_header_names)
#       #print(unique_k[i])
#       #print(col_header_names)
#       # this format isn't super helpful - need to get in data frame format
#       # i.e. nrow = n_iter, column for mu11, mu12, mu21, m22, etc...
#       
#       # kgroup_means[[i]] = array(data = unlist(ex1$means[k_index]), 
#       #                           dim = c(nrow(ex1$means[k_index][[1]]),  
#       #                                   # dimension of problem - i.e. = 2 if bivariate normal
#       #                                   ncol(ex1$means[k_index][[1]]), # number of groups (i.e. k)
#       #                                   length(ex1$means[k_index]))  # number of draws in index
#       #           )
#     }
#     
#   } else{
#     
#     # if relabel == TRUE
#     if(1 %in% unique_k){
#       # unique_k = unique_k[-1] # drop the 1 before correcting -- no permutations
#       i_init = 2  # start after k=1
#       permutation = c(1, permutation) # put placeholder in for k=1 group (not 
#       # included in stephens results bc no label switching) -- avoids error with
#       # indexing later on in for loop otherwise index doesn't match in means & perm
#       
#       # put the k=1 group names in here 
#       k_index = which(num_means == 1)
#       mean_list_by_k[[1]] = data.frame(matrix(data = unlist(mean_list[k_index]),   
#                                               ncol = npar,                   
#                                               byrow = TRUE))
#       # name columns i.e. mu23 is mean for group 2, 3rd component 
#       # (i.e. from a length 3 mean vector)
#       names(mean_list_by_k[[1]]) = paste0("mu", 1, 1:npar)
#       
#     } else{
#       
#       i_init = 1 # business as usual
#       
#     }
#     
#     for(i in i_init:length(unique_k)){
#       
#       #print(unique_k[i])
#       
#       k_index = which(num_means == unique_k[i]) # indices of all iters with k means
#       new_labs = permutation[[i]]$permutations$STEPHENS
#       
#       # fix label switching before reformatting means
#       for(j in 1:length(k_index)){
#         
#         #print(mean_list[[k_index[j]]]) # before fix
#         #print(new_labs[j,]) # new order
#         
#         mean_list[[k_index[j]]] =  mean_list[[k_index[j]]][,new_labs[j,]]
#         
#         # print(mean_list[[k_index[j]]]) # after fix
#       }
#       
#       # reformat means now that order has been corrected
#       mean_list_by_k[[i]] = data.frame(matrix(data = unlist(mean_list[k_index]), 
#                                               ncol = npar*unique_k[i], 
#                                               byrow = TRUE))
#       # name columns i.e. mu23 is mean for group 2, 3rd component 
#       # (i.e. from a length 3 mean vector)
#       col_header_names = unlist(lapply(X = 1:npar, 
#                                        FUN = function(x){
#                                          paste0("mu", 1:unique_k[i], x)
#                                        }))
#       names(mean_list_by_k[[i]]) = sort(col_header_names)
#       #print(unique_k[i])
#       #print(col_header_names)
#       # this format isn't super helpful - need to get in data frame format
#       # i.e. nrow = n_iter, column for mu11, mu12, mu21, m22, etc...
#       
#       # kgroup_means[[i]] = array(data = unlist(ex1$means[k_index]), 
#       #                           dim = c(nrow(ex1$means[k_index][[1]]),  
#       #                                   # dimension of problem - i.e. = 2 if bivariate normal
#       #                                   ncol(ex1$means[k_index][[1]]), # number of groups (i.e. k)
#       #                                   length(ex1$means[k_index]))  # number of draws in index
#       #           )
#     }
#     
#   }
#   
#   
#   return(mean_list_by_k)
# }


list_params_by_k <- function(draws, iter_list, k_vec, off_diag = FALSE, dont_drop = FALSE, # burn_in = 50, iter_threshold = 0, 
                             relabel = FALSE, permutation = NULL, param_type, equal_var = FALSE){
  ## Function takes a list of length no. MCMC iterations and outputs a list of 
  ## data frames where each list element contains the draws for iterations
  ## That found a specified number of groups k. columns in output are each parameter.
  
  # draws in a list of length # of MCMC iterations where each list element
  # contains a matrix of MCMC output with dimension (no. obs)*(no. groups) = n*k
  # burn-in is the number of initial MCMC iterations to discard
  # relabel is logical, if TRUE then use permutation to relabel means to solve
  # the label switching problems
  # permutation is a list sorted by k (# groups) resulting from the label.switching
  # package
  # param_type is a string, either mean or var
  # equal_var is a logical argument of whether the equal variance assumption is being made
  # it means there will be only a single variance param for each group
  
  # k_vec_bi = k_vec[(burn_in+1):length(k_vec)]
  # if(any(table(k_vec_bi == 1))){
  #   singleton_labs = as.numeric(names(which(table(k_vec_bi) == 1)))
  #   singleton_iters = which(k_vec_bi %in% singleton_labs)
  #   #print(singleton_iters)
  # } else{
  #   singleton_iters = NULL
  # }
  # 
  # 
  # # filter out any k with less than the threshold number of iterations
  # k_count = as.numeric(table(k_vec)) #as.numeric(table(k_vec_bi))
  # if(any(k_count < iter_threshold)){
  #   if(sum(k_count < iter_threshold) == length(k_count)){
  #     # if all k have n_iter < iter_threshold
  #     stop("Error: iter_threshold has not been met by any component k. Choose a 
  #          lower threshold or set iter_threshold = 0 for results.")
  #   } else{
  #     
  #     threshold_labs = as.numeric(names(which(table(k_vec_bi) < iter_threshold)))
  #     threshold_iters = which(k_vec_bi %in% threshold_labs)
  #     # final_k = k_vec_bi[-threshold_iters]
  #     # original_index = original_index[-threshold_iters]
  #   }
  # } else{
  #   threshold_iters = NULL
  # }
  # # all iterations to drop, including burn-in and singletons
  # # include unique function bc some singletons may come in during burn-in iterations
  # 
  # drop_iters = unique(c(1:burn_in, singleton_iters, threshold_iters))
  
  keep_iters = unlist(iter_list)
  if(param_type == "Mean"){
    
    # drop burn-in AND any singleton iterations before proceeding
    param_list = draws[keep_iters]
    param_symbol = "mu"
    k_vec = k_vec[keep_iters]
    # num_params = k_vec[keep_iters]
    num_params = sapply(X = 1:length(param_list),
                        FUN = function(x){
                          ncol(param_list[[x]])
                        })
    
    #print(table(num_params)) # debugging -- have we eliminated singletons?
    
  } else if(param_type == "Var"){
    
    # drop burn-in AND any singleton iterations before proceeding

    if(equal_var == TRUE & dont_drop == TRUE){
      temp_param_list = draws # don't drop iters
      param_symbol = "sigma"
      num_params = rep(1, length(temp_param_list))
    } else if(equal_var == TRUE){
      temp_param_list = draws[keep_iters]
      k_vec = k_vec[keep_iters]
      param_symbol = "sigma"
      num_params = rep(1, length(temp_param_list))
    } else{
      temp_param_list = draws[keep_iters]
      k_vec = k_vec[keep_iters]
      param_symbol = "sigma"
      # num_params = k_vec[keep_iters]
      num_params = sapply(X = 1:length(temp_param_list),
                          FUN = function(x){
                            length(temp_param_list[[x]])
                          })
    }
    

    
    # get var diagonal components into similar format as means
    param_list = vector(mode = "list", length = length(num_params))
    for(i in 1:length(num_params)){
      
      param_list[[i]] = sapply(X = 1:num_params[i], 
                               FUN = function(x){ # if a scalar -- diag will cause weird results
                                 if(is.null(dim(temp_param_list[[i]][[x]])) == TRUE){
                                   temp_param_list[[i]][[x]]
                                 } else{ # if a matrix
                                   diag(temp_param_list[[i]][[x]])
                                 }
                               })
      # print(param_list[[i]])
      
    }
    
  } else if(param_type == "Covar"){
    # Option to include off-diagonal elements and diagonal elements
    # drop burn-in AND any singleton iterations before proceeding
    temp_param_list = draws[keep_iters]
    k_vec = k_vec[keep_iters]
    param_symbol = "sigma"
    if(equal_var == TRUE){
      
      stop("Incorrect combination of args. Cannot have param_type = Covar and equal_var = TRUE simultaneously.")
      
    } else{
      # num_params = k_vec[keep_iters]
      # num params is really just the number of cov matrices -- used to pull out
      # values, not used to index sigma_11, sigma_12, for example
      num_params = sapply(X = 1:length(temp_param_list),
                          FUN = function(x){
                            length(temp_param_list[[x]])
                          })
    }
    
    
    
    # get var diagonal components into similar format as means
    param_list = vector(mode = "list", length = length(num_params))
    for(i in 1:length(num_params)){
      
      param_list[[i]] = sapply(X = 1:num_params[i], 
                               FUN = function(x){ # if a scalar -- diag will cause weird results
                                 if(is.null(dim(temp_param_list[[i]][[x]])) == TRUE){
                                   temp_param_list[[i]][[x]]
                                 } else{ # if a matrix
                                   c(diag(temp_param_list[[i]][[x]]), 
                                     temp_param_list[[i]][[x]][upper.tri(temp_param_list[[i]][[x]])])
                                 }
                               })
      # each column is a group k, each row is a parameter sigma (either variance or covariance)
      # that has been pulled out of the original covariance matrix and been reorganized
      # print(param_list[[i]])
      
    }
    
  } else{
    
    stop("param_type not recognized. Please specify Mean, Var, or Covar only.")
    
  }
  

    if(equal_var == TRUE & param_type == "Var"){
      
      unique_k = 1
      
    } else{
      
      unique_k = sort(unique(k_vec))
      
    }
    
 
 
  # cat("\n unique_k:", unique_k)

  # number of parameters per group, does not change with k
    if(param_type == "Covar"){
      # dimension of original covariance matrix  
      p = ifelse(is.null(nrow(param_list[[1]])) == TRUE, 1, 
                 nrow(temp_param_list[[1]][[1]]))
      npar = ifelse(is.null(nrow(param_list[[1]])) == TRUE, 1, 
                    nrow(param_list[[1]]))
      
    } else{
      
      npar = ifelse(is.null(nrow(param_list[[1]])) == TRUE, 1, nrow(param_list[[1]]))
      
    }
    
  
  
  # print(npar)
  # need to create separate object for each # of groups k
    if(equal_var == TRUE & param_type == "Var"){
      
      param_list_by_k = vector(mode = "list", length = 1) # not actually going to 
      # sort by k since variance param is pooled
      
    } else{
      
      param_list_by_k = vector(mode = "list", length = length(unique_k))
      
    }
  
  
  # if label switching solution is given, also reorder params to address this
  if(relabel == FALSE){

    for(i in 1:length(unique_k)){
      
      k_index = which(k_vec == unique_k[i]) # indices of all iters with k params
      
      if(equal_var == TRUE & param_type == "Var"){

        # pooled variance, single param
        # don't sort or do anything by k_list
        param_list_by_k[[i]] = data.frame(matrix(data = unlist(param_list), 
                                                 ncol = 1,
                                                 byrow = TRUE))
        # name columns i.e. mu23 is param for group 2, 3rd component 
        # (i.e. from a length 3 mean vector)
        col_header_names = c("sigma_1_1")
        names(param_list_by_k[[i]]) =  col_header_names
        
      } else if(param_type == "Covar"){
        
        param_list_by_k[[i]] = data.frame(matrix(data = unlist(param_list[k_index]), 
                                                 ncol = npar*unique_k[i], 
                                                 byrow = TRUE))
        
        param_mat = t(sapply(X = 1:p, 
                             FUN = function(x){
                               i = 1:p
                               paste0(x,"_",i)
                             }))
        
        param_lab = c(diag(param_mat), param_mat[upper.tri(param_mat)])
        
        col_header_names = unlist(lapply(X = 1:unique_k[i], 
                                         FUN = function(x){
                                           paste0(param_symbol, "_", x, "_", param_lab)
                                         }))
        # names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
        names(param_list_by_k[[i]]) = col_header_names # dont sort this one, 
        # the apply statement was written to ensure correct ordering
        
      } else{
        
        param_list_by_k[[i]] = data.frame(matrix(data = unlist(param_list[k_index]), 
                                                 ncol = npar*unique_k[i],
                                                 byrow = TRUE))
        # name columns i.e. mu23 is param for group 2, 3rd component 
        # (i.e. from a length 3 mean vector)
        col_header_names = unlist(lapply(X = 1:npar, 
                                         FUN = function(x){
                                           paste0(param_symbol, "_", 1:unique_k[i], "_", x)
                                         }))
        names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
        
        # if(off_diag == TRUE){
        #   
        #   stop("This function not yet available. Please contact your manager for assistance. ")
        # }
        
      }
      

      #print(unique_k[i])
      #print(col_header_names)
      # this format isn't super helpful - need to get in data frame format
      # i.e. nrow = n_iter, column for mu11, mu12, mu21, m22, etc...
      
      # kgroup_means[[i]] = array(data = unlist(ex1$means[k_index]), 
      #                           dim = c(nrow(ex1$means[k_index][[1]]),  
      #                                   # dimension of problem - i.e. = 2 if bivariate normal
      #                                   ncol(ex1$means[k_index][[1]]), # number of groups (i.e. k)
      #                                   length(ex1$means[k_index]))  # number of draws in index
      #           )
    }
    
  } else{
    
    # if relabel == TRUE
    
    if(equal_var == TRUE){
      
      i_init = 1
      
    } else{
      
      if(1 %in% unique_k){
        # unique_k = unique_k[-1] # drop the 1 before correcting -- no permutations
        i_init = 2  # start after k=1
        permutation = c(1, permutation) # put placeholder in for k=1 group (not 
        # included in stephens results bc no label switching) -- avoids error with
        # indexing later on in for loop otherwise index doesn't match in params & perm
        
        # put the k=1 group names in here 
        k_index = which(num_params == 1)
        
        if(equal_var == TRUE & param_type == "Var"){
          
          param_list_by_k[[1]] = data.frame(matrix(data = unlist(param_list[k_index]),   
                                                   ncol = 1,                   
                                                   byrow = TRUE))
          
          # name columns i.e. mu23 is param for group 2, 3rd component 
          # (i.e. from a length 3 mean vector)
          col_header_names = c("sigma_1_1")
          names(param_list_by_k[[1]]) =  col_header_names
          
        } else if(param_type == "Covar"){
          
          param_list_by_k[[1]] = data.frame(matrix(data = unlist(param_list[k_index]), 
                                                   ncol = npar, 
                                                   byrow = TRUE))
          
          param_mat = t(sapply(X = 1:p, 
                               FUN = function(x){
                                 i = 1:p
                                 paste0(x,"_",i)
                               }))
          
          param_lab = c(diag(param_mat), param_mat[upper.tri(param_mat)])
          
          col_header_names = paste0(param_symbol, "_", 1, "_", param_lab)
          
          names(param_list_by_k[[1]]) =  col_header_names
          
        } else{
          # any other params - means, vars
          
          param_list_by_k[[1]] = data.frame(matrix(data = unlist(param_list[k_index]),   
                                                   ncol = npar,                   
                                                   byrow = TRUE))
          
          
          # name columns i.e. mu23 is mean for group 2, 3rd component 
          # (i.e. from a length 3 mean vector)
          
          names(param_list_by_k[[1]]) = paste0(param_symbol, "_", 1, "_", 1:npar)
          
        }
        

        
      } else{
        
        i_init = 1 # business as usual
        
      }
      
    }
    
    
    
    for(i in i_init:length(unique_k)){
      
      # cat("\n i:", i)
      # cat("\n unique_k[i]:", unique_k[i])
      # cat("\n unique_k:", unique_k)
      # print(head(permutation[[i]]))
      
      k_index = which(k_vec == unique_k[i]) # indices of all iters with k params
      
      if(equal_var == FALSE){
        
        new_labs = permutation[[i]]$permutations$STEPHENS # moving inside if 
        # statement to avoid errors when equal_var == TRUE and we don't alter 
        # permutation object
        
        # fix label switching before reformatting params
        for(j in 1:length(k_index)){
          
          # cat("\n k_index[j]", k_index[j], "\n")
          # cat("param_vals", "\n")
          # print(param_list[[k_index[j]]])
          # cat("\n")
          # cat("new_labs", new_labs[j,], "\n")
          
          if(is.null(dim(param_list[[k_index[j]]])) == TRUE){
            # if scalar -- will return an array not a matrix
            
            # cat("\n j=", j)
            # cat("\n k_index[j]", k_index[j])
            # cat("\n param =", param_list[[k_index[j]]])
            
            param_list[[k_index[j]]] = param_list[[k_index[j]]][new_labs[j,]]
            
          } else{
            # matrix
            
            param_list[[k_index[j]]] =  param_list[[k_index[j]]][,new_labs[j,]]
            
          }
          
          
        }
        
      } else{
        
        # proceed, no need to permute assignments
        
      }
      
      
      # reformat params now that order has been corrected
      if(npar == 1){
        
        if(equal_var == TRUE & param_type == "Var"){
          
          param_list_by_k[[i]] = data.frame(matrix(data = unlist(param_list), 
                                                   ncol = 1, 
                                                   byrow = TRUE))
          
          # name columns i.e. mu23 is mean for group 2, 3rd component 
          # (i.e. from a length 3 mean vector)
          col_header_names = paste0(param_symbol, "_", 1, "_", 1)
          names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
          
        } else if(param_type == "Covar"){
          
          # covar shouldn't come up here, npar = 1 and covar aren't really compatible
          # devise a new error check for this combination or leave out entirely
          
          stop("npar = 1 and param_type = Covar not compatible.")
          
        } else{
          
          param_list_by_k[[i]] = data.frame(matrix(data = unlist(param_list[k_index]), 
                                                   ncol = unique_k[i], 
                                                   byrow = TRUE))
          
          # name columns i.e. mu23 is mean for group 2, 3rd component 
          # (i.e. from a length 3 mean vector)
          col_header_names = paste0(param_symbol, "_", 1:unique_k[i], "_1")
          names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
          
        }
        
        
      } else{
        
          param_list_by_k[[i]] = data.frame(matrix(data = unlist(param_list[k_index]), 
                                                   ncol = npar*unique_k[i], 
                                                   byrow = TRUE))
          
          # name columns i.e. mu23 is mean for group 2, 3rd component 
          # (i.e. from a length 3 mean vector)
          if(param_type == "Var"){
            
            col_header_names = unlist(lapply(X = 1:npar, 
                                             FUN = function(x){
                                               paste0(param_symbol, "_", 1:unique_k[i], "_", x)
                                             }))
            names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
            
          } else if(param_type == "Covar"){
            
            param_mat = t(sapply(X = 1:p, 
                                FUN = function(x){
                                  i = 1:p
                                  paste0(x,"_",i)
                                }))
            
            param_lab = c(diag(param_mat), param_mat[upper.tri(param_mat)])
            
            col_header_names = unlist(lapply(X = 1:unique_k[i], 
                                             FUN = function(x){
                                               paste0(param_symbol, "_", x, "_", param_lab)
                                             }))
            # names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
            names(param_list_by_k[[i]]) = col_header_names # dont sort this one, 
            # the apply statement was written to ensure correct ordering
            
          } else{
            
            col_header_names = unlist(lapply(X = 1:npar, 
                                             FUN = function(x){
                                               paste0(param_symbol, "_", 1:unique_k[i], "_", x)
                                             }))
            names(param_list_by_k[[i]]) = gtools::mixedsort(col_header_names)
            
          }

          
      }
        
      #print(unique_k[i])
      #print(col_header_names)
      # this format isn't super helpful - need to get in data frame format
      # i.e. nrow = n_iter, column for mu11, mu12, mu21, m22, etc...
      
      # kgroup_params[[i]] = array(data = unlist(ex1$params[k_index]), 
      #                           dim = c(nrow(ex1$params[k_index][[1]]),  
      #                                   # dimension of problem - i.e. = 2 if bivariate normal
      #                                   ncol(ex1$params[k_index][[1]]), # number of groups (i.e. k)
      #                                   length(ex1$params[k_index]))  # number of draws in index
      #           )
    }
    
  }
  
  
  return(list(original_draws = draws, param_list_by_k = param_list_by_k))
}

###################### POSTERIOR SUMMARY STATISTICS ############################

make_postsum <- function(mcmc_df, q = 0.95, digits = 2, caption = NULL){
  # takes a data frame of MCMC draws and outputs posterior summaries
  # mean, variance, and q% posterior credible intervals
  # digits is the argument for how many decimal places to include when reporting results
  
  # burn_in in the number of burn-in iterations to discard before performing inference
  # dont do this here - do it before you sort data by k
  # mcmc_df = mcmc_df[(burn_in+1):nrow(mcmc_df),] # discard burn-in -- already doing this in early fxns
  
  lower_q = (1 - q)/2
  upper_q  = (1 + q)/2
  
  postmean = sapply(X = 1:ncol(mcmc_df), 
                    FUN = function(x){mean(mcmc_df[,x], na.rm = TRUE)})
  
  postmed = sapply(X = 1:ncol(mcmc_df), 
                   FUN = function(x){median(mcmc_df[,x], na.rm = TRUE)})
  
  postvar = sapply(X = 1:ncol(mcmc_df), 
                   FUN = function(x){var(mcmc_df[,x], na.rm = TRUE)})
  
  postlower = sapply(X = 1:ncol(mcmc_df), 
                     FUN = function(x){quantile(mcmc_df[,x], probs = lower_q, na.rm = TRUE)})
  
  postupper = sapply(X = 1:ncol(mcmc_df), 
                     FUN = function(x){quantile(mcmc_df[,x], probs = upper_q, na.rm = TRUE)})
  
  # postlower = sapply(X = mcmc_df, FUN = quantile, probs = lower_q)
  # postupper = sapply(X = mcmc_df, FUN = quantile, probs = upper_q)
  
  
  
  postsum_df = round(data.frame(
    Mean = postmean,
    Median = postmed,
    Var = postvar,
    Lower = postlower,
    Upper = postupper
  ), digits = digits)
  
  rownames(postsum_df) = names(mcmc_df)
  names(postsum_df) = c("Mean", "Median", "Empirical SE", 
                        paste0(100*lower_q, "%"),
                        paste0(100*upper_q, "%"))
  
  if(is.null(caption) == TRUE){
    return(postsum_df)
  } else{
    
  }
  return(knitr::kable(x = postsum_df, caption = caption))
}


pairwise_prob_mat <- function(burn_in = NULL, group_assign, probs = NULL, diag_weights = FALSE){
  
  # burn-in is whether to eliminate the first M iterations
  # group_assign is a S*n matrix of group assignments
  # probs is a list of length S of n*k matrices of group assignment probs from 
  # each mcmc iteration
  # diag_weights indicates whether or not diagonal elements of pairwaise prob matrix should contain
  # weights Pr(X1 in c1)^2 + Pr(X1 in c2)^2 + ....  instead of 1's
  
  if(is.null(burn_in) == TRUE){
    
    n_iter = nrow(group_assign)
    n_obs = ncol(group_assign)
    
  } else{
    
    n_iter = nrow(group_assign) - burn_in
    n_obs = ncol(group_assign)
    group_assign = group_assign[-(1:burn_in),]
    
  }
  
  
  # preallocate list for results
  adj_matrix_list = vector(mode = "list", length = n_iter)
  # cycle through each iteration and create adjacency matrix based on 
  # final group assignments for that iteration
  # result is a symmetric matrix where each element is either 0 or 1
  for(iter in 1:n_iter){
    
    # create n*n matrix of 0 for adjacency of each observation
    adj_matrix_list[[iter]] = matrix(data = 0, 
                                     nrow = n_obs, 
                                     ncol = n_obs)
    
    for(obs in 1:n_obs){
      # check whether each observation in a given iteration is in the same group
      # as any other observation -- indicator 0/1
      adj_matrix_list[[iter]][obs,] = (as.numeric(group_assign[iter,obs] ==
                                                    group_assign[iter,]))
    }
    
    # put in option to make diagonals weights instead of 1
    if(diag_weights == TRUE){
      diag_elements = rowSums(probs[[iter]]^2)
      diag(adj_matrix_list[[iter]]) = diag_elements
    }
    
  }
  
  # find average
  avg_adj = (1/n_iter)*Reduce(f = "+", x = adj_matrix_list)
  
  return(list(avg_adj = avg_adj, adj_by_iter = adj_matrix_list))
  
}




################################ TRACEPLOTS ####################################

make_k_traceplot <- function(k, group_assign, burn_in = NULL, show_min_obs = FALSE, max_yaxis = NULL){
  
  # function to make traceplot for number of groups at each mcmc iteration
  
  if(is.null(burn_in) == TRUE){
    
    n_iter = length(k)
    
  } else{
    
    n_iter = length(k) - burn_in
    # group_assign = group_assign[-c(1:burn_in),]
    k = k[-c(1:burn_in)]
    
  }
  
  if(is.null(max_yaxis) == TRUE){
    
    max_yaxis = max(k)+1
    
  } else{
    # else use provided value unless not ideal
    if(max_yaxis < (max(k)+1)){
      stop(cat("Given maximum y axis value is less than observed max +1. Input a value of ", 
               max(k)+1, " or larger, or allow function defaults."))
    } # else proceed, all good!
  }
  
  min_obs = sapply(X = 1:nrow(group_assign), 
                   FUN = function(x){
                     min(table(group_assign[x,]))
                   })
  
  # ref: https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales
  
  ylim.prim = c(0, max(k)+1)  
  
  if(show_min_obs == TRUE){
    
    ylim.sec = c(0, ncol(group_assign))  
    
    b = diff(ylim.prim)/diff(ylim.sec)
    a = ylim.prim[1] - b*ylim.sec[1]
    
    ggplot(mapping = aes(x = 1:n_iter)) + 
      geom_col(aes(y = a + min_obs*b), color = "lightblue") +
      geom_line(aes(y = k)) +
      ggtitle("no. groups K") +
      xlab("iter") + 
      #ylab("K") + 
      theme_classic() + 
      scale_y_continuous(
        # Features of the first axis
        name = "K",
        breaks = seq(0, max_yaxis, by = 1), limits = c(0,max_yaxis),
        # Add a second axis and specify its features
        sec.axis = sec_axis(~ (. - a)/b, name="min # obs.", ))
    
  } else{
    
    ggplot(mapping = aes(x = 1:n_iter, y = k)) + 
      geom_line() +
      ggtitle("no. groups K") +
      scale_y_continuous(breaks = seq(0, max_yaxis, by = 1), 
                         limits = c(0,max_yaxis)) + 
      xlab("iter") + 
      ylab("K") + 
      theme_classic()
  }
  
  
}



make_traceplot <- function(param_list_by_k, k, component_no, param_type, p, # burn_in = 0,
                           legend = FALSE, log_scale = FALSE, title_text = NULL){
  # takes a data frame of MCMC draws and outputs labeled traceplots
  # generic version of make_mean_traceplot
  
  # param_list_by_k is the product of function list_params_by_k
  # k is the number of groups to plot (i.e. index corresponds to 
  # which number of groups we want to see a traceplot for)
  # p is the dimension of the problem (i.e. a bivariate normal will have p=2)
  # component_no is the dimension to plot (i.e. for a length 3 mean can choose
  # 1,2, or 3 (mu1, mu2, mu3))
  # title_note is a string to include at the end of the plot title
  # param_type is a string to include in the plot title, i.e. Mean,  Var, etc.
  
  num_params = sapply(X = 1:length(param_list_by_k),
                      FUN = function(x){
                        ncol(param_list_by_k[[x]][1,])
                      })
  
  k_index = which(num_params/p == k) # which list element contains the 
  # number of groups that we wish to plot
  # print(num_params)
  # print(k_index)
  
  # get data frame into normal format the ggplot can use by pivoting cols
  mcmc_df = param_list_by_k[[k_index]] #[(burn_in+1):nrow(param_list_by_k[[k_index]]),] # discard burn-in
  mcmc_df$iter = 1:nrow(mcmc_df)
  # identify columns containing draws for the component of interest
  component_cols = grep(pattern = paste0(component_no, "$"), x = names(mcmc_df))
  plot_df = tidyr::pivot_longer(data = mcmc_df %>% 
                                  dplyr::select(iter, all_of(component_cols)), 
                                cols = names(mcmc_df)[component_cols],
                                names_to = "component",
                                values_to = "param")
  
  # make color coded traceplots
  if(is.null(title_text) == TRUE){
    title_text = paste0(param_type, " Component ", component_no, " (K=", k, ")")
  } else{
    # use provided title text
    title_text = paste0(title_text, " (K=", k, ")")
  }
  
  
  if(log_scale == TRUE){
    # tol = .001
    # plot_df$param = ifelse(abs(plot_df$param) <tol, plot_df$param+1, plot_df$param) 
    # check no params  = 0
    if(sum(plot_df$param == 0)){
      stop("Log scale not defined. Parameter takes on values that are equal to 0.")
    }
    
    ggplot(data = plot_df) +
      geom_line(mapping = aes(x = iter, y = log(param), color = component)) +
      ggtitle(title_text) +
      # ylab() +
      theme_classic() +
      ggplot2::theme(plot.title = element_text(size=12), legend.position = ifelse(legend == TRUE, "right","none"))
    
  } else{
    
    ggplot(data = plot_df) +
      geom_line(mapping = aes(x = iter, y = param, color = component)) +
      ggtitle(title_text) +
      # ylab() +
      theme_classic() +
      ggplot2::theme(plot.title = element_text(size=12), legend.position = ifelse(legend == TRUE, "right","none"))
    
  }
  
  
}

make_mean_traceplot <- function(mean_list_by_k, k_index, component_no, title_note){
  # mean_list_by_k is the product of function list_means_by_k
  # k_index is the element of the list to plot (i.e. index corresponds to 
  # which number of groups)
  # component_no is the dimension to plot (i.e. for a length 3 mean can choose
  # 1,2, or 3 (mu1, mu2, mu3))
  # title_note is a string to include at the end of the plot title
  
  # get data frame into normal format the ggplot can use by pivoting cols
  mcmc_df = mean_list_by_k[[k_index]]
  mcmc_df$iter = 1:nrow(mcmc_df)
  # identify columns containing draws for the component of interest
  component_cols = grep(pattern = paste0(component_no, "$"), x = names(mcmc_df))
  plot_df = tidyr::pivot_longer(data = mcmc_df %>% 
                                  dplyr::select(iter, all_of(component_cols)), 
                                cols = names(mcmc_df)[component_cols],
                                names_to = "component",
                                values_to = "mu")
  
  # make color coded traceplots
  ggplot(data = plot_df) +
    geom_line(mapping = aes(x = iter, y = mu, color = component)) +
    ggtitle(paste0("Traceplot of Mean Component ", component_no, " (", title_note, ")")) +
    #ylab("mu") +
    theme_classic()
  
}

mixprop_traceplot <- function(group_assign, burn_in = 100){
  # function takes S*n matrix of group assignments as an input, calculates the 
  # mixing proportions pi_1,pi_K at each iteration S, and plots them
  # if a group k does not exist at a particular iteration S, then pi_k is set to 0
  
  # define maximum number of groups in data
  max_k_byiter = sapply(X = 1:nrow(group_assign), 
                        FUN = function(x){
                          length(unique(group_assign[x,]))
                        })
  
  # preallocate matrix for mixing props...prefill with zeros so any unused
  # columns will set pi_k equal to implied value of zero since group does not exist
  # at that iteration
  mix_prop = matrix(data = 0, nrow = nrow(group_assign), ncol = max(max_k_byiter))
  
  # compute mixing proportions, ordering by group size to keep them consistent
  # (not addressing label switching here)
  for(i in 1:nrow(group_assign)){
    mix_prop[i,1:max_k_byiter[i]] = sort(as.numeric(table(group_assign[i,]))/ncol(group_assign), 
                                         decreasing = TRUE)
  }
  
  # format as data frame and put into long format so that ggplot can handle
  if(burn_in > 0){
    mix_prop = as.data.frame(mix_prop)[(burn_in+1):nrow(mix_prop),]
  } else{
    mix_prop = as.data.frame(mix_prop)
  }
  names(mix_prop) = paste0("pi_", 1:max(max_k_byiter))
  mix_prop$iter = 1:nrow(mix_prop) # add iteration column for plotting
  

  mix_prop = tidyr::pivot_longer(data = mix_prop, 
                                 cols = starts_with("pi_"),
                                 names_to = "Component",
                                 # names_prefix = "pi_",
                                 values_to = "pi")
  
  # return plot
  ggplot2::ggplot(data = mix_prop, aes(x = iter, y = pi, color = Component)) +
    ggplot2::geom_line() +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0("Traceplot of Mixing Proportions \u03C0", ", max(K)=",max(max_k_byiter))) +
    ggplot2::theme(legend.position = "none")
}






######################### MULTI-CHAIN DIAGNOSTICS ##############################

# future addition -- 
# perform model diagnostics - gelman rubin?
# how to do multichain diagnostics when the number of draws for each k
# varies from chain to chain (i.e. the data won't be the same length!)



################################# METRICS ######################################

correct_group_assign <- function(group_assign_list_by_k, stephens_result){
  # need to update group labels to address label switching so that correct
  # mean estimates and variance estimates are referenced when calculating KL
  # divergence and other metrics
  group_assign_k_raw = group_assign_list_by_k
  group_assign_k_corr = vector(mode = "list", length = length(group_assign_k_raw))
  
  num_params = sapply(X = 1:length(group_assign_list_by_k),
                      FUN = function(x){
                        length(unique(group_assign_list_by_k[[x]][1,]))
                      })
  
  if(1 %in% num_params){
    start_apply = 2 # stephens function will fail if you input clustering w/ 1 group 
    group_assign_k_corr[[1]] = group_assign_k_raw[[1]] 
    
    for(list_el in start_apply:length(group_assign_k_raw)){
      
      # preallocate matrix for kth element
      group_assign_k_corr[[list_el]] = matrix(data = NA, 
                                              nrow = nrow(group_assign_list_by_k[[list_el]]), 
                                              ncol = ncol(group_assign_list_by_k[[list_el]]))
      
      for(iter in 1:nrow(group_assign_k_raw[[list_el]])){
        
        # cat("\n list_el:", list_el, "\n")
        # cat("\n iter:", iter, "\n")
        
        permutation = stephens_result[[(list_el-1)]]$permutations$STEPHENS[iter,]
        assign_old = group_assign_k_raw[[list_el]][iter,]
        
        for(perm_el in 1:length(permutation)){
          
          old_index = which(assign_old == perm_el) # which is ith element
          group_assign_k_corr[[list_el]][iter,old_index] = permutation[perm_el] 
          # assign new label
          
        }
      }
    }
    
  } else{
    
    start_apply = 1
    
    for(list_el in start_apply:length(group_assign_k_raw)){
      
      # preallocate matrix for kth element
      group_assign_k_corr[[list_el]] = matrix(data = NA, 
                                              nrow = nrow(group_assign_list_by_k[[list_el]]), 
                                              ncol = ncol(group_assign_list_by_k[[list_el]]))
      
      for(iter in 1:nrow(group_assign_k_raw[[list_el]])){
        
        # cat("\n list_el:", list_el, "\n")
        # cat("\n iter:", iter, "\n")
        
        permutation = stephens_result[[list_el]]$permutations$STEPHENS[iter,]
        assign_old = group_assign_k_raw[[list_el]][iter,]
        
        for(perm_el in 1:length(permutation)){
          
          old_index = which(assign_old == perm_el) # which is ith element
          group_assign_k_corr[[list_el]][iter,old_index] = permutation[perm_el] 
          # assign new label
          
        }
      }
    }
  }
    

  return(group_assign_k_corr)
}


calc_KL_diverg <- function(y, mu_est, Sigma_est, group_assign, true_assign, 
                           mu_true, Sigma_true, 
                           # equal_var_datagen = TRUE, consider this arg in the future
                           equal_var_assump = FALSE, off_diag = FALSE){
  
  # function to calculate KL divergence between truth and posterior mean for
  # a particular clustering (e.g. k=3) after running MCMC
  
  # y, mu_est, Sigma_est are lists - mu and sigma are MCMC output that has been
  # if equal var assumption is true -- Sigma_est and Sigma_true are matrices
  # filtered by k and processed to correct label switching so that they are aligned
  # mu_true, and Sigma_true are the ground truth from which data were simulated
  # truth is a vector of true group assignments
  # group_assign is a K-length list of S*n matrix of estimated group assignments
  # after correction using the correct_group_assign fxn
  # equal_var_assump is a logical that determines whether we assume the groups
  # had unique variances estimated as part of the MCMC
  
  # assumes MV normal likelihood
  
  if(equal_var_assump == FALSE){
    # each group has its own estimated variance
    
    # calculate density of truth

    
    true_dens = sapply(X = 1:length(y), 
                       FUN = function(x){
                         # debugging
                         # cat("\n x:", x, "true_assign[x]:", true_assign[x], "\n")
                         # cat("\n mu_true: \n")
                         # print(mu_true[[true_assign[x]]])
                         # cat("\n Sigma_true: \n")
                         # print(Sigma_true[[true_assign[x]]])
                         mvtnorm::dmvnorm(x = y[[x]][,1], 
                                          mean = mu_true[[true_assign[x]]], 
                                          sigma = Sigma_true[[true_assign[x]]])
                         
                       })

    
      
    # calculate density of estimates
    est_dens = vector(mode = "list", length = length(mu_est))
    for(k in 1:length(mu_est)){
    cat("\n k=",k,"\n")
      est_dens_k = matrix(data = NA, nrow = nrow(mu_est[[k]]), ncol = length(y))
      for(iter in 1:nrow(mu_est[[k]])){
        # cat("\n iter=",iter,"\n")
        # cat("\n Sigma_est[[k]]=")
        # print(Sigma_est[[k]][iter,])
        if(off_diag == TRUE){
          
          est_dens_k[iter,] = sapply(X = 1:length(y), 
                                     FUN = function(x){
                                       mean_ind = grep(
                                         pattern = paste0("mu_", group_assign[[k]][iter,x], "_"),
                                         x = colnames(mu_est[[k]]))
                                       p = length(mean_ind)
                                       var_ind = grep(
                                         pattern = paste0("sigma_", group_assign[[k]][iter,x]),
                                         x = colnames(Sigma_est[[k]]))
                                       cov_ind = var_ind[(p+1):length(var_ind)]
                                       est_var_k = diag( 
                                         # first diag variance part
                                         Sigma_est[[k]][iter, var_ind[1:p]]
                                         ) # then off diag part
                                       
                                       offdiag_index = lapply(X = stringr::str_split(string = stringr::str_remove_all(
                                         string = unlist(stringr::str_extract_all(
                                           string = colnames(Sigma_est[[k]])[var_ind][(p+1):length(var_ind)],
                                           pattern = "_[:digit:]+_[:digit:]$")),
                                         pattern = "^_"), pattern = "_"), FUN = as.numeric)
                                       
                                       for(entry_index in 1:length(cov_ind)){
                                         entry = cov_ind[entry_index]
                                         entry_row = offdiag_index[[entry_index]][1]
                                         entry_col = offdiag_index[[entry_index]][2]
                                         est_var_k[entry_row, entry_col] = Sigma_est[[k]][iter,entry] # pull correct element
                                         # from var_by_k results and use to reconstruct covariance matrix
                                         est_var_k[entry_col, entry_row] = Sigma_est[[k]][iter,entry] # also do reverse for
                                         # symmetry of covariance matrix
                                       }
                                       
                                       # cat("\n mu_est")
                                       # print(as.numeric(mu_est[[k]][iter,mean_ind]))
                                       # cat("\n Sigma_est")
                                       # print(Sigma_est[[k]][iter,var_ind])
                                       mvtnorm::dmvnorm(x = y[[x]][,1], 
                                                        mean = as.numeric(mu_est[[k]][iter,mean_ind]), 
                                                        sigma = est_var_k)
                                       
                                     })
          
          
        } else{
          # off_diag == FALSE, do not need to consider off-diagonal variance components
          
          est_dens_k[iter,] = sapply(X = 1:length(y), 
                                     FUN = function(x){
                                       # cat("\n x ", x)
                                       # cat("\n iter ", iter)
                                       # cat("\n group_assign[iter,x] ", group_assign[[k]][iter,x])
                                       mean_ind = grep(
                                         pattern = paste0("mu_", group_assign[[k]][iter,x], "_"),
                                         x = colnames(mu_est[[k]]))
                                       # cat("\n mean_ind ", mean_ind)
                                       var_ind = grep(
                                         pattern = paste0("sigma_", group_assign[[k]][iter,x]),
                                         x = colnames(Sigma_est[[k]]))
                                       #cat("\n", dim(y[[x]][,1]))
                                       #cat("\n mu_est ", length(as.numeric(mu_est[[k]][iter,mean_ind])))
                                       #print(as.numeric(mu_est[[k]][iter,mean_ind]))
                                       # cat("\n Sigma_est")
                                       # print(Sigma_est[[k]][iter,var_ind])
                                       mvtnorm::dmvnorm(x = y[[x]][,1], 
                                                        mean = as.numeric(mu_est[[k]][iter,mean_ind]), 
                                                        sigma = diag(as.numeric(Sigma_est[[k]][iter,var_ind]), 
                                                                     length(y[[x]][,1]))
                                       )
                                     })
        }


        
      }

      est_dens[[k]] = est_dens_k # save results for kth element to list
      
    }
    
    # clean up
    est_dens_combined = est_dens[[1]] 
    if(length(est_dens) > 1){
      for(k in 2:length(est_dens)){
        est_dens_combined = rbind(est_dens_combined, est_dens[[k]])
        # print(dim(est_dens[[k]]))
      }
    }
    
    # average over all iterations
    final_est_dens = colMeans(est_dens_combined)
    
  } else{
    # pooled variance, assumed equal across groups
    
    # calculate density of truth
    true_dens = sapply(X = 1:length(y), 
                       FUN = function(x){
                         mvtnorm::dmvnorm(x = y[[x]][,1], 
                                          mean = mu_true[[true_assign[x]]], 
                                          sigma = Sigma_true)
                         
                       })
  
    # calculate density of estimates
    est_dens = vector(mode = "list", length = length(mu_est))
    for(k in 1:length(mu_est)){
      # cat("\n k:", k, "\n")
      est_dens_k = matrix(data = NA, nrow = nrow(mu_est[[k]]), ncol = length(y))
      for(iter in 1:nrow(mu_est[[k]])){
        
        # cat("\n iter:", iter, "\n")
        est_dens_k[iter,] = sapply(X = 1:length(y), 
                                   FUN = function(x){
                                     # cat("\n x ", x)
                                     # cat("\n iter ", iter)
                                     # cat("\n group_assign[iter,x] ", group_assign[[k]][iter,x])
                                     mean_ind = grep(
                                       pattern = paste0("mu_", group_assign[[k]][iter,x], "_"), 
                                       x = colnames(mu_est[[k]]))
                                     # var_ind = grep(
                                     #   pattern = paste0("sigma", group_assign[[k]][iter,x]), 
                                     #   x = names(Sigma_est[[k]]))
                                     # cat("\n mean_ind", mean_ind, "\n")
                                     # cat("\n obs", x, "\n")
                                     # cat("\n y[[x]] \n")
                                     # print(y[[x]][,1])
                                     # cat("\n mu \n")
                                     # print(as.numeric(mu_est[[k]][iter,mean_ind]))
                                     # cat("\n sigma \n")
                                     # print(diag(as.numeric(Sigma_est[[k]][iter,1]), 
                                     #            length(y[[x]][,1])))
                                     mvtnorm::dmvnorm(x = y[[x]][,1], 
                                                      mean = as.numeric(mu_est[[k]][iter,mean_ind]), 
                                                      # equal var assumption so Sigma_set should only have 1 column
                                                      # and there should be only 1 list element after most recent update
                                                      sigma = diag(as.numeric(Sigma_est[[1]][iter,1]), 
                                                                   length(y[[x]][,1]))
                                                      )
                                   })
        
      }
      
      est_dens[[k]] = est_dens_k # save results for kth element to list
      
    }
    
    # clean up
    est_dens_combined = est_dens[[1]] 
    if(length(est_dens) > 1){
      for(k in 2:length(est_dens)){
        est_dens_combined = rbind(est_dens_combined, est_dens[[k]])
      }
    }
    
    # average over all iterations
    final_est_dens = colMeans(est_dens_combined)
  
  }
  
  # now do we average over densities in estimates then calc KL div, or  calculate KL divergence
  # at each iteration and then average over all iterations?
  
  # option 1 - average over densities, then calculate KL - YES
  # option 2 - calculate KL divergence at each iteration - NO, look at overall
  # estimated density vs truth
  # print(final_est_dens)
  # print(true_dens)
  kl_div = LaplacesDemon::KLD(px = final_est_dens, py = true_dens)
  
  return(kl_div)
}


# unit test for KL div - 2d with 2 groups
# perfect agreement -- should have KL div of 0
# Sigma_est = list(matrix(c(10,10,0.5,5,5,0.25,
#                           10,10,0.5,5,5,0.25), ncol=6, byrow=TRUE))
# colnames(Sigma_est[[1]]) = c("sigma_1_1_1","sigma_1_2_2", "sigma_1_1_2", "sigma_2_1_1","sigma_2_2_2", "sigma_2_1_2")
# mu_est = list(matrix(c(0,0,10,10,0,0,10,10), ncol=4, byrow=TRUE))
# colnames(mu_est[[1]]) = c("mu_1_1", "mu_1_2", "mu_2_1", "mu_2_2")
# calc_KL_diverg(y = list(matrix(c(0,0), ncol=1),
#                         matrix(c(10.5,10.5), ncol=1),
#                         matrix(c(0.5,0.5), ncol=1), 
#                         matrix(c(8.5,8.5), ncol=1)), 
#                 mu_est = mu_est, 
#                 Sigma_est = Sigma_est,  
#                 group_assign = list(matrix(c(1,2,1,2,1,2,1,2), ncol = 4, byrow=TRUE)), 
#                 true_assign = matrix(c(1,2,1,2), ncol = 4, byrow=TRUE), 
#                 mu_true = list(matrix(c(0,0), ncol=1), matrix(c(10,10), ncol=1)), 
#                 Sigma_true = list(matrix(c(10,0.5,0.5,10), ncol=2, byrow=TRUE),
#                                   matrix(c(5,0.25,0.25,5), ncol=2, byrow=TRUE)),
#                 equal_var_assump = FALSE, 
#                 off_diag = TRUE)




get_adjmat_groupassign = function(group_assign){
  
  # function to generate list of adjacency matrices from group_assign at each iteration
  # helpful when you only have output of dpmm_summary and not full model output but want
  # to do diagnostics using genMCMCdiag
  
  adjmat_list = lapply(X = 1:(nrow(group_assign)), 
                       FUN = function(x){
                         matrix(as.numeric(outer(
                           X = group_assign[x,], 
                           Y = group_assign[x,], 
                           FUN = "==")), nrow = ncol(group_assign))
                       })
  
  return(adjmat_list)
  
}

