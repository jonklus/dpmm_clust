# getting variances into proper format so that they match mean output -- assuming
# diagonal structure

vars = ex3$vars

lapply(X = vars[[501]], FUN = diag)

format_vars <- function(vars, num_vars){
  
  lapply(X = 1:length(num_vars), 
         FUN = function(x){
           for(i in 1:num_vars[x]){
             diag(vars[[x]][[i]])
           }
         })
  
}
# 
# test = lapply(X = 1:length(num_vars), 
#               FUN = function(x){
#                 sapply(X = 1:num_vars[x], FUN = diag(ex3$vars[[]])
#                 for(i in 1:num_vars[x]){
#                   diag(ex3$vars[[x]][[i]])
#                 }
#               })

temp_param_list = ex3$vars
param_list = vector(mode = "list", length = length(num_vars))
for(i in 1:length(num_vars)){
  
  param_list[[i]] = sapply(X = 1:num_vars[i], 
                           FUN = function(x){
                             diag(ex3$vars[[i]][[x]])
                             })
  
}


# list_means_by_k 

means = ex3$means
burn_in = 100
relabel = TRUE
permutation = stephens_result3
  ## Function takes a list of length no. MCMC iterations and outputs a list of 
  ## data frames where each list element is the contains the draws for iterations
  ## That found a specified number of groups k. columns in output are each parameter.
  
  # means in a list of length # of MCMC iterations where each list element
  # contains a matrix of MCMC output with dimension (no. obs)*(no. groups) = n*k
  # burn-in is the number of initial MCMC iterations to discard
  # relabel is logical, if TRUE then use permutation to relabel means to solve
  # the label switching problems
  # permutation is a list sorted by k (# groups) resulting from the label.switching
  # package
  
  mean_list = means[(burn_in+1):length(means)]
  num_means = sapply(X = 1:length(mean_list),
                     FUN = function(x){
                       ncol(mean_list[[x]])
                     })
  unique_k = sort(unique(num_means))
  #print(unique_k)
  npar = nrow(mean_list[[1]]) # number of parameters per group, does not change with k
  #print(npar)
  # need to create separate object for each # of groups k 
  mean_list_by_k = vector(mode = "list", length = length(unique_k))
  
  # if label switching solution is given, also reorder means to address this
  if(relabel == FALSE){
    
    for(i in 1:length(unique_k)){
      k_index = which(num_means == unique_k[i]) # indices of all iters with k means
      mean_list_by_k[[i]] = data.frame(matrix(data = unlist(mean_list[k_index]), 
                                              ncol = npar*unique_k[i],
                                              byrow = TRUE))
      # name columns i.e. mu23 is mean for group 2, 3rd component 
      # (i.e. from a length 3 mean vector)
      col_header_names = unlist(lapply(X = 1:npar, 
                                       FUN = function(x){
                                         paste0("mu", 1:unique_k[i], x)
                                       }))
      names(mean_list_by_k[[i]]) = sort(col_header_names)
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
    if(1 %in% unique_k){
      # unique_k = unique_k[-1] # drop the 1 before correcting -- no permutations
      i_init = 2  # start after k=1
      permutation = c(1, permutation) # put placeholder in for k=1 group (not 
      # included in stephens results bc no label switching) -- avoids error with
      # indexing later on in for loop otherwise index doesn't match in means & perm

      # put the k=1 group names in here 
      k_index = which(num_means == 1)
      mean_list_by_k[[1]] = data.frame(matrix(data = unlist(mean_list[k_index]),   
                                              ncol = npar,                   
                                              byrow = TRUE))
      # name columns i.e. mu23 is mean for group 2, 3rd component 
      # (i.e. from a length 3 mean vector)
      names(mean_list_by_k[[1]]) = paste0("mu", 1, 1:npar)
      
    } else{
      
      i_init = 1 # business as usual
      
    }
    
    for(i in i_init:length(unique_k)){
      
      #print(unique_k[i])
      
      k_index = which(num_means == unique_k[i]) # indices of all iters with k means
      new_labs = permutation[[i]]$permutations$STEPHENS
      
      # fix label switching before reformatting means
      for(j in 1:length(k_index)){
        
        #print(mean_list[[k_index[j]]]) # before fix
        #print(new_labs[j,]) # new order
        
        mean_list[[k_index[j]]] =  mean_list[[k_index[j]]][,new_labs[j,]]
        
        # print(mean_list[[k_index[j]]]) # after fix
      }

      # reformat means now that order has been corrected
      mean_list_by_k[[i]] = data.frame(matrix(data = unlist(mean_list[k_index]), 
                                              ncol = npar*unique_k[i], 
                                              byrow = TRUE))
      # name columns i.e. mu23 is mean for group 2, 3rd component 
      # (i.e. from a length 3 mean vector)
      col_header_names = unlist(lapply(X = 1:npar, 
                                       FUN = function(x){
                                         paste0("mu", 1:unique_k[i], x)
                                       }))
      names(mean_list_by_k[[i]]) = sort(col_header_names)
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
    
  }