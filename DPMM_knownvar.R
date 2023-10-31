##############################################################
################# SIMPLE DPMM MODEL ##########################
#################     VERSiON 2     ##########################
##############################################################
## Author: Jonathan Klus
## Date: 18 January 2023
## Description: Sampler for a mixture of normal densities using Algorithm 2
## from Neal (2000). We assume the variance is known so that we have conjugacy
## between the normal likelihood and dirichlet process prior. 

###################### load packages and set seed #############################
# set.seed(516)
# library(ggplot2)

################## simulate data ##############################################
# w = c(0.5, 0.3, 0.2)
# means = c(-20, 20, 0)
# var = rep(5, 3)
# 
# assign = sample(x = 1:3, size = 20, replace = TRUE, prob = w)
# data = data.frame(
#   assign = assign,
#   y = sapply(X = assign, 
#              FUN = function(x){
#                rnorm(n = 1, mean = means[x], sd = sqrt(var[x]))
#              })
# )
# 
# ggplot(data = data, aes(x = y)) +
#   geom_histogram(binwidth = 2) +
#   theme_classic()
# 
# # drop params used to generate data, except for var which is assumed known
# rm(w, means)

########################### Gibbs sampler #####################################

# helper functions

## calculate group membership probabilities

group_prob_calc <- function(k, n, n_j, alpha, y_i, mu, sigma2, tau2, mu0, singleton = 0, curr_group_assign = NULL, curr_labels = NULL){
  # k is the number of existing groups
  # n is total number of observations
  # n_j is a vector of length k with the total number of observations in each group
  # alpha is the scalar concentration parameter from the dirichlet process
  # y_i is the single data point we are considering
  # mu is a vector of length k
  # mu0 and tau2 are scalar priors
  # sigma2 is a scalar
  
  # calculate unnormalized probabilities of group membership for obs i
  
  #### probability of joining a current group
  
  if(singleton == 1){
    
    pr_curr = sapply(X = 1:k, 
                     FUN = function(x){
                       c1 = n_j[x]/(n-1+alpha)
                       c2 = dnorm(x = y_i, mean = mu[x], sd = sqrt(sigma2))
                       return(c1*c2)
                     })
    
  } else{
    
    pr_curr = sapply(X = 1:k, 
                     FUN = function(x){
                       #### make sure you're calculating n_{-i, c}
                       if(curr_group_assign == curr_labels[x]){
                         c1 = (n_j[x]-1)/(n-1+alpha)
                         c2 = dnorm(x = y_i, mean = mu[x], sd = sqrt(sigma2))
                       } else{
                         c1 = n_j[x]/(n-1+alpha)
                         c2 = dnorm(x = y_i, mean = mu[x], sd = sqrt(sigma2))
                       }
                       return(c1*c2)
                     })
    
  }
  

  #### probability of creating a new group
  pr_new = alpha/(n-1+alpha)*(1/sqrt(2*pi*(sigma2 + tau2)))*exp(-(1/2)*(1/(sigma2+tau2))*(y_i-mu0)^2) 
  

  
  # normalize probs to account for "b"
  pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
  
  return(pr_c)
  
}

# split_merge_MH <- function(move_type, n_1, n_2, posterior){
#   # function to calculate the acceptance probability of a split or merge
#   
#   q_ratio = (1/2)^(n_1 + n_2 - 2)
#   if(move_type == "split"){
#     q_ratio = 1/q_ratio
#   } else{
#     q_ratio = q_ratio
#   }
# }


fit_DPMM_knownvar <- function(y, alpha = 1, mu0 = 0, tau2, S = 10^3, k = 2, verbose = TRUE){
  # y is data
  # alpha is concentration parameter of DP prior
  # mu0 and tau2 are prior mean and variance
  # S is the number of MCMC iterations
  # k is the initial number of groups
  # verbose = TRUE will print summary at every 100th iteration
  
  # initialize group assignments and means
  n = length(y)
  mu = mu0
  sigma2 = var[1]
  k = 2 # number of groups
  
  # preallocate memory
  label_record = vector(mode = "list", length = S)
  group_assign = matrix(data = NA, nrow = S, ncol = n)
  # group_assign[1, ] = sample(x = 1:k, size = length(y), replace = TRUE, prob = rep(1/k, k))
  # try different group assign initialization
  group_assign[1, ] = ifelse(y > mean(y), k, k-1)
  label_record[[1]] = c(k, k-1)
  means = matrix(data = NA, nrow = S, ncol = n)
  means[1, 1:k] = sapply(X = 1:k, FUN = function(x){mean(y[group_assign[1,] == x])})
  num_groups = matrix(data = NA, nrow = S, ncol = 1)
  mu = means[1, 1:k]
  curr_labels = 1:k
  avail_labels = c(1:n)[-curr_labels]
  
  # iterate 1:S
  for(s in 2:S){
    
    #print("A", mu)
    
    ## initialize group assignments for current iteration using ending state from prior iteration
    group_assign[s, ] = group_assign[s-1, ]
    
    ## consider split/merge step for conjugate case (Jain & Neal 2004)
    
    # ### select two observations uniformly at random 
    # splitmerge_cand = sample(x = group_assign[s,], size = 2, replace = FALSE, prob = rep(1/n, n))
    # 
    # ### if selected observations belong to the same mixture component: PROPOSE SPLIT OF THEIR GROUPS
    # ### if selected observations do not belong to the same mixture component: PROPOSE MERGE OF THEIR GROUPS
    # if(splitmerge_cand[1] == splitmerge_cand[2]){ 
    #   #### SPLIT PROPOSAL
    #   c1 = splitmerge_cand[1]
    #   c2 = avail_labels[1] # new group
    #   y_c = y[group_assign[s,] %in% splitmerge_cand]
    #   
    #   
    # } else{ 
    #   #### MERGE PROPOSAL
    # 
    #   
    # }
    
    
    
    ## iterate through observations 1:n 
    for(i in 1:n){
      
      ### check if observation i is a singleton
      count_assign = as.numeric(table(group_assign[s,]))
      label_assign = as.numeric(names(table(group_assign[s,])))
      singletons = label_assign[which(count_assign == 1)]
      
      ### debugging
      # print(count_assign)
      # print(label_assign)
      # print(curr_labels)
      # print(i)
      
      ### if observation i is a singleton, remove its mean from the current state of system
      if(group_assign[s,i] %in% singletons){
        
        #### only drop observation i if it is a singleton...DO NOT drop other singleton
        #### observations at this point!!!
        singleton_index = which(label_assign == group_assign[s,i]) 
        mu = mu[-singleton_index]
        count_assign = count_assign[-singleton_index]
        label_assign = label_assign[-singleton_index]
        avail_labels = c(group_assign[s,i], avail_labels)
        curr_labels = curr_labels[-which(curr_labels %in% group_assign[s,i])]
        k = length(curr_labels) # update k 
        
        #*** also need to functionalize this part after debugging
        ### for any observation i, draw a group assignment
        pr_c = group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[i], mu = mu, sigma2 = sigma2, tau2 = tau2,
                               mu0 = mu0, singleton = 1)
        group_assign[s,i] = sample(x = c(curr_labels, avail_labels[1]), size = 1, prob = pr_c)
        
        #### if new group selected
        if(group_assign[s,i] == avail_labels[1]){
          
          #### handle bookkeeping
          curr_labels = c(curr_labels, avail_labels[1])
          avail_labels = avail_labels[-1]
          k = length(curr_labels)
          
          #### draw a mean for newly created group
          b = y[i]/sigma2 + mu0/tau2
          a = 1/sigma2 + 1/tau2
          mu_k = rnorm(n = 1, mean = b/a, sd = sqrt(1/a)) # make this th kth mean
          mu = c(mu, mu_k)
        }
      } else{
        #### if obs i is not presently a singleton
        
        
        pr_c = group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[i], mu = mu, sigma2 = sigma2, tau2 = tau2,
                               mu0 = mu0, singleton = 0, curr_group_assign = group_assign[s,i], 
                               curr_labels = curr_labels)
        group_assign[s,i] = sample(x = c(curr_labels, avail_labels[1]), size = 1, prob = pr_c)
        
        #### if new group selected
        if(group_assign[s,i] == avail_labels[1]){
          
          #### handle bookkeeping
          curr_labels = c(curr_labels, avail_labels[1])
          avail_labels = avail_labels[-1]
          k = length(curr_labels)
          
          #### draw a mean for newly created group
          b = y[i]/sigma2 + mu0/tau2
          a = 1/sigma2 + 1/tau2
          mu_k = rnorm(n = 1, mean = b/a, sd = sqrt(1/a)) # make this th kth mean
          mu = c(mu, mu_k)
          
        }
        
        
      }
      
      ### iterate through all y_i
      
    } ### end iterations from i=1:n
    
    #print("B", mu)
    
    # final update of counts after a sweep
    label_record[[s]] = curr_labels
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    num_groups[s,] = k
    
    # print results after each sweep
    
    # draw group means for all k groups
    sum_y_i = sapply(X = 1:k, FUN = function(x){sum(y[group_assign[s,] == label_assign[x]])})
    b = sum_y_i/sigma2 + mu0/tau2
    a = count_assign/sigma2 + 1/tau2
    mean_param = b/a
    var_param = 1/a
    mu = rnorm(n = k, mean = mean_param, sd = sqrt(var_param))
    means[s,1:k] = mu
    
    # print progress
    if((verbose == TRUE) & (s %% 100 == 0) & (s > 100)){
      
      print(paste("iter = ", s))
      print(paste("Current k = ", k))
      print(table(group_assign[s,]))
      print(mu)
      # print(c("Current labels: ", curr_labels))
      # print(avail_labels)
      
    }
    
    
  } ### end MCMC iterations from s=1:S
  
  ## check & print no. of groups over all MCMC iterations
  print(table(sapply(X = 1:S, FUN = function(x){length(as.numeric(table(group_assign[x,])))})))
  
  # return results
  MCMC_res = list(
    group_assign = group_assign, # matrix of group assignments from each iteration
    means = means, 
    label_record = label_record
  )
  
  return(MCMC_res)
  
}
  
  

########################### Post processing #####################################
post_process_univar <- function(group_assign, means, S, label_record){
  # group_assign is the the S*n matrix of mcmc group assignments
  # label_record is the S length list of labels used at each iteration
  # means is the S*n matrix of group means, with NAs to keep dimensions for groups that don't exist
  # S is the number of MCMC iterations
  
  relabeled_assign = group_assign
  reordered_means = vector(mode = "list", length = S)
  n = ncol(group_assign)
  for(s in 1:S){
    iter_s_order = order(means[s,is.na(means[s,])==FALSE])
    reordered_means[[s]] = means[s,iter_s_order]
    labels_s = label_record[[s]] #as.numeric(names(table(group_assign[s,])))
    # now relabel based on new order
    for(x in iter_s_order){
      for(i in 1:n){
        if(group_assign[s,i] == labels_s[x]){
          relabeled_assign[s,i] = x
        }
      }
    }
  }
  
  processed_results = list(
    relabeled_assign = relabeled_assign, 
    reordered_means = reordered_means
  )
  
  return(processed_results)
  
}


############################## Test functions ##################################

# modfit = fit_DPMM_knownvar(y = data$y, alpha = 1, mu0 = 0, tau2 = 20, S = 10^3, k = 2, verbose = TRUE)
# postprocess_mod = post_process_univar(modfit$group_assign, modfit$means, S = 10^3, modfit$label_record)
  

# table(sapply(X = 1:S, FUN = function(x){length(as.numeric(table(test[x,])))}))
# 
# 
# 
# p1 = ggplot(data = data.frame(y=y, label = factor(group_assign[100,])), aes(x = y)) +
#   geom_histogram(binwidth = 2, aes(fill = label)) +
#   #geom_text(aes(x = reordered_means[[100]], y = 2, label = reordered_means[[100]])) + 
#   theme_classic()
# 
# p2 = ggplot(data = data.frame(y=y, label = factor(relabeled_assign[100,])), aes(x = y)) +
#   geom_histogram(binwidth = 2, aes(fill = label)) +
#   theme_classic()
# 
# #gridExtra::arrangeGrob(p1, p2)
# p1
# p2
# 
# 
