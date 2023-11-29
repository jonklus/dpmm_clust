##############################################################
################# MULTIVARIATE DPMM MODEL ####################
#################  VAR UNKNOWN - ALG 2    ####################
################# INDEP IG PRIORS - DEV   ####################
##############################################################

## Author: Jonathan Klus
## Date: 6 June 2023
## Description: Sampler for a mixture of multivariate normal densities using Algorithm 2
## from Neal (2000). We assume the variance is unknown, but we impose conjugacy
## by assuming the variance of the mean and likelihood differ only be a multiplicative
## parameter r. 

###################### load packages and set seed ##############################
set.seed(516)
library(ggplot2)
library(LaplacesDemon)

############################### HELPER FUNCTIONS ###############################

## calculate group membership probabilities

group_prob_calc <- function(k, n, n_j, alpha, y_i, mu, sigma2, r, a, b, mu0, 
                            singleton = 0, curr_group_assign = NULL, curr_labels = NULL){
  # k is the number of existing groups
  # n is total number of observations
  # n_j is a vector of length k with the total number of observations in each group
  # alpha is the scalar concentration parameter from the dirichlet process
  # y_i is the single data vector of dimension p*1 that we are considering
  # mu is a matrix with k columns, contains k group mean vectors of dimension p*1
  # sigma2 is either a scalar or a list of within-group variances, assumption dependent
  # mu0 is a p*1 prior mean vector, tau2 is a scalar prior
  # sigma0 is a scalar, corresponds to a multiplier on the variance ie. mu ~ N(mu0, sigma0*Sigma)
  # r is a scalar, corresponds to prior variance
  
  # calculate unnormalized probabilities of group membership for obs i
  
  #### probability of joining a current group
  
  p = length(mu[,1])
  
  if(length(sigma2) == 1){
    # EEE case, variance assumed equal across groups
    
    if(singleton == 1){
      
      pr_curr = sapply(X = 1:k, 
                       FUN = function(x){
                         # print("Step1")
                         # print(matrix(mu[,x], nrow = p))
                         # print(y_i)
                         c1 = n_j[x]/(n-1+alpha)
                         c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                               mean = mu[,x], 
                                               sigma = diag(sigma2, p)) 
                         return(c1*c2)
                       })
      
    } else{
      
      pr_curr = sapply(X = 1:k, 
                       FUN = function(x){
                         #### make sure you're calculating n_{-i, c}
                         if(curr_group_assign == curr_labels[x]){
                           # print("Step2")
                           # print(matrix(mu[,x], nrow = p))
                           # print(y_i)
                           c1 = (n_j[x]-1)/(n-1+alpha)
                           c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                                 mean = mu[,x], 
                                                 sigma = diag(sigma2, p)) 
                         } else{
                           # print("Step3")
                           # print(matrix(mu[,x], nrow = p))
                           # print(y_i)
                           c1 = n_j[x]/(n-1+alpha)
                           c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                                 mean = mu[,x],
                                                 sigma = diag(sigma2, p)) 
                         }
                         return(c1*c2)
                       })
      
    }
    
  } else{
    
    # VEE case, variance differs by group but still diagonal within group
    
    if(singleton == 1){
      
      pr_curr = sapply(X = 1:k, 
                       FUN = function(x){
                         # print("Step1")
                         # print(matrix(mu[,x], nrow = p))
                         # print(y_i)
                         c1 = n_j[x]/(n-1+alpha)
                         c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                               mean = mu[,x], 
                                               sigma = diag(sigma2[[x]], p)) 
                         return(c1*c2)
                       })
      
    } else{
      
      pr_curr = sapply(X = 1:k, 
                       FUN = function(x){
                         #### make sure you're calculating n_{-i, c}
                         if(curr_group_assign == curr_labels[x]){
                           # print("Step2")
                           # print(matrix(mu[,x], nrow = p))
                           # print(y_i)
                           c1 = (n_j[x]-1)/(n-1+alpha)
                           c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                                 mean = mu[,x], 
                                                 sigma = diag(sigma2[[x]], p)) 
                         } else{
                           # print("Step3")
                           # print(matrix(mu[,x], nrow = p))
                           # print(y_i)
                           c1 = n_j[x]/(n-1+alpha)
                           c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                                 mean = mu[,x],
                                                 sigma = diag(sigma2[[x]], p)) 
                         }
                         return(c1*c2)
                       })
      
    }
    
  }
  
  
  
  
  #### probability of creating a new group
  #### now a multivariate t distribution --- need to fix this and add required
  #### elements from prior on mean, variance
  pr_new = alpha/(n-1+alpha)*LaplacesDemon::dmvt(x = y_i[,1], 
                                                 mu = mu0[,1], 
                                                 S = diag((b/a)*(r+1), p),
                                                 df = 2*a)
  
  #### normalize probs to account for "b"
  pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
  
  return(pr_c)
  
}

prior_pred_NinvGa <- function(y_i, mu0, r, a, b){
  
  dens = LaplacesDemon::dmvt(x = y_i[,1], 
                             mu = mu0[,1], 
                             S = diag((r+1)*(b/a), length(y_i[,1])),
                             df = 2*a)
  
  return(dens)
}



post_pred_DEV <- function(obs, which_group, group_assign, split_labs, r, sm_counts, a, b, y, ybar, loss_ybar, mu0){
  
  # see BDA 3/e p. 73 for MVN unknown mean and var
  # or PDF in references folder
  
  # restricted gibbs sampling scans
  
  sum_ysq = lapply(X = 1:2, 
                     FUN = function(x){
                       
                       col_ind = x  # from outer apply
                       group_ind = which(group_assign == split_labs[x])
                       if(obs %in% group_ind){
                         obs_ind = which(obs == group_ind)
                         group_ind = group_ind[-obs_ind]
                       } # else continue
                       
                       Reduce(f = "+", 
                              x = lapply(X = group_ind, FUN = function(x){
                                t(y[[x]])%*%y[[x]]}))
                       
                     })
  
  loss_mu0 = (ybar[[which_group]] - mu0)%*%t(ybar[[which_group]] - mu0)
  
  mu_n = ((1/r)*mu0 + sm_counts[which_group]*ybar[[which_group]])/((1/r) + sm_counts[which_group])
  
  # cat("obs:", obs)
  # cat("\n")
  # cat("mu_n:", mu_n)
  # cat("\n")
  # cat("dim(mu):", dim(mu_n))
  # cat("\n")
  # cat("y:", y[[obs]])
  # cat("\n")
  # cat("dim(y):", dim(y[[obs]]))
  # cat("\n")
  # cat(r, sm_counts[which_group])
  # cat("\n")
  # print(ybar[[which_group]])
  # cat("\n")
  
  # need to recalculate the posterior predictives yourself...not super confident
  # in PDF this came from
  
  a_n = a + sm_counts[which_group]/2
  
  b_n = b + (t(mu0)%*%mu0/r + sum_ysq[[which_group]] - (1/r + sm_counts[which_group])*t(mu_n)%*%mu_n)/2

  lambda_n = diag(b_n[,1]*(1+(1/r + sm_counts[which_group])^(-1))/a_n, length(mu_n[,1]))
  # take first "column" of scalar b_n snce R was still recognizing as a matrix
  
  # print(mu_n)
  # print(lambda_n)
  # print(nu_n)
  
  val = sm_counts[which_group]*LaplacesDemon::dmvt(x = y[[obs]][,1], 
                                                   mu = mu_n[,1], 
                                                   S = lambda_n, 
                                                   df = 2*a_n)
  
  return(val)
  
}

final_post_pred_UVV <- function(y_i, r, y, mu0, a, b){
  
  # use to calculate values after final gibbs scan 
  ## y_i is observation of interest, y is all data being considered for posterior,
  ## exclusive of the current observation y_i under consideration
  
  ybar = Reduce(f = "+", x = y)/length(y)
  
  loss_ybar = Reduce(f = "+", 
                     x = lapply(X = 1:length(y), 
                                FUN = function(x){
                                  (y[[x]] - ybar)%*%t(y[[x]] - ybar)
                                })
  )
  
  sm_counts = length(y)
  
  loss_mu0 = (ybar - mu0)%*%t(ybar - mu0)
  
  sum_ysq = Reduce(f = "+", x = lapply(X = 1:length(y), 
                                       FUN = function(x){
                                         t(y[[x]])%*%y[[x]]
                                         }))
                     
  loss_mu0 = (ybar - mu0)%*%t(ybar - mu0)
  
  mu_n = ((1/r)*mu0 + sm_counts*ybar)/((1/r) + sm_counts)
  a_n = a + sm_counts/2
  b_n = b + (t(mu0)%*%mu0/r + sum_ysq - (1/r + sm_counts)*t(mu_n)%*%mu_n)/2
  
  lambda_n = diag(b_n[,1]*(1+(1/r + sm_counts)^(-1))/a_n, length(mu_n[,1]))

  val = sm_counts*LaplacesDemon::dmvt(x = y_i, 
                                      mu = mu_n[,1], 
                                      S = lambda_n, 
                                      df = 2*a_n)
  
  return(val)
  
  
}

ll_components <- function(subset_index, obs_ind, y, mu0, lambda0, r, nu){
  # function to calculate the components of the likelihood ratio for the MH
  # acceptance probability after n_iter split merge restricted Gibbs scans
  
  if(obs_ind == 1){
    # first observation --- prior predictive
    val = prior_pred_NinvGa(y_i = y[[subset_index[obs_ind]]], 
                           mu0 = mu0, r = r, 
                           a = a, b = b)
  } else{
    # posterior predictive
    # need to calculate for all obs, in same group
    subset_yvals = y[subset_index[1:obs_ind]] # what y values to include at each iter
    val = final_post_pred_DEV(y_i = y[[subset_index[obs_ind]]], r = r, 
                              y = y[subset_index[1:(obs_ind-1)]], 
                              mu0 = mu0, a = a, b = b)
    

  }
  
  return(val)  
  
}

# test line
# post_pred_UVV(obs=2,which_group=1, r=10, sm_counts=c(10,11), nu=2, y=y, ybar=ybar, 
#               loss_ybar=loss_ybar, mu0=matrix(data=0,nrow=2), lambda0=diag(10,2))
######################

split_merge_prob_DEV <- function(obs, split_labs, group_assign, r, a, b, y, mu0){
  # split_labs is an array of length 2 indicating which entries in counts correspond
  # to the groups that are part of the split/merge
  # which_group is a scalar valued 1 or 2 indicating which of the groups we are considering
  # obs is the index for the observation being considered at this point
  # group_assign is an array of length n corresponding to group assignments for each obs
  # r, a, and b are scalar hyperparameters
  # counts is an array of the number of obs assigned to each group
  # y is the data
  # mu0 is the prior mean 
  
  sm_counts = sapply(X = split_labs, FUN = function(x){sum(group_assign[-obs] == x)})
  
  # cat("\n")
  # cat("sm_counts:", sm_counts)
  # cat("\n")
  # print(group_assign)
  
  ybar = lapply(X = split_labs, 
                FUN = function(x){
                  group_ind = which(group_assign == x)
                  if(obs %in% group_ind){
                    obs_ind = which(obs == group_ind)
                    group_ind = group_ind[-obs_ind]
                  } # else continue
                  
                  ysum = Reduce(f = "+", 
                                x = lapply(X = group_ind, 
                                           FUN = function(x){(y[[x]])}))
                  
                  return(ysum/length(group_ind))
                })
  
  # cat("ybar:")
  # print(ybar)
  # cat("\n")
  
  loss_ybar = lapply(X = 1:2, 
                     FUN = function(x){
                       
                       col_ind = x  # from outer apply
                       group_ind = which(group_assign == split_labs[x])
                       if(obs %in% group_ind){
                         obs_ind = which(obs == group_ind)
                         group_ind = group_ind[-obs_ind]
                       } # else continue
                       
                       Reduce(f = "+", 
                              x = lapply(X = group_ind, FUN = function(x){
                                (y[[x]] - ybar[[col_ind]])%*%t(y[[x]] - ybar[[col_ind]])}))
                       
                     })
  
  # cat("loss_ybar:")
  # print(loss_ybar)
  # cat("\n")
  
  
  
  
  # ratio = sapply(X = 1:2,
  #                FUN = function(x){
  #                  
  #                  xx = ifelse(x==1,2,1)
  #                  num = post_pred_UVV(obs = obs, which_group = x, r = r, # which group is which????
  #                                      sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
  #                                      loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
  #                  denom = num + post_pred_UVV(obs = obs, which_group = xx, r = r,
  #                                              sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
  #                                              loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
  #                  return(num/denom)
  #                  
  #                })
  
  if(0 %in% sm_counts){
    # if there is a singleton, take action to prevent issues with ybar and loss_ybar
    # need to use prior predictive instead of posterior predictive for this group
    
    which_zero = which(sm_counts == 0)
    
    if(which_zero == 1){
      
      num = prior_pred_NinvGa(y_i = y[[obs]], mu0 = mu0, r = r, 
                             a = a, b = b)
      
      denom = num + post_pred_DEV(obs = obs, which_group = 2, r = r, group_assign = group_assign,
                                  sm_counts = sm_counts, y = y, ybar = ybar, 
                                  loss_ybar = loss_ybar, mu0 = mu0, a = a, b = b)
      
      ratio = c(num/denom, 1-(num/denom))
      
    } else{ 
      # which_zero == 2
      
      num = post_pred_DEV(obs = obs, which_group = 1, r = r, 
                          group_assign = group_assign,
                          sm_counts = sm_counts, y = y, ybar = ybar, 
                          loss_ybar = loss_ybar, mu0 = mu0, a = a, b = b)
      
      denom = num + prior_pred_NinvGa(y_i = y[[obs]], mu0 = mu0, r = r, 
                                      a = a, b = b)
      
      ratio = c(num/denom, 1-(num/denom))
      
    }
    
    
    
  } else{
    # proceed as usual 
    
    num = post_pred_DEV(obs = obs, which_group = 1, r = r, 
                        group_assign = group_assign,
                        sm_counts = sm_counts, y = y, ybar = ybar, 
                        loss_ybar = loss_ybar, mu0 = mu0, a = a, b = b)
    
    denom = num + post_pred_DEV(obs = obs, which_group = 2, r = r,
                                group_assign = group_assign,
                                sm_counts = sm_counts, y = y, ybar = ybar, 
                                loss_ybar = loss_ybar, mu0 = mu0, a = a, b = b)
    
    ratio = c(num/denom, 1-(num/denom))
    
  }
  
  return(ratio)
  
}

############################ INDEPENDENT IG PRIORS ############################# 

MVN_CRP_sampler_DEV <- function(S = 10^3, seed = 516, y, r = 2, alpha = 1, a = 1/2, b = 10, mu0, k_init = 2,
                                d = 1, f = 1, g = 1, h = 1, sigma_hyperprior = TRUE, fix_r = FALSE,
                                diag_weights = FALSE, verbose = TRUE, print_iter = 100){
  
  # S is number of MCMC iterations
  # y is a list of data of length n
  # Sigma0 is the prior variance
  # alpha is the dirichlet process concentration parameter
  # a,b are prior variance hyperparameter on the IG distribution of sigma2. b is initial
  # value for b if sigma_hyperprior option is used
  # d, f are the hyperprior parameters on the Gamma dist placed on b, assume a known
  # sigma_hyperprior is a logical argument of whether to use the hyperprior or assume indep.
  # mu0 is the prior mean - must be of dimension p*1
  # K_init is the initial number of groups
  # diag_weights is an argument to the Laplacian matrix - whether the diagonal should be 1 or not
  # verbose & printmod tells the function if and how often to print a progress summary
  
  set.seed(seed = seed)
  
  # data is y - a list of length n
  n = length(y)
  p = length(y[[1]]) # dimensionality of MVN
  k = k_init # initial number of groups
  
  # preallocate memory and set initial values
  accept_ind = matrix(data = NA, nrow = S, ncol = n)
  group_assign = matrix(data = NA, nrow = S, ncol = n)
  
  group_assign[1, ] = sample(x = 1:k, size = length(y), replace = TRUE, prob = rep(1/k, k))
  # try different group assign initialization
  # group_assign[1, ] = ifelse(y > mean(y), k, k-1)  doesn't work for MVN, try kmeans?
  
  means = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  vars = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  probs = vector(mode = "list", length = S)  # group assignment probabiltiies
  
  emp_means = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  emp_vars = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  
  if(sigma_hyperprior == TRUE & fix_r == FALSE){
    extra_params = matrix(data = NA, nrow = S, ncol = 2)
    colnames(extra_params) = c("b", "r")
    extra_params[1,1] = b # note that here b has dimension 1, diagonal var
    extra_params[1,2] = r
  } else if(sigma_hyperprior == FALSE & fix_r == FALSE){
    extra_params = matrix(data = NA, nrow = S, ncol = 1)
    colnames(extra_params) = c("r")
    extra_params[1,1] = r
  } else if(sigma_hyperprior == TRUE & fix_r == TRUE){
    extra_params = matrix(data = NA, nrow = S, ncol = 1)
    colnames(extra_params) = c("b")
    extra_params[1,1] = b
  } 
  
  
  # need to find p*1 vector of means based on this list of observed p*1 y_i values
  means[[1]] = sapply(X = 1:k, 
                      FUN = function(x){
                        rowMeans(matrix(unlist(y[group_assign[1,] == x]), nrow = p))
                        # unravel list of p*1 observations, put in matrix, find mean
                      })
  
  vars[[1]] = lapply(X = 1:k,
                     FUN = function(x){
                       mean(apply(X = matrix(unlist(y[group_assign[1,] == x]), nrow = p),
                                  MARGIN = 1, FUN = var))
                       # collapse Sigma vector initially - start by assuming diagonal & equal
                       # for convenience
                       # if you want matrix instead, replace "mean" with "diag"
                       # unravel list of p*1 observations, put in matrix, find variance
                       # only works for diagonal
                       # for off-diagonal elements use:
                       # cov(t(matrix(unlist(y[group_assign[1,] == x]), nrow = p)))
                       # but note that initial off-diag elements may be extreme
                     })
  
  # initial allocation probs
  if(k == 1){
    probs[[1]] = matrix(data = 0, nrow = n, ncol = 1) 
    probs[[1]][,1] = rep(1, n)
  } else{
    probs[[1]] = matrix(data = 0, nrow = n, ncol = k) 
    probs[[1]][,1:k] = rep(1/k, k)
  }
  
  num_groups = matrix(data = NA, nrow = S, ncol = 1)
  num_groups[1,1] = k_init
  mu = means[[1]]
  sigma2 = vars[[1]]
  curr_labels = 1:k
  avail_labels = c(1:n)[-curr_labels]
  
  # iterate 1:S
  for(s in 2:S){
    
    ## initialize group assignments for current iteration using ending state from prior iteration
    group_assign[s, ] = group_assign[s-1, ]
    
    ## iterate through observations 1:n 
    for(i in 1:n){
      
      #cat("S =", S, " i =", i, "\n")
      
      ### check if observation i is a singleton
      count_assign = as.numeric(table(group_assign[s,]))
      label_assign = as.numeric(names(table(group_assign[s,])))
      singletons = label_assign[which(count_assign == 1)]
      
      ### record current group assignment & parameters for observation i
      ### used to compute acceptance prob in later steps
      curr_assign = which(label_assign == group_assign[s,i]) 
      mu_curr = mu[,curr_assign]
      sigma2_curr = sigma2[[curr_assign]]
      
      # print("Current state")
      # print(curr_assign)
      # print(mu_curr)
      # print(Sigma_curr)
      
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
        # check if there are singletons identified, otherwise do not touch mu
        
        # save current values for acceptance prob
        # mu_curr = mu[,singleton_index]
        # Sigma_curr = Sigma[[singleton_index]]
        
        # remove current values for singleton from index
        mu = matrix(mu[,-singleton_index], nrow = p)
        sigma2 = sigma2[-singleton_index]
        
        count_assign = count_assign[-singleton_index]
        label_assign = label_assign[-singleton_index]
        
        avail_labels = c(group_assign[s,i], avail_labels)
        curr_labels = curr_labels[-which(curr_labels %in% group_assign[s,i])]
        k = length(curr_labels) # update k 
        
        #### calculate proposal distribution for group assignment
        ### for any observation i, calculate group membership probabilities
        pr_c = group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[[i]], mu = mu, sigma2 = sigma2, r = r, 
                               a = a, b = b, mu0 = mu0, singleton = 1)
        
        
      } else{
        
        #### calculate proposal distribution for group assignment
        #### if obs i is not presently a singleton
        pr_c = group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[[i]], mu = mu, sigma2 = sigma2, r = r, 
                               a = a, b = b, mu0 = mu0, singleton = 0, 
                               curr_group_assign = group_assign[s,i], 
                               curr_labels = curr_labels)
        
      }
      
      ### draw a group assignment conditional on group membership probs
      group_assign[s,i] = sample(x = c(curr_labels, avail_labels[1]), 
                                 size = 1, prob = pr_c)
      
      #### if new group selected
      if(group_assign[s,i] == avail_labels[1]){
        
        #### handle bookkeeping
        curr_labels = c(curr_labels, avail_labels[1])
        avail_labels = avail_labels[-1]
        k = length(curr_labels)
        
        
        ### using only the ith observation:
        
        #### draw variance for newly created group from FC posterior of sigma2
        #### according to algo, but this is conditional on mean so do prior for now
        sigma2_k = 1/rgamma(n = 1, shape = a, rate = b)
        sigma2 = c(sigma2, sigma2_k)
        
        #### draw a mean for newly created group from FC posterior of mu 
        mu_mean = (y[[i]] + (1/r)*mu0)/(1/r + 1)
        mu_cov = diag(sigma2_k/(1/r + 1), p)
        mu_k = matrix(mvtnorm::rmvnorm(n = 1, mean = mu_mean, sigma = mu_cov), nrow = p) # kth mean
        mu = cbind(mu, mu_k)
        
      }
      
      
      # print(dim(y[[i]]))
      # print(dim(Sigma_k))
      # print(dim(Sigma_curr))
      # print(c("Proposed group", prop_group_assign, prop_group_assign_index))
      # print(mu)
      # print(Sigma)
      
      
      ### iterate through all y_i
      
      # cat("Labels", curr_labels)
      # cat("Probs", pr_c)
      #### save allocation probabiltiies after each observation i 
      # probs[[s]][i,curr_labels] = pr_c
      
    } ### end iterations from i=1:n
    
    if((s %% print_iter == 0) & (s > print_iter) & (verbose == TRUE)){
      print("End of CRP step") # just create a new line for separation
      print(paste("iter = ", s))
      print(paste("Current k = ", k))
      print(table(group_assign[s,]))
      print(mu)
      print(sigma2)
    }
    
    # final update of counts after a sweep
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    num_groups[s,] = k
    
    # proceed to split-merge step if true
    
    # print results after each sweep
    
    # Proceed to Gibbs step
    
    # draw group means for all K groups
    # Sigma = diag(x = sigma2, nrow = p)
    
    sum_y_i = sapply(X = 1:k, 
                     FUN = function(x){
                       rowSums(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p))
                       # unravel list of p*1 observations, put in matrix, find sum
                     })
    
    mu_cov = lapply(X = 1:k, 
                    FUN = function(x){diag(sigma2[[x]]/(1/r + count_assign[x]), p)}) 
    
    mu_mean = lapply(X = 1:k, 
                     FUN = function(x){(sum_y_i[,x] + mu0/r)/(1/r + count_assign[x])})
    
    mu_list = lapply(X = 1:k, 
                     FUN = function(x){
                       t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                          mean = mu_mean[[x]], 
                                          sigma = mu_cov[[x]]))
                     }) 
    
    mu = matrix(data = unlist(x = mu_list), nrow = p) # put draws of mu back into same
    # p*K matrix format we were using before
    
    # draw group variances for all K groups
    loss_y_i = sapply(X = 1:k, 
                      FUN = function(x){
                        rowSums((matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p) - mu[,x])^2)
                        # unravel list of p*1 observations, put in matrix, find sum
                      })
    
    
    loss_mu_k = sapply(X = 1:k, 
                       FUN = function(x){
                         t(mu0 - matrix(mu[,x], nrow = p))%*%(mu0 - matrix(mu[,x], nrow = p))
                       })
    
    # Sigma = lapply(X = 1:k, 
    #                FUN = function(x){
    #                  diag(1/rgamma(n = p, 
    #                                shape = rep(count_assign[x]/2 + a, p),
    #                                rate = loss_y_i[,x]/2 + b))
    #                })
    sigma2 = lapply(X = 1:k, 
                    FUN = function(x){
                      1/rgamma(n = 1, 
                               shape = (p*(count_assign[x]+1) + 2*a)/2,
                               rate = sum(loss_y_i[,x])/2 + loss_mu_k[x]/(2*r) + b)
                    })
    
    # draw r parameter for variance of means
    if(fix_r == FALSE){
      
      loss_mu_k = sapply(X = 1:k, 
                         FUN = function(x){
                           t(matrix(mu[,x], nrow = p) - mu0)%*%(matrix(mu[,x], nrow = p) - mu0)/sigma2[[x]]
                         })
      
      r = 1/rgamma(n = 1, 
                   shape = (p*k + 2*g)/2, 
                   rate = sum(loss_mu_k)/2 + h)
      
      extra_params[s,"r"] = r
      
    }
    
    
    # draw b parameter if indicated
    if(sigma_hyperprior == TRUE){
      
      # find sum of precision for each variance element across k groups
      # sum_prec_k = Reduce(f = "+", 
      #                    x = lapply(X = 1:k,
      #                               FUN = function(x){
      #                                 1/diag(Sigma[[x]])
      #                               })
      #                    )
      
      inv_vars = lapply(X = 1:k,
                        FUN = function(x){
                          1/sigma2[[x]]
                        })
      
      b = rgamma(n = 1, 
                 shape = k*a + d, 
                 rate = Reduce(f = "+", x = inv_vars) + f)
      
      extra_params[s,"b"] = b
      
    } # else continue as usual
    
    # save draws of mu and Sigma
    means[[s]] = mu
    vars[[s]] = sigma2
    
    # save empirical mean and variance
    
    emp_means[[s]] = sapply(X = 1:k, 
                            FUN = function(x){
                              rowMeans(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p))
                              # unravel list of p*1 observations, put in matrix, find empirical mean
                            })
    
    emp_vars[[s]] = lapply(X = 1:k, 
                           FUN = function(x){
                             var(t(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p)))
                             # unravel list of p*1 observations, put in matrix, find emp cov matrix
                           })
    
    # print progress
    if((s %% print_iter == 0) & (s > print_iter) & (verbose == TRUE)){
      print("After Gibbs step:") # just create a new line for separate
      # print(paste("iter = ", s))
      # print(paste("Current k = ", k))
      # print(table(group_assign[s,]))
      cat("mu",mu)
      cat("sigma2",as.numeric(paste(sigma2)))
      cat("r",r)
      cat("b",b)
      # print(c("Current labels: ", curr_labels))
      # print(avail_labels)
      
      if(nrow(mu0) == 2){
        # if this is a 2D problem, can make scatterplot of group assign
        yvals = matrix(data = unlist(y), ncol = nrow(mu0), byrow = TRUE)
        plot_y = data.frame(
          y1 = yvals[,1],
          y2 = yvals[,2],
          curr_assign = group_assign[s,]
        )
        print(plot_y$curr_assign)
        prog_plot = ggplot(data = plot_y, aes(x = y1, y = y2, label = rownames(plot_y))) +
          #geom_point(color = assign) +
          #geom_text(size = 3, hjust = 0, nudge_x = 0.5, color = assign) +
          geom_text(size = 3, color = plot_y$curr_assign) +
          ggtitle(paste0("Group Assignments at Iteration s=", s, ", k=", k)) + 
          theme_classic()
        print(prog_plot)
      }
      
      
    }
    
    # at the end of each iteration, recalculate group probs based on final k
    # and save for label switching fix
    pr_c = sapply(X = 1:length(y), FUN = function(x){
      group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, a = a, b = b, 
                      y_i = y[[x]], mu = mu, sigma2 = sigma2, r = r, mu0 = mu0,
                      singleton = 0, curr_group_assign = group_assign[s,x], 
                      curr_labels = curr_labels)
    }) # result is a k*n matrix
    
    #print(dim(pr_c))
    if(k == 1){
      
      probs[[s]] = matrix(pr_c, ncol = 1)  # fixing issue with dimensions of prob matrix
      # when only one group is found
      
    } else{
      
      probs[[s]] = t(pr_c) # save probabilities for each obs, group
      # transpose so output is a n*k matrix - chop off pr_new -- not needed here
      
    }
    
    
    
    ### something is happening above --- output is sometimes 1d instead of n*k
    ### look at this
    
  } ### end MCMC iterations from s=1:S
  
  
  ## calculate Laplacian and average adjacency matrices across all MCMC runs
  pairwise_mats = pairwise_prob_mat(group_assign = group_assign, probs = probs, diag_weights = diag_weights)
  
  if(sigma_hyperprior == TRUE | fix_r == FALSE){
    
    settings = list(S = S, y = y, alpha = alpha, a = a, b = b, mu0 = mu0, 
                    k_init = k_init, d = d, f = f, g = g, h = h, r = r)
    
    return_list = list(settings = settings,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       emp_means = emp_means,
                       emp_vars = emp_vars,
                       extra_params = extra_params,
                       accept = accept_ind,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  } else{
    
    settings = list(S = S, seed = seed, y = y, alpha = alpha,
                    a = a, b = b, mu0 = mu0, k_init = k_init, 
                    d = d, f = f, r = r)
    
    return_list = list(settings = settings,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       emp_means = emp_means,
                       emp_vars = emp_vars,
                       #extra_params = extra_params,
                       accept = accept_ind,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  }
  
  return(return_list)
  
} ### end MCMC function



