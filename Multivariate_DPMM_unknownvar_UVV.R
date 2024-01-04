##############################################################
################# MULTIVARIATE DPMM MODEL ####################
#################  VAR UNKNOWN - ALG 2    ####################
################# INV WISH PRIOR - UVV    ####################
##############################################################

## Author: Jonathan Klus
## Date: 6 June 2023
## Description: Sampler for a mixture of multivariate normal densities using Algorithm 2
## from Neal (2000). We assume the variance is unknown, but we impose conjugacy
## by assuming the variance of the mean and likelihood differ only be a multiplicative
## parameter r. 

###################### load packages and set seed #############################
set.seed(516)
library(ggplot2)
library(LaplacesDemon)

########################### Gibbs sampler #####################################

########################### HELPER FUNCTIONS #################################

## calculate group membership probabilities

prior_pred_NinvW <- function(y_i, mu0, r, lambda0, nu){
  
  dens = LaplacesDemon::dmvt(x = y_i[,1], 
                             mu = mu0[,1], 
                             S = (r+1)*lambda0,
                             df = nu)
}

group_prob_calc_UVV <- function(k, n, n_j, alpha, y_i, mu, Sigma, r, mu0, lambda0, nu,
                                singleton = 0, curr_group_assign = NULL, curr_labels = NULL){
  
  # k is the number of existing groups
  # n is total number of observations
  # n_j is a vector of length k with the total number of observations in each group
  # alpha is the scalar concentration parameter from the dirichlet process
  # y_i is the single data vector of dimension p*1 that we are considering
  # mu is a matrix with k columns, contains k group mean vectors of dimension p*1
  # Sigma is a list of p*p covariance matrices from the Gibbs step for FC posteriors
  # mu0 is a p*1 prior mean vector, tau2 is a scalar prior
  # lambda0 is a matrix, corresponds to the prior covariance matrix
  # r is a scalar, corresponds to prior variance
  
  # calculate unnormalized probabilities of group membership for obs i
  
  p = length(mu[,1])
  
  #### probability of joining a current group
  
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
                                             sigma = Sigma[[x]]) 
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
                                               sigma = Sigma[[x]]) 
                       } else{
                         # print("Step3")
                         # print(matrix(mu[,x], nrow = p))
                         # print(y_i)
                         c1 = n_j[x]/(n-1+alpha)
                         c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                                               mean = mu[,x],
                                               sigma = Sigma[[x]]) 
                       }
                       return(c1*c2)
                     })
    
  }
  
  
  #### probability of creating a new group
  pr_new = alpha/(n-1+alpha)*prior_pred_NinvW(y_i = y_i, mu0 = mu0,
                                              r = r, lambda0 = lambda0, 
                                              nu = nu)
  
  #### normalize probs to account for "b"
  pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
  
  return(pr_c)
  
}

post_pred_UVV <- function(obs, which_group, r, sm_counts, nu, y, ybar, loss_ybar, mu0, lambda0){
  
  # restricted gibbs sampling scans
  n_minus = sm_counts[which_group] - 1
  loss_mu0 = (ybar[[which_group]] - mu0)%*%t(ybar[[which_group]] - mu0)
  mu_n = ((1/r)*mu0 + n_minus*ybar[[which_group]])/((1/r) + n_minus)
  
  # print(mu_n)
  # print(sm_counts)
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
  # 
  nu_n = nu + n_minus - nrow(mu0) + 1
  
  k_n = (r+n_minus+1)/((r+n_minus)*(nu_n))
  
  lambda_n = k_n*(lambda0 + loss_ybar[[which_group]] + 
                    loss_mu0*((n_minus/r)/(1/r + n_minus)))
  
  # print(mu_n)
  # print(lambda_n)
  # print(nu_n)
  
  val = n_minus*LaplacesDemon::dmvt(x = y[[obs]][,1], 
                                     mu = mu_n[,1], 
                                     S = lambda_n, 
                                     df = nu_n)
  
  return(val)
  
}

final_post_pred_UVV <- function(y_i, r, nu, y, mu0, lambda0){
  
  # use to calculate values after final gibbs scan 
  ## y_i is observation of interest, y is all data being considered for posterior,
  ## exlusive of the current observation y_i under consideration
  
  ybar = Reduce(f = "+", x = y)/length(y)
  
  loss_ybar = Reduce(f = "+", 
                     x = lapply(X = 1:length(y), 
                                FUN = function(x){
                                  (y[[x]] - ybar)%*%t(y[[x]] - ybar)
                                })
  )
  
  sm_counts = length(y)
  
  loss_mu0 = (ybar - mu0)%*%t(ybar - mu0)
  
  mu_n = ((1/r)*mu0 + sm_counts*ybar)/((1/r) + sm_counts)
  
  nu_n = nu + sm_counts - nrow(mu0) + 1
  
  k_n = (r+sm_counts+1)/((r+sm_counts)*(nu_n))
  
  lambda_n = k_n*(lambda0 + loss_ybar + 
                    loss_mu0*((sm_counts/r)/(1/r + sm_counts)))
  
  # print(mu_n)
  # print(lambda_n)
  # print(nu_n)
  
  val = LaplacesDemon::dmvt(x = y_i[,1], 
                            mu = mu_n[,1], 
                            S = lambda_n, 
                            df = nu_n)
  
  
}

ll_components_UVV <- function(subset_index, obs_ind, y, mu0, lambda0, r, nu){
  # function to calculate the components of the likelihood ratio for the MH
  # acceptance probability after n_iter split merge restricted Gibbs scans
  
  if(obs_ind == 1){
    # first observation --- prior predictive
    val = prior_pred_NinvW(y_i = y[[subset_index[obs_ind]]], 
                           mu0 = mu0, r = r, 
                           lambda0 = lambda0, nu = nu)
  } else{
    # posterior predictive
    # need to calculate for all obs, in same group
    subset_yvals = y[subset_index[1:obs_ind]] # what y values to include at each iter
    val = final_post_pred_UVV(y_i = y[[subset_index[obs_ind]]], r = r, nu = nu, 
                              y = y[subset_index[1:(obs_ind-1)]], 
                              mu0 = mu0, lambda0 = lambda0)
    
  }
  
  return(val)  
  
}

# test line
# post_pred_UVV(obs=2,which_group=1, r=10, sm_counts=c(10,11), nu=2, y=y, ybar=ybar, 
#               loss_ybar=loss_ybar, mu0=matrix(data=0,nrow=2), lambda0=diag(10,2))
######################

split_merge_prob_UVV <- function(obs, split_labs, group_assign, r, nu, y, mu0, lambda0){
  # split_labs is an array of length 2 indicating which entries in counts correspond
  # to the groups that are part of the split/merge
  # which_group is a scalar valued 1 or 2 indicating which of the groups we are considering
  # obs is the index for the observation being considered at this point
  # group_assign is an array of length n corresponding to group assignments for each obs
  # r and nu are scalar hyperparameters
  # counts is an array of the number of obs assigned to each group
  # y is the data
  # mu0 and lambda0 are the prior mean and covariance matrix
  
  # sm_counts = sapply(X = split_labs, FUN = function(x){sum(group_assign[-obs] == x)})
  # need to make up your mind -- drop obs or leave in here??
  sm_counts = sapply(X = split_labs, FUN = function(x){sum(group_assign == x)})

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
  
  cat("ybar:")
  print(ybar)
  cat("\n")
  
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
  
  cat("loss_ybar:")
  print(loss_ybar)
  cat("\n")
  
  
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
  
  if(1 %in% sm_counts){
    # if there is a singleton, take action to prevent issues with ybar and loss_ybar
    # need to use prior predictive instead of posterior predictive for this group

    which_one = which(sm_counts == 1)

    # check length of which_zero
    if(length(which_one) == 1){
      # proceed as usual
      # cat("1 singleton", "\n")
      if(which_one == 1){
        
        num = prior_pred_NinvW(y_i = y[[obs]], mu0 = mu0, r = r, 
                               lambda0 = lambda0, nu = nu)
        
        denom = num + post_pred_UVV(obs = obs, which_group = 2, r = r, 
                                    sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                                    loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
        
        ratio = c(num/denom, 1-(num/denom))
        
      } else{ 
        # which_one == 2
        
        num = post_pred_UVV(obs = obs, which_group = 1, r = r, 
                            sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                            loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
        
        denom = num + prior_pred_NinvW(y_i = y[[obs]], mu0 = mu0, r = r, 
                                       lambda0 = lambda0, nu = nu)
        
        ratio = c(num/denom, 1-(num/denom))
        
      }
      
    } else{
      # two singletons being considered
      # cat("2 singleton", "\n")
      num = prior_pred_NinvW(y_i = y[[obs]], mu0 = mu0, r = r, 
                             lambda0 = lambda0, nu = nu)
      
      denom = num + prior_pred_NinvW(y_i = y[[obs]], mu0 = mu0, r = r, 
                                     lambda0 = lambda0, nu = nu)
      
      ratio = c(num/denom, 1-(num/denom))
      
    }
    
  } else if(0 %in% sm_counts){
    # only comes up at very end when calculating merge transition prob
    which_one = which(sm_counts == 0)
    # cat("1 singleton", "\n")
    if(which_one == 1){
      
      num = prior_pred_NinvW(y_i = y[[obs]], mu0 = mu0, r = r, 
                             lambda0 = lambda0, nu = nu)
      
      denom = num + post_pred_UVV(obs = obs, which_group = 2, r = r, 
                                  sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                                  loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
      
      ratio = c(num/denom, 1-(num/denom))
      
    } else{ 
      # which_one == 2
      
      num = post_pred_UVV(obs = obs, which_group = 1, r = r, 
                          sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                          loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
      
      denom = num + prior_pred_NinvW(y_i = y[[obs]], mu0 = mu0, r = r, 
                                     lambda0 = lambda0, nu = nu)
      
      ratio = c(num/denom, 1-(num/denom))
      
    }
    
  } else{
    # proceed as usual 
    # cat("normal", "\n")
    num = post_pred_UVV(obs = obs, which_group = 1, r = r, 
                        sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                        loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
    
    denom = num + post_pred_UVV(obs = obs, which_group = 2, r = r,
                                sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                                loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
    
    ratio = c(num/denom, 1-(num/denom))
    
  }
  
  return(ratio)
  
}

MVN_CRP_sampler_UVV <- function(S = 10^3, seed = 516, y, r = 2, alpha = 1, lambda0, mu0, k_init = 2,
                                g = 1, h = 1, nu = 2, nu_hyperprior = FALSE, fix_r = FALSE, split_merge = FALSE,
                                sm_iter = 5, diag_weights = FALSE, verbose = TRUE, print_iter = 100){
  
  # S is number of MCMC iterations
  # y is a list of data of length n
  # lambda0 is the prior variance
  # alpha is the dirichlet process concentration parameter
  # nu is a prior hyperparameter for the invWish density. must be >1 or will cause error
  # d, f are the hyperprior parameters on the Gamma dist placed on b, assume a known
  # nu_hyperprior is a logical argument of whether to use the hyperprior or not
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
  
  # split merge step - only used if needed
  sm_results = matrix(data = NA, nrow = 1, ncol = 7)
  colnames(sm_results) = c("s", "sm_iter", "move_type","accept", "prob", "k_start", "k_end")
  
  if(nu_hyperprior == TRUE & fix_r == FALSE){
    extra_params = matrix(data = NA, nrow = S, ncol = 2)
    colnames(extra_params) = c("nu", "r")
    extra_params[1,1] = nu # note that here nu has dimension 1
    extra_params[1,2] = r
  } else if(nu_hyperprior == FALSE & fix_r == FALSE){
    extra_params = matrix(data = NA, nrow = S, ncol = 1)
    colnames(extra_params) = c("r")
    extra_params[1,1] = r
  } else if(nu_hyperprior == TRUE & fix_r == TRUE){
    extra_params = matrix(data = NA, nrow = S, ncol = 1)
    colnames(extra_params) = c("nu")
    extra_params[1,1] = nu
  } 
  
  
  # need to find p*1 vector of means based on this list of observed p*1 y_i values
  means[[1]] = sapply(X = 1:k, 
                      FUN = function(x){
                        rowMeans(matrix(unlist(y[group_assign[1,] == x]), nrow = p))
                        # unravel list of p*1 observations, put in matrix, find mean
                      })
  
  vars[[1]] = lapply(X = 1:k,
                     FUN = function(x){
                       diag(apply(X = matrix(unlist(y[group_assign[1,] == x]), nrow = p),
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
  Sigma = vars[[1]]
  curr_labels = 1:k
  avail_labels = c(1:n)[-curr_labels]
  
  # iterate 1:S
  for(s in 2:S){
    
    # print progress
    if((verbose == TRUE) & (s %% print_iter == 0)){
      
      cat("\n\n")
      cat("*******************************************************************")
      cat("\n")
      cat("Starting iter: ", s)
      cat("\n")
      cat("*******************************************************************")
      
    }    
    
    
    
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
      Sigma_curr = Sigma[[curr_assign]]
      
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
        
        # cat("\n")
        # cat("Singleton")
        
        #### only drop observation i if it is a singleton...DO NOT drop other singleton
        #### observations at this point!!!
        singleton_index = which(label_assign == group_assign[s,i]) 
        # check if there are singletons identified, otherwise do not touch mu
        
        # save current values for acceptance prob
        # mu_curr = mu[,singleton_index]
        # Sigma_curr = Sigma[[singleton_index]]
        
        # remove current values for singleton from index
        mu = matrix(mu[,-singleton_index], nrow = p)
        Sigma = Sigma[-singleton_index]
        
        count_assign = count_assign[-singleton_index]
        label_assign = label_assign[-singleton_index]
        
        avail_labels = c(group_assign[s,i], avail_labels)
        curr_labels = curr_labels[-which(curr_labels %in% group_assign[s,i])]
        k = length(curr_labels) # update k 
        
        #### calculate proposal distribution for group assignment
        ### for any observation i, calculate group membership probabilities
        pr_c = group_prob_calc_UVV(k = k, n = n, n_j = count_assign, alpha = alpha, 
                                   y_i = y[[i]], mu = mu, Sigma = Sigma, r = r, nu = nu,
                                   lambda0 = lambda0, mu0 = mu0, singleton = 1)
        
        
      } else{
        
        # cat("\n")
        # cat("non-singleton")
        # cat("\n")
        # cat("count_assign:", count_assign)
        # cat("\n")
        # cat("group_assign", group_assign[s,])
        # cat("\n")
        # cat("mu")
        # print(mu)
        # cat("\n")
        # print(Sigma)
        # cat("\n")
        #### calculate proposal distribution for group assignment
        #### if obs i is not presently a singleton
        pr_c = group_prob_calc_UVV(k = k, n = n, n_j = count_assign, alpha = alpha, 
                                   y_i = y[[i]], mu = mu, Sigma = Sigma, r = r, nu = nu,
                                   lambda0 = lambda0, mu0 = mu0, singleton = 0, 
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
        # cat("\n")
        # cat("New curr labels array after bookeeping:", curr_labels)
        avail_labels = avail_labels[-1]
        # cat("\n")
        # cat("New avail labels array after bookeeping:", head(avail_labels))
        k = length(curr_labels)
        
        ### using only the ith observation:
        
        #### draw variance for newly created group from FC posterior of sigma2
        #### according to algo, but this is conditional on mean so do prior for now
        Sigma_k = LaplacesDemon::rinvwishart(nu = nu, S = lambda0)
        length_Sigma = length(Sigma)
        Sigma[[length_Sigma+1]] = Sigma_k
        
        #### draw a mean for newly created group from FC posterior of mu 
        mu_mean = (y[[i]] + (1/r)*mu0)/(1/r + 1)
        mu_cov = Sigma_k/(1/r + 1)
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
    
    if((s %% print_iter == 0) & (s >= print_iter) & (verbose == TRUE)){
      cat("\n")
      cat("End of CRP step") # just create a new line for separation
      cat("\n")
      # print(paste("iter = ", s))
      # cat("\n")
      print(paste("Current k = ", k))
      cat("\n")
      print(group_assign[s,])
      cat("\n")
      print(table(group_assign[s,]))
      cat("\n")
      print(mu)
      cat("\n")
      print(Sigma)
      cat("\n")
    }
    
    # final update of counts after a sweep
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    num_groups[s,] = k
    
    # print results after each sweep
    
    # Split-Merge step --- every 10 iterations
    
    if((split_merge == TRUE) & (s %% sm_iter == 0)){

      k_start = k
      
      temp_group_assign = matrix(data = NA, nrow = sm_iter + 1, ncol = length(group_assign[s,]))
      sm_probs = matrix(data = NA, nrow = sm_iter + 1, ncol = length(y))
      temp_group_assign[1,] = group_assign[s,]
      
      # randomly select two observed data points y
      sampled_obs = sample(x = 1:length(y), size = 2, replace = FALSE)
      # cat("sampled_obs:", sampled_obs)
      # cat("\n")
      
      # check if in same group - if yes SPLIT
      # if no, MERGE
      lab1 = temp_group_assign[1, sampled_obs[1]]
      lab2 = temp_group_assign[1, sampled_obs[2]]
      move_type = ifelse(lab1 == lab2, "SPLIT", "MERGE")
      
      cat("move_type:", move_type)
      cat("\n")
      cat("sampled_obs:", sampled_obs)
      cat("\n")
      cat("group_labs:", c(lab1, lab2))
      cat("\n")
      
      # bookkeeping - group labels
      subset_index = which(temp_group_assign[1,] %in% c(lab1, lab2)) 
      # anchor_obs_index = which(subset_index %in% sampled_obs)
      # cat("anchor_obs_index", anchor_obs_index)
      # cat("\n")
      # check whether subset_index_minus contains only 2 observations
      # if(length(subset_index[-anchor_obs_index]) == 0){
      #   # when a single observation is in each group only
      #   # edge case - will return numeric zero
      #   subset_index_minus = NA
      #   
      # } else{
      #   # proceed as usual
      #   subset_index_minus = subset_index[-anchor_obs_index] # except sampled observations i and j
      # }
      
      existing_group_index = which(label_assign == lab1) 
      
      # perform restricted Gibbs scans
      
      # if SPLIT
      if(move_type == "SPLIT"){
        
        # specify random launch state
        split_lab = c(lab1, avail_labels[1]) # keep original label, new one for 2nd group
        
        cat("split_labs:", split_lab)
        cat("\n")
        # cat("subset_index:", subset_index)
        # cat("\n")

        for(scan in 1:(sm_iter+1)){
          
          cat("Scan #", scan)
          cat("\n")
          
          # initialize next scan with result of previous scan
          if(scan == 1){
            
            # intialize anchors so that there are no "zeros" for split groups
            # under consideration when starting below
            temp_group_assign[scan,sampled_obs[1]] = split_lab[1] 
            temp_group_assign[scan,sampled_obs[2]] = split_lab[2] 
            
          } else{
            
            temp_group_assign[scan,] = temp_group_assign[(scan-1),] # initialize
            # moved this up, was previously nested below, which means this step would
            # be *repeated* for each observation, thus washing out any split/merge that
            # was occurring as a result of the algorithm at each step
          }
          

          # print(table(temp_group_assign[scan,]))
          
          for(obs in subset_index){
            
            # cat("Obs:", obs)
            
              if(scan == 1){
                # specify random launch state
                if(obs == sampled_obs[2]){
                  temp_group_assign[scan,obs] = split_lab[2] 
                } else if(obs == sampled_obs[1]){
                  temp_group_assign[scan,obs] = split_lab[1] 
                } else{
                  # random launch state
                  temp_group_assign[scan,obs] = sample(x = split_lab, size = 1)
                }

                # cat(" Assign:",temp_group_assign[scan,obs])
                # cat("\n")
                
              } else{ # for remaining scans after random launch state set
                
                if(obs %in% sampled_obs){
                  # cat(" Anchor, ")
                  # dont sample if anchor obs -- assignment cannot change
                  # specify new group label to 2nd anchor point as well
                  if(obs == sampled_obs[2]){
                    temp_group_assign[scan,obs] = split_lab[2] 
                  } else if(obs == sampled_obs[1]){
                    temp_group_assign[scan,obs] = split_lab[1] 
                  }
                  
                  # cat(" Assign:",temp_group_assign[scan,obs])
                  # cat("\n")
                  
                  split_assign_prob = split_merge_prob_UVV(
                    obs = obs, split_labs = split_lab, r=r, 
                    group_assign = temp_group_assign[scan,], nu = nu, 
                    y = y, mu0 = mu0, lambda0= lambda0)
                  
                  # dont sample --- but do record prob of group anchor obs in in!
                  which_split_lab_anchor = which(split_lab == temp_group_assign[scan,obs])
                  # temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                  # dont need to assign again - already initialized since anchor
                  sm_probs[scan,obs] = split_assign_prob[which_split_lab_anchor]
                  
                } else{
                  
                  sm_count_assign = as.numeric(table(temp_group_assign[scan,]))
                  sm_label_assign = as.numeric(names(table(temp_group_assign[scan,])))
                  sm_singletons = sm_label_assign[which(sm_count_assign == 1)]
                  
                  if(temp_group_assign[scan,obs] %in% sm_singletons){ # changing this to if singleton
                    # maybe borrow some code from CRP stuff???
                    # still need to deal with anchor obs because group assignment
                    # is deterministic for these observations
                    
                    # ad hoc fix for now!
                    sm_counts = table(temp_group_assign[scan,]) # dont drop obs from count when singleton
                    split_group_count_index = which(as.numeric(names(sm_counts)) %in% split_lab)
                    #current_obs_index = which(temp_group_assign[scan,] == obs)
                    #split_group_lab_index1 = which(temp_group_assign[scan,] == split_lab[1])
                    #split_group_lab_index2 = which(temp_group_assign[scan,] == split_lab[2])
                    
                    # current observation under consideration cannot be included in posterior
                    # when singleton, revert to prior
                    
                    split_assign_prob = split_merge_prob_UVV(
                      obs = obs, split_labs = split_lab, r=r, 
                      group_assign = temp_group_assign[scan,], nu = nu, 
                      y = y, mu0 = mu0, lambda0= lambda0)
                    
                    
                    # proceed as usual
                    sm_prop_index = sample(x = 1:2, size = 1, 
                                           prob = split_assign_prob)
                    
                    temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                    sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
                    
                    
                    # cat("Assign:",temp_group_assign[scan,obs])
                    # cat("\n")
                    
                    
                  } else{
                    
                    sm_counts = table(temp_group_assign[scan,-obs])
                    
                    #current_obs_index = which(temp_group_assign[scan,] == obs)
                    #split_group_lab_index1 = which(temp_group_assign[scan,] == split_lab[1])
                    #split_group_lab_index2 = which(temp_group_assign[scan,] == split_lab[2])
                    
                    # current observation under consideration cannot be included here
                    
                    split_assign_prob = split_merge_prob_UVV(
                      obs = obs, split_labs = split_lab, r=r, 
                      group_assign = temp_group_assign[scan,], nu = nu, 
                      y = y, mu0 = mu0, lambda0= lambda0)
                    
                    sm_prop_index = sample(x = 1:2, size = 1, 
                                           prob = split_assign_prob)
                    
                    temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                    sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
                    
                    # cat("Assign:",temp_group_assign[scan,obs])
                    # cat("\n")
                  }
                  
                }
              
            } 
              

 
          } # iterate through all observations in the two split groups under consideration
        } # scans 1:(sm_iter+1)
        
        print(temp_group_assign)
        cat("\n")
        
        # calculate & evaluate acceptance prob
        sm_counts = table(temp_group_assign[sm_iter+1,]) # update counts after scans
        split_group_count_index = which(as.numeric(names(sm_counts)) %in% split_lab)
        
        cat("Final counts & labels", "\n")
        print(sm_counts)
        cat("\n")
        ## proposal probability

        #print(sm_probs[sm_iter+1,])
        # print(sm_counts)
        
        prob1 = -Reduce(f = "+", x = log(sm_probs[sm_iter+1,subset_index])) # log1 - sum(logs)
        #print(sm_probs)
        ## prior ratio
        
        cat("split_group_count_index", split_group_count_index, "\n")
        cat("sm_counts", sm_counts, "\n")
        # print(temp_group_assign[scan,])
        
        prob2_num = factorial(sm_counts[[split_group_count_index[1]]] -1)*factorial(sm_counts[[split_group_count_index[2]]] -1)
        prob2_denom = factorial(sm_counts[[split_group_count_index[1]]]+sm_counts[[split_group_count_index[2]]]-1)
        prob2 = log(alpha) + (log(prob2_num) - log(prob2_denom))
        
        ## likelihood ratio
        subset_index_grp1 = which(temp_group_assign[sm_iter+1,] %in% split_lab[1]) 
        subset_index_grp2 = which(temp_group_assign[sm_iter+1,] %in% split_lab[2]) 
        
        ### component 1 - numerator I (group 1)
        prob3_num1 = 0
        for(obs_ind in 1:length(subset_index_grp1)){
          val = ll_components_UVV(subset_index = subset_index_grp1, obs_ind = obs_ind, y = y, 
                              mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
          prob3_num1 = prob3_num1 + log(val)
        }
        
        ### component 2 - numerator II (group 2)
        prob3_num2 = 0
        for(obs_ind in 1:length(subset_index_grp2)){
          val = ll_components_UVV(subset_index = subset_index_grp2, obs_ind = obs_ind, y = y, 
                              mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
          prob3_num2 = prob3_num2 + log(val)
        }
        
        
        ### component 3 - denominator (all in original group)
        prob3_denom = 0
        for(obs_ind in 1:length(subset_index)){
          val = ll_components_UVV(subset_index = subset_index, obs_ind = obs_ind, y = y, 
                              mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
          prob3_denom = prob3_denom + log(val)
        }
        
        ## evaluate acceptance prob
        prob3 = prob3_num1 + prob3_num2 - prob3_denom
        accept_prob = min(1, exp(prob1 + prob2 + prob3))
        u = runif(n = 1)
        if(accept_prob > u){
          
          # accept
          accept = 1
          group_assign[s,] = temp_group_assign[sm_iter+1,]
          
          # print(table(group_assign[s,]))
          
          # take new group lab out of reserve
          curr_labels = c(curr_labels, avail_labels[1])
          avail_labels = avail_labels[-1]
          k = length(curr_labels)
          # update labels, etc
          count_assign = as.numeric(table(group_assign[s,]))
          label_assign = as.numeric(names(table(group_assign[s,])))
          which_split_labs = which(label_assign %in% split_lab) 
          num_groups[s,] = k
          
          # if new group created by split, give it a mean and variance
          
          ## draw variances for both groups - use marginal posterior variance 
          
          ybar = lapply(X = split_lab, 
                        FUN = function(x){
                          group_ind = which(group_assign[s,] == x)
                          if(obs %in% group_ind){
                            obs_ind = which(obs == group_ind)
                            # group_ind = group_ind[-obs_ind] # include all obs this time
                          } # else continue
                          
                          ysum = Reduce(f = "+", 
                                        x = lapply(X = group_ind, 
                                                   FUN = function(x){(y[[x]])}))
                          
                          return(ysum/length(group_ind))
                        })

          emp_loss = lapply(X = 1:2, 
                           FUN = function(x){
                             
                             col_ind = x  # from outer apply
                             group_ind = which(group_assign[s,] == split_lab[x])
                             if(obs %in% group_ind){
                               obs_ind = which(obs == group_ind)
                               #group_ind = group_ind[-obs_ind] # include all obs this time
                             } # else continue
                             
                             emp_loss = Reduce(f = "+", 
                                               x = lapply(X = group_ind, FUN = function(x){
                                                 (y[[x]] - mu0)%*%t(y[[x]] - mu0)}))
                             
                             return(emp_loss)
                             
                           })
          
          Sigma_temp = lapply(X = 1:2, 
                         FUN = function(x){
                           n_k = count_assign[which_split_labs[x]]
                           LaplacesDemon::rinvwishart(nu = nu + n_k + 1, 
                                                      S = emp_loss[[x]]/(1+r) + nu*lambda0)
                         })
          
          ## draw means for both groups conditional on empirical variances...
          # cat("\n")
          # cat("split_lab", split_lab)
          # cat("\n")
          # cat("grp_assgn", group_assign[s,])
          sum_y_i = sapply(X = split_lab, 
                           FUN = function(x){
                             rowSums(matrix(unlist(y[group_assign[s,] == x]), nrow = p))
                             # unravel list of p*1 observations, put in matrix, find sum
                           })
          
          mu_cov = lapply(X = 1:2, 
                          FUN = function(x){
                            n_k = count_assign[which_split_labs[x]]
                            return(Sigma_temp[[x]]/(1/r + n_k))
                            }) 
          
          mu_mean = lapply(X = 1:2, 
                           FUN = function(x){
                             n_k = count_assign[which_split_labs[x]]
                             return((sum_y_i[,x] + mu0/r)/(1/r + n_k))
                             })
          
          mu_list = lapply(X = 1:2, 
                           FUN = function(x){
                             t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                                mean = mu_mean[[x]], 
                                                sigma = mu_cov[[x]]))
                           }) 
          
          ## add new means and variances to relevant vectors/lists
          length_Sigma = length(Sigma)
          Sigma[[which_split_labs[1]]] = Sigma_temp[[1]]
          Sigma[[length_Sigma+1]] = Sigma_temp[[2]]
          rm(Sigma_temp)
          
          mu[,which_split_labs[1]] = mu_list[[1]]
          mu = cbind(mu, mu_list[[2]])
          
        } else{
          # reject
          accept = 0
          # group assign remains unchanged
        }
        
        
        
        
        # if MERGE    
      } else if(move_type == "MERGE"){
        
        # specify random launch state
        split_lab = c(lab1, lab2) # original labels, assign merged obs to lab1
        
        count_assign = as.numeric(table(temp_group_assign[1,]))
        label_assign = as.numeric(names(table(temp_group_assign[1,])))
        which_split_labs = which(label_assign %in% split_lab) 
        
        scan = sm_iter + 1 # skip to last row of table
        temp_group_assign[scan,] = temp_group_assign[1,] # initialize
        
        sm_counts_before = table(temp_group_assign[scan,])
        split_group_count_index = which(as.numeric(names(sm_counts_before)) %in% split_lab)
        
        temp_group_assign[scan,sampled_obs] = lab1 # anchor observations that aren't
        # in subset-index-minus...
        
        for(obs in subset_index){
          
          temp_group_assign[scan,obs] = lab1 # keep 1st group label
          
          # sm_counts_before = table(temp_group_assign[scan,-obs])
          # split_group_count_index_before = which(as.numeric(names(sm_counts)) %in% split_lab)
          
          # current observation under consideration cannot be included here
          
          split_assign_prob = split_merge_prob_UVV(
            obs = obs, split_labs = split_lab, r=r, 
            group_assign = temp_group_assign[scan,], nu = nu, 
            y = y, mu0 = mu0, lambda0= lambda0)
          
          sm_prop_index = sample(x = 1:2, size = 1,
                                 prob = split_assign_prob)
          # 
          # temp_group_assign[scan,obs] = split_lab[sm_prop_index]
          sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
          
        }  # iterate through all observations in the two split groups under consideration
        

        # calculate & evaluate acceptance prob
        sm_counts = table(temp_group_assign[sm_iter+1,]) # update counts after scans
        ## proposal probability
        prob1 = Reduce(f = "+", x = log(sm_probs[sm_iter+1,subset_index])) # log1 - sum(logs)
        # need to index by subset since NAs for observation not part of split/merge -- as well as the
        # two anchor observations. shoudl you be calculating a prob for those as well though???
        
        ## prior ratio
        
        # cat("\n")
        # cat("sm_counts_before:")
        # print(sm_counts_before)
        # cat("\n")
        # cat("split_group_count_index:", split_group_count_index)
        
        prob2_num = factorial(sm_counts_before[[split_group_count_index[1]]] -1)*factorial(sm_counts_before[[split_group_count_index[2]]] -1)
        prob2_denom = factorial(sm_counts_before[[split_group_count_index[1]]]+sm_counts_before[[split_group_count_index[2]]]-1)
        prob2 = -(log(alpha) + (log(prob2_num) - log(prob2_denom)))
        
        ## likelihood ratio
        subset_index_grp1 = which(temp_group_assign[1,] %in% split_lab[1]) 
        subset_index_grp2 = which(temp_group_assign[1,] %in% split_lab[2]) 
        
        ### component 1 - numerator I (group 1)
        prob3_num1 = 0
        for(obs_ind in 1:length(subset_index_grp1)){
          val = ll_components_UVV(subset_index = subset_index_grp1, obs_ind = obs_ind, y = y, 
                              mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
          prob3_num1 = prob3_num1 + log(val)
        }
        
        ### component 2 - numerator II (group 2)
        prob3_num2 = 0
        for(obs_ind in 1:length(subset_index_grp2)){
          val = ll_components_UVV(subset_index = subset_index_grp2, obs_ind = obs_ind, y = y, 
                              mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
          prob3_num2 = prob3_num2 + log(val)
        }
        
        
        ### component 3 - denominator (all in original group)
        prob3_denom = 0
        for(obs_ind in 1:length(subset_index)){
          val = ll_components_UVV(subset_index = subset_index, obs_ind = obs_ind, y = y, 
                              mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
          prob3_denom = prob3_denom + log(val)
        }
        
        ## evaluate acceptance prob
        prob3 = -(prob3_num1 + prob3_num2 - prob3_denom)
        accept_prob = min(1, exp(prob1 + prob2 + prob3))
        u = runif(n = 1)
        if(accept_prob > u){
          # accept
          accept = 1
          group_assign[s,] = temp_group_assign[sm_iter+1,]
          
          # bookeeping if merge step ACCEPTED only, keep first group -- remove 2nd 
          mu = matrix(mu[,-which_split_labs[2]], nrow = p)
          Sigma = Sigma[-which_split_labs[2]]
          
          # count_assign = count_assign[-which_split_labs[2]]
          # label_assign = label_assign[-which_split_labs[2]]
          # avail_labels = c(temp_group_assign[1, sampled_obs[2]], avail_labels)
          # i think the line above is duplicating bookeeping step of adding to availabe
          # labels...comment out for now and test...
          
          # put old group lab out of reserve
          curr_label_del = which(curr_labels == split_lab[2])
          avail_labels = c(curr_labels[curr_label_del], avail_labels)
          curr_labels = curr_labels[-curr_label_del]
          
          k = length(curr_labels)
          # update labels, etc
          count_assign = as.numeric(table(group_assign[s,]))
          label_assign = as.numeric(names(table(group_assign[s,])))
          which_split_lab = which(label_assign == split_lab[1]) 
          num_groups[s,] = k
          
          # save info about accept/reject & probs
          
          # merged group...give it a mean and variance
          
          ## draw variances for new groups - use empirical variance since mean not known yet
          ## may want to change this in the future if you come up with a better idea...
          
          # ybar = lapply(X = split_lab[1], 
          #               FUN = function(x){
          #                 group_ind = which(group_assign[s,] == x)
          #                 if(obs %in% group_ind){
          #                   obs_ind = which(obs == group_ind)
          #                   # group_ind = group_ind[-obs_ind] # include all obs this time
          #                 } # else continue
          #                 
          #                 ysum = Reduce(f = "+", 
          #                               x = lapply(X = group_ind, 
          #                                          FUN = function(x){(y[[x]])}))
          #                 
          #                 return(ysum/length(group_ind))
          #               })
          
          group_ind = which(group_assign[s,] == split_lab[1])
          ybar = Reduce(f = "+", x = y[group_ind])/length(group_ind) 
          
          emp_loss = Reduce(f = "+", 
                            x = lapply(X = group_ind, FUN = function(x){
                              (y[[x]] - mu0)%*%t(y[[x]] - mu0)}))
          
          # print(emp_loss/(1+r))
          # print(nu*lambda0)
          
          Sigma_temp = LaplacesDemon::rinvwishart(nu = nu + count_assign[which_split_lab[1]] + 1, 
                                                  S = emp_loss/(1+r) + nu*lambda0)

          
          ## draw means for both groups conditional on drawn variances...
          
          sum_y_i = rowSums(matrix(unlist(y[group_assign[s,] == split_lab[1]]), nrow = p))
          # unravel list of p*1 observations, put in matrix, find sum
          
          mu_cov = Sigma_temp/(1/r + count_assign[which_split_lab[1]])
          
          mu_mean = (sum_y_i + mu0/r)/(1/r + count_assign[which_split_lab[1]])
          
          mu_list = t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                       mean = mu_mean, 
                                       sigma = mu_cov))
          
          ## add new means and variances to relevant vectors/lists
          length_Sigma = length(Sigma)
          Sigma[[which_split_lab]] = Sigma_temp
          
          mu[,which_split_lab] = mu_list
          
          
        } else{
          # reject
          accept = 0
          # group assign remains unchanged
        }
        
        
        
        
        
      }
      
      k_end = k
      sm_results = rbind(sm_results, c(s, sm_iter, move_type, accept, accept_prob, k_start, k_end))
      
      # cat("\n")
      # cat("SM Step complete:")
      # cat("\n")
      # cat("iter", s)
      # cat("\n")
      # cat("move_type", move_type)
      # cat("\n")
      # cat("accept", accept)
      # cat("\n")
      # cat("prob", accept_prob)
      # cat("\n")
      # print(table(group_assign[s,]))
      
    }
    
    
    
    # Proceed to Gibbs step
    
    # draw group means for all K groups
    # Sigma = diag(x = sigma2, nrow = p)
    
    sum_y_i = sapply(X = 1:k, 
                     FUN = function(x){
                       rowSums(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p))
                       # unravel list of p*1 observations, put in matrix, find sum
                     })
    
    mu_cov = lapply(X = 1:k, 
                    FUN = function(x){Sigma[[x]]/(1/r + count_assign[x])}) 
    
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
    
    #### need to redo loss fxns for cov matrix --- X %*% X' now --- p*p matrix!!!
    loss_y_i = lapply(X = 1:k, 
                      FUN = function(x){
                        
                        col_ind = x  # from outer apply
                        group_ind = which(group_assign[s,] == label_assign[x])
                        
                        Reduce(f = "+", 
                               x = lapply(X = group_ind, FUN = function(x){
                                 (y[[x]] - mu[,col_ind])%*%t(y[[x]] - mu[,col_ind])}))
                        
                      })
    
    
    loss_mu_k = lapply(X = 1:k, 
                       FUN = function(x){
                         (mu0 - matrix(mu[,x], nrow = p))%*%t(mu0 - matrix(mu[,x], nrow = p))
                       })
    
    # Sigma = lapply(X = 1:k, 
    #                FUN = function(x){
    #                  diag(1/rgamma(n = p, 
    #                                shape = rep(count_assign[x]/2 + a, p),
    #                                rate = loss_y_i[,x]/2 + b))
    #                })
    Sigma = lapply(X = 1:k, 
                   FUN = function(x){
                     LaplacesDemon::rinvwishart(nu = nu + count_assign[x] + 1, 
                                                S = loss_y_i[[x]] + loss_mu_k[[x]]/r + nu*lambda0)
                   })
    
    # draw r parameter for variance of means
    if(fix_r == FALSE){
      
      loss_mu_k = sapply(X = 1:k, 
                         FUN = function(x){
                           t(matrix(mu[,x], nrow = p) - mu0)%*%solve(Sigma[[x]])%*%(matrix(mu[,x], nrow = p) - mu0)
                         })
      
      r = 1/rgamma(n = 1, 
                   shape = g + p/2, 
                   rate = sum(loss_mu_k)/2 + h)
      
      extra_params[s,"r"] = r
      
    }
    
    
    # draw b parameter if indicated
    if(nu_hyperprior == TRUE){
      
      # find sum of precision for each variance element across k groups
      # sum_prec_k = Reduce(f = "+", 
      #                    x = lapply(X = 1:k,
      #                               FUN = function(x){
      #                                 1/diag(Sigma[[x]])
      #                               })
      #                    )
      
      ##
      ## NEED TO REWRITE THIS PIECE IN TERMS OF NU
      ##
      ##
      
      # inv_vars = lapply(X = 1:k,
      #                   FUN = function(x){
      #                     1/sigma2[[x]]
      #                   })
      # 
      # b = rgamma(n = 1, 
      #            shape = k*a + d, 
      #            rate = Reduce(f = "+", x = inv_vars) + f)
      # 
      # extra_params[s,"b"] = b
      
    } # else continue as usual
    
    # save draws of mu and Sigma
    means[[s]] = mu
    vars[[s]] = Sigma
    
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
    if((s %% print_iter == 0) & (s >= print_iter) & (verbose == TRUE)){
      cat("After Gibbs step:") # just create a new line for separate
      cat("\n")
      cat("mu")
      cat("\n")
      print(mu)
      cat("\n")
      cat("Sigma")
      cat("\n")
      print(Sigma)
      cat("\n")
      cat("r",r)
      cat("\n")
      cat("nu",nu)
      cat("\n")
      
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
      group_prob_calc_UVV(k = k, n = n, n_j = count_assign, alpha = alpha, nu = nu, 
                          y_i = y[[x]], mu = mu, Sigma = Sigma, r = r, mu0 = mu0, 
                          lambda0 = lambda0, curr_group_assign = group_assign[s,x], 
                          singleton = 0, curr_labels = curr_labels)
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
  
  if(nu_hyperprior == TRUE | fix_r == FALSE){
    
    settings = list(S = S, alpha = alpha, mu0 = mu0, lambda0 = lambda0,
                    k_init = k_init, g = g, h = h, r = r) # d = d, f = f,
    
    return_list = list(settings = settings,
                       data = y,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       emp_means = emp_means,
                       emp_vars = emp_vars,
                       extra_params = extra_params,
                       accept = accept_ind,
                       sm_results = sm_results,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  } else{
    
    settings = list(S = S, seed = seed, alpha = alpha, lambda0 = lambda0,
                    nu = nu, mu0 = mu0, k_init = k_init, g = g, h = h, r = r)
    
    return_list = list(settings = settings,
                       data = y, 
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       emp_means = emp_means,
                       emp_vars = emp_vars,
                       #extra_params = extra_params,
                       accept = accept_ind,
                       sm_results = sm_results,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  }
  
  return(return_list)
  
} ### end MCMC function
