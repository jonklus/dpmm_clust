##############################################################
################# MULTIVARIATE DPMM MODEL ####################
###############  NO GAPS NON CONJ - ALG 8 ####################
################# INDEP IG PRIORS - DEV   ####################
##############################################################

## Author: Jonathan Klus
## Date: 6 June 2023
## Description: Sampler for a mixture of multivariate normal densities using Algorithm 8
## from Neal (2000). We assume the variance is unknown, We assume the variance is unknown, 
## and no longer assume conjugacy.

############################# load packages ####################################
# set.seed(516)
library(ggplot2)
library(LaplacesDemon)

############################### HELPER FUNCTIONS ###############################

## calculate group membership probabilities

group_prob_calc_DEV <- function(k, n, n_j, alpha, y_i, mu, sigma2, a, b, mu0, sigma0,
                                singleton = 0, curr_group_assign = NULL, 
                                curr_labels = NULL){
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
  
  # DEV case, variance differs by group but still diagonal within group
  
  k_minus = k-1 # number of distinct groups excluding the jth group
  
  if(singleton == 1){
    
    
    # label at k_minus + 1
    
    ## draw params for new group from G_0 in step below -- will be same
    ## procedure as non-singleton
    
    ## calculate prob of current groups
    ## not sure that this is correct way to handle singletons here === come back
    ## looks like we're already dropping ith observation if singleton before
    ## sending to function so it's kosher
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
  
  #### probability of creating a new group
  
  ##### draw new values of phi from G_0
  mu_new = t(mvtnorm::rmvnorm(n = 1, mean = mu0[,1], sigma = diag(sigma0, p))) # column vec
  sigma2_new = 1/rgamma(n = 1, shape = a, rate = b)
  
  ##### calculate probs
  pr_new = ((alpha/1)/(n-1+alpha))*mvtnorm::dmvnorm(x = y_i[,1], 
                                                    mean = mu_new[,1],
                                                    sigma = diag(sigma2_new, p)) 
  #### normalize probs to account for "b"
  pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
  
  return(list(pr_c = pr_c, sigma2_new = sigma2_new))
  
}



update_phi_DEV <- function(curr_label, group_assign, count_assign, y, 
                           mu, mu0, Sigma, Sigma0, a, b){
  # function to sample from full conditional posteriors of DEV model
  
  # group_index is the current split or merge group for which new parameters are 
  # being drawn
  
  # distinguish between split, where mu and Sigma will have two entries, and merge
  # where they will have just 1???? do we need to do this?
  
  if(is.null(nrow(mu)) == TRUE){
    p = length(mu)
  } else{
    p = nrow(mu)
  }
  
  sigma2 = diag(Sigma)[1]
  sigma0 = diag(Sigma0)[1]
  
  
  
  # cat("\n p=", p)
  # cat("\n curr_label=", curr_label)
  # print(unlist(y[group_assign == curr_label]))
  # print(matrix(unlist(y[group_assign == curr_label]), nrow = p))
  
  # draw group variance
  loss_y_i = rowSums((matrix(unlist(y[group_assign == curr_label]), nrow = p) - c(mu))^2)
  
  loss_mu_k = t(mu0 - matrix(mu, nrow = p))%*%(mu0 - matrix(mu, nrow = p))
  
  # Sigma = lapply(X = 1:k, 
  #                FUN = function(x){
  #                  diag(1/rgamma(n = p, 
  #                                shape = rep(count_assign[x]/2 + a, p),
  #                                rate = loss_y_i[,x]/2 + b))
  #                })
  Sigma = diag(1/rgamma(n = 1, 
                        shape = (p*count_assign + 2*a)/2,
                        rate = sum(loss_y_i)/2 + b), p)
  
  
  # draw group mean
  sum_y_i = rowSums(matrix(unlist(y[group_assign == curr_label]), nrow = p))
  
  mu_cov = 1/(count_assign/sigma2 + 1/sigma0)
  
  mu_mean = (sum_y_i/sigma2 + mu0/sigma0)*mu_cov
  
  mu = t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                          mean = mu_mean, 
                          sigma = diag(mu_cov,p)))
  
  
  

  
  return(list(mu = mu, Sigma = Sigma))
  
}

nonconj_prior_dens_DEV <- function(mu, mu0, Sigma, Sigma0, a, b){
  
  # this is the function g(\phi_c) defined in Jain & Neal (2007) - used to 
  # evaluate group-specific parameters \phi_c under the joint prior density placed
  # on these parameters
  
  # mu and mu0 are p-dimensional arrays
  # Sigma Sigma0 are p*p dimensional matrices
  # a and b are scalar hyperparameters for the IG prior
  
  sigma2 = diag(Sigma)[1]
  # sigma0 = diag(Sigma0)[1]
  # p = nrow(mu)
  
  dens = log(mvtnorm::dmvnorm(x = c(mu), mean = c(mu0), sigma = Sigma0)) + 
    log(LaplacesDemon::dinvgamma(x = sigma2, shape = a, scale = b))
  
  return(exp(dens))
}


nonconj_phi_prob_DEV <- function(curr_label, group_assign, count_assign, y, 
                                 mu, mu_L, mu0, Sigma, Sigma_L, Sigma0, a, b){
  # This is P_GS(phi*|phi^L,...) from Jain & Neal 2007 
  
  
  # for the kth component under a DEV assumption
  # mu and mu0 are p-dimensional arrays
  # Sigma Sigma0 are p*p dimensional matrices
  # mu_L represents the launch state, mu is the current state
  # a and b are scalar hyperparameters for the IG prior
  
    if(is.null(nrow(mu)) == TRUE){
      p = length(mu)
    } else{
      p = nrow(mu)
    }
  
    sigma2 = diag(Sigma)[1]
    sigma2_L = diag(Sigma_L)[1]
    sigma0 = diag(Sigma0)[1]
    
    # mean - posterior density of newly drawn mu* given launch state
    sum_y_i = rowSums(matrix(unlist(y[group_assign == curr_label]), nrow = p))
    
    mu_cov = 1/(count_assign/sigma2_L + 1/sigma0)
    
    mu_mean = (sum_y_i/sigma2_L + mu0/sigma0)*mu_cov
    
    ldens_mu = mvtnorm::dmvnorm(
      x = c(mu),
      mean = c(mu_mean), 
      sigma = diag(mu_cov,p))
  
  # variance - posterior density of newly drawn sigma2* given launch state
  loss_y_i = rowSums((matrix(unlist(y[group_assign == curr_label]), nrow = p) - c(mu_L))^2)
  
  loss_mu_k = t(mu0 - matrix(mu_L, nrow = p))%*%(mu0 - matrix(mu_L, nrow = p))
  
  ldens_sigma2 = LaplacesDemon::dinvgamma(
    x = sigma2,
    shape = (p*count_assign + 2*a)/2,
    scale = sum(loss_y_i)/2 + b, 
    log = TRUE)
  
  return(ldens_mu + ldens_sigma2)
}

ll_components_DEV <- function(obs_ind, y, mu, Sigma){
  # function to calculate the components of the likelihood ratio for the MH
  # acceptance probability after n_iter split merge restricted Gibbs scans
  
  val = mvtnorm::dmvnorm(x = y[[obs_ind]][,1], 
                         mean = c(mu), 
                         sigma = Sigma)
  
  return(val)  
  
}

nonconj_component_prob_c <- function(obs, split_labs, group_assign, y, mu, Sigma,
                                     proposal_calc = FALSE){
  # This is P_GS(c*|c^L,...) from Jain & Neal 2007 
  
  # split_labs is an array of length 2 indicating which entries in counts correspond
  # to the groups that are part of the split/merge
  # which_group is a scalar valued 1 or 2 indicating which of the groups we are considering
  # obs is the index for the observation being considered at this point
  # group_assign is an array of length n corresponding to group assignments for each obs
  # counts is an array of the number of obs assigned to each group
  # y is the data
  # mu and Sigma are lists of the the current mean and covariance matrix for each 
  # component 
  
  # sm_counts = sapply(X = split_labs, FUN = function(x){sum(group_assign[-obs] == x)})
  # need to make up your mind -- drop obs or leave in here??
  sm_counts = sapply(X = split_labs, FUN = function(x){sum(group_assign == x)})
  cat("\n sm_counts:", sm_counts, "\n")
  which_group_k = which(split_labs == group_assign[obs])
  sm_counts[which_group_k] = sm_counts[which_group_k] - 1
  cat("\n sm_counts[which_group_k]:", sm_counts[which_group_k], "\n")
  
  if(0 %in% sm_counts){
    
    which_one = which(sm_counts == 1)
    
    # check length of which_zero
    if(length(which_one) == 1){
      # proceed as usual
      # cat("1 singleton", "\n")
      if(which_one == 1){
        
        num = (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]][,1], 
                                              mean = mu[[2]][,1], 
                                              sigma = Sigma[[2]]) 
        
        denom = num + (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]][,1], 
                                                      mean = mu[[1]][,1], 
                                                      sigma = Sigma[[1]])
        
        # will just be (1,0)... seems extreme
        if(is.nan(num/denom)){
          cat("\n Warning: NaN detected in SM group assignment calculation. ")
          ratio = rep(1/2,2)
        } else{
          ratio =  ratio = c(1-(num/denom), num/denom)
        }
        
      } else{ 
        # which_one == 2
        
        num = (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]][,1], 
                                              mean = mu[[1]][,1], 
                                              sigma = Sigma[[1]]) 
        
        denom = num + (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]][,1], 
                                                      mean = mu[[2]][,1], 
                                                      sigma = Sigma[[2]])
        
        # will just be (1,0)... seems extreme
        if(is.nan(num/denom)){
          cat("\n Warning: NaN detected in SM group assignment calculation. ")
          ratio = rep(1/2,2)
        } else{
          ratio = c(num/denom, 1-(num/denom))
        }
        
      }
      
    } else{
      # two singletons being considered
      # cat("2 singleton", "\n")
      # num = sm_counts[1]*mvtnorm::dmvnorm(x = y[[obs]], 
      #                                     mean = mu[[1]], 
      #                                     sigma = Sigma[[1]]) 
      # 
      # denom = num + sm_counts[2]*mvtnorm::dmvnorm(x = y[[obs]], 
      #                                             mean = mu[[2]], 
      #                                             sigma = Sigma[[2]])
      
      ratio = c(0.5, 0.5)
      
    }
    
    
    # } else if(0 %in% sm_counts){
    #   
    #   which_one = which(sm_counts == 0)
    #   
    #   if(which_one == 1){
    #     
    #     num = (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]], 
    #                                           mean = mu[[2]], 
    #                                           sigma = Sigma[[2]]) 
    #     
    #     denom = num + (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]], 
    #                                                   mean = mu[[1]], 
    #                                                   sigma = Sigma[[1]])
    #     
    #     # will just be (1,0)... seems extreme
    #     ratio = c(1-(num/denom), num/denom)
    
    # } else{ 
    #   # which_one == 2
    #   
    #   num = (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]], 
    #                                         mean = mu[[1]], 
    #                                         sigma = Sigma[[1]]) 
    #   
    #   denom = num + (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]], 
    #                                                 mean = mu[[2]], 
    #                                                 sigma = Sigma[[2]])
    #   
    #   # will just be (1,0)... seems extreme
    #   ratio = c(num/denom, 1-(num/denom))
    #   
    # }
    
  } else{
    
    num = (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]][,1], 
                                          mean = mu[[1]][,1], 
                                          sigma = Sigma[[1]]) 
    # cat("\n num:",num,"\n")
    # print(sm_counts[1])
    # print(y[[obs]][,1])
    # print(mu[[1]][,1])
    # print(Sigma[[1]])
    denom = num + (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]][,1], 
                                                  mean = mu[[2]][,1], 
                                                  sigma = Sigma[[2]])
    # cat("\n denom:",denom,"\n")
    if(is.nan(num/denom)){
      cat("\n Warning: NaN detected in SM group assignment calculation. ")
      ratio = rep(1/2,2)
    } else{
      ratio = c(num/denom, 1-(num/denom))
    }
    
    
  }
  
  
  if(proposal_calc == TRUE){
    
    return(ratio[which_group_k])
    
  } else{
    
    return(ratio)
    
  }
  
  
  
}






############################ INDEPENDENT IG PRIORS ############################# 

MVN_CRP_nonconj_DEV <- function(S = 10^3, seed = 516, y, alpha = 1, 
                                a = 1, b = 10, mu0, sigma0, 
                                k_init = 3, init_method = "kmeans",
                                # d = 1, f = 1, 
                                sigma_hyperprior = FALSE, standardize_y = FALSE,
                                split_merge = FALSE, sm_iter = 5, truth = NA,
                                diag_weights = FALSE, 
                                verbose = TRUE, print_iter = 100, print_start = 0){
  
  # S is number of MCMC iterations
  # y is a list of data of length n
  # Sigma0 is the prior variance
  # alpha is the dirichlet process concentration parameter
  # a,b are prior variance hyperparameter on the IG distribution of sigma2. b is initial
  # value for b if sigma_hyperprior option is used
  # d, f are the hyperprior parameters on the Gamma dist placed on b, assume a known
  # sigma_hyperprior is a logical argument of whether to use the hyperprior or assume indep.
  # standardize_y is a logical argument for whether to center & scale data y before fitting model
  # mu0 is the prior mean - must be of dimension p*1
  # K_init is the initial number of groups
  # diag_weights is an argument to the Laplacian matrix - whether the diagonal should be 1 or not
  # verbose & printmod tells the function if and how often to print a progress summary
  # truth is an optional list argument with the true parameters used to generate simulated data
  
  start = Sys.time()
  
  set.seed(seed = seed)
  
  # data is y - a list of length n
  n = length(y)
  p = length(y[[1]]) # dimensionality of MVN
  k = k_init # initial number of groups
  y_matrix = matrix(data = unlist(y), ncol = p, byrow = TRUE)
  
  # center and scale data if standardize_y == TRUE
  if(standardize_y == TRUE){
    std_y_matrix = scale(y_matrix)
    y = lapply(X = 1:nrow(std_y_matrix), 
               FUN = function(x){matrix(std_y_matrix[x,], ncol=1)})
  }
  
  # preallocate memory and set initial values
  accept_ind = matrix(data = NA, nrow = S, ncol = n)
  group_assign = matrix(data = NA, nrow = S, ncol = n)
  
  if(init_method == "kmeans"){
    group_assign[1, ] = kmeans(x = y_matrix, centers = k_init, iter.max = 10)$cluster
  } else{
    # random
    group_assign[1, ] = sample(x = 1:k_init, size = length(y), 
                               replace = TRUE, prob = rep(1/k_init, k_init))
    # try different group assign initialization
    # group_assign[1, ] = ifelse(y > mean(y), k, k-1)  doesn't work for MVN, try kmeans?
  }
  
  
  means = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  vars = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  probs = vector(mode = "list", length = S)  # group assignment probabiltiies
  
  # split merge step - only used if needed
  sm_results = matrix(data = NA, nrow = 1, ncol = 7)
  colnames(sm_results) = c("s", "sm_iter", "move_type","accept", "prob", "k_start", "k_end")
  
  if(sigma_hyperprior == TRUE){
    extra_params = matrix(data = NA, nrow = S, ncol = 2)
    colnames(extra_params) = c("b")
    extra_params[1,1] = b # note that here b has dimension 1, diagonal var
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
  
  if(verbose == TRUE){
    cat("\n Initial means: \n")
    print(means[[1]])
    cat("\n Initial vars: \n")
    print(vars[[1]])
  }
  
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
  Sigma0 = diag(sigma0, p)
  curr_labels = 1:k
  avail_labels = c(1:n)[-curr_labels]
  
  # iterate 1:S
  for(s in 2:S){
    
    # print progress
    if((s %% print_iter == 0) & (s >= print_start) & (verbose == TRUE)){
      
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
        pr_res = group_prob_calc_DEV(k = k, n = n, n_j = count_assign, alpha = alpha, 
                                     y_i = y[[i]], mu = mu, sigma2 = sigma2,
                                     a = a, b = b, mu0 = mu0, sigma0 = sigma0, 
                                     singleton = 1)
        pr_c = pr_res$pr_c
        
      } else{
        
        #### calculate proposal distribution for group assignment
        #### if obs i is not presently a singleton
        pr_res = group_prob_calc_DEV(k = k, n = n, n_j = count_assign, alpha = alpha, 
                                     y_i = y[[i]], mu = mu, sigma2 = sigma2, 
                                     a = a, b = b, mu0 = mu0, sigma0 = sigma0, singleton = 0, 
                                     curr_group_assign = group_assign[s,i], 
                                     curr_labels = curr_labels)
        
        pr_c = pr_res$pr_c
        
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
        sigma2_k = pr_res$sigma2_new #1/rgamma(n = 1, shape = a, rate = b)
        sigma2 = c(sigma2, sigma2_k)
        
        #### draw a mean for newly created group from FC posterior of mu 
        mu_cov_k = 1/(1/sigma2_k + 1/sigma0)
        mu_mean_k = (y[[i]]/sigma2_k + mu0/sigma0)*mu_cov_k
        mu_k = matrix(mvtnorm::rmvnorm(n = 1, 
                                       mean = mu_mean_k, 
                                       sigma = diag(mu_cov_k, p)), 
                      nrow = p) # kth mean
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
    
    if((s %% print_iter == 0) & (s >= print_start) & (verbose == TRUE)){
      cat("\n End of CRP step, iter: ", s, "\n") # just create a new line for separation
      # print(paste("iter = ", s))
      print(paste("Current k = ", k))
      cat("\n")
      print(group_assign[s,])
      cat("\n")
      print(table(group_assign[s,]))
      cat("\n Curr labels: \n")
      print(curr_labels)
      cat("\n")
      print(mu)
      cat("\n")
      print(sigma2)
      cat("\n")
    }
    
    # final update of counts after a sweep
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    num_groups[s,] = k
    
    # proceed to split-merge step if true
    
    # Split-Merge step --- every 10 iterations
    
    if(split_merge == TRUE & (s %% sm_iter) == 0){
      
      k_start = k
      
      split_temp_group_assign = matrix(data = NA, nrow = (sm_iter + 1), ncol = length(group_assign[s,]))
      merge_temp_group_assign = matrix(data = NA, nrow = (sm_iter + 1), ncol = length(group_assign[s,]))
      
      split_sm_probs = matrix(data = NA, nrow = (sm_iter + 1), ncol = length(y))
      merge_sm_probs = matrix(data = NA, nrow = (sm_iter + 1), ncol = length(y))
      
      split_temp_group_assign[1,] = group_assign[s,]
      merge_temp_group_assign[1,] = group_assign[s,]
      
      split_means = vector(mode = "list", length = (sm_iter + 1)) 
      split_vars = vector(mode = "list", length = (sm_iter + 1)) 
      
      merge_means = vector(mode = "list", length = (sm_iter + 1)) 
      merge_vars = vector(mode = "list", length = (sm_iter + 1)) 
      
      # randomly select two observed data points y
      sampled_obs = sample(x = 1:length(y), size = 2, replace = FALSE)
      # cat("sampled_obs:", sampled_obs)
      # cat("\n")
      
      # check if in same group - if yes SPLIT
      # if no, MERGE
      lab1 = split_temp_group_assign[1, sampled_obs[1]]
      lab2 = split_temp_group_assign[1, sampled_obs[2]]
      move_type = ifelse(lab1 == lab2, "SPLIT", "MERGE")
      
      # cat("move_type:", move_type)
      # cat("\n")
      # cat("sampled_obs:", sampled_obs)
      # cat("\n")
      # cat("group_labs:", c(lab1, lab2))
      # cat("\n")
      # cat("curr_labels:", curr_labels)
      # cat("\n")
      
      # bookkeeping - group labels
      subset_index = which(split_temp_group_assign[1,] %in% c(lab1, lab2)) 
      existing_group_index = which(label_assign == lab1) 
      
      # identify original mean and variance for merge proposal
      original_param_index1 = which(curr_labels == lab1)
      original_mu1 = mu[,original_param_index1]
      original_sigma1 = sigma2[[original_param_index1]]
      
      original_param_index2 = which(curr_labels == lab2)
      original_mu2 = mu[,original_param_index2]
      original_sigma2 = sigma2[[original_param_index2]]
      
      # cat("split_labs:", split_lab)
      # cat("\n")
      # cat("subset_index:", subset_index)
      # cat("\n")
      
      
      # perform restricted Gibbs scans
      
      # if SPLIT
      if(move_type == "SPLIT"){
        
        # cat("\n Curr labels (at top of split): ", curr_labels, "\n")
        
        # need to set random launch states for both split and merge in non conj algo
        # specify random launch state for split
        split_lab = c(lab1, avail_labels[1]) # keep original label, new one for 2nd group
        
        # set launch states for merge
        merge_lab = lab1
        
        for(scan in 1:(sm_iter+1)){
          
          # initialize next scan with result of previous scan
          if(scan == 1){
            
            # intialize anchors so that there are no "zeros" for split groups
            # under consideration when starting below
            
            # split start
            split_temp_group_assign[scan,sampled_obs[1]] = split_lab[1] 
            split_temp_group_assign[scan,sampled_obs[2]] = split_lab[2] 
            
            # merge start
            merge_temp_group_assign[scan,sampled_obs[1]] = merge_lab
            merge_temp_group_assign[scan,sampled_obs[2]] = merge_lab
            
            # draw params from prior - random launch state for split proposal
            split_means[[1]] = lapply(X = 1:2, 
                                      FUN = function(x){
                                        t(mvtnorm::rmvnorm(n = 1, 
                                                           mean = mu0,
                                                           sigma = Sigma0))
                                      })
            
            split_vars[[1]] = lapply(X = 1:2, 
                                     FUN = function(x){
                                       diag(1/rgamma(n = 1, shape = a, rate = b), p)
                                       # for UVV
                                       # LaplacesDemon::rinvwishart(nu = nu, 
                                       #                            S = lambda0)
                                     })
            
            # cat("\n Split means (init): \n")
            # print(split_means[[1]])
            # cat("\n Split vars: \n")
            # print(split_vars[[1]])
            
            # draw params from prior - random launch state for merge proposal
            merge_means[[1]] = t(mvtnorm::rmvnorm(n = 1, mean = mu0, sigma = Sigma0))
            merge_vars[[1]] = diag(1/rgamma(n = 1, shape = a, rate = b), p)  
            # merge_vars[[1]] = LaplacesDemon::rinvwishart(nu = nu, S = lambda0)
            
            # cat("\n Merge means (init): \n")
            # print(merge_means[[1]])
            # cat("\n Merge vars: \n")
            # print(merge_vars[[1]])
            # 
          } else{
            
            # initialize current restricted Gibbs scan iteration with previous result
            split_temp_group_assign[scan,] = split_temp_group_assign[(scan-1),] 
            merge_temp_group_assign[scan,] = merge_temp_group_assign[(scan-1),] 
            
            split_means[[scan]] = split_means[[scan-1]]
            split_vars[[scan]] = split_vars[[scan-1]] 
            
            merge_means[[scan]] = merge_means[[scan-1]]
            merge_vars[[scan]] = merge_vars[[scan-1]] 
            
            # moved this up, was previously nested below, which means this step would
            # be *repeated* for each observation, thus washing out any split/merge that
            # was occurring as a result of the algorithm at each step
          }
          
          # print(table(temp_group_assign[scan,]))
          
          for(obs in subset_index){
            
            # cat("Obs:", obs)
            
            if(scan == 1){
              
              # specify component assignment - random launch state for split proposal
              # yes this is redundant, but only way to make using subset_index instead
              # of subset_index_minus in for loop work properly 
              if(obs == sampled_obs[2]){
                split_temp_group_assign[scan,obs] = split_lab[2] 
                merge_temp_group_assign[scan,obs] = merge_lab
              } else if(obs == sampled_obs[1]){
                split_temp_group_assign[scan,obs] = split_lab[1] 
                merge_temp_group_assign[scan,obs] = merge_lab
              } else{
                # random launch state
                split_temp_group_assign[scan,obs] = sample(x = split_lab, size = 1)
                # no need to sample for merge...only one option
                merge_temp_group_assign[scan,obs] = merge_lab
              }
              
            } else{ # for remaining scans after random launch state set
              # perform restricted gibbs scans for both component assign and params
              
              if(obs %in% sampled_obs){
                # cat(" Anchor, ")
                # dont sample if anchor obs -- assignment cannot change
                # specify new group label to 2nd anchor point as well
                if(obs == sampled_obs[2]){
                  split_temp_group_assign[scan,obs] = split_lab[2] 
                  merge_temp_group_assign[scan,obs] = merge_lab
                } else if(obs == sampled_obs[1]){
                  split_temp_group_assign[scan,obs] = split_lab[1] 
                  merge_temp_group_assign[scan,obs] = merge_lab
                }
                
                # update parameter vector
                
                # cat(" Assign:",temp_group_assign[scan,obs])
                # cat("\n")
                
                # update group assignment probabilities 
                
                split_assign_prob = nonconj_component_prob_c(
                  obs = obs, split_labs = split_lab,
                  group_assign = split_temp_group_assign[scan,], 
                  y = y, mu = split_means[[scan]], Sigma = split_vars[[scan]])
                
                ### merge_assign prob by definition always 1 when proposing merge
                ### launch state from c_i = c_j, no other way to permute assignment
                merge_assign_prob = 1
                
                # dont sample --- but do record prob of group anchor obs in in!
                which_split_labs_anchor = which(split_lab == split_temp_group_assign[scan,obs])
                # how to do this best for merge????
                which_merge_lab_anchor = 1  # why is this 1?? need to figure out what this should be 
                
                # temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                # dont need to assign again - already initialized since anchor
                split_sm_probs[scan,obs] = split_assign_prob[which_split_labs_anchor]
                merge_sm_probs[scan,obs] = merge_assign_prob[which_merge_lab_anchor]
                
              } else{
                
                split_count_assign = as.numeric(table(split_temp_group_assign[scan,]))
                split_label_assign = as.numeric(names(table(split_temp_group_assign[scan,])))
                # split_singletons = split_label_assign[which(split_count_assign == 1)]
                
                merge_count_assign = as.numeric(table(merge_temp_group_assign[scan,]))
                merge_label_assign = as.numeric(names(table(merge_temp_group_assign[scan,])))
                # merge_singletons = merge_label_assign[which(merge_count_assign == 1)]
                # should never be a merge singleton --- need at least 2 anchor observations
                
                # not subtracted 1 for kth observation here -- do within function
                split_counts = table(split_temp_group_assign[scan,])
                merge_counts = table(merge_temp_group_assign[scan,])
                
                #current_obs_index = which(temp_group_assign[scan,] == obs)
                #split_group_lab_index1 = which(temp_group_assign[scan,] == split_lab[1])
                #split_group_lab_index2 = which(temp_group_assign[scan,] == split_lab[2])
                
                # current observation under consideration cannot be included here
                
                split_assign_prob = nonconj_component_prob_c(
                  obs = obs, split_labs = split_lab,
                  group_assign = split_temp_group_assign[scan,], 
                  y = y, mu = split_means[[scan]], Sigma = split_vars[[scan]])
                
                sm_prop_index = sample(x = 1:2, size = 1, 
                                       prob = split_assign_prob)
                
                split_temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                split_sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
                
                ### merge_assign prob by definition always 1 when proposing merge
                ### launch state from c_i = c_j, no other way to permute assignment
                merge_temp_group_assign[scan,obs] = merge_lab
                merge_sm_probs[scan,obs] = 1
                
                
              }
              
            } 
            
          } # iterate through all observations in the two split groups under consideration
          
          # update phi -- after each scan
          
          # update counts after scan
          split_count_assign = as.numeric(table(split_temp_group_assign[scan,]))
          split_label_assign = as.numeric(names(table(split_temp_group_assign[scan,])))
          # split_singletons = split_label_assign[which(split_count_assign == 1)]
          
          merge_count_assign = as.numeric(table(merge_temp_group_assign[scan,]))
          merge_label_assign = as.numeric(names(table(merge_temp_group_assign[scan,])))
          # merge_singletons = merge_label_assign[which(merge_count_assign == 1)]
          # should never be a merge singleton --- need at least 2 anchor observations
          
          split_counts = table(split_temp_group_assign[scan,])
          merge_counts = table(merge_temp_group_assign[scan,])
          
          split_group_count_index = which(as.numeric(names(split_counts)) %in% split_lab)
          merge_group_count_index = which(as.numeric(names(merge_counts)) %in% merge_lab)
          
          split_phi = lapply(X = 1:2, #split_lab, 
                             FUN = function(x){
                               update_phi_DEV(curr_label = split_lab[x], 
                                              group_assign = split_temp_group_assign[scan,], 
                                              count_assign = split_count_assign[split_group_count_index][x], 
                                              y = y, 
                                              mu = split_means[[scan]][[x]], 
                                              mu0 = mu0, 
                                              Sigma = split_vars[[scan]][[x]], 
                                              Sigma0 = Sigma0, a = a, b = b)
                             })
          
          split_means[[scan]] = list(split_phi[[1]]$mu, split_phi[[2]]$mu)
          split_vars[[scan]] = list(split_phi[[1]]$Sigma, split_phi[[2]]$Sigma)
          
          if(scan == sm_iter + 1){
            # stop, record launch state (previous scan result), no need to do extra scan
            
            merge_means[[scan]] = merge_phi$mu
            merge_vars[[scan]] = merge_phi$Sigma
            
          } else{
            
            merge_phi = update_phi_DEV(curr_label = merge_lab, 
                                       group_assign = merge_temp_group_assign[scan,], 
                                       count_assign = merge_count_assign[merge_group_count_index], 
                                       y = y, 
                                       mu = merge_means[[scan]], 
                                       mu0 = mu0, 
                                       Sigma = merge_vars[[scan]], 
                                       Sigma0 = Sigma0, a = a, b = b)
            
            merge_means[[scan]] = merge_phi$mu
            merge_vars[[scan]] = merge_phi$Sigma
            
          }
          
          
          # cat("\n Updated split means: \n")
          # print(split_means[[scan]])
          # cat("\n Updated split vars: \n")
          # print(split_vars[[scan]])
          # 
          # cat("\n Updated merge means: \n")
          # print(merge_means[[scan]])
          # cat("\n Updated merge vars: \n")
          # print(merge_vars[[scan]])
          
        } # scans 1:(sm_iter+1)
        
        
        # calculate & evaluate acceptance prob
        
        # update counts after scans
        # split_count_assign = as.numeric(table(split_temp_group_assign[sm_iter+1,]))
        # split_label_assign = as.numeric(names(table(split_temp_group_assign[sm_iter+1,])))
        # # split_singletons = split_label_assign[which(split_count_assign == 1)]
        # 
        # merge_count_assign = as.numeric(table(merge_temp_group_assign[sm_iter+1,]))
        # merge_label_assign = as.numeric(names(table(merge_temp_group_assign[sm_iter+1,])))
        # # merge_singletons = merge_label_assign[which(merge_count_assign == 1)]
        # # should never be a merge singleton --- need at least 2 anchor observations
        # 
        split_counts = table(split_temp_group_assign[sm_iter+1,])
        merge_counts = table(merge_temp_group_assign[sm_iter+1,])
        
        split_lab_assign = as.numeric(names(split_counts))
        merge_lab_assign = as.numeric(names(merge_counts))
        
        split_count_assign = as.numeric(split_counts)
        merge_count_assign = as.numeric(merge_counts)
        
        split_group_count_index = which(split_lab_assign %in% split_lab)
        merge_group_count_index = which(merge_lab_assign %in% merge_lab)
        
        ## proposal probability
        
        # compute P_GS(phi) from launch state to final scan for both split and merge proposals
        # cat("\n split_lab", split_lab, "\n")
        # cat("\n index", split_group_count_index, "\n")
        # cat("\n count", split_count_assign, "\n")
        # cat("\n group_assign", split_temp_group_assign[sm_iter+1,], "\n")
        # cat("\n split_means", "\n")
        # print(split_means[[scan]])
        # cat("\n split_vars", "\n")
        # print(split_vars[[scan]])
        
        split_phi_prob = sapply(X = 1:2, 
                                FUN = function(x){
                                  nonconj_phi_prob_DEV(
                                    curr_label = split_lab[x],
                                    group_assign = split_temp_group_assign[sm_iter+1,], 
                                    count_assign = split_count_assign[split_group_count_index][x], 
                                    y = y, 
                                    mu_L = split_means[[scan-1]][[x]], # launch state
                                    mu = split_means[[scan]][[x]], # final proposal
                                    mu0 = mu0, 
                                    Sigma_L = split_vars[[scan-1]][[x]], # launch state
                                    Sigma = split_vars[[scan]][[x]], # final proposal
                                    Sigma0 = Sigma0, a = a, b = b)
                                })
        
        # cat("\n merge_lab", merge_lab, "\n")
        # cat("\n count_assign", merge_count_assign, "\n")
        # cat("\n index", merge_group_count_index, "\n")
        # cat("\n group_assign", merge_temp_group_assign[sm_iter+1,], "\n")
        # cat("\n merge_means", "\n")
        # print(merge_means[[scan]])
        # cat("\n merge_vars", "\n")
        # print(merge_vars[[scan]])
        
        # when doing a SPLIT proposal, this is the \phi component of 
        # q(\gamma|\gamma^{split}) - i.e. the Gibbs sampling transition kernel
        # from the merge launch state to the *original merged configuration*
        merge_phi_prob = nonconj_phi_prob_DEV(curr_label = merge_lab, 
                                              group_assign = merge_temp_group_assign[sm_iter+1,], 
                                              count_assign = merge_count_assign[merge_group_count_index], 
                                              y = y, 
                                              mu_L = merge_means[[scan-1]], # launch state
                                              mu = original_mu1,
                                              mu0 = mu0, 
                                              Sigma_L = merge_vars[[scan-1]], # launch state
                                              Sigma = diag(original_sigma1,p),
                                              Sigma0 = Sigma0, a = a, b = b)
        
        # for reverse split proposal, only one way to arrange --- merge, so overall prob = 1
        prob1_c_num = log(1) # Reduce(f = "+", x = log(merge_sm_probs[sm_iter+1,subset_index]))
        # check to prevent numeric overflow from small densities (happens when proposal is two large
        # groups that are not compatible)
        prob1_phi_num = ifelse(merge_phi_prob < 10^(-300), log(10^(-300)), merge_phi_prob)  
        
        prob1_c_denom = Reduce(f = "+", x = log(split_sm_probs[sm_iter+1,subset_index]))
        
        # check to prevent numeric overflow from small densities (happens when proposal is two large
        # groups that are not compatible)
        if(any(split_phi_prob < 10^(-300)) == TRUE){
          
          which_below_tol = which(split_phi_prob < 10^(-300))
          split_phi_prob[which_below_tol] = 10^-300 # if below tol, set equal to tol
          
        } 
        
        # else, business as usual
        prob1_phi_denom = Reduce(f = "+", x = split_phi_prob) # only calculated at end
        # so no need to index
        
        prob1 = (prob1_c_num + prob1_phi_num) - (prob1_c_denom + prob1_phi_denom)
        
        ## prior ratio
        ## don't log nonconj_prior_dens again, already taking log in fxn!
        prob2_num = sum(log(1:(split_counts[[split_group_count_index[1]]]-1))) + 
          sum(log(1:(split_counts[[split_group_count_index[2]]]-1))) + 
          nonconj_prior_dens_DEV(mu = split_means[[scan]][[1]], mu0 = mu0, 
                                 Sigma = split_vars[[scan]][[1]], 
                                 Sigma0 = Sigma0, a = a, b = b) + 
          nonconj_prior_dens_DEV(mu = split_means[[scan]][[2]], mu0 = mu0, 
                                 Sigma = split_vars[[scan]][[2]], 
                                 Sigma0 = Sigma0, a = a, b = b)
        
        prob2_denom = sum(log(1:(split_counts[[split_group_count_index[1]]] + 
                                   split_counts[[split_group_count_index[2]]]-1))) +
          nonconj_prior_dens_DEV(mu = original_mu1, mu0 = mu0, 
                                 Sigma = diag(original_sigma1,p), 
                                 Sigma0 = Sigma0, a = a, b = b)
        
        prob2 = log(alpha) + prob2_num - prob2_denom
        
        ## likelihood ratio
        subset_index_grp1 = which(split_temp_group_assign[sm_iter+1,] %in% split_lab[1]) 
        subset_index_grp2 = which(split_temp_group_assign[sm_iter+1,] %in% split_lab[2]) 
        
        ### component 1 - numerator I (group 1 - split proposal)
        prob3_num1 = 0
        for(obs_ind in 1:length(subset_index_grp1)){
          val = ll_components_DEV(obs_ind = subset_index_grp1[obs_ind], y = y, 
                                  mu = split_means[[sm_iter+1]][[1]], 
                                  Sigma = split_vars[[sm_iter+1]][[1]])
          prob3_num1 = prob3_num1 + log(val)
        }
        
        ### component 2 - numerator II (group 2 - split proposal)
        prob3_num2 = 0
        for(obs_ind in 1:length(subset_index_grp2)){
          val = ll_components_DEV(obs_ind = subset_index_grp2[obs_ind], y = y, 
                                  mu = split_means[[sm_iter+1]][[2]], 
                                  Sigma = split_vars[[sm_iter+1]][[2]])
          prob3_num2 = prob3_num2 + log(val)
        }
        
        
        ### component 3 - denominator (all in original group w/ merge proposal params)
        prob3_denom = 0
        for(obs_ind in 1:length(subset_index)){
          val = ll_components_DEV(obs_ind = subset_index[obs_ind], y = y, 
                                  mu = original_mu1, 
                                  Sigma = diag(original_sigma1, p))
          prob3_denom = prob3_denom + log(val)
        }
        
        prob3 = prob3_num1 + prob3_num2 - prob3_denom
        
        ## evaluate acceptance prob
        prob_components = c(prob1, prob2, prob3)
        
        ## check for indeterminate forms in accept prob - will cause errors
        if(is.nan(exp(prob1 + prob2 + prob3)) == TRUE){
          cat("move_type:", move_type)
          cat("\n")
          cat("\n split_counts: \n")
          print(split_counts)
          cat("\n log prob_components = ", c(prob1, prob2, prob3), "\n")
          cat("\n exp accept_prob = ", exp(prob1 + prob2 + prob3))
          
          prob_sign = sign(c(prob1, prob2, prob3))
          inf_index = which(is.infinite(prob_components))
          
          # if NaN but no inf components detected - throw error
          if(length(inf_index) == 0){
            stop("NaN detected in SM acceptance probability. Unable to resolve.")
          }
          
          # if inf detected, resolve indeterminate forms by setting Inf equal to 
          # some arbitrarily large number
          for(prob_i in 1:length(inf_index)){
            prob_components[inf_index[prob_i]] = prob_sign[inf_index[prob_i]]*5000
          }
        }
        
        
        accept_prob = min(1, exp(sum(prob_components)))
        # cat("\n accept prob:", round(accept_prob,5), "\n")
        u = runif(n = 1)
        if(accept_prob > u){
          
          # accept
          accept = 1
          group_assign[s,] = split_temp_group_assign[sm_iter+1,]
          
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
          
          # if new group created by split, update mean and variance
          ## add new means and variances from final Gibbs scan to relevant vectors/lists
          length_sigma2 = length(sigma2)
          sigma2[[which_split_labs[1]]] = split_vars[[sm_iter+1]][[1]][1] # scalar in DEV case
          sigma2[[length_sigma2+1]] = split_vars[[sm_iter+1]][[2]][1] # scalar in DEV case
          
          mu[,which_split_labs[1]] = split_means[[sm_iter+1]][[1]]
          mu = cbind(mu, split_means[[sm_iter+1]][[2]])
          
        } else{
          # reject
          accept = 0
          # group assign remains unchanged
        }
        
        
        # cat("\n accept = ", accept, "\n")
        
        # if MERGE    
      } else if(move_type == "MERGE"){
        
        #cat("\n Curr labels (at top of merge): ", curr_labels, "\n")
        
        # need to set random launch states for both split and merge in non conj algo
        # specify random launch state for split
        split_lab = c(lab1, lab2) # keep original labels for both groups for split 
        which_split_labs = which(label_assign %in% split_lab) # check before updating
        # launch state in merge proposal
        
        # set launch states for merge
        merge_lab = lab1
        
        for(scan in 1:(sm_iter+1)){
          
          # initialize next scan with result of previous scan
          if(scan == 1){
            
            # intialize anchors so that there are no "zeros" for split groups
            # under consideration when starting below
            
            # split start
            split_temp_group_assign[scan,sampled_obs[1]] = split_lab[1] 
            split_temp_group_assign[scan,sampled_obs[2]] = split_lab[2] 
            
            # merge start
            merge_temp_group_assign[scan,sampled_obs[1]] = merge_lab
            merge_temp_group_assign[scan,sampled_obs[2]] = merge_lab
            
            # draw params from prior - random launch state for split proposal
            split_means[[1]] = lapply(X = 1:2, 
                                      FUN = function(x){
                                        t(mvtnorm::rmvnorm(n = 1, 
                                                           mean = mu0,
                                                           sigma = Sigma0))
                                      })
            
            split_vars[[1]] = lapply(X = 1:2, 
                                     FUN = function(x){
                                       diag(1/rgamma(n = 1, shape = a, rate = b), p)
                                       # for UVV
                                       # LaplacesDemon::rinvwishart(nu = nu, 
                                       #                            S = lambda0)
                                     })
            
            # draw params from prior - random launch state for merge proposal
            merge_means[[1]] = t(mvtnorm::rmvnorm(n = 1, mean = mu0, sigma = Sigma0))
            merge_vars[[1]] = diag(1/rgamma(n = 1, shape = a, rate = b), p)  
            # merge_vars[[1]] = LaplacesDemon::rinvwishart(nu = nu, S = lambda0)
            
          } else{
            
            # initialize current restricted Gibbs scan iteration with previous result
            split_temp_group_assign[scan,] = split_temp_group_assign[(scan-1),] 
            merge_temp_group_assign[scan,] = merge_temp_group_assign[(scan-1),] 
            
            split_means[[scan]] = split_means[[scan-1]]
            split_vars[[scan]] = split_vars[[scan-1]] 
            
            merge_means[[scan]] = merge_means[[scan-1]]
            merge_vars[[scan]] = merge_vars[[scan-1]] 
            
            # moved this up, was previously nested below, which means this step would
            # be *repeated* for each observation, thus washing out any split/merge that
            # was occurring as a result of the algorithm at each step
          }
          
          # print(table(temp_group_assign[scan,]))
          
          for(obs in subset_index){
            
            # cat("Obs:", obs)
            
            if(scan == 1){
              
              # specify component assignment - random launch state for split proposal
              # yes this is redundant, but only way to make using subset_index instead
              # of subset_index_minus in for loop work properly 
              if(obs == sampled_obs[2]){
                split_temp_group_assign[scan,obs] = split_lab[2] 
                merge_temp_group_assign[scan,obs] = merge_lab
              } else if(obs == sampled_obs[1]){
                split_temp_group_assign[scan,obs] = split_lab[1] 
                merge_temp_group_assign[scan,obs] = merge_lab
              } else{
                # random launch state
                split_temp_group_assign[scan,obs] = sample(x = split_lab, size = 1)
                # no need to sample for merge...only one option
                merge_temp_group_assign[scan,obs] = merge_lab
              }
              
            } else{ # for remaining scans after random launch state set
              # perform restricted gibbs scans for both component assign and params
              
              if(obs %in% sampled_obs){
                # cat(" Anchor, ")
                # dont sample if anchor obs -- assignment cannot change
                # specify new group label to 2nd anchor point as well
                if(obs == sampled_obs[2]){
                  split_temp_group_assign[scan,obs] = split_lab[2] 
                  merge_temp_group_assign[scan,obs] = merge_lab
                } else if(obs == sampled_obs[1]){
                  split_temp_group_assign[scan,obs] = split_lab[1] 
                  merge_temp_group_assign[scan,obs] = merge_lab
                }
                
                # update parameter vector
                
                # cat(" Assign:",temp_group_assign[scan,obs])
                # cat("\n")
                
                # update group assignment probabilities 
                
                split_assign_prob = nonconj_component_prob_c(
                  obs = obs, split_labs = split_lab,
                  group_assign = split_temp_group_assign[scan,], 
                  y = y, mu = split_means[[scan]], Sigma = split_vars[[scan]])
                
                ### merge_assign prob by definition always 1 when proposing merge
                ### launch state from c_i = c_j, no other way to permute assignment
                merge_assign_prob = 1
                
                # dont sample --- but do record prob of group anchor obs in in!
                which_split_labs_anchor = which(split_lab == split_temp_group_assign[scan,obs])
                # how to do this best for merge????
                which_merge_lab_anchor = 1  # why is this 1?? need to figure out what this should be 
                
                # temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                # dont need to assign again - already initialized since anchor
                split_sm_probs[scan,obs] = split_assign_prob[which_split_labs_anchor]
                merge_sm_probs[scan,obs] = merge_assign_prob[which_merge_lab_anchor]
                
              } else{
                
                split_count_assign = as.numeric(table(split_temp_group_assign[scan,]))
                split_label_assign = as.numeric(names(table(split_temp_group_assign[scan,])))
                # split_singletons = split_label_assign[which(split_count_assign == 1)]
                
                merge_count_assign = as.numeric(table(merge_temp_group_assign[scan,]))
                merge_label_assign = as.numeric(names(table(merge_temp_group_assign[scan,])))
                # merge_singletons = merge_label_assign[which(merge_count_assign == 1)]
                # should never be a merge singleton --- need at least 2 anchor observations
                
                # not subtracted 1 for kth observation here -- do within function
                split_counts = table(split_temp_group_assign[scan,])
                merge_counts = table(merge_temp_group_assign[scan,])
                
                #current_obs_index = which(temp_group_assign[scan,] == obs)
                #split_group_lab_index1 = which(temp_group_assign[scan,] == split_lab[1])
                #split_group_lab_index2 = which(temp_group_assign[scan,] == split_lab[2])
                
                # current observation under consideration cannot be included here
                
                if(scan == sm_iter + 1){
                  # final scan is hypothetical transition from split launch state
                  # back to original split configuration (use \phi from launch
                  # state but calculate prob of going back to original group assigns)
                  
                  split_assign_prob = nonconj_component_prob_c(
                    obs = obs, split_labs = split_lab, y = y, 
                    group_assign = group_assign[s,], # original group assign 
                    mu = split_means[[scan]], Sigma = split_vars[[scan]],
                    proposal_calc = TRUE)
                  
                  # sm_prop_index = sample(x = 1:2, size = 1, 
                  #                        prob = split_assign_prob)
                  
                  # dont actually change assignments from launch state
                  split_temp_group_assign[scan,obs] =  split_temp_group_assign[scan-1,obs]
                  
                  # record probs of returning to original configuration
                  split_sm_probs[scan,obs] = split_assign_prob # function should
                  # only return a single prob when proposal_calc == TRUE
                  
                } else{
                  
                  split_assign_prob = nonconj_component_prob_c(
                    obs = obs, split_labs = split_lab,
                    group_assign = split_temp_group_assign[scan,], 
                    y = y, mu = split_means[[scan]], Sigma = split_vars[[scan]])
                  
                  sm_prop_index = sample(x = 1:2, size = 1, 
                                         prob = split_assign_prob)
                  
                  split_temp_group_assign[scan,obs] = split_lab[sm_prop_index]
                  split_sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
                  
                }

                
                ### merge_assign prob by definition always 1 when proposing merge
                ### launch state from c_i = c_j, no other way to permute assignment
                merge_temp_group_assign[scan,obs] = merge_lab
                merge_sm_probs[scan,obs] = 1
                
                
              }
              
            } 
            
          } # iterate through all observations in the two split groups under consideration
          
          # update phi -- after each scan
          
          # update counts after scan
          split_count_assign = as.numeric(table(split_temp_group_assign[scan,]))
          split_label_assign = as.numeric(names(table(split_temp_group_assign[scan,])))
          # split_singletons = split_label_assign[which(split_count_assign == 1)]
          
          merge_count_assign = as.numeric(table(merge_temp_group_assign[scan,]))
          merge_label_assign = as.numeric(names(table(merge_temp_group_assign[scan,])))
          # merge_singletons = merge_label_assign[which(merge_count_assign == 1)]
          # should never be a merge singleton --- need at least 2 anchor observations
          
          split_counts = table(split_temp_group_assign[scan,])
          merge_counts = table(merge_temp_group_assign[scan,])
          
          split_group_count_index = which(as.numeric(names(split_counts)) %in% split_lab)
          merge_group_count_index = which(as.numeric(names(merge_counts)) %in% merge_lab)
          
          if(scan == sm_iter + 1){
            # don't do final scan - stop at split launch state
            split_means[[scan]] = split_means[[scan-1]] 
            split_vars[[scan]] = split_vars[[scan-1]] 
            
          } else{
            
            split_phi = lapply(X = 1:2, #split_lab, 
                               FUN = function(x){
                                 update_phi_DEV(curr_label = split_lab[x], 
                                                group_assign = split_temp_group_assign[scan,], 
                                                count_assign = split_count_assign[split_group_count_index][x], 
                                                y = y, 
                                                mu = split_means[[scan]][[x]], 
                                                mu0 = mu0, 
                                                Sigma = split_vars[[scan]][[x]], 
                                                Sigma0 = Sigma0, a = a, b = b)
                               })
            
            split_means[[scan]] = list(split_phi[[1]]$mu, split_phi[[2]]$mu)
            split_vars[[scan]] = list(split_phi[[1]]$Sigma, split_phi[[2]]$Sigma)
            
          }
          
          
          merge_phi = update_phi_DEV(curr_label = merge_lab, 
                                     group_assign = merge_temp_group_assign[scan,], 
                                     count_assign = merge_count_assign[merge_group_count_index], 
                                     y = y, 
                                     mu = merge_means[[scan]], 
                                     mu0 = mu0, 
                                     Sigma = merge_vars[[scan]], 
                                     Sigma0 = Sigma0, a = a, b = b)
          
          merge_means[[scan]] = merge_phi$mu
          merge_vars[[scan]] = merge_phi$Sigma
          
        } # scans 1:(sm_iter+1)
        
        
        # calculate & evaluate acceptance prob
        
        # update counts after scans
        split_counts = table(split_temp_group_assign[sm_iter+1,])
        merge_counts = table(merge_temp_group_assign[sm_iter+1,])
        
        split_lab_assign = as.numeric(names(split_counts))
        merge_lab_assign = as.numeric(names(merge_counts))
        
        split_count_assign = as.numeric(split_counts)
        merge_count_assign = as.numeric(merge_counts)
        
        split_group_count_index = which(split_lab_assign %in% split_lab)
        merge_group_count_index = which(merge_lab_assign %in% merge_lab)
        
        ## proposal probability
        
        # compute P_GS(phi) from launch state to final scan for both split and merge proposals
        
        # cat("\n split_lab", split_lab, "\n")
        # cat("\n index", split_group_count_index, "\n")
        # cat("\n count", split_count_assign, "\n")
        # cat("\n group_assign", split_temp_group_assign[sm_iter+1,], "\n")
        # cat("\n split_means", "\n")
        # print(split_means[[scan]])
        # cat("\n split_vars", "\n")
        # print(split_vars[[scan]])
        
        ## proposal probability
        split_phi_prob = sapply(X = 1:2, 
                                FUN = function(x){
                                  nonconj_phi_prob_DEV(
                                    curr_label = split_lab[x],
                                    group_assign = split_temp_group_assign[sm_iter+1,], 
                                    count_assign = split_count_assign[split_group_count_index][x], 
                                    y = y, 
                                    mu_L = split_means[[scan-1]][[x]], 
                                    mu = list(original_mu1, original_mu2)[[x]],
                                    mu0 = mu0, 
                                    Sigma_L = split_vars[[scan-1]][[x]], 
                                    Sigma = list(diag(original_sigma1, p),
                                                 diag(original_sigma2, p))[[x]],
                                    Sigma0 = Sigma0, a = a, b = b)
                                })
        
        # cat("\n merge_lab", merge_lab, "\n")
        # cat("\n index", merge_group_count_index, "\n")
        # cat("\n count", merge_count_assign, "\n")
        # cat("\n group_assign", merge_temp_group_assign[sm_iter+1,], "\n")
        # cat("\n merge_means", "\n")
        # print(merge_means[[scan]])
        # cat("\n merge_vars", "\n")
        # print(merge_vars[[scan]])
        
        ## proposal probability
        
        merge_phi_prob = nonconj_phi_prob_DEV(curr_label = merge_lab, 
                                              group_assign = merge_temp_group_assign[sm_iter+1,], 
                                              count_assign = merge_count_assign[merge_group_count_index], 
                                              y = y, 
                                              mu_L = merge_means[[scan-1]], 
                                              mu = merge_means[[scan]], 
                                              mu0 = mu0, 
                                              Sigma_L = merge_vars[[scan-1]], 
                                              Sigma = merge_vars[[scan]], 
                                              Sigma0 = Sigma0, a = a, b = b)
        #cat("\n merge_phi_prob: ", merge_phi_prob, "\n")
        
        prob1_c_denom = log(1) # Reduce(f = "+", x = log(merge_sm_probs[sm_iter+1,subset_index]))
        prob1_phi_denom = ifelse(merge_phi_prob < 10^(-300), log(10^(-300)), merge_phi_prob)  
        
        prob1_c_num = Reduce(f = "+", x = log(split_sm_probs[sm_iter+1,subset_index]))
        
        if(any(split_phi_prob < 10^(-300)) == TRUE){
          
          which_below_tol = which(split_phi_prob < 10^(-300))
          split_phi_prob[which_below_tol] = 10^-300 # if below tol, set equal to tol
          
        } 
        
        prob1_phi_num = Reduce(f = "+", x = split_phi_prob)
        # cat("\n split_phi_prob: ", split_phi_prob, "\n")
        
        prob1 = (prob1_c_num + prob1_phi_num) - (prob1_c_denom + prob1_phi_denom)
        
        ## prior ratio
        ## dont log prior density again here -- already taking log in fxn!
        prob2_num = sum(log(1:(split_counts[[split_group_count_index[1]]] + 
                                 split_counts[[split_group_count_index[2]]]-1))) +
          nonconj_prior_dens_DEV(mu = merge_means[[scan]], mu0 = mu0, 
                                 Sigma = merge_vars[[scan]], 
                                 Sigma0 = Sigma0, a = a, b = b)
        
        # cat("\n prob2_num", prob2_num, "\n")
        # cat("\n merge factorial: ", sum(log(1:(split_counts[[split_group_count_index[1]]] + 
        #                                          split_counts[[split_group_count_index[2]]]-1))), "\n")
        # cat("\n num dens: ", log(nonconj_prior_dens_DEV(mu = merge_means[[scan]], mu0 = mu0, 
        #                                                 Sigma = merge_vars[[scan]], 
        #                                                 Sigma0 = Sigma0, a = a, b = b)), "\n")
        
        prob2_denom = sum(log(1:(split_counts[[split_group_count_index[1]]]-1))) + 
          sum(log(1:(split_counts[[split_group_count_index[2]]]-1))) +
          nonconj_prior_dens_DEV(mu = original_mu1, mu0 = mu0, 
                                 Sigma = diag(original_sigma1,p), 
                                 Sigma0 = Sigma0, a = a, b = b) +
          nonconj_prior_dens_DEV(mu = original_mu2, mu0 = mu0, 
                                 Sigma = diag(original_sigma2,p), 
                                 Sigma0 = Sigma0, a = a, b = b)
        
        # cat("\n prob2_denom", prob2_denom)
        # cat("\n split factorial: ", sum(log(1:(split_counts[[split_group_count_index[1]]]-1))) + 
        #       sum(log(1:(split_counts[[split_group_count_index[2]]]-1))), "\n")
        # cat("\n num dens 1: ", nonconj_prior_dens_DEV(mu = original_mu1, mu0 = mu0,
        #                                               Sigma = diag(original_sigma1,p),
        #                                               Sigma0 = Sigma0, a = a, b = b), "\n")
        # cat("\n num dens 2: ", nonconj_prior_dens_DEV(mu = original_mu2, mu0 = mu0,
        #                                               Sigma = diag(original_sigma2,p),
        #                                               Sigma0 = Sigma0, a = a, b = b), "\n")
        
        prob2 = log(alpha) + prob2_num - prob2_denom
        
        ## likelihood ratio
        subset_index_grp1 = which(split_temp_group_assign[sm_iter+1,] %in% split_lab[1]) 
        subset_index_grp2 = which(split_temp_group_assign[sm_iter+1,] %in% split_lab[2]) 
        
        ### component 1 - numerator I (group 1 - split proposal)
        prob3_num1 = 0
        for(obs_ind in 1:length(subset_index_grp1)){
          val = ll_components_DEV(obs_ind = subset_index_grp1[obs_ind], y = y, 
                                  mu = original_mu1, 
                                  Sigma = diag(original_sigma1,p))
          prob3_num1 = prob3_num1 + log(val)
        }
        
        ### component 2 - numerator II (group 2 - split proposal)
        prob3_num2 = 0
        for(obs_ind in 1:length(subset_index_grp2)){
          val = ll_components_DEV(obs_ind = subset_index_grp2[obs_ind], y = y, 
                                  mu = original_mu2, # change, needs tot be merge mu
                                  Sigma = diag(original_sigma2,p))
          prob3_num2 = prob3_num2 + log(val)
        }
        
        
        ### component 3 - denominator (all in original group w/ merge proposal params)
        prob3_denom = 0
        for(obs_ind in 1:length(subset_index)){
          val = ll_components_DEV(obs_ind = subset_index[obs_ind], y = y, 
                                  mu = merge_means[[sm_iter+1]], 
                                  Sigma = merge_vars[[sm_iter+1]])
          prob3_denom = prob3_denom + log(val)
        }
        
        # cat("\n")
        # cat("merge_vars[[sm_iter+1]]", merge_vars[[sm_iter+1]])
        # cat("\n")
        # cat("merge_means[[sm_iter+1]]", merge_means[[sm_iter+1]])
        
        # flip this for merge step
        prob3 = prob3_denom - (prob3_num1 + prob3_num2)
        
        ## evaluate acceptance prob
        prob_components = c(prob1, prob2, prob3)
        
        ## check for indeterminate forms in accept prob - will cause errors
        if(is.nan(exp(prob1 + prob2 + prob3)) == TRUE){
          cat("move_type:", move_type)
          cat("\n")
          cat("\n split_counts: \n")
          print(split_counts)
          cat("\n log prob_components = ", c(prob1, prob2, prob3), "\n")
          cat("\n exp accept_prob = ", exp(prob1 + prob2 + prob3))
          
          prob_sign = sign(c(prob1, prob2, prob3))
          inf_index = which(is.infinite(prob_components))
          
          # if NaN but no inf components detected - throw error
          if(length(inf_index) == 0){
            stop("NaN detected in SM acceptance probability. Unable to resolve.")
          }
          
          # if inf detected, resolve indeterminate forms by setting Inf equal to 
          # some arbitrarily large number
          for(prob_i in 1:length(inf_index)){
            prob_components[inf_index[prob_i]] = prob_sign[inf_index[prob_i]]*5000
          }
        }
        
        
        accept_prob = min(1, exp(sum(prob_components)))
        u = runif(n = 1)
        if(accept_prob > u){
          
          # accept
          accept = 1
          group_assign[s,] = merge_temp_group_assign[sm_iter+1,]
          
          # print(table(group_assign[s,]))
          
          # put old group lab out of reserve
          # cat("\n Curr labels (in accept if statement): ", curr_labels, "\n")
          
          curr_label_del = which(curr_labels == split_lab[2])
          avail_labels = c(curr_labels[curr_label_del], avail_labels)
          curr_labels = curr_labels[-curr_label_del]
          
          k = length(curr_labels)
          
          # update labels, etc
          count_assign = as.numeric(table(group_assign[s,]))
          label_assign = as.numeric(names(table(group_assign[s,])))
          num_groups[s,] = k
          
          # if new group created by merge, update mean and variance
          ## add new means and variances from final Gibbs scan to relevant vectors/lists
          length_sigma2 = length(sigma2)
          sigma2[[which_split_labs[1]]] = merge_vars[[sm_iter+1]][1] # scalar in DEV
          sigma2 = sigma2[-which_split_labs[2]]
          
          mu[,which_split_labs[1]] = merge_means[[sm_iter+1]]
          mu = mu[,-which_split_labs[2]]
          
        } else{
          # reject
          accept = 0
          # group assign remains unchanged
        }
        
        # cat("\n accept = ", accept, "\n")
        
      } # end merge proposal
      
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
      
      if((s %% print_iter == 0) & (s >= print_start) & (verbose == TRUE)){
        cat("\n End of split/merge step, iter: ", s, "\n") # just create a new line for separation
        cat("\n PROPOSED MOVE: ", move_type, " Accepted = ", accept, "\n")
        print(paste("Current k = ", k))
        cat("\n")
        print(group_assign[s,])
        cat("\n")
        print(table(group_assign[s,]))
        cat("\n Curr labels: \n")
        print(curr_labels)
        cat("\n")
        print(mu)
        cat("\n")
        print(sigma2)
        cat("\n")
        
        # print plot
        if(move_type == "SPLIT"){
          proposed_assign = split_temp_group_assign[sm_iter+1,]
        } else{
          # merge
          proposed_assign = merge_temp_group_assign[sm_iter+1,]
        }
        
        
        if(nrow(mu0) == 2){
          # if this is a 2D problem, can make scatterplot of group assign
          yvals = matrix(data = unlist(y), ncol = nrow(mu0), byrow = TRUE)
          plot_y = data.frame(
            y1 = yvals[,1],
            y2 = yvals[,2],
            curr_assign = proposed_assign
          )
          
          prog_plot = ggplot(data = plot_y, aes(x = y1, y = y2, label = rownames(plot_y))) +
            #geom_point(color = assign) +
            #geom_text(size = 3, hjust = 0, nudge_x = 0.5, color = assign) +
            geom_text(size = 3, color = plot_y$curr_assign) +
            ggtitle(paste0("Proposed ", move_type, " s=", s, ", k=", k, " Accept=", accept)) + 
            theme_classic()
          print(prog_plot)
          
        } else if(nrow(mu0) == 3){
          # if this is a 3D problem, can make scatterplot of group assign
          yvals = matrix(data = unlist(y), ncol = nrow(mu0), byrow = TRUE)
          plot_y = data.frame(
            y1 = yvals[,1],
            y2 = yvals[,2],
            y3 = yvals[,3],
            curr_assign = proposed_assign
          )
          
          split_obs_col = rep(1, nrow(plot_y))
          split_obs_col[sampled_obs] = 3 # color SM candidates green
          
          prog_plot = scatterplot3d(x = plot_y$y1, y = plot_y$y2, z = plot_y$y3, 
                                    color = plot_y$curr_assign, angle = -45, cex.symbols = 0.5, 
                                    xlab = "y1", ylab = "y2", zlab = "y3", pch = 20,
                                    main = paste0("Proposed ", move_type, " s=", s, ", k=", k, ", Accept=", accept))
          
          text(prog_plot$xyz.convert(plot_y[,1:3]), labels = rownames(plot_y), 
               pos = 4, cex = 0.75, col = split_obs_col)
          # print(prog_plot)
        }
        
        
      }
      
    }
    
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
                    # FUN = function(x){diag(1/(count_assign[x]/sigma2[[x]] + 1/sigma0), p)}) 
                    FUN = function(x){1/(count_assign[x]/sigma2[[x]] + 1/sigma0)}) 
    
    mu_mean = lapply(X = 1:k, 
                     FUN = function(x){(sum_y_i[,x]/sigma2[[x]] + mu0/sigma0)*mu_cov[[x]]})
    
    mu_list = lapply(X = 1:k, 
                     FUN = function(x){
                       t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                          mean = mu_mean[[x]], 
                                          sigma = diag(mu_cov[[x]],p)))
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
                               shape = (p*count_assign[x] + 2*a)/2,
                               rate = sum(loss_y_i[,x])/2 + b)
                    })
    
    
    # draw b parameter if indicated
    if(sigma_hyperprior == TRUE){
      
      # not implemented in this model at this time
      print("Sigma hyperprior not implemented in this model at this time.")
      
      # find sum of precision for each variance element across k groups
      # sum_prec_k = Reduce(f = "+", 
      #                    x = lapply(X = 1:k,
      #                               FUN = function(x){
      #                                 1/diag(Sigma[[x]])
      #                               })
      #                    )
      
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
    vars[[s]] = sigma2
    
    
    # print progress
    if((s %% print_iter == 0) & (s >= print_start) & (verbose == TRUE)){
      cat("After Gibbs step:") # just create a new line for separate
      cat("\n")
      cat("mu")
      cat("\n")
      print(mu)
      cat("\n")
      cat("sigma2")
      cat("\n")
      print(sigma2)
      cat("\n")
      cat("b",b)
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
        prog_plot = ggplot(data = plot_y, aes(x = y1, y = y2, label = curr_assign)) +
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
      group_prob_calc_DEV(k = k, n = n, n_j = count_assign, alpha = alpha, a = a, b = b, 
                          y_i = y[[x]], mu = mu, sigma2 = sigma2, mu0 = mu0, sigma0 = sigma0,
                          singleton = 0, curr_group_assign = group_assign[s,x], 
                          curr_labels = curr_labels)$pr_c
    }) # result is a k*n matrix
    
    # cat("\n Final prob calc dimension dim:", dim(pr_c), " length:", length(pr_c), "\n")
    
    #print(dim(pr_c))
    if(k == 1){
      
      probs[[s]] = matrix(pr_c, nrow=length(y))  # fixing issue with dimensions of prob matrix
      # when only one group is found
      
      # previously were fixing number of cols but this was an issue when the pr_new
      # was included -- - R coerced number of rows to be 2*n when ncol set to 1 for k=1
      
    } else{
      
      probs[[s]] = t(pr_c) # save probabilities for each obs, group
      # transpose so output is a n*k matrix - chop off pr_new -- not needed here
      
    }
    
  } ### end MCMC iterations from s=1:S
  
  
  ## calculate Laplacian and average adjacency matrices across all MCMC runs
  pairwise_mats = pairwise_prob_mat(group_assign = group_assign, probs = probs, diag_weights = diag_weights)
  
  end = Sys.time()
  
  if(sigma_hyperprior == TRUE){
    
    settings = list(S = S, alpha = alpha, a = a, b = b, 
                    mu0 = mu0, sigma0 = sigma0,
                    k_init = k_init, init_method = init_method,
                    mod_type = "nonconjDEV", 
                    split_merge = split_merge, sm_iter = sm_iter)
    
    return_list = list(settings = settings,
                       runtime = difftime(end, start, units = "m"),
                       data = y,
                       truth = truth,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       extra_params = extra_params,
                       accept = accept_ind,
                       sm_results = sm_results,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  } else{
    
    settings = list(S = S, seed = seed, alpha = alpha,
                    a = a, b = b, mu0 = mu0, sigma0 = sigma0, 
                    k_init = k_init, init_method = init_method,
                    #d = d, f = f,
                    mod_type = "nonconjDEV", 
                    split_merge = split_merge, sm_iter = sm_iter)
    
    return_list = list(settings = settings,
                       runtime = difftime(end, start, units = "m"),
                       data = y,
                       truth = truth,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       #extra_params = extra_params,
                       accept = accept_ind,
                       sm_results = sm_results,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  }
  
  return(return_list)
  
} ### end MCMC function



