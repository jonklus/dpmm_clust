##############################################################
################# MULTIVARIATE DPMM MODEL ####################
#################     VAR UNKNOWN         ####################
##############################################################

## Author: Jonathan Klus
## Date: 6 June 2023
## Description: Sampler for a mixture of multivariate normal densities using Algorithm 4
## from Neal (2000). We assume the variance is unknown, so we no longer have conjugacy
## between the normal likelihood and dirichlet process prior. 

###################### load packages and set seed #############################
set.seed(516)
library(ggplot2)

########################### Gibbs sampler #####################################

# helper functions

## calculate group membership probabilities

prop_calc <- function(k, n, n_j, alpha, p, singleton = 0, curr_group_assign = NULL, curr_labels = NULL){
  # k is the number of existing groups
  # n is total number of observations
  # n_j is a vector of length k with the total number of observations in each group
  # alpha is the scalar concentration parameter from the dirichlet process
  # mu is a matrix with k columns, contains k group mean vectors of dimension p*1
  
  # calculate conditional priors of group membership for obs i
  
  #### probability of joining a current group
  

  
  if(singleton == 1){
    
    # if singleton, probability of remaining is zero
    
    pr_curr = sapply(X = 1:k, 
                     FUN = function(x){
                       # print("Step1")
                       # print(matrix(mu[,x], nrow = p))
                       # print(y_i)
                       c1 = n_j[x]/(n-1+alpha)
                       return(c1)
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
                         # c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                         #                       mean = mu[,x], 
                         #                       sigma = diag(sigma2, p)) 
                       } else{
                         # print("Step3")
                         # print(matrix(mu[,x], nrow = p))
                         # print(y_i)
                         c1 = n_j[x]/(n-1+alpha)
                         # c2 = mvtnorm::dmvnorm(x = y_i[,1], 
                         #                       mean = mu[,x],
                         #                       sigma = diag(sigma2, p)) 
                       }
                       return(c1)
                     })
    
  }
  
  
  #### probability of creating a new group
  pr_new = alpha/(n-1+alpha)
  
  #### normalize probs to account for "b"
  pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
  
  return(pr_c)
  
}

## calculate group membership probabilities after CRP step to use later to fix
## label switching problem

group_prob_calc <- function(k, n, n_j, alpha, y_i, mu, Sigma, singleton = 0, curr_group_assign = NULL, curr_labels = NULL){
  
  ## THIS FUNCTION IS NOT TO BE USED TO CALCULATE GROUP INCLUSION PROBS FOR THE 
  ## CRP, ESPECIALLY IN THE NON-CONJUGATE SETTING. ONLY USED TO CALCULATE PROBS
  ## **AFTER** GROUP ASSIGNMENT HAS ALREADY OCCURRED 
  
  # k is the number of existing groups
  # n is total number of observations
  # n_j is a vector of length k with the total number of observations in each group
  # alpha is the scalar concentration parameter from the dirichlet process
  # y_i is the single data vector of dimension p*1 that we are considering
  # mu is a matrix with k columns, contains k group mean vectors of dimension p*1
  # mu0 is a p*1 prior mean vector, tau2 is a scalar prior
  # sigma2 is a scalar, corresponds to diagonal covariance matrix
  # r is a scalar, corresponds to prior variance
  
  # calculate unnormalized probabilities of group membership for obs i
  
  #### probability of joining a current group
  
  p = length(mu[,1])
  
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

  #### normalize probs to account for "b"
  pr_c = pr_curr/sum(pr_curr)
  
  return(pr_c)
  
}


MVN_CRP_sampler <- function(S = 10^3, seed = 516, y, Sigma0, alpha = 1, a = 1/2, b = 10, mu0, k_init = 2,
                            d = 1, f = 1, sigma_hyperprior = TRUE, diag_weights = FALSE, verbose = TRUE, print_iter = 100){
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
  
  if(sigma_hyperprior == TRUE){
    extra_params = matrix(data = NA, nrow = S, ncol = p)
    colnames(extra_params) = paste0("b_", 1:p)
    b = rep(b, p)
    extra_params[1,1:p] = b # note that here b has dimension p, one for each dim
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
        pr_c = prop_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                         p = p, singleton = 1)
        
        
      } else{
        
        #### calculate proposal distribution for group assignment
        #### if obs i is not presently a singleton
        pr_c = prop_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                         p = p, singleton = 0, curr_group_assign = group_assign[s,i], 
                         curr_labels = curr_labels)
        
      }
      
      ### draw a group assignment conditional on group membership probs
      prop_group_assign = sample(x = c(curr_labels, avail_labels[1]), size = 1, prob = pr_c)
      
      #### if new group proposed - draw parameters for that group from base distn
      if(prop_group_assign == avail_labels[1]){
        
        #### draw means and covariance matrix for newly created group from G_0
        ##### draw variance components
        Sigma_k = diag(1/rgamma(n = p, shape = a, rate = b))
        
        # print("new var")
        # print(Sigma_k)
        
        ##### draw mean components conditional on variance components 
        ## should this also account for the value of y_i somehow???
        mu_k = matrix(mvtnorm::rmvnorm(n = 1, mean = c(mu0), sigma = Sigma0), nrow = p) # kth mean
        # print("new mu")
        # print(mu_k)
        #mu = cbind(mu, mu_k)
        
      } else{
        
        # if existing group chosen, identify its parameters
        prop_group_assign_index = which(label_assign == prop_group_assign)
        mu_k = matrix(mu[, prop_group_assign_index], nrow = p)
        Sigma_k = Sigma[[prop_group_assign_index]]
        
        # print("existing mu")
        # print(mu_k)
        
      }
      
      # print(dim(y[[i]]))
      # print(dim(Sigma_k))
      # print(dim(Sigma_curr))
      # print(c("Proposed group", prop_group_assign, prop_group_assign_index))
      # print(mu)
      # print(Sigma)
      
      if(group_assign[s,i] %in% singletons){
        
        # if singleton, we ALWAYS accept (either new singleton group or join existing)
        accept = 1
        
      } else{
        
        
        # compute acceptance probability
        prop_dens = mvtnorm::dmvnorm(x = c(y[[i]]), mean = c(mu_k), sigma = Sigma_k)
        curr_dens = mvtnorm::dmvnorm(x = c(y[[i]]), 
                                     mean = c(mu_curr), 
                                     sigma = Sigma_curr)
        accept_prob = min(1, prop_dens/curr_dens)
        
        #print(c(prop_dens, curr_dens))
        #print(c(mu_k))
        #print(c(mu_curr))
        
        u = runif(n = 1) # draw uniform random number on 0,1
        if(accept_prob >= u){
          accept = 1
        } else{
          accept = 0
        }
        
        # cannot do this way -- not numerically stable if accept_prob too large or small
        # accept = sample(x = c(1,0), prob = c(accept_prob, 1-accept_prob))
        
      }
      
      #print(accept)
      
      # if accepted --- update 
      if(accept == 1){
        
        accept_ind[s,i] = 1
        
        
        # if new group has been created 
        if(prop_group_assign == avail_labels[1]){
          
          # print(prop_group_assign)
          # print(avail_labels[1])
          
          #### update parameter vectors/lists
          mu = cbind(mu, mu_k)
          length_Sigma = length(Sigma)
          Sigma[[length_Sigma+1]] = Sigma_k
          
          #### handle bookkeeping
          curr_labels = c(curr_labels, avail_labels[1])
          avail_labels = avail_labels[-1]
          k = length(curr_labels)
          
          group_assign[s,i] = prop_group_assign
          
          # print("after update")
          # print(mu)
          # print(Sigma)
          
          
          
        } else{
          # if new group has not been created
          
          group_assign[s,i] = prop_group_assign
          # all else remains the same
          
        }
        
        
        
      } else{
        # if not accepted, remain the same
        accept_ind[s,i] = 0
        # group_assign[s,i] remains unchanged
        # parameter lists/vectors remain unchanged
        
      }
      
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
      print(Sigma)
    }
    
    # final update of counts after a sweep
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    num_groups[s,] = k
    
    # print results after each sweep
    
    # draw group means for all K groups
    sum_y_i = sapply(X = 1:k, 
                     FUN = function(x){
                       rowSums(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p))
                       # unravel list of p*1 observations, put in matrix, find sum
                     })
    
    mu_cov = lapply(X = 1:k, 
                    FUN = function(x){solve(count_assign[x]*solve(Sigma[[x]]) + solve(Sigma0))}) 
    
    mu_mean = lapply(X = 1:k, 
                     FUN = function(x){(t(sum_y_i[,x])%*%solve(Sigma[[x]]) + t(mu0)%*%solve(Sigma0))%*%mu_cov[[x]]})
    
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
    
    Sigma = lapply(X = 1:k, 
                   FUN = function(x){
                     diag(1/rgamma(n = p, 
                                   shape = rep(count_assign[x]/2 + a, p),
                                   rate = loss_y_i[,x]/2 + b))
                   })
    
    
    # draw b parameter if indicated
    if(sigma_hyperprior == TRUE){
      
      # find sum of precision for each variance element across k groups
      sum_prec_k = Reduce(f = "+", 
                         x = lapply(X = 1:k,
                                    FUN = function(x){
                                      1/diag(Sigma[[x]])
                                    })
                         )
        
      b = rgamma(n = p, shape = k*a + d, rate = sum_prec_k + f)
      
      extra_params[s,1:p] = b
      
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
    if((s %% print_iter == 0) & (s > print_iter) & (verbose == TRUE)){
      print("After Gibbs step:") # just create a new line for separate
      # print(paste("iter = ", s))
      # print(paste("Current k = ", k))
      # print(table(group_assign[s,]))
      print(mu)
      print(Sigma)
      print(b)
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
      group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                      y_i = y[[x]], mu = mu, Sigma = Sigma,
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
  
  if(sigma_hyperprior == TRUE){
    
    settings = list(S = S, y = y, Sigma0 = Sigma0, alpha = alpha, a = a, b = b, mu0 = mu0, k_init = k_init,
                    d = d, f = f)
    
    return_list = list(settings = settings,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       emp_means = emp_means,
                       emp_vars = emp_vars,
                       b = extra_params,
                       accept = accept_ind,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  } else{
    
    settings = list(S = S, seed = seed, y = y, Sigma0 = Sigma0, alpha = alpha, 
                    a = a, b = b, mu0 = mu0, k_init = k_init, d = d, f = f)
    
    return_list = list(settings = settings,
                       k = num_groups, 
                       means = means,
                       vars = vars,
                       emp_means = emp_means,
                       emp_vars = emp_vars,
                       b = extra_params,
                       accept = accept_ind,
                       group_probs = probs,
                       group_assign = group_assign,
                       pairwise_mats = pairwise_mats)
    
  }
  
  return(return_list)
  
} ### end MCMC function

## run function
# test1 = MVN_CRP_sampler(S = 10^3, y = y, sigma2 = var[1], alpha = 1, r = 1, 
#                         mu0 = matrix(data = rep(0, 3), nrow = 3), 
#                         k_init = 2)
# 
# ## check labels
# table(sapply(X = 1:nrow(test1$group_assign), FUN = function(x){length(as.numeric(table(test1$group_assign[x,])))}))

## write post-processing function to deal with label switching


