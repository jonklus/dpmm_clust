##############################################################
################# MULTIVARIATE DPMM MODEL ####################
###############  NO GAPS NON CONJ - ALG 8 ####################
################# INDEP IG PRIORS - DEV   ####################
##############################################################

## Author: Jonathan Klus
## Date: 6 June 2023
## Description: Sampler for a mixture of multivariate normal densities using Algorithm 2
## from Neal (2000). We assume the variance is unknown, but we impose conjugacy
## by assuming the variance of the mean and likelihood differ only be a multiplicative
## parameter r. 

############################# load packages ####################################
# set.seed(516)
library(ggplot2)
library(LaplacesDemon)
library(mvtnorm)

############################### HELPER FUNCTIONS ###############################

## calculate group membership probabilities

group_prob_calc_diag <- function(k, n, n_j, alpha, y_i, mu, Sigma, mu0, Sigma0, 
                                 nu, Lambda0, singleton = 0, curr_group_assign = NULL, 
                                  curr_labels = NULL){
  # k is the number of existing groups
  # n is total number of observations
  # n_j is a vector of length k with the total number of observations in each group
  # alpha is the scalar concentration parameter from the dirichlet process
  # y_i is the single data vector of dimension p*1 that we are considering
  # mu is a matrix with k columns, contains k group mean vectors of dimension p*1

  
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
    
    ##### draw new values of phi from G_0
    mu_new = t(mvtnorm::rmvnorm(n = 1, mean = mu0[,1], sigma = Sigma0)) # column vec
    Sigma_new = LaplacesDemon::rinvwishart(nu = nu, S = Lambda0)
    
    ##### calculate probs
    pr_new = ((alpha/1)/(n-1+alpha))*mvtnorm::dmvnorm(x = y_i[,1], 
                                                mean = mu_new[,1],
                                                sigma = Sigma_new) 
    #### normalize probs to account for "b"
    pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
    
    return(list(pr_c = pr_c, Sigma_new = Sigma_new))
    
  }
  
  
  
  

  



############################ INDEPENDENT IG PRIORS ############################# 

MVN_CRP_nonconj_UVV <- function(S = 10^3, seed = 516, y, alpha = 1, k_init = 2,
                                mu0, Sigma0, nu, Lambda0,
                                split_merge = FALSE, sm_iter = 5, truth = NA,
                                diag_weights = FALSE, verbose = TRUE, print_iter = 100){
  
  # S is number of MCMC iterations
  # y is a list of data of length n
  # mu0 is the prior mean - must be of dimension p*1
  # Sigma0 is the p*p prior variance
  # nu is the prior df for the inverse Wishart prior on Sigma
  # Lambda0 is the prior covariance matrix for teh inverse Wishart
  # alpha is the dirichlet process concentration parameter
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
  
  # if(sigma_hyperprior == TRUE){
  #   extra_params = matrix(data = NA, nrow = S, ncol = 2)
  #   colnames(extra_params) = c("b")
  #   extra_params[1,1] = b # note that here b has dimension 1, diagonal var
  # } 
  
  # need to find p*1 vector of means based on this list of observed p*1 y_i values
  means[[1]] = sapply(X = 1:k, 
                      FUN = function(x){
                        rowMeans(matrix(unlist(y[group_assign[1,] == x]), nrow = p))
                        # unravel list of p*1 observations, put in matrix, find mean
                      })
  
  vars[[1]] = lapply(X = 1:k,
                     FUN = function(x){
                       cov(t(matrix(unlist(y[group_assign[1,] == x]), nrow = p)))
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
    if((s %% print_iter == 0) & (verbose == TRUE)){
      
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
        pr_res = group_prob_calc_diag(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[[i]], mu = mu, Sigma = Sigma,
                               a = a, b = b, mu0 = mu0, Sigma0 = Sigma0,
                               nu = nu, Lambda0 = Lambda0, 
                               singleton = 1)
        pr_c = pr_res$pr_c
        
      } else{
        
        #### calculate proposal distribution for group assignment
        #### if obs i is not presently a singleton
        pr_res = group_prob_calc_diag(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[[i]], mu = mu, Sigma = Sigma, 
                               mu0 = mu0, Sigma0 = Sigma0, nu = nu, Lambda0 = Lambda0, 
                               singleton = 0, curr_group_assign = group_assign[s,i], 
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
        
        #### draw variance for newly created group from FC posterior of Sigma
        #### according to algo, but this is conditional on mean so do prior for now
        Sigma_k = pr_res$Sigma_new #1/rgamma(n = 1, shape = a, rate = b)
        length_Sigma = length(Sigma)
        Sigma[[length_Sigma+1]] = Sigma_k
        
        #### draw a mean for newly created group from FC posterior of mu
        mu_cov_k = solve(solve(Sigma_k) + solve(Sigma0))
        mu_mean_k = (t(y[[i]])%*%solve(Sigma_k) + t(mu0)%*%solve(Sigma0))%*%mu_cov_k
        mu_k = matrix(mvtnorm::rmvnorm(n = 1, 
                                       mean = mu_mean_k, 
                                       sigma = mu_cov_k), 
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
    
    # proceed to split-merge step if true
    
    # Split-Merge step --- every 10 iterations
    
    if(split_merge == TRUE & (s %% sm_iter) == 0){
      
      ##  NONCONJ SPLIT MERGE GOES HERE ONCE COMPLETE
 
    }
    
    # print results after each sweep
    
    # Proceed to Gibbs step
    
    # draw group means for all K groups
    # Sigma = diag(x = Sigma, nrow = p)
    
    sum_y_i = sapply(X = 1:k, 
                     FUN = function(x){
                       rowSums(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p))
                       # unravel list of p*1 observations, put in matrix, find sum
                     })
    
    mu_cov = lapply(X = 1:k, 
                    # FUN = function(x){diag(1/(count_assign[x]/Sigma[[x]] + 1/Sigma0), p)}) 
                    FUN = function(x){
                      solve(count_assign[x]*solve(Sigma[[x]]) + solve(Sigma0))
                      }) 
    
    mu_mean = lapply(X = 1:k, 
                     FUN = function(x){
                       (t(sum_y_i[,x])%*%solve(Sigma[[x]]) + 
                          t(mu0)%*%solve(Sigma0))%*%mu_cov[[x]]
                       })
    
    mu_list = lapply(X = 1:k, 
                     FUN = function(x){
                       t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                          mean = mu_mean[[x]], 
                                          sigma = mu_cov[[x]]))
                     }) 
    
    mu = matrix(data = unlist(x = mu_list), nrow = p) # put draws of mu back into same
    # p*K matrix format we were using before
    
    # draw group variances for all K groups
    loss_y_i = lapply(X = 1:k, 
                      FUN = function(x){
                        
                        col_ind = x  # from outer apply
                        group_ind = which(group_assign[s,] == label_assign[x])
                        
                        Reduce(f = "+", 
                               x = lapply(X = group_ind, FUN = function(x){
                                 (y[[x]] - mu[,col_ind])%*%t(y[[x]] - mu[,col_ind])}))
                        
                      })

    
    Sigma = lapply(X = 1:k, 
                    FUN = function(x){
                      LaplacesDemon::rinvwishart(
                        nu = nu + count_assign[x], 
                        S = loss_y_i[[x]] + Lambda0)
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
      #                     1/Sigma[[x]]
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
      group_prob_calc_diag(k = k, n = n, n_j = count_assign, alpha = alpha, 
                      y_i = y[[x]], mu = mu, Sigma = Sigma, mu0 = mu0, 
                      Sigma0 = Sigma0, nu = nu, Lambda0 = Lambda0,
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
    
    
    
    ### something is happening above --- output is sometimes 1d instead of n*k
    ### look at this
    
  } ### end MCMC iterations from s=1:S
  
  
  ## calculate Laplacian and average adjacency matrices across all MCMC runs
  pairwise_mats = pairwise_prob_mat(group_assign = group_assign, probs = probs, diag_weights = diag_weights)
  
  end = Sys.time()
  
  # if(sigma_hyperprior == TRUE){
  #   
  #   settings = list(S = S, alpha = alpha, a = a, b = b, mu0 = mu0, Sigma0 = Sigma0,
  #                   k_init = k_init, #d = d, f = f, g = g, h = h,
  #                   mod_type = "conjDEV", 
  #                   split_merge = split_merge, sm_iter = sm_iter)
  #   
  #   return_list = list(settings = settings,
  #                      runtime = difftime(end, start, units = "m"),
  #                      data = y,
  #                      truth = truth,
  #                      k = num_groups, 
  #                      means = means,
  #                      vars = vars,
  #                      emp_means = emp_means,
  #                      emp_vars = emp_vars,
  #                      extra_params = extra_params,
  #                      accept = accept_ind,
  #                      sm_results = sm_results,
  #                      group_probs = probs,
  #                      group_assign = group_assign,
  #                      pairwise_mats = pairwise_mats)
  #   
  # } else{
  #   
    settings = list(S = S, seed = seed, alpha = alpha,
                    mu0 = mu0, Sigma0 = Sigma0, 
                    nu = nu, Lambda0 = Lambda0,
                    k_init = k_init, 
                    mod_type = "nonconjUVV", 
                    split_merge = split_merge, sm_iter = sm_iter)
    
    return_list = list(settings = settings,
                       runtime = difftime(end, start, units = "m"),
                       data = y,
                       truth = truth,
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
    
 # }
  
  return(return_list)
  
} ### end MCMC function



