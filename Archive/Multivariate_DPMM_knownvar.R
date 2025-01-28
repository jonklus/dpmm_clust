##############################################################
################# MULTIVARIATE DPMM MODEL ####################
#################   VAR KNOWN & EQUAL    #####################
##############################################################

## Author: Jonathan Klus
## Date: 15 March 2023
## Description: Sampler for a mixture of multivariate normal densities using Algorithm 2
## from Neal (2000). We assume the variance is known and equal so that we have conjugacy
## between the normal likelihood and dirichlet process prior. 

###################### load packages and set seed #############################
set.seed(516)
library(ggplot2)

################## simulate data ##############################################
# w = c(0.4, 0.3, 0.3)
# means = list(
#   c(-20, 20, 0),
#   c(20, -20, 10),
#   c(0, 0, -20)
# )
# 
# var = diag(5, length(means[[1]])) # variances known, diagonal, and equal
# 
# assign = sample(x = 1:length(means), size = 20, replace = TRUE, prob = w)
# y = lapply(X = assign,
#            FUN = function(x){
#              t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
#            })

# ggplot(data = data, aes(x = y)) +
#   geom_histogram(binwidth = 2) +
#   theme_classic()

# drop params used to generate data, except for var which is assumed known
#rm(w, means)

########################### Gibbs sampler #####################################

# helper functions

## calculate group membership probabilities

group_prob_calc <- function(k, n, n_j, alpha, y_i, mu, sigma2, r, mu0, singleton = 0, curr_group_assign = NULL, curr_labels = NULL){
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
  
  
  #### probability of creating a new group
  pr_new = alpha/(n-1+alpha)*mvtnorm::dmvnorm(x = y_i[,1], 
                                              mean = mu0,
                                              sigma = diag(sigma2+r*sigma2, p))

  #### normalize probs to account for "b"
  pr_c = c(pr_curr, pr_new)/(sum(pr_curr) + pr_new)
  
  return(pr_c)
  
}


MVN_CRP_sampler <- function(S = 10^3, y, sigma2, alpha = 1, r = 1, mu0, k_init = 2,
                            verbose = TRUE, print_iter = 100){
  # S is number of MCMC iterations
  # y is a list of data of length n
  # sigma2 is the known variance
  # alpha is the dirichlet process concentration parameter
  # r is the prior variance hyperparamter on the distribution of mu
  # mu0 is the prior mean - must be of dimension p*1
  # K_init is the initial number of groups
  # verbose & printmod tells the function if and how often to print a progress summary
  
  # data is y - a list of length n
  n = length(y)
  p = length(y[[1]]) # dimensionality of MVN
  k = k_init # initial number of groups
  
  # preallocate memory and set initial values
  group_assign = matrix(data = NA, nrow = S, ncol = n)
  
  group_assign[1, ] = sample(x = 1:k, size = length(y), replace = TRUE, prob = rep(1/k, k))
  # try different group assign initialization
  # group_assign[1, ] = ifelse(y > mean(y), k, k-1)  doesn't work for MVN, try kmeans?
  
  means = vector(mode = "list", length = S) #matrix(data = NA, nrow = S, ncol = n)
  probs = vector(mode = "list", length = S)  # group assignment probabiltiies
  # need to find p*1 vector of means based on this list of observed p*1 y_i values
  means[[1]] = sapply(X = 1:k, 
                      FUN = function(x){
                        rowMeans(matrix(unlist(y[group_assign[1,] == x]), nrow = p))
                        # unravel list of p*1 observations, put in matrix, find mean
                      })
  
  # initial allocation probs
  if(k == 1){
    probs[[1]] = matrix(data = 0, nrow = n, ncol = n) 
    probs[[1]][,1] = rep(1, n)
  } else{
    probs[[1]] = matrix(data = 0, nrow = n, ncol = n) 
    probs[[1]][,1:k] = rep(1/k, k)
  }

  num_groups = matrix(data = NA, nrow = S, ncol = 1)
  num_groups[1,1] = k_init
  mu = means[[1]]
  curr_labels = 1:k
  avail_labels = c(1:n)[-curr_labels]
  
  # iterate 1:S
  for(s in 2:S){
    
    ## initialize group assignments for current iteration using ending state from prior iteration
    group_assign[s, ] = group_assign[s-1, ]
    #probs[[s]] = matrix(data = 0, nrow = n, ncol = n) # rows will be individuals i=1:n
    # columns will be groups 1:k, there are max n groups
    
    ## consider split/merge step for conjugate case (Jain & Neal 2004) ?? 
    
    
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
        # check if there are singletons identified, otherwise do not touch mu
        mu = matrix(mu[,-singleton_index], nrow = p)
        count_assign = count_assign[-singleton_index]
        label_assign = label_assign[-singleton_index]
        avail_labels = c(group_assign[s,i], avail_labels)
        curr_labels = curr_labels[-which(curr_labels %in% group_assign[s,i])]
        k = length(curr_labels) # update k 
        
        ### for any observation i, calculate group membership probabilities
        pr_c = group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[[i]], mu = mu, sigma2 = sigma2, r = r,
                               mu0 = mu0, singleton = 1)
        

        ### draw a group assignment conditional on calculate group membership probs
        group_assign[s,i] = sample(x = c(curr_labels, avail_labels[1]), size = 1, prob = pr_c)
        
        #### if new group selected
        if(group_assign[s,i] == avail_labels[1]){
          
          #### handle bookkeeping
          curr_labels = c(curr_labels, avail_labels[1])
          avail_labels = avail_labels[-1]
          k = length(curr_labels)
          
          #### draw a mean for newly created group from FC posterior of mu using
          #### only the ith observation
          mu_mean = (y[[i]] + (1/r)*mu0)/(1/r + 1)
          mu_cov = diag(sigma2/(1/r + 1), p)
          mu_k = matrix(mvtnorm::rmvnorm(n = 1, mean = mu_mean, sigma = mu_cov), nrow = p) # kth mean
          mu = cbind(mu, mu_k)
        }
      } else{
        
        #### if obs i is not presently a singleton
        
        
        pr_c = group_prob_calc(k = k, n = n, n_j = count_assign, alpha = alpha, 
                               y_i = y[[i]], mu = mu, sigma2 = sigma2, r = r,
                               mu0 = mu0, singleton = 0, curr_group_assign = group_assign[s,i], 
                               curr_labels = curr_labels)
        

        group_assign[s,i] = sample(x = c(curr_labels, avail_labels[1]), size = 1, prob = pr_c)
        
        #### if new group selected
        if(group_assign[s,i] == avail_labels[1]){
          
          #### handle bookkeeping
          curr_labels = c(curr_labels, avail_labels[1])
          avail_labels = avail_labels[-1]
          k = length(curr_labels)
          
          #### draw a mean for newly created group from FC posterior of mu using
          #### only the ith observation
          mu_mean = (y[[i]] + (1/r)*mu0)/(1/r + 1)
          mu_cov = diag(sigma2/(1/r + 1), p)
          mu_k = matrix(mvtnorm::rmvnorm(n = 1, mean = mu_mean, sigma = mu_cov), nrow = p) # kth mean
          mu = cbind(mu, mu_k)
          
        }
        
      }
      
      ### iterate through all y_i
      
      # cat("Labels", curr_labels)
      # cat("Probs", pr_c)
      #### save allocation probabiltiies after each observation i 
      # probs[[s]][i,curr_labels] = pr_c
      
    } ### end iterations from i=1:n
    
    #print("B", mu)
    
    # final update of counts after a sweep
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    num_groups[s,] = k
    
    # print results after each sweep
    
    # draw group means for all k groups
    sum_y_i = sapply(X = 1:k, 
                     FUN = function(x){
                       rowSums(matrix(unlist(y[group_assign[s,] == label_assign[x]]), nrow = p))
                       # unravel list of p*1 observations, put in matrix, find sum
                     })
    
    mu_mean = lapply(X = 1:k, FUN = function(x){(sum_y_i[,x] + mu0/k)/(1/r + count_assign[x])})
    mu_cov = lapply(X = 1:k, FUN = function(x){diag(sigma2/(1/r + count_assign[x]), p)})   
    mu_k = lapply(X = 1:k, 
                  FUN = function(x){
                    mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                     mean = mu_mean[[x]], 
                                     sigma = mu_cov[[x]])
                  }) 
    
    mu = matrix(data = unlist(x = mu_k), nrow = p) # put draws of mu back into same
    # format we were using before
    
    # save draws of mu
    means[[s]] = mu
    
    # print progress
    if((s %% print_iter == 0) & (s > print_iter) & (verbose == TRUE)){
      print(" ") # just create a new line for separate
      print(paste("iter = ", s))
      print(paste("Current k = ", k))
      print(table(group_assign[s,]))
      print(mu)
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
        #print(plot_y)
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
                      y_i = y[[x]], mu = mu, sigma2 = sigma2, r = r,
                      mu0 = mu0, singleton = 0, curr_group_assign = group_assign[s,x], 
                      curr_labels = curr_labels)
    }) # result is a (k+1)*n matrix -- still including new group
    
    probs[[s]] = t(pr_c)[,1:k] # save probabilities for each obs, group
    # transpose so output is a n*k matrix - chop off pr_new -- not needed here
    
  } ### end MCMC iterations from s=1:S
  
  return(list(k = num_groups, 
              means = means,
              group_probs = probs,
              group_assign = group_assign))
  
} ### end MCMC function

## run function
# test1 = MVN_CRP_sampler(S = 10^3, y = y, sigma2 = var[1], alpha = 1, r = 1, 
#                         mu0 = matrix(data = rep(0, 3), nrow = 3), 
#                         k_init = 2)
# 
# ## check labels
# table(sapply(X = 1:nrow(test1$group_assign), FUN = function(x){length(as.numeric(table(test1$group_assign[x,])))}))

## write post-processing function to deal with label switching


