##############################################################
################# MULTIVARIATE DPMM MODEL ####################
#################  NON-CONJ SPLIT MERGE   ####################
##############################################################

## Author: Jonathan Klus
## Date: 30 January 2024
## Description: Split-merge algorithm for a mixture of multivariate normal densities 
# using algorithm from Jain & Neal (2007). 

###################### load packages and set seed ##############################
library(ggplot2)
library(mvtnorm)
library(LaplacesDemon)


############################### HELPER FUNCTIONS ###############################

update_phi_DEV <- function(curr_label, group_assign, count_assign, y, 
                           mu, mu0, Sigma, Sigma0, a, b){
  # function to sample from full conditional posteriors of DEV model
  
  # group_index is the current split or merge group for which new parameters are 
  # being drawn
  p = nrow(mu)
  sigma2 = diag(Sigma)[1]
  sigma0 = diag(Sigma0)[1]
  
  # draw group mean
  sum_y_i = rowSums(matrix(unlist(y[group_assign == curr_label]), nrow = p))
  
  mu_cov = 1/(count_assign/sigma2 + 1/sigma0)
  
  mu_mean = (sum_y_i/sigma2 + mu0/sigma0)*mu_cov
  
  mu_list = t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                              mean = mu_mean, 
                              sigma = diag(mu_cov,p)))
  
  
  
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
  
  dens = log(mvtnorm::dmvnorm(x = mu, mean = mu0, sigma = Sigma0)) + 
    log(LaplacesDemon::dinvgamma(x = sigma2, shape = a, scale = b))
  
  return(exp(dens))
}


nonconj_phi_prob_DEV <- function(curr_label, group_assign, count_assign, y, 
                             mu, mu0, Sigma, Sigma0, a, b){
  # This is P_GS(phi*|phi^L,...) from Jain & Neal 2007 
  
  
  # for the kth component under a DEV assumption
  sigma2 = diag(Sigma)[1]
  sigma0 = diag(Sigma0)[1]
  p = nrow(mu)
  loss_y_i = rowSums((matrix(unlist(y[group_assign == curr_label]), nrow = p) - mu)^2)
  loss_mu_k = t(mu0 - matrix(mu, nrow = p))%*%(mu0 - matrix(mu, nrow = p))
  # density of posterior up to a constant...
  dens = (sigma2^(-(p/2+a-1)))*exp(-0.5*(loss_y_i/sigma2 + 
                                               2*b/sigma2 + loss_mu_k/sigma0))
  
  # for the kth component under a UVV assumption
  
  return(dens)
}

ll_components_DEV <- function(obs_ind, y, mu, Sigma){
  # function to calculate the components of the likelihood ratio for the MH
  # acceptance probability after n_iter split merge restricted Gibbs scans

    val = mvtnorm::dmvnorm(x = y[[obs_ind]][,1], 
                           mean = mu[,1], 
                           sigma = Sigma)
    
   return(val)  
  
}

nonconj_component_prob_c <- function(obs, split_labs, group_assign, y, mu, Sigma){
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
  which_group_k = which(split_labs == group_assign[obs])
  sm_counts[which_group_k] = sm_counts[which_group_k] - 1
  
  if(0 %in% sm_counts){
    
    which_one = which(sm_counts == 1)
    
    # check length of which_zero
    if(length(which_one) == 1){
      # proceed as usual
      # cat("1 singleton", "\n")
      if(which_one == 1){
        
        num = (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]], 
                                            mean = mu[[2]], 
                                            sigma = Sigma[[2]]) 
        
        denom = num + (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]], 
                                                    mean = mu[[1]], 
                                                    sigma = Sigma[[1]])
        
        # will just be (1,0)... seems extreme
        ratio = c(1-(num/denom), num/denom)
        
      } else{ 
        # which_one == 2
        
        num = (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]], 
                                              mean = mu[[1]], 
                                              sigma = Sigma[[1]]) 
        
        denom = num + (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]], 
                                                        mean = mu[[2]], 
                                                        sigma = Sigma[[2]])
        
        # will just be (1,0)... seems extreme
        ratio = c(num/denom, 1-(num/denom))
        
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
    
    num = (sm_counts[1])*mvtnorm::dmvnorm(x = y[[obs]], 
                                        mean = mu[[1]], 
                                        sigma = Sigma[[1]]) 
    
    denom = num + (sm_counts[2])*mvtnorm::dmvnorm(x = y[[obs]], 
                                                mean = mu[[2]], 
                                                sigma = Sigma[[2]])
    
    ratio = c(num/denom, 1-(num/denom))
    
  }

    
  
  return(ratio)
  
}


################################ MAIN FUNCTIONS ################################

# inputs
p = 2
mu = matrix(rep(0, p), ncol=1)
mu0 = matrix(rep(0, p), ncol=1)
Sigma0 = diag(10, p)
Sigma = diag(10, p)
sm_iter = 5

# Split-Merge step --- every 10 iterations

# final update of counts after previous sweep
count_assign = as.numeric(table(group_assign[s,]))
label_assign = as.numeric(names(table(group_assign[s,])))
num_groups[s,] = k

if((split_merge == TRUE) & (s %% sm_iter == 0)){
  
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
  
  # bookkeeping - group labels
  subset_index = which(split_temp_group_assign[1,] %in% c(lab1, lab2)) 
  existing_group_index = which(label_assign == lab1) 
  original_param_index = which(curr_labels == lab1)
  original_mu = means[[original_param_index]]
  original_sigma = vars[[original_param_index]]
  
  # cat("split_labs:", split_lab)
  # cat("\n")
  # cat("subset_index:", subset_index)
  # cat("\n")
  

  # perform restricted Gibbs scans
  
  # if SPLIT
  if(move_type == "SPLIT"){
    
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
                                    diag(rgamma(n = 1, shape = a, rate = b), p)
                                    # for UVV
                                    # LaplacesDemon::rinvwishart(nu = nu, 
                                    #                            S = lambda0)
                                  })
        
        # draw params from prior - random launch state for merge proposal
        merge_means[[1]] = t(mvtnorm::rmvnorm(n = 1, mean = mu0, sigma = Sigma0))
        merge_vars[[1]] = diag(rgamma(n = 1, shape = a, rate = b), p)  
        # merge_vars[[1]] = LaplacesDemon::rinvwishart(nu = nu, S = lambda0)
        
      } else{
        
        # initialize current restricted Gibbs scan iteration with previous result
        split_temp_group_assign[scan,] = split_temp_group_assign[(scan-1),] 
        merge_temp_group_assign[scan,] = merge_temp_group_assign[(scan-1),] 
        split_means[[scan]] = split_means[[scan-1]]
        split_vars[[scan]] = split_vars[[scan-1]] 
        
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
            
            split_assign_prob = nonconj_component_prob(
                obs = obs, split_labs = split_lab,
                group_assign = split_temp_group_assign[scan,], 
                y = y, mu = split_means[[scan]], Sigma = split_vars[[scan]])
            
            ### merge_assign prob by definition always 1 when proposing merge
            ### launch state from c_i = c_j, no other way to permute assignment
            merge_assign_prob = 1
            
            # dont sample --- but do record prob of group anchor obs in in!
            which_split_lab_anchor = which(split_lab == split_temp_group_assign[scan,obs])
            # how to do this best for merge????
            which_merge_lab_anchor = 1  # why is this 1?? need to figure out what this should be 
            
            # temp_group_assign[scan,obs] = split_lab[sm_prop_index]
            # dont need to assign again - already initialized since anchor
            split_sm_probs[scan,obs] = split_assign_prob[which_split_lab_anchor]
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
            
            split_assign_prob = nonconj_component_prob(
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
      
      merge_phi = update_phi_DEV(curr_label = merge_lab, 
                                 group_assign = merge_temp_group_assign[scan,], 
                                 count_assign = merge_count_assign[merge_group_count_index][x], 
                                 y = y, 
                                 mu = merge_means[[scan]], 
                                 mu0 = mu0, 
                                 Sigma = merge_vars[[scan]], 
                                 Sigma0 = Sigma0, a = a, b = b)

      split_means[[scan]] = list(split_phi[[1]]$mu, split_phi[[2]]$mu)
      split_vars[[scan]] = list(split_phi[[1]]$Sigma, split_phi[[2]]$Sigma)
      
      merge_means[[scan]] = merge_phi$mu
      merge_vars[[scan]] = merge_phi$Sigma
      
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
    # split_counts = table(split_temp_group_assign[sm_iter+1,])
    # merge_counts = table(split_temp_group_assign[sm_iter+1,])
    # 
    # split_group_count_index = which(as.numeric(names(split_counts)) %in% split_lab)
    
    ## proposal probability
    
    # compute P_GS(phi) from launch state to final scan for both split and merge proposals
    split_phi_prob = nonconj_phi_prob_DEV(group_assign = split_temp_group_assign[sm_iter+1,], 
                                      count_assign = split_count_assign, y = y, 
                                      mu = split_means[[scan]], 
                                      mu0 = mu0, 
                                      Sigma = split_vars[[scan]], 
                                      Sigma0 = Sigma0, a = a, b = b)
    
    merge_phi_prob = nonconj_phi_prob_DEV(group_assign = merge_temp_group_assign[sm_iter+1,], 
                                      count_assign = merge_count_assign, y = y, 
                                      mu = merge_means[[scan]], 
                                      mu0 = mu0, 
                                      Sigma = merge_vars[[scan]], 
                                      Sigma0 = Sigma0, a = a, b = b)
    
    prob1_c_num = Reduce(f = "+", x = log(split_sm_prob[sm_iter+1,subset_index]))
    prob1_phi_num = Reduce(f = "+", x = log(split_phi_prob[sm_iter+1,subset_index]))
    
    prob1_c_denom = Reduce(f = "+", x = log(merge_sm_prob[sm_iter+1,subset_index]))
    prob1_phi_denom = Reduce(f = "+", x = log(merge_phi_prob[sm_iter+1,subset_index]))
    
    prob1 = (prob1_c_num + prob1_phi_num) - (prob1_c_denom + prob1_phi_denom)
    
    ## prior ratio
    prob2_num = factorial(sm_counts[[split_group_count_index[1]]] -1)*
      factorial(sm_counts[[split_group_count_index[2]]] -1)*
      nonconj_prior_dens_DEV(mu = split_means[[scan]][[1]], mu0 = mu0, 
                             Sigma = split_vars[[scan]][[1]], 
                             Sigma0 = Sigma0, a = a, b = b)*
      nonconj_prior_dens_DEV(mu = merge_means[[scan]][[2]], mu0 = mu0, 
                             Sigma = merge_vars[[scan]][[2]], 
                             Sigma0 = Sigma0, a = a, b = b)
      
    prob2_denom = factorial(sm_counts[[split_group_count_index[1]]] + 
                              sm_counts[[split_group_count_index[2]]]-1)*
      nonconj_prior_dens_DEV(mu = merge_means[[scan]], mu0 = mu0, 
                             Sigma = merge_vars[[scan]], 
                             Sigma0 = Sigma0, a = a, b = b)

    prob2 = log(alpha) + (log(prob2_num) - log(prob2_denom))
    
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
                              mu = original_mu, 
                              Sigma = original_sigma)
      prob3_denom = prob3_denom + log(val)
    }
    
    prob3 = prob3_num1 + prob3_num2 - prob3_denom
    
    ## evaluate acceptance prob
    accept_prob = min(1, exp(prob1 + prob2 + prob3))
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
      ## add new means and variances from finall Gibbs scan to relevant vectors/lists
      
      ### TEST THIS PIECE, NOT SURE IF CORRECT!
      length_Sigma = length(Sigma)
      Sigma[[which_split_labs[1]]] = merge_vars[[sm_iter+1]][[1]]
      Sigma[[length_Sigma+1]] = merge_vars[[sm_iter+1]][[2]]
      
      mu[,which_split_labs[1]] = merge_means[[sm_iter+1]][[1]]
      mu = cbind(mu, merge_vars[[sm_iter+1]][[2]])
      
    } else{
      # reject
      accept = 0
      # group assign remains unchanged
    }
    
    
    
    
    # if MERGE    
  } else if(move_type == "MERGE"){
    
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
                                    mvtnorm::rmvnorm(n = 1, 
                                                     mean = mu0,
                                                     sigma = Sigma0)
                                  })
        
        split_vars[[1]] = lapply(X = 1:2, 
                                 FUN = function(x){
                                   diag(rgamma(n = 1, shape = a, rate = b), p)
                                   # for UVV
                                   # LaplacesDemon::rinvwishart(nu = nu, 
                                   #                            S = lambda0)
                                 })
        
        # draw params from prior - random launch state for merge proposal
        merge_means[[1]] = mvtnorm::rmvnorm(n = 1, mean = mu0, sigma = Sigma0)
        merge_vars[[1]] = diag(rgamma(n = 1, shape = a, rate = b), p)
        # merge_vars[[1]] = LaplacesDemon::rinvwishart(nu = nu, S = lambda0)
        
      } else{
        
        # initialize current restricted Gibbs scan iteration with previous result
        split_temp_group_assign[scan,] = split_temp_group_assign[(scan-1),] 
        merge_temp_group_assign[scan,] = merge_temp_group_assign[(scan-1),] 
        
        split_means[[scan]] = split_means[[scan-1]]
        split_vars[[scan]] = split_vars[[scan-1]] 
        
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
            
            split_assign_prob = nonconj_component_prob(
              obs = obs, split_labs = split_lab,
              group_assign = split_temp_group_assign[scan,], 
              y = y, mu = split_means[[scan]], Sigma = split_vars[[scan]])
            
            ### merge_assign prob by definition always 1 when proposing merge
            ### launch state from c_i = c_j, no other way to permute assignment
            merge_assign_prob = 1
            
            # dont sample --- but do record prob of group anchor obs in in!
            which_split_lab_anchor = which(split_lab == split_temp_group_assign[scan,obs])
            # how to do this best for merge????
            which_merge_lab_anchor = 1  # why is this 1?? need to figure out what this should be 
            
            # temp_group_assign[scan,obs] = split_lab[sm_prop_index]
            # dont need to assign again - already initialized since anchor
            split_sm_probs[scan,obs] = split_assign_prob[which_split_lab_anchor]
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
            merge_counts = table(split_temp_group_assign[scan,])
            
            #current_obs_index = which(temp_group_assign[scan,] == obs)
            #split_group_lab_index1 = which(temp_group_assign[scan,] == split_lab[1])
            #split_group_lab_index2 = which(temp_group_assign[scan,] == split_lab[2])
            
            # current observation under consideration cannot be included here
            
            split_assign_prob = nonconj_component_prob(
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
      
      # update phi -- should this be done after each scan or after each observation?
      split_phi = lapply(X = split_lab, 
                         FUN = function(x){
                           update_phi_DEV(curr_label = x, 
                                          group_assign = split_temp_group_assign[scan,], 
                                          count_assign = split_count_assign, y = y, 
                                          mu = split_means[[scan]][[x]], 
                                          mu0 = mu0, 
                                          Sigma = split_vars[[scan]][[x]], 
                                          Sigma0 = Sigma0, a = a, b = b)
                         })
      
      merge_phi = update_phi_DEV(curr_label = merge_lab, 
                                 group_assign = merge_temp_group_assign[scan,], 
                                 count_assign = merge_count_assign, y = y, 
                                 mu = merge_means[[scan]], 
                                 mu0 = mu0, 
                                 Sigma = merge_vars[[scan]], 
                                 Sigma0 = Sigma0, a = a, b = b)
      
      split_means[[scan]] = list(split_phi[[1]]$mu, split_phi[[2]]$mu)
      split_vars[[scan]] = list(split_phi[[1]]$Sigma, split_phi[[2]]$Sigma)
      
      merge_means[[scan]] = merge_phi$mu
      merge_vars[[scan]] = merge_phi$Sigma
      
    } # scans 1:(sm_iter+1)
    
    
    # calculate & evaluate acceptance prob
    
    # update counts after scans
    split_count_assign = as.numeric(table(split_temp_group_assign[sm_iter+1,]))
    split_label_assign = as.numeric(names(table(split_temp_group_assign[sm_iter+1,])))
    # split_singletons = split_label_assign[which(split_count_assign == 1)]
    
    merge_count_assign = as.numeric(table(merge_temp_group_assign[sm_iter+1,]))
    merge_label_assign = as.numeric(names(table(merge_temp_group_assign[sm_iter+1,])))
    # merge_singletons = merge_label_assign[which(merge_count_assign == 1)]
    # should never be a merge singleton --- need at least 2 anchor observations
    
    split_counts = table(split_temp_group_assign[sm_iter+1,])
    merge_counts = table(split_temp_group_assign[sm_iter+1,])
    
    split_group_count_index = which(as.numeric(names(split_counts)) %in% split_lab)
    
    ## proposal probability
    
    # compute P_GS(phi) from launch state to final scan for both split and merge proposals
    split_phi_prob = nonconj_phi_prob_DEV(group_assign = split_temp_group_assign[sm_iter+1,], 
                                          count_assign = split_count_assign, y = y, 
                                          mu = split_means[[scan]], 
                                          mu0 = mu0, 
                                          Sigma = split_vars[[scan]], 
                                          Sigma0 = Sigma0, a = a, b = b)
    
    merge_phi_prob = nonconj_phi_prob_DEV(group_assign = merge_temp_group_assign[sm_iter+1,], 
                                          count_assign = merge_count_assign, y = y, 
                                          mu = merge_means[[scan]], 
                                          mu0 = mu0, 
                                          Sigma = merge_vars[[scan]], 
                                          Sigma0 = Sigma0, a = a, b = b)
    
    prob1_c_num = Reduce(f = "+", x = log(split_sm_prob[sm_iter+1,subset_index]))
    prob1_phi_num = Reduce(f = "+", x = log(split_phi_prob[sm_iter+1,subset_index]))
    
    prob1_c_denom = Reduce(f = "+", x = log(merge_sm_prob[sm_iter+1,subset_index]))
    prob1_phi_denom = Reduce(f = "+", x = log(merge_phi_prob[sm_iter+1,subset_index]))
    
    prob1 = (prob1_c_num + prob1_phi_num) - (prob1_c_denom + prob1_phi_denom)
    
    ## prior ratio
    prob2_num = factorial(sm_counts[[split_group_count_index[1]]] -1)*
      factorial(sm_counts[[split_group_count_index[2]]] -1)*
      nonconj_prior_dens_DEV(mu = split_means[[scan]][[1]], mu0 = mu0, 
                             Sigma = split_vars[[scan]][[1]], 
                             Sigma0 = Sigma0, a = a, b = b)*
      nonconj_prior_dens_DEV(mu = merge_means[[scan]][[2]], mu0 = mu0, 
                             Sigma = merge_vars[[scan]][[2]], 
                             Sigma0 = Sigma0, a = a, b = b)
    
    prob2_denom = factorial(sm_counts[[split_group_count_index[1]]] + 
                              sm_counts[[split_group_count_index[2]]]-1)*
      nonconj_prior_dens_DEV(mu = merge_means[[scan]], mu0 = mu0, 
                             Sigma = merge_vars[[scan]], 
                             Sigma0 = Sigma0, a = a, b = b)
    
    prob2 = log(alpha) + (log(prob2_num) - log(prob2_denom))
    
    ## likelihood ratio
    subset_index_grp1 = which(temp_group_assign[sm_iter+1,] %in% split_lab[1]) 
    subset_index_grp2 = which(temp_group_assign[sm_iter+1,] %in% split_lab[2]) 
    
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
                              mu = merge_means[[sm_iter+1]][[1]], 
                              Sigma = merge_vars[[sm_iter+1]][[1]])
      prob3_denom = prob3_denom + log(val)
    }
    
    prob3 = prob3_num1 + prob3_num2 - prob3_denom
    
    ## evaluate acceptance prob
    accept_prob = min(1, exp(prob1 + prob2 + prob3))
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
      ## add new means and variances from finall Gibbs scan to relevant vectors/lists
      length_Sigma = length(Sigma)
      Sigma[[which_split_labs[1]]] = merge_vars[[sm_iter+1]][[1]]
      Sigma[[length_Sigma+1]] = merge_vars[[sm_iter+1]][[2]]
      
      mu[,which_split_labs[1]] = merge_means[[sm_iter+1]][[1]]
      mu = cbind(mu, merge_vars[[sm_iter+1]][[2]])
      
    } else{
      # reject
      accept = 0
      # group assign remains unchanged
    }
    
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
  
}
