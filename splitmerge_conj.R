# SPLIT MERGE CODE - CONJUGATE CASE
# FROM JAIN & NEAL 2004 

post_pred_EVV <- function(obs, which_group, r, sm_counts, nu, y, ybar, loss_ybar, mu0, lambda0){
  
    # restricted gibbs sampling scans
    
    loss_mu0 = (ybar[[which_group]] - mu0)%*%t(ybar[[which_group]] - mu0)
    
    mu_n = ((1/r)*mu0[,1] + sm_counts[which_group]*ybar[[which_group]])/((1/r) + sm_counts[which_group])
    
    nu_n = nu + sm_counts[which_group] - nrow(mu0) + 1
    
    k_n = (r+sm_counts[which_group]+1)/((r+sm_counts[which_group])*(nu_n))
    
    lambda_n = k_n*(lambda0 + loss_ybar[[which_group]] + 
                      loss_mu0*((sm_counts[which_group]/r)/(1/r + sm_counts[which_group])))
    
    # print(mu_n)
    # print(lambda_n)
    # print(nu_n)
    
    val = sm_counts[which_group]*LaplacesDemon::dmvt(x = y[[obs]][,1], 
                                                     mu = mu_n[,1], 
                                                     S = lambda_n, 
                                                     df = nu_n)
    
  return(val)
  
}

final_post_pred_EVV <- function(y_i, r, nu, y, mu0, lambda0){
  
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
  
  mu_n = ((1/r)*mu0[,1] + sm_counts*ybar)/((1/r) + sm_counts)
  
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

ll_components <- function(subset_index, obs_ind, y, mu0, lambda0, r, nu){
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
    val = final_post_pred_EVV(y_i = y[[subset_index[obs_ind]]], r = r, nu = nu, 
                              y = y[subset_index[1:(obs_ind-1)]], 
                              mu0 = mu0, lambda0 = lambda0)
    
  }
    
  return(val)  
  
}

# test line
# post_pred_EVV(obs=2,which_group=1, r=10, sm_counts=c(10,11), nu=2, y=y, ybar=ybar, 
#               loss_ybar=loss_ybar, mu0=matrix(data=0,nrow=2), lambda0=diag(10,2))
######################

split_merge_prob_EVV <- function(obs, split_labs, group_assign, r, nu, y, mu0, lambda0){
  # split_labs is an array of length 2 indicating which entries in counts correspond
  # to the groups that are part of the split/merge
  # which_group is a scalar valued 1 or 2 indicating which of the groups we are considering
  # obs is the index for the observation being considered at this point
  # group_assign is an array of length n corresponding to group assignments for each obs
  # r and nu are scalar hyperparameters
  # counts is an array of the number of obs assigned to each group
  # y is the data
  # mu0 and lambda0 are the prior mean and covariance matrix
  
  sm_counts = sapply(X = split_labs, FUN = function(x){sum(group_assign[-obs] == x)})
  
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
  
  
  
  
  ratio = sapply(X = 1:2,
                 FUN = function(x){
                   
                   xx = ifelse(x==1,2,1)
                   num = post_pred_EVV(obs = obs, which_group = x, r = r, # which group is which????
                                       sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                                       loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
                   denom = num + post_pred_EVV(obs = obs, which_group = xx, r = r,
                                               sm_counts = sm_counts, nu = nu, y = y, ybar = ybar, 
                                               loss_ybar = loss_ybar, mu0 = mu0, lambda0 = lambda0)
                   return(num/denom)
                   
                 })
  
  return(ratio)
  
}

# needs to go at top of CRP script, do not want overwritten each time SM step run
sm_results = matrix(data = NA, nrow = 1, ncol = 5)
colnames(sm_results) = c("s", "sm_iter", "move_type","accept", "prob")

temp_group_assign = matrix(data = NA, nrow = sm_iter + 1, ncol = length(group_assign[s,]))
sm_probs = matrix(data = NA, nrow = sm_iter + 1, ncol = length(y))
temp_group_assign[1,] = group_assign[s,]

# randomly select two observed data points y
sampled_obs = sample(x = 1:length(y), size = 2, replace = FALSE)

# check if in same group - if yes SPLIT
# if no, MERGE
lab1 = temp_group_assign[1, sampled_obs[1]]
lab2 = temp_group_assign[1, sampled_obs[2]]
move_type = ifelse(lab1 == lab2, "SPLIT", "MERGE")

# bookkeeping - group labels
subset_index = which(temp_group_assign[1,] %in% c(lab1, lab2)) 
anchor_obs_index = which(subset_index %in% sampled_obs)
subset_index_minus = subset_index[-anchor_obs_index] # except sampled observations i and j
existing_group_index = which(label_assign == lab1) 

# perform restricted Gibbs scans

# if SPLIT
if(move_type == "SPLIT"){
  
  # specify random launch state
  split_lab = c(lab1, avail_labels[1]) # keep original label, new one for 2nd group
  
  
  for(scan in 1:(sm_iter+1)){
    
    for(obs in subset_index_minus){
      
      if(scan == 1){
        
        # specify random launch state
        temp_group_assign[scan,obs] = sample(x = split_lab, size = 1)
        # specify new group label to 2nd anchor point as well
        temp_group_assign[scan,anchor_obs_index[2]] = split_lab[2] 
        
      } else{ # for remaining scans after random launch state set
        
        temp_group_assign[scan,] = temp_group_assign[(scan-1),] # initialize
        sm_counts = table(temp_group_assign[scan,-obs])
        split_group_count_index = which(as.numeric(names(sm_counts)) %in% split_lab)
        #current_obs_index = which(temp_group_assign[scan,] == obs)
        #split_group_lab_index1 = which(temp_group_assign[scan,] == split_lab[1])
        #split_group_lab_index2 = which(temp_group_assign[scan,] == split_lab[2])
        
        # current observation under consideration cannot be included here
          
        split_assign_prob = split_merge_prob_EVV(
          obs = obs, split_labs = split_lab, r=r, 
          group_assign = temp_group_assign[scan,], nu = nu, 
          y = y, mu0 = mu0, lambda0= lambda0)

        sm_prop_index = sample(x = 1:2, size = 1, 
                               prob = split_assign_prob)
        
        temp_group_assign[scan,obs] = split_lab[sm_prop_index]
        sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
  
      }  
    } # iterate through all observations in the two split groups under consideration
  } # scans 1:(sm_iter+1)
  
  # calculate & evaluate acceptance prob
  sm_counts = table(temp_group_assign[sm_iter+1,]) # update counts after scans
  ## proposal probability
  prob1 = -Reduce(f = "+", x = log(sm_probs[sm_iter+1,])) # log1 - sum(logs)
  ## prior ratio
  prob2_num = factorial(sm_counts[[split_group_count_index[1]]] -1)*factorial(sm_counts[[split_group_count_index[2]]] -1)
  prob2_denom = factorial(sm_counts[[split_group_count_index[1]]]+sm_counts[[split_group_count_index[2]]]-1)
  prob2 = log(alpha) + (log(prob2_num) - log(prob2_denom))
  
  ## likelihood ratio
  subset_index_grp1 = which(temp_group_assign[sm_iter+1,] %in% split_lab[1]) 
  subset_index_grp2 = which(temp_group_assign[sm_iter+1,] %in% split_lab[2]) 
  
  ### component 1 - numerator I (group 1)
  prob3_num1 = 0
  for(obs_ind in 1:length(subset_index_grp1)){
    val = ll_components(subset_index = subset_index_grp1, obs_ind = obs_ind, y = y, 
                        mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
    prob3_num1 = prob3_num1 + log(val)
  }

  ### component 2 - numerator II (group 2)
  prob3_num2 = 0
  for(obs_ind in 1:length(subset_index_grp2)){
    val = ll_components(subset_index = subset_index_grp2, obs_ind = obs_ind, y = y, 
                        mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
    prob3_num2 = prob3_num2 + log(val)
  }
  
  
  ### component 3 - denominator (all in original group)
  prob3_denom = 0
  for(obs_ind in 1:length(subset_index)){
    val = ll_components(subset_index = subset_index, obs_ind = obs_ind, y = y, 
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
    # take new group lab out of reserve
    curr_labels = c(curr_labels, avail_labels[1])
    avail_labels = avail_labels[-1]
    k = length(curr_labels)
    # update labels, etc
    count_assign = as.numeric(table(group_assign[s,]))
    label_assign = as.numeric(names(table(group_assign[s,])))
    which_split_labs = which(label_assign == split_lab) 
    
  } else{
    # reject
    accept = 0
    # group assign remains unchanged
  }
  
  
  # if new group created by split, give it a mean and variance
  
  ## draw variances for both groups - use empirical variance since mean not known yet
  ## may want to change this in the future if you come up with a better idea...
  
  ybar = lapply(X = split_labs, 
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
  
  emp_var = lapply(X = 1:2, 
                   FUN = function(x){
                     
                     col_ind = x  # from outer apply
                     group_ind = which(group_assign[s,] == split_lab[x])
                     if(obs %in% group_ind){
                       obs_ind = which(obs == group_ind)
                       #group_ind = group_ind[-obs_ind] # include all obs this time
                     } # else continue
                     
                     emp_loss = Reduce(f = "+", 
                                       x = lapply(X = group_ind, FUN = function(x){
                                         (y[[x]] - ybar[[col_ind]])%*%t(y[[x]] - ybar[[col_ind]])}))
                     
                     return(emp_loss/length(group_ind))
                     
                   })
  
  ## draw means for both groups conditional on empirical variances...
  
  sum_y_i = sapply(X = split_lab, 
                   FUN = function(x){
                     rowSums(matrix(unlist(y[group_assign[s,] == x]), nrow = p))
                     # unravel list of p*1 observations, put in matrix, find sum
                   })
  
  mu_cov = lapply(X = 1:2, 
                  FUN = function(x){emp_var[[x]]/(1/r + count_assign[x])}) 
  
  mu_mean = lapply(X = 1:2, 
                   FUN = function(x){(sum_y_i[,x] + mu0/r)/(1/r + count_assign[x])})
  
  mu_list = lapply(X = 1:2, 
                   FUN = function(x){
                     t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                        mean = mu_mean[[x]], 
                                        sigma = mu_cov[[x]]))
                   }) 

  ## add new means and variances to relevant vectors/lists
  length_Sigma = length(Sigma)
  Sigma[[which_split_labs[1]]] = emp_var[[1]]
  Sigma[[length_Sigma+1]] = emp_var[[2]]
  
  mu[,which_split_labs[1]] = mu_list[[1]]
  mu = cbind(mu, mu_list[[2]])

  # if MERGE    
} else if(move_type == "MERGE"){
  
  # specify random launch state
  split_lab = c(lab1, lab2) # original labels, assign merged obs to lab1
  temp_group_assign[scan,] = temp_group_assign[1,] # initialize
  
  count_assign = as.numeric(table(temp_group_assign[1,]))
  label_assign = as.numeric(names(table(temp_group_assign[1,])))
  which_split_labs = which(label_assign %in% split_lab) 
  
  # specify merge group
  scan = sm_iter + 1 # skip to last row of table
  temp_group_assign[scan,sampled_obs] = lab1 # anchor observations that arent
  # in subset-index-minus...
  
  sm_counts_before = table(temp_group_assign[scan,-obs])
  split_group_count_index = which(as.numeric(names(sm_counts_before)) %in% split_lab)
  
  for(obs in subset_index_minus){
      
      temp_group_assign[scan,obs] = lab1 # keep 1st group label
      
      sm_counts_before = table(temp_group_assign[scan,-obs])
      split_group_count_index_before = which(as.numeric(names(sm_counts)) %in% split_lab)
      
      # current observation under consideration cannot be included here
      
      split_assign_prob = split_merge_prob_EVV(
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
    prob1 = Reduce(f = "+", x = log(sm_probs[sm_iter+1,subset_index_minus])) # log1 - sum(logs)
    # need to index by subset since NAs for observation not part of split/merge -- as well as the
    # two anchor observations. shoudl you be calculating a prob for those as well though???
    
    ## prior ratio
    prob2_num = factorial(sm_counts_before[[split_group_count_index[1]]] -1)*factorial(sm_counts_before[[split_group_count_index[2]]] -1)
    prob2_denom = factorial(sm_counts_before[[split_group_count_index[1]]]+sm_counts_before[[split_group_count_index[2]]]-1)
    prob2 = -(log(alpha) + (log(prob2_num) - log(prob2_denom)))
    
    ## likelihood ratio
    subset_index_grp1 = which(temp_group_assign[1,] %in% split_lab[1]) 
    subset_index_grp2 = which(temp_group_assign[1,] %in% split_lab[2]) 
    
    ### component 1 - numerator I (group 1)
    prob3_num1 = 0
    for(obs_ind in 1:length(subset_index_grp1)){
      val = ll_components(subset_index = subset_index_grp1, obs_ind = obs_ind, y = y, 
                          mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
      prob3_num1 = prob3_num1 + log(val)
    }
    
    ### component 2 - numerator II (group 2)
    prob3_num2 = 0
    for(obs_ind in 1:length(subset_index_grp2)){
      val = ll_components(subset_index = subset_index_grp2, obs_ind = obs_ind, y = y, 
                          mu0 = mu0, lambda0 = lambda0, r = r, nu = nu)
      prob3_num2 = prob3_num2 + log(val)
    }
    
    
    ### component 3 - denominator (all in original group)
    prob3_denom = 0
    for(obs_ind in 1:length(subset_index)){
      val = ll_components(subset_index = subset_index, obs_ind = obs_ind, y = y, 
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
            avail_labels = c(temp_group_assign[1, sampled_obs[2]], avail_labels)
      
      # put old group lab out of reserve
      curr_label_del = which(curr_labels == split_lab[2])
      avail_labels = c(curr_labels[curr_label_del], avail_labels)
      curr_labels = curr_labels[-curr_label_del]

      k = length(curr_labels)
      # update labels, etc
      count_assign = as.numeric(table(group_assign[s,]))
      label_assign = as.numeric(names(table(group_assign[s,])))
      which_split_lab = which(label_assign == split_lab[1]) 
      
      # save info about accept/reject & probs
      
      # merged group...give it a mean and variance
      
      ## draw variances for new groups - use empirical variance since mean not known yet
      ## may want to change this in the future if you come up with a better idea...
      
      ybar = lapply(X = split_lab[1], 
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
      
      emp_var = lapply(X = 1, 
                       FUN = function(x){
                         
                         col_ind = x  # from outer apply
                         group_ind = which(group_assign[s,] == split_lab[x])
                         if(obs %in% group_ind){
                           obs_ind = which(obs == group_ind)
                           #group_ind = group_ind[-obs_ind] # include all obs this time
                         } # else continue
                         
                         emp_loss = Reduce(f = "+", 
                                           x = lapply(X = group_ind, FUN = function(x){
                                             (y[[x]] - ybar[[col_ind]])%*%t(y[[x]] - ybar[[col_ind]])}))
                         
                         return(emp_loss/length(group_ind))
                         
                       })
      
      ## draw means for both groups conditional on empirical variances...
      
      sum_y_i = sapply(X = split_lab[1], 
                       FUN = function(x){
                         rowSums(matrix(unlist(y[group_assign[s,] == x]), nrow = p))
                         # unravel list of p*1 observations, put in matrix, find sum
                       })
      
      mu_cov = lapply(X = which_split_lab, 
                      FUN = function(x){emp_var[[x]]/(1/r + count_assign[x])}) 
      
      mu_mean = lapply(X = which_split_lab, 
                       FUN = function(x){(sum_y_i[,x] + mu0/r)/(1/r + count_assign[x])})
      
      mu_list = lapply(X = which_split_lab, 
                       FUN = function(x){
                         t(mvtnorm::rmvnorm(n = 1, # make this the kth mean
                                            mean = mu_mean[[x]], 
                                            sigma = mu_cov[[x]]))
                       }) 
      
      ## add new means and variances to relevant vectors/lists
      length_Sigma = length(Sigma)
      Sigma[[which_split_lab]] = emp_var[[1]]
      
      mu[,which_split_lab] = mu_list[[1]]

      
    } else{
      # reject
      accept = 0
      # group assign remains unchanged
    }
    
    
    
    

}
  
sm_results = rbind(sm_results, c(s, sm_iter, move_type,accept, prob))


# calculate proposal probability

# calculate & evaluate acceptance prob

# if new group created by merge, give it a mean and variance

# final bookkeeping 


}