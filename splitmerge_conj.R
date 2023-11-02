# SPLIT MERGE CODE - CONJUGATE CASE
# FROM JAIN & NEAL 2004 

post_pred_EVV <- function(obs, which_group, r, sm_counts, nu, y, ybar, loss_ybar, mu0, lambda0, type){
  
  if(type == "scans"){
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
    
    
  } else if(type = "ll_ratio"){
    # components of likelihood ratio for final acceptance prob calc
    
    loss_mu0 = (ybar - mu0)%*%t(ybar - mu0)
    
    mu_n = ((1/r)*mu0[,1] + sm_counts*ybar)/((1/r) + sm_counts)
    
    nu_n = nu + sm_counts - nrow(mu0) + 1
    
    k_n = (r+sm_counts+1)/((r+sm_counts)*(nu_n))
    
    lambda_n = k_n*(lambda0 + loss_ybar + 
                      loss_mu0*((sm_counts/r)/(1/r + sm_counts)))
    
    # print(mu_n)
    # print(lambda_n)
    # print(nu_n)
    
    val = sm_counts[which_group]*LaplacesDemon::dmvt(x = y[[obs]][,1], 
                                                     mu = mu_n[,1], 
                                                     S = lambda_n, 
                                                     df = nu_n)
    
    
  }
  
  
  
  
  return(val)
  
}

# test line
# post_pred_EVV(obs=2,which_group=1, r=10, sm_counts=c(10,11), nu=2, y=y, ybar=ybar, 
#               loss_ybar=loss_ybar, mu0=matrix(data=0,nrow=2), lambda0=diag(10,2))
######################

split_merge_prop_prob_EVV <- function(obs, split_labs, group_assign, r, nu, y, mu0, lambda0){
  # split_labs is an array of length 2 indicating which entries in counts correspond
  # to the groups that are part of the split/merge
  # which_group is a scalar valued 1 or 2 indicating which of the groups we are considering
  # obs is the index for the observation being considered at this point
  # group_assign is an array of length n corresponding to group assignments for each obs
  # r and nu are scalar hyperparameters
  # counts is an array of the number of obs assigned to each group
  # y is the data
  # mu0 and lambda0 are the prior mean and covariance matrix
  
  
  # which(as.numeric(names(test)) %in% c(1,3))
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


sm_iter = 10

temp_group_assign = matrix(data = NA, nrow = sm_iter + 1, ncol = length(group_assign[s,]))
sm_probs = matrix(data = NA, nrow = sm_iter, ncol = length(y))
temp_group_assign[1,] = group_assign[s,]

# randomly select two observed data points y
sampled_obs = sample(x = 1:length(y), size = 2, replace = FALSE)

# check if in same group - if yes SPLIT
# if no, MERGE
lab1 = temp_group_assign[1, sampled_obs[1]]
lab2 = temp_group_assign[1, sampled_obs[2]]
move_type = ifelse(lab1 == lab2, "SPLIT", "MERGE")

# bookkeeping - group labels
subset_index = which(temp_group_assign %in% c(lab1, lab2)) 
anchor_obs_index = which(subset_index %in% sampled_obs)
subset_index_minus = subset_index[-anchor_obs_index] # except sampled observations i and j

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
        
      } else{ # for remaining scans after random launch state set
        
        sm_counts = table(temp_group_assign[scan,-obs])
        split_group_count_index = which(as.numeric(names(sm_counts)) %in% split_lab)
        current_obs_index = which(temp_group_assign == obs)
        split_group_lab_index1 = which(temp_group_assign[(sm_iter+1),] == split_lab[1])
        split_group_lab_index2 = which(temp_group_assign[(sm_iter+1),] == split_lab[2])
        
        # current observation under consideration cannot be included here
        if(obs %in% split_group_lab_index1)
          
          split_assign_prob = split_merge_prop_prob_EVV(
            obs = obs, split_labs = split_lab, r=r, 
            group_assign = temp_group_assign[scan,], nu = nu, 
            y = y, mu0 = mu0, lambda0= lambda0)
        
        
        
        sm_prop_index = sample(x = 1:2, size = 1, 
                               prob = split_assign_prob)
        
        temp_group_assign[scan,obs] = split_lab[sm_prop_index]
        sm_probs[scan,obs] = split_assign_prob[sm_prop_index]
        
      }
      
    }
  }
  
  # calculate & evaluate acceptance prob
  sm_counts = table(temp_group_assign[sm_iter+1,]) # update counts after scans
  # proposal probability
  prob1 = -Reduce(f = "+", x = log(sm_probs[sm_iter+1,])) # log1 - sum(logs)
  # prior ratio
  prob2_num = factorial(sm_counts[[split_group_count_index[1]]] -1)*factorial(sm_counts[[split_group_count_index[2]]] -1)
  prob2_denom = factorial(sm_counts[[split_group_count_index[1]]]+sm_counts[[split_group_count_index[2]]]-1)
  prob2 = log(alpha) + (log(prob2_num) - log(prob2_denom))
  
  # likelihood ratio
  
  ## component 1 - numerator I (group 1)
  prob3_num1 = 0
  for(obs_ind in 1:length(subset_index)){
    if(obs_ind == 1){
      # first observation --- prior predictive
      val = prior_pred_NinvW(y_i = y_i, mu0 = mu0,
                             r = r, lambda0 = lambda0, 
                             nu = nu)
      prob3_num1 = prob3_num1 + log(val)
      
    } else{
      # posterior predictive
      val = post_pred_EVV()
      
    }
  }
  
  
  ## component 2 - numerator II (group 2)
  
  ## component 3 - denominator (all in original group)
  
  # if new group created by split, give it a mean and variance
  
  # final bookkeeping 
  
  
  
  # if MERGE    
} else if(move_type == "MERGE"){
  
}



# calculate proposal probability

# calculate & evaluate acceptance prob

# if new group created by split, give it a mean and variance

# final bookkeeping 


}