post_pred_EVV <- function(obs, which_group, r, sm_counts, nu, y, ybar, loss_ybar, mu0, lambda0){
  
  loss_mu0 = (ybar[[which_group]] - mu0)%*%t(ybar[[which_group]] - mu0)
  
  mu_n = ((1/r)*mu0[,1] + sm_counts[which_group]*ybar[[which_group]])/((1/r) + sm_counts[which_group])
  
  nu_n = nu + sm_counts[which_group] - nrow(mu0) + 1
  
  k_n = (r+sm_counts[which_group]+1)/((r+sm_counts[which_group])*(nu_n))
  
  lambda_n = k_n*(lambda0 + loss_ybar[[which_group]] + 
    loss_mu0*((sm_counts[which_group]/r)/(1/r + sm_counts[which_group])))
  
  # print(mu_n)
  # print(lambda_n)
  # print(nu_n)
  
  val = LaplacesDemon::dmvt(x = y[[obs]][,1], 
                             mu = mu_n[,1], 
                             S = lambda_n, 
                             df = nu_n)

  return(sm_counts[which_group]*val)
  
}

# test line
post_pred_EVV(obs=2,which_group=1, r=10, sm_counts=c(10,11), nu=2, y=y, ybar=ybar, 
              loss_ybar=loss_ybar, mu0=matrix(data=0,nrow=2), lambda0=diag(10,2))
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

# test line
split_merge_prop_prob_EVV(obs=2, split_labs = c(2,3), r=10, 
                          group_assign = group_assign, nu=2, y=y,
                          mu0=matrix(data=0,nrow=2), lambda0=diag(10,2))

