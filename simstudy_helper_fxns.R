################################################################################
################### SIMULATION STUDY HELPER FUNCTIONS ##########################
##################               VERSION 1            ##########################
################################################################################
## AUTHOR: JONATHAN KLUS
## DATE: 11 JULY 2023
## DESCRIPTION: CONVENIENCE FUNCTIONS FOR PERFORMING SIMULATION STUDIES USING
## DPMM FUNCTIONS AVAILABLE IN OTHER SCRIPTS IN THIS REPO. 

####################### LOAD ANY REQUIRED PACKAGES #############################
library(ggplot2)
library(tidyr)
library(parallel)
library(LaplacesDemon)

########################### GENERATE DATA ######################################

########################### SIMULATION CODE ####################################

# check parallelization settings
# n_cores = parallel::detectCores()
# if(n_cores > 10){
#   n_cores = 10 # on server, usually 12 cores 
# } else{
#   n_cores = round(n_cores/2)
# }

# grid of sample sizes
n_grid = c(10,30, 100, 200)

# grid of priors



# other sim settings
n_reps = 10^0  # number of times to repeat each simulation
n_iter = 10^3  # number of iterations

# grid of seeds
seeds = sample(x = 1:10^5, 
               size = n_reps)

# test combining grids
grid_names = c("alpha", "mu0", "Sigma0","a", "b", "d", "f", "n")
## on concentration param
alpha_grid = NULL 
## on mean
mu0_grid = NULL
Sigma0_grid = NULL
## on variance
a_grid = c(0.5, 1, 2)
b_grid = c(1, 2, 5, 10)
d_grid = NULL
f_grid = NULL
# on sample size
n_grid = c(10,20,30)

grid_list = list(alpha_grid, mu0_grid, Sigma0_grid, a_grid, b_grid, d_grid, f_grid, n_grid)
grid_include = sapply(grid_list, is.null) # check whether
# a,b or mu0, Sigma0 are the ones being varied -- i.e. if their grid args are NULL or not
master_grid = expand.grid(grid_list[grid_include == FALSE]) # include if not null, i.e. FALSE
colnames(master_grid) = grid_names[grid_include == FALSE]

# results = vector(mode = "list", length = nrow(master_grid))



for(i in 1:nrow(master_grid)){

      
      # names(results[count]) = paste0("n=", n_grid[k], " a=", alpha_grid[i], "b=", beta_grid[j])
      # results[[i]] = vector(mode = "list", length = 2)
      # names(results[[i]]) = c("Settings", "Results") # c("Settings", "Mean", "Variance", "k", "Traceplots")
      # results[[i]]$Settings = master_grid[i,]
      

      # need to add some sort of "try" function so all results aren't lost in case of failure
      # test out with a negative SD in normal 
      
      # use mclapply or some other parallelization here to cut down on compute time across reps
      
      results = parallel::mclapply(
        
        X = 1:n_reps,
        
        mc.preschedule = FALSE, 
        
        mc.cores = n_cores, 
        
        FUN = function(x){

 
           
           # generate data to use in each rep
           
           ################## simulate data ####################################
           
           #### need to turn this part into a helper function that will take arguments of seed, n, etc
           #### and output data set.... BUT ALSO RECALL THAT n_obs is a setting of the function so this
           #### step actually needs to occur within each rep
           set.seed(seeds[x])
          
            w = c(.35, .25, .4) # should this be varied as well as part of the sim study??? make this optional
            means = list(
              #c(10, 10),
              c(0, -5),
              c(0, 10),
              c(12,2)
            )
           
           var = diag(9, length(means[[1]])) # variances known, diagonal, and equal
           
           assign = sample(x = 1:length(means), size = 30, replace = TRUE, prob = w)
           y = lapply(X = assign,
                      FUN = function(x){
                        t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
                      })
           y_matrix = matrix(data = unlist(y), ncol = 2, byrow = TRUE)
           
           data = data.frame(
             assign = assign,
             y1 = y_matrix[,1],
             y2 = y_matrix[,2]
           )
           
           settings = list(S = S, n_reps = n_reps, seed = seeds[x], y = y, Sigma0 = Sigma0, alpha = alpha, 
                           a = a, b = b, mu0 = mu0, k_init = k_init, d = d, f = f)
           
           mcmc_res = try(MVN_CRP_sampler(S = n_iter, seed = seeds[x], y = y, alpha = 1, Sigma0 = Sigma0,
                               mu0 = mu0, a = 1, b = 10, d = 2, f = 1, sigma_hyperprior = TRUE,
                               k_init = 1, verbose = TRUE, print_iter = 1000)) 
           
           return(list(mcmc = mcmc_res, data = y, settings = settings))
           
         })


      # save results
      
      SLURM_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # slurm job id to use 
      obj_name = paste0("dpmm_sim_", paste0(colnames(master_grid), master_grid[1,], collapse ="_"), 
                        "_" , Sys.Date(), "_", SLURM_ID, ".rds")
      saveRDS(object = results, file = obj_name)
      

    }



################################# METRICS ######################################

calc_KL_diverg <- function(y, mu_est, Sigma_est, group_assign, true_assign, mu_true, 
                           Sigma_true, equal_var_assump = FALSE){
  
  # function to calculate KL divergence between truth and posterior mean for
  # a particular clustering (e.g. k=3) after running MCMC
  
  # y, mu_est, Sigma_est are lists - mu and sigma are MCMC output that has been
  # if equal var assumption is true -- Sigma_est and Sigma_true are matrices
  # filtered by k and processed to correct label switching so that they are aligned
  # mu_true, and Sigma_true are the ground truth from which data were simulated
  # truth is a vector of true group assignments
  # group_assign is a S*n matrix of estimated group assignments
  # equal_var_assump is a logical that determines whether we assume the groups
  # had unique variances estimated as part of the MCMC
  
  # assumes MV normal likelihood
  
    if(equal_var_assump == FALSE){
      # each group has its own estimated variance
      
      # calculate density of truth
      true_dens = sapply(X = 1:length(y), 
                         FUN = function(x){
                           mvtnorm::dmvnorm(x = y[[x]][,1], 
                                            mean = mu_true[[truth[x]]][,1], 
                                            sigma = Sigma_true[[truth[x]]])
                           
                         })
        
        
      
      # calculate density of estimates
      est_dens_i = matrix(data = NA, nrow = length(mu), ncol = length(y))
      for(iter in 1:length(mu)){
        
        est_dens_i[iter,] = sapply(X = 1:length(y), 
                                 FUN = function(x){
                                   mvtnorm::dmvnorm(x = y[[x]][,1], 
                                                    mean = mu_est[[group_assign[iter,x]]][,1], 
                                                    sigma = Sigma_est[[group_assign[iter,x]]])
                                   
                                 })
        
      }
      
      est_dens = colMeans(est_dens_i) # average over all iterations
    
  } else{
    # pooled variance, assumed equal across groups
    
    # calculate density of truth
    true_dens = sapply(X = 1:length(y), 
                       FUN = function(x){
                         mvtnorm::dmvnorm(x = y[[x]][,1], 
                                          mean = mu_true[[truth[x]]][,1], 
                                          sigma = Sigma_true)
                         
                       })
    
    
    
    # calculate density of estimates
    est_dens_i = matrix(data = NA, nrow = length(mu), ncol = length(y))
    for(iter in 1:length(mu)){
      
      est_dens_i[iter,] = sapply(X = 1:length(y), 
                                 FUN = function(x){
                                   mvtnorm::dmvnorm(x = y[[x]][,1], 
                                                    mean = mu_est[[group_assign[iter,x]]][,1], 
                                                    sigma = Sigma_est)
                                   
                                 })
      
    }
    
    est_dens = colMeans(est_dens_i) # average over all iterations
    
  }
  
  # now do we average over densities in estimates then calc KL div, or  calculate KL divergence
  # at each iteration and then average over all iterations?
  
  # option 1 - average over densities, then calculate KL
  
  
  # option 2 - calculate KL divergence at each iteration
  
  kl_div = LaplacesDemon::KLD(px = true_dens, py = est_dens)
  
  return(kl_div)
}

