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



