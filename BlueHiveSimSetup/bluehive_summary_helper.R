##############################################################
############### BLUEHIVE SIM SUMMARY SCRIPT ##################
######################  DPM MODELS ###########################
##############################################################

## Author: Jonathan Klus
## Date: 24 February 2024
## Description: Call to helper function to summarize simulation study results
## run on the BlueHive server

###################### load packages and set seed ##############################
library(parallel)

source("./post_processing_inf.R")

############################### HELPER FUNCTIONS ###############################

# select results file
n=30
output_path = paste0("conjDEV_wellsep_n", n, "_noSM_sim_results_2024_02_20") 
SLURM_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

dir_name = paste0("./Summary/MODSUM_", output_path)

if(dir.exists(dir_name) == FALSE){
  dir.create(dir_name)
} 

# define truth  & run summary 
# mod_sum = mclapply(X = 1:100, mc.cores = n_cores, 
#                 FUN = function(x){
#                   output = readRDS(file = paste0("./", output_path, "/output_", x, ".rds"))
#                   sum1 = dpmm_summary(output = output, 
#                                print_phi_sum = TRUE,
#                                print_k_sum = TRUE, 
#                                make_traceplot = FALSE,
#                                burn_in = 2000, t_hold = 250, 
#                                num_dims = 2, 
#                                calc_perf = TRUE, 
#                                mu_true = output$truth$mu_true, 
#                                var_true = lapply(X = 1:length(output$truth$mu_true), 
#                                                  FUN = function(x){output$truth$var_true}), 
#                                assign_true = output$truth$assign_true, 
#                                equal_var = FALSE)
#                   return(sum1)
#                   })


output = readRDS(file = paste0("./", output_path, "/output_", SLURM_ID, ".rds"))
mod_sum = dpmm_summary(output = output, 
                       print_phi_sum = TRUE,
                       print_k_sum = TRUE, 
                       make_traceplot = FALSE,
                       burn_in = 2000, t_hold = 250, 
                       num_dims = 2, 
                       calc_perf = TRUE, 
                       mu_true = output$truth$mu_true, 
                       var_true = lapply(X = 1:length(output$truth$mu_true), 
                                         FUN = function(x){output$truth$var_true}), 
                       assign_true = output$truth$assign_true, 
                       equal_var = FALSE)



# save summary
saveRDS(object = mod_sum, 
        file = paste0("./Summary/MODSUM_", output_path,"/sum_",  SLURM_ID,".rds"))