################################################################################
############ POST PROCESSING AND FORMAL INFERENCE FUNCTIONS ####################
##################               VERSION 1            ##########################
################################################################################
## AUTHOR: JONATHAN KLUS
## DATE: 10 JANUARY 2024
## DESCRIPTION: CONVENIENCE FUNCTIONS FOR FINDING POSTERIOR SUMMARIES FOR INFERENCE
## AND DIAGNOSTICS FROM MCMC DRAWS IN A SIMULATION STUDY WITH MULTIPLE DATA SETS. 
## CAN PROVIDE OVERALL SUMMARY OR BY DATA SET. 

####################### LOAD ANY REQUIRED PACKAGES #############################

# homemade functions
source("./posterior_helper_fxns.R")
source("./post_processing_inf.R")

# load R libraries
library(ggplot2)
library(gridExtra)
library(plotly)
library(label.switching)
library(LaplacesDemon)
library(parallel)
library(stringr)
library(dplyr)
library(mclust)

########################### LOAD R DATA FILES ##################################
path = "~/../../projects/jklus/JKSTproj/BlueHive_Sim_Results/"
sim = "conjDEV_close_n300_noSM_sim_results_2024_02_19/output_"
x = 10
output = readRDS(file = paste0(path, sim, x, ".rds"))

########################### DO POST PROCESSING #################################

sum_output = dpmm_summary(output = output, print_phi_sum = TRUE,
                          print_k_sum = TRUE, make_traceplot = FALSE,
                          burn_in = 2000, t_hold = 250, num_dims = 2, 
                          calc_perf = TRUE, 
                          mu_true = output$truth$mu_true, 
                          var_true = lapply(X = 1:length(output$truth$mu_true), 
                                            FUN = function(x){output$truth$var_true}), 
                          assign_true = output$truth$assign_true, 
                          equal_var = FALSE)

