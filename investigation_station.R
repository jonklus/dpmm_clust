# checking DEV results

####################### LOAD ANY REQUIRED PACKAGES #############################

# homemade functions
source("./posterior_helper_fxns.R")

# load R libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(knitr) # kable
library(kableExtra) # use for complex tables http://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf 

######################## DEFINE HELPER FUNCTIONS################################
# running locally
# data_path = "//GENE.bst.rochester.edu/Projects/JKSTproj/BlueHive_Sim_Results/"

# running on server

# model summary
sum_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/Summary"
mod_sum = readRDS(file = paste0(sum_path, "/MODSUM_conjDEV_close_n30_withSM_sim_results_2024_02_20/sum_15.rds"))

# raw output
output_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results"
output = readRDS(file = paste0(output_path, "/conjDEV_close_n30_withSM_sim_results_2024_02_20/output_15.rds"))