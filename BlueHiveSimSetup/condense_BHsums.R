################################################################################
##################  CONDENSE POST PROCESSING RESULTS   #########################
##################               VERSION 1            ##########################
################################################################################
## AUTHOR: JONATHAN KLUS
## DATE: 26 FEBRUARY 2024
## DESCRIPTION: CONDENSE OUTPUT FROM SUMMARY FUNCTION TO MAKE RESULT TABLES

####################### LOAD ANY REQUIRED PACKAGES #############################

# homemade functions
source("./posterior_helper_fxns.R")

# load R libraries
library(ggplot2)
library(dplyr)
library(stringr)

######################## DEFINE HELPER FUNCTIONS################################
# running locally
# data_path = "//GENE.bst.rochester.edu/Projects/JKSTproj/BlueHive_Sim_Results/"

# running on server
data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/"

# list file extensions
file_ext = c(
  "MODSUM_conjUVV_wellsep_n30_noSM_sim_results_2024_02_18",
  "MODSUM_conjUVV_wellsep_n100_noSM_sim_results_2024_02_18",
  "MODSUM_conjUVV_wellsep_n300_noSM_sim_results_2024_02_18"
)


# for each extension, make a new element of a data table with simulation attributes
summary_table = data.frame(Model = NA, Scenario = NA, SM = NA, n_obs = NA, n_datasets = NA,
                           ARI = NA, KL = NA, Time = NA) # should add , SumTime = NA for mod summary compute time

for(ext_index in 1:length(file_ext)){
  # list files
  file_list = list.files(path = paste0(data_path, file_ext[ext_index]))
  
  # parse out attributes and save
  parsed_dirname = unlist(stringr::str_extract_all(string = file_ext[ext_index], pattern = "_[:alnum:]+"))
  summary_table[ext_index, "Model"] = stringr::str_remove(string = parsed_dirname[1], pattern = "_") 
  summary_table[ext_index, "Scenario"] = stringr::str_remove(string = parsed_dirname[2], pattern = "_") 
  summary_table[ext_index, "n_obs"] = stringr::str_remove(string = parsed_dirname[3], pattern = "_n") 
  summary_table[ext_index, "SM"] = stringr::str_remove(string = parsed_dirname[4], pattern = "_") 
  summary_table[ext_index, "n_datasets"] = length(file_list) 
  
  # loop through all files in f
  sum_bydataset = sapply(X = 1:length(file_list), 
         FUN = function(x){
           mod_sum = readRDS(paste0(data_path, file_ext[ext_index], "/", file_list[x]))
           c(mod_sum$mean_ARI, mod_sum$kl_div, as.numeric(mod_sum$fit_runtime))
         })
  avg_sum = colMeans(t(sum_bydataset))
  summary_table[ext_index, "ARI"] = avg_sum[1]
  summary_table[ext_index, "KL"] = avg_sum[2]
  summary_table[ext_index, "Time"] = avg_sum[3]
  
}

