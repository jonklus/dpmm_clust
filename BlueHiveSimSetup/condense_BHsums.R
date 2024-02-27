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
library(knitr) # kable
library(kableExtra) # use for complex tables http://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf 

######################## DEFINE HELPER FUNCTIONS################################
# running locally
# data_path = "//GENE.bst.rochester.edu/Projects/JKSTproj/BlueHive_Sim_Results/"

# running on server
data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/Summary"

# list file extensions
file_ext = list.dirs(data_path)

# for each extension, make a new element of a data table with simulation attributes
summary_table = data.frame(Model = NA, Scenario = NA, SM = NA, n_obs = NA, n_datasets = NA, dataset_no = NA,
                           ARI = NA, KL = NA, Time = NA) # should add , SumTime = NA for mod summary compute time

for(ext_index in 2:length(file_ext)){ # skip first one, listing main directory
  # list files
  file_list = list.files(path = file_ext[ext_index])
  if(length(file_list) > 0){
    
    temp_summary_table = data.frame(Model = NA, Scenario = NA, SM = NA, n_obs = NA, n_datasets = NA, dataset_no = NA,
                                    ARI = NA, KL = NA, Time = NA) # should add , SumTime = NA for mod summary compute time
    # parse out attributes and save
    short_file_ext = stringr::str_split(string = file_ext[ext_index], pattern = "/MODSUM")[[1]][2]
    parsed_dirname = unlist(stringr::str_extract_all(string = file_ext[ext_index], pattern = "_[:alnum:]+"))
    temp_summary_table[1:length(file_list), "Model"] = stringr::str_remove(string = parsed_dirname[3], pattern = "_") 
    temp_summary_table[1:length(file_list), "Scenario"] = stringr::str_remove(string = parsed_dirname[4], pattern = "_") 
    temp_summary_table[1:length(file_list), "n_obs"] = stringr::str_remove(string = parsed_dirname[5], pattern = "_n") 
    temp_summary_table[1:length(file_list), "SM"] = stringr::str_remove(string = parsed_dirname[6], pattern = "_") 
    temp_summary_table[1:length(file_list), "n_datasets"] = length(file_list) 
    temp_summary_table[, "dataset_no"] = as.numeric(stringr::str_remove_all(string = file_list, pattern = "sum_|.rds"))
    
    # loop through all files in f
    sum_bydataset = t(sapply(X = 1:length(file_list), 
                             FUN = function(x){
                               mod_sum = readRDS(paste0(file_ext[ext_index], "/", file_list[x]))
                               c(mod_sum$mean_ARI, mod_sum$kl_div, as.numeric(mod_sum$fit_runtime))
                             }))
    temp_summary_table[, "ARI"] = sum_bydataset[,1]
    temp_summary_table[, "KL"] = sum_bydataset[,2]
    temp_summary_table[, "Time"] = sum_bydataset[,3]
    summary_table = rbind(summary_table, temp_summary_table)
    
  }
  
}




# make individual summary tables for ARI, KL, Time
summary_table[-1,] %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(n())

summary_table[-1,] %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(ARI = mean(ARI))
