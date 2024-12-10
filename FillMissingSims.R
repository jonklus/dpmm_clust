# script to identify and list missing simulation runs with their seed numbers

source("./Samplers/posterior_helper_fxns.R")
source("./Samplers/post_processing_inf.R")

# load R libraries
library(ggplot2)
library(mclust)
library(LaplacesDemon)
library(mvtnorm)
library(stringr)
library(gridExtra)
library(dplyr)
library(knitr) # kable
library(kableExtra) # use for complex tables http://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf


######################## DEFINE HELPER FUNCTIONS################################
# running locally
# data_path = "//GENE.bst.rochester.edu/Projects/JKSTproj/BlueHive_Sim_Results/"

# running on server
# data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/SummaryManuscript"
directory_path = "/scratch/jklus/dpmSims/SummaryLargePriorSS"

# list file extensions
file_ext = list.dirs(directory_path)

# for each extension, make a new element of a data table with simulation attributes
summary_list = vector(mode = "list", length = length(file_ext))

for(ext_index in 2:length(file_ext)){ # skip first one, listing main directory
  # skip second one, Archive directory
  # list files
  file_list = list.files(path = file_ext[ext_index])
  if(length(file_list) > 0){
    
    temp_summary_table = data.frame(Path = NA, Model = NA, Scenario = NA, SM = NA, n_obs = NA, dataset_no = NA) 
    # parse out attributes and save
    short_file_ext = stringr::str_split(string = file_ext[ext_index], pattern = "/MODSUM")[[1]][2]
    parsed_dirname = unlist(stringr::str_extract_all(string = file_ext[ext_index], pattern = "_[:alnum:]+"))
    temp_summary_table[1:length(file_list), "Path"] = file_ext[ext_index]
    temp_summary_table[1:length(file_list), "Model"] = stringr::str_remove(string = parsed_dirname[1], pattern = "_")
    temp_summary_table[1:length(file_list), "Scenario"] = stringr::str_remove(string = parsed_dirname[2], pattern = "_")
    temp_summary_table[1:length(file_list), "n_obs"] = stringr::str_remove(string = parsed_dirname[3], pattern = "_n")
    temp_summary_table[1:length(file_list), "SM"] = stringr::str_remove(string = parsed_dirname[4], pattern = "_")
    # temp_summary_table[1:length(file_list), "Tag"] = stringr::str_remove(string = parsed_dirname[length(parsed_dirname)], pattern = "_")
    temp_summary_table[, "dataset_no"] = as.numeric(stringr::str_remove_all(string = file_list, pattern = "sum_|.rds"))
    
    summary_list[[ext_index-1]] = temp_summary_table
    
  }
  
}

head(summary_list[[1]], 10)