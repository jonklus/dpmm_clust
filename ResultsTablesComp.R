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
data_path = "/scratch/jklus/dpmSims/SummaryKMeansInit"

# dev_withSM = readRDS(paste0(data_path, "/MODSUM_conjDEV_3close_n30_withSM_sim_results_ENAR/sum_1.rds"))

# list file extensions
file_ext = list.dirs(data_path)

# for each extension, make a new element of a data table with simulation attributes
summary_table = data.frame(Model = NA, Scenario = NA, SM = NA, Tag = NA, n_obs = NA, dataset_no = NA,
                           ARI = NA, KL = NA, Time = NA, MAP_K = NA) # should add , SumTime = NA for mod summary compute time

for(ext_index in 2:length(file_ext)){ # skip first one, listing main directory
  # skip second one, Archive directory
  # list files
  file_list = list.files(path = file_ext[ext_index])
  if(length(file_list) > 0){
    
    temp_summary_table = data.frame(Model = NA, Scenario = NA, SM = NA, Tag = NA, n_obs = NA, dataset_no = NA,
                                    ARI = NA, KL = NA, Time = NA, MAP_K = NA) # should add , SumTime = NA for mod summary compute time
    # parse out attributes and save
    short_file_ext = stringr::str_split(string = file_ext[ext_index], pattern = "/MODSUM")[[1]][2]
    parsed_dirname = unlist(stringr::str_extract_all(string = file_ext[ext_index], pattern = "_[:alnum:]+"))
    temp_summary_table[1:length(file_list), "Model"] = stringr::str_remove(string = parsed_dirname[1], pattern = "_")
    temp_summary_table[1:length(file_list), "Scenario"] = stringr::str_remove(string = parsed_dirname[2], pattern = "_")
    temp_summary_table[1:length(file_list), "n_obs"] = stringr::str_remove(string = parsed_dirname[3], pattern = "_n")
    temp_summary_table[1:length(file_list), "SM"] = stringr::str_remove(string = parsed_dirname[4], pattern = "_")
    # temp_summary_table[1:length(file_list), "Tag"] = stringr::str_remove(string = parsed_dirname[length(parsed_dirname)], pattern = "_")
    temp_summary_table[, "dataset_no"] = as.numeric(stringr::str_remove_all(string = file_list, pattern = "sum_|.rds"))
    
    # loop through all files in directory
    sum_bydataset = t(sapply(X = 1:length(file_list), 
                             FUN = function(x){
                               mod_sum = readRDS(paste0(file_ext[ext_index], "/", file_list[x]))
                               c(mod_sum$mean_ARI, 
                                 mod_sum$kl_div, 
                                 as.numeric(mod_sum$fit_runtime),
                                 as.numeric(names(mod_sum$k_freqtab)[which.max(as.numeric(mod_sum$k_freqtab))])
                               )
                             }))
    
    temp_summary_table[, "ARI"] = sum_bydataset[,1]
    temp_summary_table[, "KL"] = sum_bydataset[,2]
    temp_summary_table[, "Time"] = sum_bydataset[,3]
    temp_summary_table[, "MAP_K"] = sum_bydataset[,4]
    summary_table = rbind(summary_table, temp_summary_table)
    
  }
  
}

summary_table = summary_table[-1,] # get rid of first row of NAs 
#vdplyr::filter(Scenario == "5grp3d")

head(summary_table, 10)

summary_table %>%
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(Count = n()) %>%
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = Count) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2) %>%
  knitr::kable(caption = "Number of data sets completed")

# time
summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(Time = round(mean(Time),2)) %>%
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = Time) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2) %>%
  knitr::kable(caption = "Average computation time")

# ARI and KL are typically skewed, so use the median (min, max) to summarize
# results across the 100 data sets!
summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(mean_ARI = paste0(round(median(ARI),2), " (", round(min(ARI),2), ",", round(max(ARI),2), ")")) %>% 
                   # mean_ARI = round(mean(ARI),2), 
                   # min_ARI = round(min(ARI),2),
                   # max_ARI = round(max(ARI),2)) %>%
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = mean_ARI) %>%
  # tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = c(mean_ARI, min_ARI, max_ARI)) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2) %>%
  knitr::kable(caption = "Adjusted RAND Index")

summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(mean_KL = paste0(round(median(KL),2), " (", round(min(KL),2), ",", round(max(KL),2), ")")) %>% 
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = mean_KL) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2) %>%
  knitr::kable(caption = "KL Divergence")

summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(MAP_K = paste0(round(median(MAP_K),2), " (", round(min(MAP_K),2), ",", round(max(MAP_K),2), ")")) %>% 
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = MAP_K) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2) %>%
  knitr::kable(caption = "MAP(K)")