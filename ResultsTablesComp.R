source("./Samplers/posterior_helper_fxns.R")
source("./Samplers/post_processing_inf.R")

# load R libraries
library(ggplot2)
library(stringr)
library(gridExtra)
library(dplyr)
library(knitr) # kable
library(kableExtra) # use for complex tables http://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf

# ARI and KL are typically skewed, so use the median (min, max) to summarize
# results across the 100 data sets!

# read in results database
data_path = "./results_database_SummaryLargePriorSS_2024_12_22.rds"
summary_table = readRDS(data_path)

####################################################################################################################

# output LATEX results to .txt file
sink(file = paste0("./ResultsCompOutput_",
                   stringr::str_extract(string = stringr::str_extract(string = data_path, pattern = "database_[:alnum:]+"), pattern = "[:alnum:]+$"),
                   "_",
                   stringr::str_extract(string = data_path, pattern = "[:digit:]{4}_[:digit:]{2}_[:digit:]{2}"),
                   ".txt"),
     type = "output")

# number of data sets completed
table1 = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(num_datasets = n()) %>%
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = num_datasets) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2)

xtable::xtable(table1, caption = "Counts")

ari_table = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(mean_ARI = paste0(round(median(ARI),2), " (", round(IQR(ARI),2), ")")) %>% 
  # dplyr::summarize(mean_ARI = paste0(round(median(ARI),2), " (", round(min(ARI),2), ",", round(max(ARI),2), ")")) %>% 
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = mean_ARI) %>%
  # tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = c(mean_ARI, min_ARI, max_ARI)) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2) 

kld_table = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(mean_KL = paste0(round(median(KL),2), " (", round(IQR(KL),2), ")")) %>% 
  # dplyr::summarize(mean_KL = paste0(round(median(KL),2), " (", round(min(KL),2), ",", round(max(KL),2), ")")) %>% 
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = mean_KL) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2)

mapk_table = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(MAP_K = paste0(round(median(MAP_K),2), " (", round(IQR(MAP_K),2),")")) %>% 
  # dplyr::summarize(MAP_K = paste0(round(median(MAP_K),2), " (", round(min(MAP_K),2), ",", round(max(MAP_K),2), ")")) %>% 
  tidyr::pivot_wider(names_from = c(Scenario, n_obs), values_from = MAP_K) %>%
  # keep first 2 columns with model info, sort remaining cols with results by n
  dplyr::select(1,2, gtools::mixedorder(names(.)[3:length(names(.))])+2)

s1_table = cbind(
  kld_table$Model, kld_table$SM,
  kld_table %>% dplyr::ungroup() %>% dplyr::select(`3wellsep_30`, `3wellsep_100`, `3wellsep_300`),
  ari_table %>% dplyr::ungroup() %>% dplyr::select(`3wellsep_30`, `3wellsep_100`, `3wellsep_300`),
  mapk_table %>% dplyr::ungroup() %>% dplyr::select(`3wellsep_30`, `3wellsep_100`, `3wellsep_300`)
)

xtable::xtable(s1_table, 
               digits = 1, 
               caption = "3wellsep", )

s2_table = cbind(
  kld_table$Model, kld_table$SM,
  kld_table %>% dplyr::ungroup() %>% dplyr::select(`3close_30`, `3close_100`, `3close_300`),
  ari_table %>% dplyr::ungroup() %>% dplyr::select(`3close_30`, `3close_100`, `3close_300`),
  mapk_table %>% dplyr::ungroup() %>% dplyr::select(`3close_30`, `3close_100`, `3close_300`)
)

xtable::xtable(s2_table, 
               digits = 1, 
               caption = "3close", )

s3_table = cbind(
  kld_table$Model, kld_table$SM,
  kld_table %>% dplyr::ungroup() %>% dplyr::select(`5grp3d_30`, `5grp3d_100`, `5grp3d_300`),
  ari_table %>% dplyr::ungroup() %>% dplyr::select(`5grp3d_30`, `5grp3d_100`, `5grp3d_300`),
  mapk_table %>% dplyr::ungroup() %>% dplyr::select(`5grp3d_30`, `5grp3d_100`, `5grp3d_300`)
)

xtable::xtable(s3_table, 
               digits = 1, 
               caption = "5grp3d", )


# important: return output to console!!
sink()



# make visualizations for results

ari_table_long = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(med_ARI = median(ARI), iqr_ARI = IQR(ARI)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_obs = as.numeric(n_obs))

kld_table_long = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(med_KL = median(KL), iqr_KL = IQR(KL)) %>%
  dplyr::mutate(n_obs = as.numeric(n_obs))

mapk_table_long = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(med_mapk = median(MAP_K), iqr_mapk = IQR(MAP_K)) %>%
  dplyr::mutate(n_obs = as.numeric(n_obs))

ggplot2::ggplot(data = ari_table_long %>% dplyr::filter(Scenario == "3wellsep"), 
                aes(x = n_obs, y = med_ARI)) +
  ggplot2::geom_line(aes(group = Model)) + # need to group by model AND split/merge!!! otherwise plot will look wierd
  ggplot2::geom_point(shape = 9, color = "blue") +
  ggplot2::theme_classic()
