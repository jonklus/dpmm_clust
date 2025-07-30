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
data_path = "./Results/results_database_SummaryLargePriorSS_2024_12_22.rds"
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

# function to get Model and SM columns formatted for paper 
col_format <- function(x, type){
  # x is a column of labels 
  # type is string, either "Model" or "SM" corresponding to x
  
  if(type == "Model"){
    
    conjT = ifelse(stringr::str_detect(string = x, pattern = "^conj"),
                   "Conjugate ", "Non-Conjugate ")
    
    label = paste0(conjT,
      stringr::str_extract(string = x, pattern = "[:alpha:]{3}$"))
    
  } else if(type == "SM"){
    
    label = ifelse(stringr::str_detect(string = x, pattern = "noSM"),
                   "FALSE", "TRUE")
  }
  
  return(label)
  
}


s1_table = cbind(
  Model = kld_table$Model, SM = kld_table$SM,
  kld_table %>% dplyr::ungroup() %>% dplyr::select(`3wellsep_30`, `3wellsep_100`, `3wellsep_300`),
  ari_table %>% dplyr::ungroup() %>% dplyr::select(`3wellsep_30`, `3wellsep_100`, `3wellsep_300`),
  mapk_table %>% dplyr::ungroup() %>% dplyr::select(`3wellsep_30`, `3wellsep_100`, `3wellsep_300`)
)

s1_table$Model = col_format(x = s1_table$Model, type = "Model")
s1_table$SM = col_format(x = s1_table$SM, type = "SM")

xtable::xtable(s1_table, 
               digits = 1, 
               caption = "3wellsep")

s2_table = cbind(
  Model = kld_table$Model, SM = kld_table$SM,
  kld_table %>% dplyr::ungroup() %>% dplyr::select(`3close_30`, `3close_100`, `3close_300`),
  ari_table %>% dplyr::ungroup() %>% dplyr::select(`3close_30`, `3close_100`, `3close_300`),
  mapk_table %>% dplyr::ungroup() %>% dplyr::select(`3close_30`, `3close_100`, `3close_300`)
)

s2_table$Model = col_format(x = s2_table$Model, type = "Model")
s2_table$SM = col_format(x = s2_table$SM, type = "SM")

xtable::xtable(s2_table, 
               digits = 1, 
               caption = "3close")

s3_table = cbind(
  Model = kld_table$Model, SM = kld_table$SM,
  kld_table %>% dplyr::ungroup() %>% dplyr::select(`5grp3d_30`, `5grp3d_100`, `5grp3d_300`),
  ari_table %>% dplyr::ungroup() %>% dplyr::select(`5grp3d_30`, `5grp3d_100`, `5grp3d_300`),
  mapk_table %>% dplyr::ungroup() %>% dplyr::select(`5grp3d_30`, `5grp3d_100`, `5grp3d_300`)
)

s3_table$Model = col_format(x = s3_table$Model, type = "Model")
s3_table$SM = col_format(x = s3_table$SM, type = "SM")

xtable::xtable(s3_table, 
               digits = 1, 
               caption = "5grp3d")


# important: return output to console!!
sink()



# make visualizations for results

ari_table_long = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(med = median(ARI), iqr = IQR(ARI)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n_obs = factor(as.numeric(n_obs)), metric = "ARI")

ari_table_long$Model = col_format(x = ari_table_long$Model, type = "Model")
ari_table_long$SM = col_format(x = ari_table_long$SM, type = "SM")

kld_table_long = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(med = median(KL), iqr = IQR(KL)) %>%
  dplyr::mutate(n_obs = factor(as.numeric(n_obs)), metric = "KLD")

kld_table_long$Model = col_format(x = kld_table_long$Model, type = "Model")
kld_table_long$SM = col_format(x = kld_table_long$SM, type = "SM")

mapk_table_long = summary_table %>% 
  # dplyr::filter(Tag == "ENAR") %>%
  dplyr::group_by(Model, Scenario, SM, n_obs) %>%
  dplyr::summarize(med = median(MAP_K), iqr = IQR(MAP_K)) %>%
  dplyr::mutate(n_obs = factor(as.numeric(n_obs)), metric = "MAPK")

mapk_table_long$Model = col_format(x = mapk_table_long$Model, type = "Model")
mapk_table_long$SM = col_format(x = mapk_table_long$SM, type = "SM")

results_plot_table_long = rbind(ari_table_long, kld_table_long) # mapk_table_long
results_plot_table_long$SM = factor(results_plot_table_long$SM)
# results_plot_table_long$Model = factor(results_plot_table_long$Model)
results_plot_table_long$metric = factor(results_plot_table_long$metric)
results_plot_table_long$Conjugate = (!stringr::str_detect(string = results_plot_table_long$Model, pattern = "Non-Conjugate"))
results_plot_table_long$Type = factor(stringr::str_extract(string = results_plot_table_long$Model, pattern = "[:alpha:]{3}$"))
# thought -- pivot below by scenario for final plots??


ggplot2::ggplot(data = results_plot_table_long %>% 
                  dplyr::filter(Scenario == "3wellsep"), 
                                # Model %in% c("Conjugate DEE", "Conjugate DEV", "Conjugate UVV")), 
                aes(x = n_obs, y = med, shape = SM, lty = Model, color = Conjugate,
                    group = interaction(Model, SM))) +
  ggplot2::facet_wrap(facets = vars(metric), ncol = 3, scales = "free") +
  ggplot2::geom_line() + 
  ggplot2::geom_point(size = 2) +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = c("blue", "black")) + 
  # ggplot2::theme(legend.position = "bottom") + 
  ggplot2::xlab("Sample Size") +
  ggplot2::ylab("Value") +
  ggplot2::ggtitle("Results for Well-Separated Scenarios")

# do with just ARI instead...
ggplot2::ggplot(data = results_plot_table_long %>% 
                  dplyr::filter(Scenario == "3wellsep", metric == "ARI"), 
                # Model %in% c("Conjugate DEE", "Conjugate DEV", "Conjugate UVV")), 
                aes(x = n_obs, y = med, shape = SM, lty = Type, # color = Conjugate,
                    group = interaction(Model, SM))) +
  ggplot2::facet_wrap(facets = vars(Conjugate), ncol = 3, scales = "free") +
  ggplot2::geom_line() + 
  ggplot2::geom_point(size = 2) +
  ggplot2::theme_classic() +
  # ggplot2::scale_color_manual(values = c("blue", "black")) + 
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::xlab("Sample Size") +
  ggplot2::ylab("Value") +
  ggplot2::ggtitle("Results for Well-Separated Scenarios")

ggsave("./3wellsep_results.png")


ggplot2::ggplot(data = results_plot_table_long %>% 
                  dplyr::filter(Scenario == "3close"), 
                                # Model %in% c("Conjugate DEE", "Conjugate DEV", "Conjugate UVV")),
                aes(x = n_obs, y = med, shape = SM, lty = Model, color = Conjugate,
                    group = interaction(Model, SM))) +
  ggplot2::facet_wrap(facets = vars(metric), ncol = 3, scales = "free") +
  ggplot2::geom_line() + 
  ggplot2::geom_point(size = 2) +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = c("blue", "black")) + 
  # ggplot2::theme(legend.position = "bottom") + 
  ggplot2::xlab("Sample Size") +
  ggplot2::ylab("Value") +
  ggplot2::ggtitle("Results for Close Together Scenario: Conjugate Models")
ggsave("./3close_results.png")

ggplot2::ggplot(data = results_plot_table_long %>% 
                  dplyr::filter(Scenario == "5grp3d"),
                                # Model %in% c("Conjugate DEE", "Conjugate DEV", "Conjugate UVV")),
                aes(x = n_obs, y = med, shape = SM, lty = Model, color = Conjugate,
                    group = interaction(Model, SM))) +
  ggplot2::facet_wrap(facets = vars(metric), ncol = 3, scales = "free") +
  ggplot2::geom_line() + 
  ggplot2::geom_point(size = 2) +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = c("blue", "black")) + 
  # ggplot2::theme(legend.position = "bottom") + 
  ggplot2::xlab("Sample Size") +
  ggplot2::ylab("Value") +
  ggplot2::ggtitle("Results for 5 group 3D Scenario: Conjugate Models")
ggsave("./5grp3d_results.png")
