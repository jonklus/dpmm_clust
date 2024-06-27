# model summary lookup

# DEE 
data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/Summary_ENAR/"
folder_name = "MODSUM_conjDEE_3close_n30_noSM_sim_results_ENAR/"
file_name = "sum_1.rds"

result = readRDS(file = paste0(data_path, folder_name, file_name))
result$settings

# DEV
data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/Summary_ENAR/"
folder_name = "MODSUM_conjDEV_3close_n300_withSM_sim_results_ENAR/"
file_name = "sum_1.rds"

result = readRDS(file = paste0(data_path, folder_name, file_name))
result$settings

# UVV 
data_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/Summary_ENAR/"
folder_name = "MODSUM_conjUVV_3close_n30_withSM_sim_results_ENAR/"
file_name = "sum_1.rds"

result = readRDS(file = paste0(data_path, folder_name, file_name))
result$settings
