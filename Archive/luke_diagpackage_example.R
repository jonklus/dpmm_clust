# Test Luke diagnostics package

# installation
# devtools::install_github("LukeDuttweiler/genMCMCDiag") 

# load package
library(genMCMCDiag)

# load simulation results
# result = readRDS("./MCMC_Runs/conjDEVsamp_minisimstudy_noSM_2024_01_09.rds")[[1]]
result_path = "/projects/jklus/JKSTproj/BlueHive_Sim_Results/Thyroid"

chain1 = readRDS(paste0(result_path, "/MODFIT__DEV_noSM_1.rds"))
chain2 = readRDS(paste0(result_path, "/MODFIT__DEV_noSM_2.rds"))
chain3 = readRDS(paste0(result_path, "/MODFIT__DEV_noSM_3.rds"))
chain4 = readRDS(paste0(result_path, "/MODFIT__DEV_noSM_4.rds"))

chain1_sm = readRDS(paste0(result_path, "/MODFIT__DEV_withSM_1.rds"))
chain2_sm = readRDS(paste0(result_path, "/MODFIT__DEV_withSM_2.rds"))
chain3_sm = readRDS(paste0(result_path, "/MODFIT__DEV_withSM_3.rds"))
chain4_sm = readRDS(paste0(result_path, "/MODFIT__DEV_withSM_4.rds"))

# mc_k = list(list(val = as.list(result$k[,1]), Posterior = NA))
# k_list = list(chain1$k, chain2$k, chain3$k, chain4$k)
mcmc_list = list(chain1$pairwise_mats$adj_by_iter, 
                 chain2$pairwise_mats$adj_by_iter, 
                 chain3$pairwise_mats$adj_by_iter, 
                 chain4$pairwise_mats$adj_by_iter)

mcmc_list_sm = list(chain1_sm$pairwise_mats$adj_by_iter, 
                 chain2_sm$pairwise_mats$adj_by_iter, 
                 chain3_sm$pairwise_mats$adj_by_iter, 
                 chain4_sm$pairwise_mats$adj_by_iter)

# make traceplots & get diagnostics
genDiagnostic(mhDraws = mcmc_list, 
              method = "lanfear", 
              distance = hammingDist)

genDiagnostic(mhDraws = mcmc_list_sm, 
              method = "lanfear", 
              distance = hammingDist)


# try dpmmDistance function -- need to have a list of chains, for each chain,
# need a list where each list element is an iteration with elements Zscore and Adj


# can also modify ggplot formatting
p1@diagnostics$traceplot + ggplot2::theme_classic()

# make traceplots for means