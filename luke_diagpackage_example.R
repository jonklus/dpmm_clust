# Test Luke diagnostics package

# installation
# devtools::install_github(‘LukeDuttweiler/genMCMCDiag’) 

# load package
library(genMCMCDiag)

# load simulation results
result = readRDS("./MCMC_Runs/conjDEVsamp_minisimstudy_noSM_2024_01_09.rds")[[1]]
mc_k = list(list(val = as.list(result$k[,1]), Posterior = NA))

# make traceplots for k
p1 = genDiagnostic(mhDraws = mc_k, 
              method = "standard", 
              diagnostics = "traceplot") 


# can also modify ggplot formatting
p1@diagnostics$traceplot + ggplot2::theme_classic()

# make traceplots for means