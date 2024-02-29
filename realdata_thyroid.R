## DPMM Real Data Example

# load necessary functions
source("./Multivariate_DPMM_unknownvar_DEV.R")
source("./Multivariate_DPMM_unknownvar_DEE.R")
source("./Multivariate_DPMM_unknownvar_UVV.R")
source("./posterior_helper_fxns.R")
source("./post_processing_inf.R")

# load R libraries
library(ggplot2)
library(gridExtra)
library(plotly)
library(label.switching)
library(LaplacesDemon)
library(parallel)
library(stringr)
library(mclust)

# visualize data
data = mclust::thyroid
true_assign = data$Diagnosis
exposure = scale(data[,2:6]) # center and scale exposure data

# put exposure data into format that model will accept
y = lapply(X = 1:nrow(exposure), 
           FUN = function(x){
             matrix(exposure[x,], nrow = length(exposure[x,]))
           })

# fit model

mod1 = MVN_CRP_sampler_DEV(
  S = 12000, seed = 516, y = y,
  alpha = 1, r = 10, g = 1, h = 50,
  sigma_hyperprior = FALSE, fix_r = FALSE,
  mu0 = matrix(rep(0, ncol(exposure)), ncol = 1),
  a = 1, b = 50,
  k_init = 1, diag_weights = FALSE,
  verbose = TRUE, split_merge = TRUE)

# summarize model fit