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

# 5d plot?? or just do multiple 2d plots ## should this data be log transformed
# before centering and scaling???
ggplot2::ggplot(data = exposure, mapping = aes(x = TSH, y = T3)) +
  ggplot2::geom_point(aes(color = data$Diagnosis))

ggplot2::ggplot(data = exposure, mapping = aes(x = TSH, y = T4)) +
  ggplot2::geom_point(aes(color = data$Diagnosis))

ggplot2::ggplot(data = exposure, mapping = aes(x = T3, y = T4)) +
  ggplot2::geom_point(aes(color = data$Diagnosis))


# 3d plot

# fig = plot_ly(data, x = ~TSH, y = ~T3, z = ~T4) #, color = ~am, colors = c('#BF382A', '#0C4B8E'))
# 
# fig = fig %>% add_markers()
# 
# fig = fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
#                                    
#                                    yaxis = list(title = 'Gross horsepower'),
#                                    
#                                    zaxis = list(title = '1/4 mile time')))
# 
# 
# fig

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
mod1_sum = dpmm_summary(output = output,
                        print_phi_sum = TRUE,
                        print_k_sum = TRUE,
                        make_traceplot = TRUE,
                        burn_in = 2000, t_hold = 250,
                        num_dims = 5,
                        calc_perf = FALSE,
                        equal_var = FALSE)
# save model summary

