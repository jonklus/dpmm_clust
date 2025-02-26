################## make traceplots for manuscript ########################

library(ggplot2)
library(gridExtra)

source("./Samplers/posterior_helper_fxns.R")

############################### DEV model #####################################

# sum_13.rds

output1 = readRDS(paste0("../MCMC_Runs/",
              "MODSUM_conjDEV_3wellsep_n100_noSM_sim_results_2024_12_07/",
              "sum_13.rds"))

# num groups
num_groups1 = sapply(X=1:nrow(output1$group_assign), 
                     FUN = function(x){
                       length(unique(output1$group_assign[x,]))
                     })[2001:12000]
ggplot2::ggplot(data = NULL, 
                aes(x=1:10^4, 
                    y = num_groups1)) +
  ggplot2::geom_line() +
  ggplot2::theme_classic() +
  ggplot2::xlab("iter") + 
  ggplot2::ylab("param") 

# means
make_traceplot(param_list_by_k = output1$mean_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Mean", p = 2, k = 3,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

make_traceplot(param_list_by_k = output1$mean_list_by_k_stephens[[2]],
               component_no = 2, param_type = "Mean", p = 2, k = 3,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

make_traceplot(param_list_by_k = output1$mean_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Mean", p = 2, k = 4,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

make_traceplot(param_list_by_k = output1$mean_list_by_k_stephens[[2]],
               component_no = 2, param_type = "Mean", p = 2, k = 4,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

# variances
make_traceplot(param_list_by_k = output1$var_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Var", p = 1, k = 3,
               legend = FALSE, log_scale = TRUE, title_text = NULL)

make_traceplot(param_list_by_k = output1$var_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Var", p = 1, k = 4,
               legend = FALSE, log_scale = TRUE, title_text = NULL)

# r
ggplot2::ggplot(data = NULL, 
                aes(x=1:10^4, 
                    y = output1$extra_params[2001:12000,"r"])) +
  ggplot2::geom_line() +
  ggplot2::theme_classic() +
  ggplot2::xlab("iter") + 
  ggplot2::ylab("param") 

# mixing proportions
mixprop_traceplot(output1$group_assign, burn_in = 2000)


# sum_15.rds 

output2 = readRDS(paste0("../MCMC_Runs/",
                         "MODSUM_conjDEV_3wellsep_n100_noSM_sim_results_2024_12_07/",
                         "sum_15.rds"))

# num groups
num_groups2 = sapply(X=1:nrow(output2$group_assign), 
                     FUN = function(x){
                       length(unique(output2$group_assign[x,]))
                     })[2001:12000]
ggplot2::ggplot(data = NULL, 
                aes(x=1:10^4, 
                    y = num_groups2)) +
  ggplot2::geom_line() +
  ggplot2::theme_classic() +
  ggplot2::xlab("iter") + 
  ggplot2::ylab("param") 

# means
make_traceplot(param_list_by_k = output2$mean_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Mean", p = 2, k = 4,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

make_traceplot(param_list_by_k = output2$mean_list_by_k_stephens[[2]],
               component_no = 2, param_type = "Mean", p = 2, k = 4,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

make_traceplot(param_list_by_k = output2$mean_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Mean", p = 2, k = 5,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

make_traceplot(param_list_by_k = output2$mean_list_by_k_stephens[[2]],
               component_no = 2, param_type = "Mean", p = 2, k = 5,
               legend = FALSE, log_scale = FALSE, title_text = NULL)

# variances
make_traceplot(param_list_by_k = output2$var_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Var", p = 1, k = 4,
               legend = FALSE, log_scale = TRUE, title_text = NULL)

make_traceplot(param_list_by_k = output2$var_list_by_k_stephens[[2]],
               component_no = 1, param_type = "Var", p = 1, k = 5,
               legend = FALSE, log_scale = TRUE, title_text = NULL)

# r
ggplot2::ggplot(data = NULL, 
                aes(x=1:10^4, 
                    y = output2$extra_params[2001:12000,"r"])) +
  ggplot2::geom_line() +
  ggplot2::theme_classic() +
  ggplot2::xlab("iter") + 
  ggplot2::ylab("param") 

# mixing proportions
mixprop_traceplot(output2$group_assign, burn_in = 2000)


