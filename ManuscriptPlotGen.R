library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(scatterplot3d)

#################### helper functions #########################################

simData <- function(scenario = "3wellsep", seed = 516, n = 30, means, var, w){
  
  set.seed(seeds[x])
  assign = sample(x = 1:length(means), size = n, replace = TRUE, prob = w)
  y = lapply(X = assign,
             FUN = function(x){
               t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
             })
  
  y_matrix = matrix(data = unlist(y), ncol = length(means[[1]]), byrow = TRUE)
  
  if(length(means[[1]]) == 2){
    
    data = data.frame(
      assign = factor(assign),
      y1 = y_matrix[,1],
      y2 = y_matrix[,2]
    )
    
  } else if(length(means[[1]]) == 3){
    
    data = data.frame(
      assign = factor(assign),
      y1 = y_matrix[,1],
      y2 = y_matrix[,2],
      y3 = y_matrix[,3]
    )
    
  }
  

  return(data)
  
}

################## simulate data ##############################################
seeds = readRDS("./BHsimseeds.rds")
x = 5 
scenario = "3close"


if(scenario == "3close"){
  
  w = c(.35, .25, .4)
  
  means = list(
    #c(10, 10),
    c(-1, -5),
    c(0, 10),
    c(10,2)
  )
  
  var = diag(4, length(means[[1]])) # variances known, diagonal, and equal
  
  y30 = simData(scenario = "3close", seed = seeds[x], n = 30, 
                means = means, var = var, w = w)
  y100 = simData(scenario = "3close", seed = seeds[x], n = 100, 
                 means = means, var = var, w = w)
  y300 = simData(scenario = "3close", seed = seeds[x], n = 300, 
                 means = means, var = var, w = w)
  
  p1 = ggplot(data = y30, aes(x = y1, y = y2, color = assign)) +
    geom_point(size = 0.75) +
    theme_classic() +
    theme(plot.title = element_text(size=12)) +
    ggtitle("n=30") +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(-10, 15)) +
    scale_y_continuous(limits = c(-10, 15)) 
  
  p2 = ggplot(data = y100, aes(x = y1, y = y2, color = assign)) + 
    geom_point(size = 0.75) +
    theme_classic() +
    theme(plot.title = element_text(size=12)) +
    ggtitle("n=100") +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(-10, 15)) +
    scale_y_continuous(limits = c(-10, 15)) 
  
  p3 = ggplot(data = y300, aes(x = y1, y = y2, color = assign)) + 
    geom_point(size = 0.75) +
    theme_classic() +
    theme(plot.title = element_text(size=12)) +
    ggtitle("n=300") +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(-10, 15)) +
    scale_y_continuous(limits = c(-10, 15)) 
  
  png(filename = "./manuscript_plots/3close_3blockplot.png", res = 300, width = 9, height = 3, units = "in")
  grid.arrange(p1, p2, p3, nrow = 1, respect = TRUE)
  dev.off()
  
  
} else if(scenario == "3wellsep"){
  
  w = c(1/3, 1/3, 1/3)
  
  means = list(
    c(-20, 20),
    c(20, -20),
    c(0, 0)
  )
  
  var = diag(10, length(means[[1]])) # variances diagonal, and equal
  
  y30 = simData(scenario = "3wellsep", seed = seeds[x], n = 30, 
                means = means, var = var, w = w)
  y100 = simData(scenario = "3wellsep", seed = seeds[x], n = 100, 
                 means = means, var = var, w = w)
  y300 = simData(scenario = "3wellsep", seed = seeds[x], n = 300, 
                 means = means, var = var, w = w)
  
  p1 = ggplot(data = y30, aes(x = y1, y = y2, color = assign)) + 
    geom_point(size = 0.75) +
    theme_classic() +
    theme(plot.title = element_text(size=12)) +
    ggtitle("n=30") +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(-30, 30), breaks = seq(-30,30,15)) +
    scale_y_continuous(limits = c(-30, 30), breaks = seq(-30,30,15))
  
  p2 = ggplot(data = y100, aes(x = y1, y = y2, color = assign)) + 
    geom_point(size = 0.75) +
    theme_classic() +
    theme(plot.title = element_text(size=12)) +
    ggtitle("n=100") +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(-30, 30), breaks = seq(-30,30,15)) +
    scale_y_continuous(limits = c(-30, 30), breaks = seq(-30,30,15)) 
  
  p3 = ggplot(data = y300, aes(x = y1, y = y2, color = assign)) + 
    geom_point(size = 0.75) +
    theme_classic() +
    theme(plot.title = element_text(size=12)) +
    ggtitle("n=300") +
    theme(legend.position="none") +
    scale_x_continuous(limits = c(-30, 30), breaks = seq(-30,30,15)) +
    scale_y_continuous(limits = c(-30, 30), breaks = seq(-30,30,15)) 
  
  png(filename = "./manuscript_plots/3wellsep_3blockplot.png", res = 300, width = 9, height = 3, units = "in")
  grid.arrange(p1, p2, p3, nrow = 1, respect = TRUE)
  dev.off()
  
} else if(scenario == "5grp3d"){
  
  w = rep(0.2, 5)
  means = list(
    c(4, 4, -2),
    c(2, 3, 2),
    c(3, 2, 5),
    c(8, 8, 0),
    c(9, 7, 3)
  )
  
  var = diag(0.25, length(means[[1]])) #+ matrix(data = 0.1, ncol=3, nrow=3) - diag(0.1,3)
  
  y30 = simData(scenario = "5grp3d", seed = seeds[x], n = 30, 
                means = means, var = var, w = w)
  y100 = simData(scenario = "5grp3d", seed = seeds[x], n = 100, 
                means = means, var = var, w = w)
  y300 = simData(scenario = "5grp3d", seed = seeds[x], n = 300, 
                means = means, var = var, w = w)
  
  # make plots
  png(filename = "./manuscript_plots/k53d_3blockplot_30.png", res = 300, width = 6, height = 6, units = "in")
  p1 = scatterplot3d(x = y30$y1, y = y30$y2, z = y30$y3, color = y30$assign, angle = -45, 
                     xlab = "", ylab = "", zlab = "",
                     pch = 20, main = "n=30", 
                     xlim = c(0,12), ylim = c(0,10), zlim = c(-4,8))
  dev.off()
  
  png(filename = "./manuscript_plots/k53d_3blockplot_100.png", res = 300, width = 6, height = 6, units = "in")
  p2 = scatterplot3d(x = y100$y1, y = y100$y2, z = y100$y3, color = y100$assign, angle = -45, 
                     xlab = "", ylab = "", zlab = "",
                     pch = 20, main = "n=100", 
                     xlim = c(0,12), ylim = c(0,10), zlim = c(-4,8))
  dev.off()
  
  
  png(filename = "./manuscript_plots/k53d_3blockplot_300.png", res = 300, width = 6, height = 6, units = "in")
  p3 = scatterplot3d(x = y300$y1, y = y300$y2, z = y300$y3, color = y300$assign, angle = -45, 
                     xlab = "", ylab = "", zlab = "",
                     pch = 20, main = "n=300", 
                     xlim = c(0,12), ylim = c(0,10), zlim = c(-4,8))
  dev.off()

}