library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(scatterplot3d)

################## simulate data ##############################################



seeds = readRDS("./BHsimseeds.rds")
nreps = 10

### well separate K=3 case

# w = c(0.4, 0.3, 0.3)
# means = list(
#   c(-20, 20),
#   c(20, -20),
#   c(0, 0)
# )
# 
# var = diag(10, length(means[[1]])) # variances known, diagonal, and equal

## close K=3 case
# w = c(.35, .25, .4)
# means = list(
#   #c(10, 10),
#   c(0, -5),
#   c(0, 8),
#   c(10,0)
# )
# 
# var = diag(5, length(means[[1]])) # variances diagonal, and equal

# 3D K=5 case
w = rep(0.2, 5)
means = list(
  c(4, 4, -2),
  c(2, 3, 2),
  c(3, 2, 5),
  c(8, 8, 0),
  c(9, 7, 3)
)

var = diag(0.25, length(means[[1]])) #+ matrix(data = 0.1, ncol=3, nrow=3) - diag(0.1,3)

yreps30 = lapply(X = 1:nreps, 
                 FUN = function(x){
                   set.seed(seeds[x])
                   assign = sample(x = 1:length(means), size = 30, replace = TRUE, prob = w)
                   y = lapply(X = assign,
                              FUN = function(x){
                                t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
                              })
                   return(list(y = y, assign = assign))
                 })

yreps100 = lapply(X = 1:nreps, 
                  FUN = function(x){
                    set.seed(seeds[x])
                    assign = sample(x = 1:length(means), size = 100, replace = TRUE, prob = w)
                    y = lapply(X = assign,
                               FUN = function(x){
                                 t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
                               })
                    return(list(y = y, assign = assign))
                  })

yreps300 = lapply(X = 1:nreps, 
                  FUN = function(x){
                    set.seed(seeds[x])
                    assign = sample(x = 1:length(means), size = 300, replace = TRUE, prob = w)
                    y = lapply(X = assign,
                               FUN = function(x){
                                 t(mvtnorm::rmvnorm(n = 1, mean = means[[x]], sigma = var))
                               })
                    return(list(y = y, assign = assign))
                  })

# make a grid of plots
pltlist = list()

rep = 5

######### n=30
###############
y_matrix = matrix(data = unlist(yreps30[[rep]]$y), ncol = 3, byrow = TRUE)
assign = yreps30[[rep]]$assign

data = data.frame(
  assign = assign,
  y1 = y_matrix[,1],
  y2 = y_matrix[,2],
  y3 = y_matrix[,3]
)

# p = ggplot(data = data, aes(x = y1, y = y2)) + #, label = rownames(data))) +
#   geom_point(color = assign, size = 0.75) +
#   #geom_text(size = 3, hjust = 0, nudge_x = 0.5, color = assign) +
#   # geom_text(size = 3, color = assign) +
#   theme_classic() +
#   theme(plot.title = element_text(size=18)) + 
#   ggtitle("n=30")

attach(data)
p = scatterplot3d(x = y1, y = y2, z = y3, color = assign, angle = -45, 
                  pch = 20, main = "n=30", xlim = c(0,12), ylim = c(0,10), zlim = c(-4,8))
ggsave(filename = "./k53d_3blockplot_30.png", device = "png", plot = p)

pltlist[[1]] = p

######### n=100
###############
y_matrix = matrix(data = unlist(yreps100[[rep]]$y), ncol = 3, byrow = TRUE)
assign = yreps100[[rep]]$assign

data = data.frame(
  assign = assign,
  y1 = y_matrix[,1],
  y2 = y_matrix[,2],
  y3 = y_matrix[,3]
)

# p = ggplot(data = data, aes(x = y1, y = y2)) + #, label = rownames(data))) +
#   geom_point(color = assign, size = 0.75) +
#   #geom_text(size = 3, hjust = 0, nudge_x = 0.5, color = assign) +
#   # geom_text(size = 3, color = assign) +
#   theme_classic() +
#   theme(plot.title = element_text(size=18)) +
#   ggtitle("n=100")

attach(data)
p = scatterplot3d(x = y1, y = y2, z = y3, color = assign, angle = -45, 
                  pch = 20, main = "n=100", xlim = c(0,12), ylim = c(0,10), zlim = c(-4,8))
ggsave(filename = "./k53d_3blockplot_100.png", device = "png", plot = p)

pltlist[[2]] = p

######### n=300
###############
y_matrix = matrix(data = unlist(yreps300[[rep]]$y), ncol = 3, byrow = TRUE)
assign = yreps300[[rep]]$assign

data = data.frame(
  assign = assign,
  y1 = y_matrix[,1],
  y2 = y_matrix[,2],
  y3 = y_matrix[,3]
)

# p = ggplot(data = data, aes(x = y1, y = y2)) + #, label = rownames(data))) +
#   geom_point(color = assign, size = 0.75) +
#   #geom_text(size = 3, hjust = 0, nudge_x = 0.5, color = assign) +
#   # geom_text(size = 3, color = assign) +
#   theme_classic() +
#   theme(plot.title = element_text(size=18)) +
#   ggtitle("n=300")

attach(data)
p = scatterplot3d(x = y1, y = y2, z = y3, color = assign, angle = -45, 
                  pch = 20, main = "n=300")
ggsave(filename = "./k53d_3blockplot_300.png", device = "png", plot = p)

pltlist[[3]] = p

# 
# for(rep in 1:nreps){
#   
#   y_matrix = matrix(data = unlist(yreps[[rep]]$y), ncol = 2, byrow = TRUE)
#   assign = yreps[[rep]]$assign
#   
#   data = data.frame(
#     assign = assign,
#     y1 = y_matrix[,1],
#     y2 = y_matrix[,2]
#   )
#   
#   p = ggplot(data = data, aes(x = y1, y = y2)) + #, label = rownames(data))) +
#     geom_point(color = assign, size = 0.75) +
#     #geom_text(size = 3, hjust = 0, nudge_x = 0.5, color = assign) +
#     # geom_text(size = 3, color = assign) +
#     theme_classic()
#   
#   pltlist[[rep]] = p
# }

# do.call(grid.arrange, pltlist)#, top = "Simulated 2D Data and True Group Assignments")
gridplot = arrangeGrob(pltlist[[1]], pltlist[[2]], pltlist[[3]], 
                       nrow = 1, 
                       respect = TRUE) #

plotly::layout(pltlist[[1]], pltlist[[2]], pltlist[[3]])
# widths = rep(2,3), 
# heights = 2)
ggsave(filename = "./k53d_3blockplot_30.png", device = "png", plot = gridplot)
