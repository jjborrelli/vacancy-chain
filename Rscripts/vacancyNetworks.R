library(igraph)
library(compiler)
library(ggplot2)
library(data.table)

# R cookbook function to plot multiple ggplot objects in one window
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Function to plot a spatially explicit network given inputs of an edgelist, xy positions of nodes, and node sizes
plot_net <- function(edges, pos, shell){
  require(ggplot2)
  
  x1 <- c()
  y1 <- c()
  x2 <- c()
  y2 <- c()
  for(i in 1:nrow(edges)){
    x1[i] <- pos$x[edges[i,1]]
    y1[i] <- pos$y[edges[i,1]]
    x2[i] <- pos$x[edges[i,2]]
    y2[i] <- pos$y[edges[i,2]]
  }
  e <- data.frame(x1, y1, x2, y2)
  
  
  p <- ggplot(pos, aes(x = x, y = y)) 
  p <- p + geom_segment(data = e, aes(x = x1, y = y1, xend = x2, yend = y2), alpha = .75)
  p <- p + geom_point(aes(size = factor(shell)), col = "darkgreen") 
  # the rest is just eliminating the background
  p <- p + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) 
  p <- p + theme(panel.background = element_blank()) + theme(legend.position="none")
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
  p <- p + theme( legend.background = element_rect(colour = NA)) 
  p <- p + theme(panel.background = element_rect(fill = "white", colour = NA)) 
  p <- p + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  print(p)
}


############################################################################################
#########                             ######################################################
#########       MODEL FUNCTIONS       ######################################################
#########                             ######################################################
############################################################################################

# Check to make sure the parameters are input correctly
test_params <- function(shellpar, limit){
  if(!length(shellpar) == 2){
    return(stop("shellpar must be a vector of length 2"))
  }
  if(!length(limit) == 2){
    return(stop("limit must be a vector of length 2"))
  }
}

# Create the size distribution of the nodes
size_distr <- function(n, shellpar, mode){
  # Shell size distribution
  if(mode == "lnorm"){
    shell <- rlnorm(n, shellpar[1], shellpar[2])
  }else if(mode == "unif"){
    shell <- runif(n, shellpar[1], shellpar[2])
  }else if(mode == "exp"){
    shell <- rexp(n)
  }else if(mode == "norm"){
    shell <- abs(rnorm(n, shellpar[1], shellpar[2]))
  }
  return(shell)
}

size_distr.c <- cmpfun(size_distr)

# Create the random spatial distribution of the nodes
spat_distr <- function(n, spatial){
  if(spatial == "unif"){
    #where is each crab?
    x <- runif(n)
    y <- runif(n)
  }
  if(spatial == "normal"){
    x <- rnorm(n)
    y = rnorm(n)
  }
  if(spatial == "lognormal"){
    x <- rlnorm(n)
    y <- rlnorm(n)
  }
  return(matrix(c(x,y), ncol = 2))
}

spat_distr.c <- cmpfun(spat_distr)

# Check if there is an edge based on threshold distance
edge_distance <-function(d, thres){
  #are they within visual range?
  edge <- which(d < thres, arr.ind = T)
  return(edge)
}

edge_distance.c <- cmpfun(edge_distance)

# Limit edges by size differences
size_distance <- function(shell, edge, limit){
  lower.lim <- limit[1] # how much bigger the shell must be
  upper.lim <- limit[2] # how much bigger the shell can be
  
  min.size <- shell[edge[,1]]*lower.lim
  max.size <- shell[edge[,1]]*upper.lim
  
  cond <- apply(cbind(shell[edge[,2]], min.size, max.size), 1, function(x){return(x[1] > x[2] & x[1] <= x[3])})
  
  return(matrix(edge[cond,], ncol = 2))
}

size_distance.c <- cmpfun(size_distance)

# Plot web
plot_web <- function(edge, shell, pos){
  plot(graph.edgelist(edge), vertex.size = sqrt(shell), layout = pos,
       vertex.label = NA, margin = c(0, 0, 0, 0), edge.arrow.size = .5)
}

graph.props <- function(ed){
  if(nrow(ed) == 0){return(c(0,0))}
  
  g <- graph.edgelist(ed)
  d <- diameter(g)
  apl <- average.path.length(g)
  return(c(d, apl))
}

graph.props.c <- cmpfun(graph.props)

web_iters <- function(iter, n, sp, t, lim){
  diff <- lim[2] - lim[1]
  res <- matrix(nrow = iter, ncol = 5)
  for(i in 1:iter){
    shellsize <- size_distr.c(n = n, shellpar = sp, mode = "lnorm")
    spatial.d <- spat_distr.c(n = n, spatial = "unif")
    edge.d <- edge_distance.c(d = as.matrix(dist(as.data.frame(spatial.d))), thres = t)
    edges <- size_distance.c(shell = shellsize, edge = edge.d, limit = lim)
    
    res[i,] <- c(graph.props.c(edges), n, t, diff)
  }
  colnames(res) <- c("diam", "avpath", "N", "Th", "diff")
  return(res)
}

# number of individuals
n <- seq(50, 500, 50)
# shell parameters of lognormal distr
spar <- c(.5, 1)
# threshold
t <- seq(.1, 1, .1)
# min/max shell size for swapping
lim <- matrix(c(rep(1,10), seq(1.1, 3, .2)), nrow = 10, ncol = 2)

pars <- expand.grid(n, t, lim[,1], lim[,2])

diff1 <- pars[,4]-pars[,3]


require(doSNOW)
require(parallel)
require(data.table)

#make the cluster
cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)

system.time(
RESULT <- foreach(i = 1:10000) %dopar% {
  source("./Rscripts/vacancyNetScript.R")
  paths <- web_iters(iter = 100, n = pars[i,1], sp = spar, t = pars[i,2], lim = c(pars[i, 3], pars[i,4]))
  write(i, file = "C:/Users/jjborrelli/Dropbox/vacancy-runs.txt", append = T)
  return(as.data.frame(paths))
}
)

stopCluster(cl)

allDAT <- rbindlist(RESULT)

ggplot(allDAT, aes(x = Th, y = avpath)) + geom_point() + facet_grid(diff~N)
ggplot(allDAT, aes(x = factor(Th), y = avpath)) + geom_boxplot() + facet_grid(diff~N, scales = "free_y") + theme_bw() + xlab("Threshold") + ylab("Average Chain Length")

ggplot(allDAT[which(allDAT$N == 50 | allDAT$N == 250 | allDAT$N == 500),], aes(x = factor(Th), y = avpath, fill = factor(diff))) + geom_boxplot() + facet_grid(diff~N) + theme_bw() + xlab("Threshold") + ylab("Average Chain Length")

ggplot(allDAT, aes(x = Th, y = diff, fill = avpath)) + geom_point() + facet_grid(.~N)

all2 <- aggregate(allDAT$avpath, by = list(allDAT$N, allDAT$Th, allDAT$diff), mean)
ggplot(all2, aes(x = Group.2, y = Group.3, col = x)) + geom_point(size = 3) + facet_grid(.~Group.1) + xlab("Threshold") + ylab("diff")

############################################################################################
#save.image(file = "C:/Users/jjborrelli/Desktop/vacChain.Rdata")
#load("./Data/vacChain.Rdata")
