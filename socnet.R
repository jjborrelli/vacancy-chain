library(ggplot2)
library(grid)
library(igraph)
library(igraph)
library(NetIndices)
library(data.table)

#plotting functions
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


makeNET <- function(n, shellpar, thres = 0.2, limit = 1.5, mode = "lnorm", plot = F){
  if(!length(shellpar) == 2){
    stop("shellpar must be a vector of length 2")
  }
  # how big is each crab's shell
  if(mode == "lnorm"){
    shell <- rlnorm(n, shellpar[1], shellpar[2])
  }else if(mode == "unif"){
    shell <- runif(n, shellpar[1], shellpar[2])
  }else if(mode == "exp"){
    shell <- rexp(n)
  }
  
  #where is each crab?
  x <- runif(n)
  y <- runif(n)
  
  #compute distance between each crab
  d <- as.matrix(dist(data.frame(x, y)))
  
  #how far can they see
  thres <- thres
  
  #are they within visual range?
  edge <- matrix(0, ncol = 2)
  for(i in 1:n){
    for(j in 1:n){
      if(d[i,j] <= thres){adj <- matrix(c(i,j), nrow = 1, ncol = 2)}else{next}
      edge <- rbind(edge, adj)
    }
  }
  edge <- edge[-1,]
  
  #is it a bigger shell
  lim <- limit # how much bigger the shell can be
  rem <- c()
  for(i in 1:nrow(edge)){
    if(shell[edge[i,1]]*lim[1] > shell[edge[i,2]]){
      rem <- c(rem, i)
    }else if(edge[i,1] == edge[i,2]){
      rem <- c(rem, i)
    }else if(shell[edge[i,2]] > shell[edge[i,1]]*lim[2] ){
      rem <- c(rem, i)
    }
  } 
  edge <- edge[-rem,]
  
  w <- c()
  for(i in 1:nrow(edge)){
    w[i] <- d[edge[i,1], edge[i,2]]
  }
  
  edges <- cbind(edge, w)
  
  if(plot){
    plot_net(edges = edge, pos = data.frame(x, y), shell = shell)
  
  }
  return(list(edge = edges, shell = shell))
}

meanPATH <- function(edge){
  require(igraph)
  require(NetIndices)
  if(ncol(edge) > 2){edge <- edge[,1:2]}
  g <- graph.edgelist(edge)
  gmat <- get.adjacency(g, sparse = F)
  tind <- TrophInd(gmat)
  apl <- average.path.length(g)
  return(list(avpath = apl, TL = tind))
}



densities <- seq(50, 1000, 50)
avpth <- c()
tpos <- c()
system.time(
for(i in length(densities):1){
  net <- makeNET(densities[i], shellpar = c(0, 1), .2, mode = "unif")
  mpth <- meanPATH(net$edge)
  avpth[i] <- unlist(mpth["avpath"])
  tpos[[i]] <- cbind(mpth$TL, shell = net$shell)
  print(densities[i])
}
)


plot(unlist(avpth)~densities)

d1 <- data.frame(x = densities, y = avpth, typ = "density")

thres <- seq(.05, .9, .025)
avpth2 <- c()
tpos2 <- c()
for(i in 1:length(thres)){
  net <- makeNET(200, shellpar = c(0, 1), thres = thres[i], mode = "unif")
  mpth <- meanPATH(net$edge)
  avpth2[i] <- unlist(mpth["avpath"])
  tpos2[[i]] <- cbind(mpth$TL, shell = net$shell)
  cat(i, thres[i], "\n")
}

plot(avpth2~thres)

d2 <- data.frame(x = thres, y = avpth2, typ = "threshold")

mu <- rep(seq(-10, 10, .5), each = 3)
avpth3 <- c()
tpos3 <- list()
for(i in 108:length(mu)){
  net <- makeNET(200, shellpar = c(mu[i], 1), .2, mode = "unif")
  mpth <- meanPATH(net$edge)
  avpth3[i] <- unlist(mpth["avpath"])
  tpos3[[i]] <- cbind(mpth$TL, shell = net$shell)
  cat(i, mu[i], "\n")
}

plot(avpth3~mu)

d3 <- data.frame(x = mu, y = avpth3, typ = "meansize")

lim <- rep(seq(1.1, 10, .1), each = 3)
avpth4 <- c()
tpos4 <- list()
for(i in 1:length(lim)){
  net <- makeNET(200, shellpar = c(0, 1), .2, limit = c(1.1, lim[i]))
  mpth <- meanPATH(net$edge)
  avpth4[i] <- unlist(mpth["avpath"])
  #tpos4[[i]] <- cbind(mpth$TL, shell = net$shell)
  cat(i, lim[i], "\n")
}

plot(avpth4~lim)

d4 <- data.frame(x = lim, y = avpth4, typ = "swapsizelimit")


#d.all <- rbind(d1,d2,d3,d4)
dim(d.all)

require(ggplot2)
g1 <- ggplot(d1, aes(x,y)) + geom_point() + xlab("Density") + ylab("Average Path")
g2 <- ggplot(d2, aes(x,y)) + geom_point() + xlab("Distance Threshold") + ylab("Average Path")
g3 <- ggplot(d3, aes(x,y)) + geom_point() + xlab("Mean Size") + ylab("Average Path")
g4 <- ggplot(d4, aes(x,y)) + geom_point() + xlab("Max Swap Size") + ylab("Average Path")


multiplot(g1, g2, g3, g4, cols = 2)

#dataTAB1 <- cbind(rbindlist(tpos),factors = rep(densities, sapply(tpos, nrow)), type = rep("density", 10500), mode = rep("lnorm", 10500))

#dataTAB2 <- cbind(rbindlist(tpos2),factors = rep(thres, sapply(tpos2, nrow)), type = rep("density", 7000), mode = rep("lnorm", 7000))

#dataTAB3 <- cbind(rbindlist(tpos3),factors = rep(mu, sapply(tpos3, nrow)), type = rep("density", 24600), mode = rep("lnorm", 24600))

#dataTAB4 <- cbind(rbindlist(tpos4),factors = rep(lim, sapply(tpos4, nrow)), type = rep("density", 54000), mode = rep("lnorm", 54000))

dataTAB5 <- cbind(rbindlist(tpos),factors = rep(densities, sapply(tpos, nrow)), type = rep("density", 10500), mode = rep("unif", 10500))

dataTAB6 <- cbind(rbindlist(tpos2),factors = rep(thres, sapply(tpos2, nrow)), type = rep("density", 7000), mode = rep("unif", 7000))

dataTAB7 <- cbind(rbindlist(tpos3),factors = rep(mu, sapply(tpos3, nrow)), type = rep("density", 24600), mode = rep("unif", 24600))

dataTAB8 <- cbind(rbindlist(tpos4),factors = rep(lim, sapply(tpos4, nrow)), type = rep("density", 54000), mode = rep("unif", 54000))

dataTABLE <- rbind(dataTAB1, dataTAB2, dataTAB3, dataTAB4)

###
### Varying all params

dens <- c(100, 200, 300, 400, 500)
thres <- c(.1, .2, .3, .6, .8)
sizes <- c(-5, -2.5, 0, 2.5, 5)
swaps <- c(1.5, 2.5, 5, 7.5, 10)

params <- expand.grid(dens, thres, sizes, swaps)

s.avpth <- c()
for(i in 1:nrow(params)){
  net <- makeNET(n = params[i,1], mean.shell = params[i,3], sd.shell = 1,
                 thres = params[i,2], limit = params[i,4])
  s.avpth[i] <- unlist(meanPATH(net)["avpath"])
  cat(i, "of", nrow(params), "\n")
}

plot(s.avpth[params$Var2 == 0.1]~params$Var1[params$Var2 == 0.1])

setwd("C:/Users/borre_000/Documents/")
save.image("socnet.RData")
load("socnet.Rdata")

t <- c()
for(i in 1:ncol(gmat)){
  for(j in 1:nrow(gmat)){
    t[i] <- gmat[i,j] == gmat[j,i]
  }
}
sum(t)
