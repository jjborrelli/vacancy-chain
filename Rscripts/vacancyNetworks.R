library(igraph)
library(compiler)
library(ggplot2)
library(data.table)
library(reshape2)

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
  plh <- path.length.hist(g, directed = T)
  pl.sd <- sd(rep(1:length(plh), each = plh))
  return(c(d, apl, pl.sd))
}

graph.props.c <- cmpfun(graph.props)

web_iters <- function(iter, n, sp, t, lim, shelldist = "lnorm", spatdist = "unif"){
  diff <- lim[2] - lim[1]
  res <- matrix(nrow = iter, ncol = 6)
  for(i in 1:iter){
    shellsize <- size_distr.c(n = n, shellpar = sp, mode = shelldist)
    spatial.d <- spat_distr.c(n = n, spatial = spatdist)
    edge.d <- edge_distance.c(d = as.matrix(dist(as.data.frame(spatial.d))), thres = t)
    edges <- size_distance.c(shell = shellsize, edge = edge.d, limit = lim)
    
    res[i,] <- c(graph.props.c(edges), n, t, diff)
  }
  colnames(res) <- c("diam", "avpath", "pathSD", "N", "Th", "diff")
  return(res)
}

# number of individuals
n <- seq(50, 350, 100)
# shell parameters of lognormal distr
spar <- c(.5, 1)
# threshold
t <- seq(.1, 1, .1)
# min/max shell size for swapping
lim <- matrix(c(rep(1,10), seq(1.1, 3, .2)), nrow = 10, ncol = 2)

#pars <- expand.grid(n, t, lim[,1], lim[,2])
pars <- rbind(expand.grid(n, t, 1.05, 1.25), expand.grid(n, t, 1.25, 1.45), expand.grid(n, t, 1.45, 1.65))

diff1 <- pars[,4]-pars[,3]


require(doSNOW)
require(parallel)
require(data.table)

#make the cluster
cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)

system.time(
RESULT3 <- foreach(i = 1:120) %dopar% {
  source("./Rscripts/vacancyNetScript.R")
  paths <- web_iters(iter = 100, n = pars[i,1], sp = spar, t = pars[i,2], lim = c(pars[i, 3], pars[i,4]))
  write(i, file = "C:/Users/jjborrelli/Dropbox/vacancy-runs.txt", append = T)
  return(as.data.frame(paths))
}
)

stopCluster(cl)

allDAT.3 <- cbind(rbindlist(RESULT3), win = rep(c("A", "B", "C"), each = 4000))

allDAT.4 <- rbind(cbind(allDAT.2, winSize = ".1"), cbind(allDAT.3, winSize = ".2"))


ggplot(allDAT, aes(x = Th, y = avpath)) + geom_point() + facet_grid(diff~N)
ggplot(allDAT.2, aes(x = factor(Th), y = avpath, fill = win)) + geom_boxplot() + facet_grid(N~., scales = "free_y") + theme_bw() + xlab("Threshold") + ylab("Average Chain Length")

plotDAT <- allDAT[which(allDAT$N == 50 | allDAT$N == 250 | allDAT$N == 500),]
plotDAT <- plotDAT[which(plotDAT$diff == "0.1" | plotDAT$diff == "0.3" | plotDAT$diff == "0.7" | plotDAT$diff == "1.3" | plotDAT$diff == "1.9"),]
ggplot(plotDAT , aes(x = factor(Th), y = avpath)) + geom_boxplot() + facet_grid(diff~N, scales = "free_y") + theme_bw() + xlab("Threshold") + ylab("Average Chain Length")

ggplot(allDAT, aes(x = Th, y = diff, fill = avpath)) + geom_point() + facet_grid(.~N)

all2 <- aggregate(allDAT$avpath, by = list(allDAT$N, allDAT$Th, allDAT$diff), mean)
ggplot(all2, aes(x = Group.2, y = Group.3, col = x)) + geom_point(size = 3) + facet_grid(.~Group.1) + xlab("Threshold") + ylab("diff")

############################################################################################
#save.image(file = "C:/Users/jjborrelli/Desktop/vacChain2.Rdata")
#load("./Data/vacChain.Rdata")
#load("C:/Users/jjborrelli/Desktop/vacChain.Rdata")

############################################################################################
#
# Size distributions : lnorm , unif , exp , norm
# Spatial distributions : unif , normal , lognormal

sizes <- c("lnorm", "unif", "exp", "norm")
spatial <- c("unif", "normal", "lognormal")


# number of individuals
n <- seq(100, 300, 100)
# shell parameters of lognormal distr
spar <- list(c(.5, 1), c(0, 1), c(NA, NA), c(.5, 1))
spar2 <- rep(spar, 3)
# threshold
t <- seq(.1, 1, .2)
# min/max shell size for swapping
lim <- matrix(c(rep(1,3), seq(1.2, 2, .4)), nrow = 3, ncol = 2)

pars <- cbind(expand.grid(n, t, lim[,2]), 1)

diff1 <- pars[,4]-pars[,3]

distros <- as.matrix(expand.grid(sizes, spatial))

require(doSNOW)
require(parallel)
require(data.table)

#make the cluster
cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)

allDAT2 <- list()
for(k in 1:nrow(distros)){
  RESULT <- foreach(i = 1:nrow(pars)) %dopar% {
    source("./Rscripts/vacancyNetScript.R")
    paths <- web_iters(iter = 100, n = pars[i,1], sp = spar2[[k]], t = pars[i,2], lim = c(pars[i, 4], pars[i,3]),
                       shelldist = distros[k,1], spatdist = distros[k,2])
    write(paste(k, i, sep = "-"), file = "C:/Users/jjborrelli/Dropbox/vacancy-runs.txt", append = T)
    return(as.data.frame(paths))
  }
  print(k)
  allDAT2[[k]] <- rbindlist(RESULT)
}

stopCluster(cl)

#allDAT <- rbindlist(RESULT)

test <- allDAT2[[1]]
for(i in 1:12){
  allDAT2[[i]] <- cbind(allDAT2[[i]], shell = distros[i,1], spat = distros[i,2])
}

rblAD <- rbindlist(allDAT2)

ggplot(rbindlist(allDAT2), aes(x = Th, y = avpath, col = factor(N), shape = factor(diff), alpha = 0.5)) + geom_point() + geom_jitter() + facet_grid(shell ~ spat)

agDAT <- aggregate(rblAD$avpath, list(rblAD$N, rblAD$Th, rblAD$diff, rblAD$shell, rblAD$spat), mean)
agDAT2 <- aggregate(rblAD$pathSD, list(rblAD$N, rblAD$Th, rblAD$diff, rblAD$shell, rblAD$spat), mean)

ggplot(agDAT, aes(x = Group.2, y = x, col = factor(Group.1), shape = factor(Group.3))) + geom_point(size = 3) + scale_color_discrete(name="N") + scale_shape_discrete(name = "Diff") + facet_grid(Group.4~Group.5) + xlab("Threshold") + ylab("Chain Length") + theme_bw()

ggplot(agDAT2, aes(x = Group.2, y = x, col = factor(Group.1), shape = factor(Group.3))) + geom_point() + scale_color_discrete(name="N") + scale_shape_discrete(name = "Diff") + facet_grid(Group.4~Group.5) + xlab("Threshold") + ylab("Chain StDev") + theme_bw()



# density test

tm <- matrix(c(runif(500, 0, 1), runif(500, 0, 1)), nrow = 500, ncol = 2)

test1 <- c()
for(i in 1:100){
  r1 <- runif(2, 0, 1)
  test1[i] <- sum((tm[,1] - r1[1])^2 + (tm[,2] - r1[2])^2 <= .1)
}

#### the average number of points within the .1 threshold is ~120

tm <- matrix(c(runif(250, 0, 1), runif(250, 0, 1)), nrow = 250, ncol = 2)

test1 <- c()
for(i in 1:100){
  r1 <- runif(2, 0, 1)
  test1[i] <- sum((tm[,1] - r1[1])^2 + (tm[,2] - r1[2])^2 <= .4)
}

mean(test1)
