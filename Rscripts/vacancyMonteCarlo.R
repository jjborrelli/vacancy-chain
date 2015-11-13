library(igraph)
library(compiler)
library(ggplot2)
library(data.table)


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


setup <- function(n, shellpar, mode, spatial){
  sizes <- size_distr.c(n, shellpar, mode)
  spat <- spat_distr.c(n, spatial) 
  
  return(matrix(c(spat, sizes), ncol = 3))
}

makeCHAIN <- function(intro, arena, threshold = .2, limits = c(1.01,1.5)){
  # intro is vector of new coordinates and shell size, , arena is result of setup funct
  all <- rbind(intro, arena)
  ed <- edge_distance.c(as.matrix(dist(all[,1:2])), thres = threshold)
  ed2 <- size_distance.c(all[,3], ed, limit = limits)
  g <- graph.edgelist(ed2)
  
  v <- 1
  i = 1
  condition <- TRUE
  while(condition){
    nb <- unlist(neighborhood(g, 1, v[i], mode = "in"))[-1]
    condition <- length(nb) != 0
    if(condition) v[i+1] <- nb[sample(1:length(nb), 1)]
    i <- i + 1
  }
  
  return(v)
}

sizes <- c("lnorm", "unif", "exp", "norm")
ll10 <- list()
for(k in 1:4){
  
  arena<- setup(n = 200, shellpar = c(.5,10), mode = sizes[k], spatial = "unif")
  q <- quantile(arena[,3], probs = seq(.1, 1, .1))
  coords <- matrix(runif(400), ncol = 2) 
  
  lens <- matrix(nrow = length(q), ncol = 200)
  for(j in 1:length(q)){
    chlen <- list()
    for(i in 1:200){
      ini <- c(coords[i,], q[j])
      chlen[[i]] <- makeCHAIN(ini, arena)
    }
    lens[j,] <- sapply(chlen, length)
  }
  
  ll10[[k]] <- lens
  print(k)
}

matplot(seq(.1, 1, .1), t(do.call(rbind, lapply(ll, apply, 1, median))), typ  = "l", ylab = "chainlength")
