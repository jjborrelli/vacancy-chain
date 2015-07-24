library(igraph)
library(compiler)
library(ggplot2)
library(data.table)

## Model functions


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

web_iters <- function(iter, n, sp, t, lim, shelldist = "lnorm", spatdist = "unif"){
  diff <- lim[2] - lim[1]
  res <- matrix(nrow = iter, ncol = 5)
  for(i in 1:iter){
    shellsize <- size_distr.c(n = n, shellpar = sp, mode = shelldist)
    spatial.d <- spat_distr.c(n = n, spatial = spatdist)
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

#diff1 <- pars[,4]-pars[,3]

#paths <- matrix(0, nrow = 100, ncol = 5)