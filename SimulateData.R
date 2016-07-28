###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.IMIFA      <- function(N = 500, G = 3, P = 25, Q = rep(5, G), nn = NULL) {
  
  if(any(N  < 0, P  < 0, Q < 0, G <= 0)) stop("N, P, and Q must be strictly non-negative and G must be strictly positive")
  if(any(Q >= P, Q >= N - 1))            stop(paste0("Cannot generate this many factors relative to N=", N, " and P=", P))
  if(length(Q) != G)                    
    if(!missing(Q))  {
      if(length(Q) == 1)  {
        Q      <- rep(Q, G)
      } else if(length(Q != G))  {       stop(paste0("Q must supplied for each of the G=", G, " groups"))  
      }                                       
    }
  if(missing(nn)) {
    pie        <- rdirichlet(1, rep(0.5, G))
    ind.mat    <- rmultinom(N, size=1, prob=pie)
    labs       <- factor(which(ind.mat != 0, arr.ind=TRUE)[,1], levels=seq_len(G))
    nn         <- tabulate(labs, nbins=G)
  } else if(any(length(nn) != G, 
                sum(nn) != N))   {       stop(paste0("nn must be a vector of length ", G, " which sums to ", N))
  }
  if(sum(nn == 0) > 0) {
    nn.ind     <- which(nn == 0)
    n.max      <- which.max(nn)
    n.empty    <- length(nn.ind)
    nn[nn.ind] <- nn[nn.ind] + 1
    nn[n.max]  <- nn[n.max] - n.empty
  }                           
                                              
  gnames        <- paste0("Group", seq_len(G))
  vnames        <- paste0("Var ", seq_len(P))
  SimData       <- matrix(0, nr=0, nc=P)
  true.mu       <- setNames(vector("list", G), gnames)
  true.l        <- setNames(vector("list", G), gnames)
  true.psi      <- setNames(vector("list", G), gnames)
  true.cov      <- setNames(vector("list", G), gnames)
  
# Simulate true parameter values
  for(g in seq_len(G)) {
    Q.g         <- Q[g]
    U.mu        <- sqrt(abs(rnorm(P, Q.g, g)) * diag(P))
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- U.mu %*% z.mu
    mu.mu       <- rnorm(P, 0, max(Q.g, 1))
    mu.true     <- as.vector(mu.mu + v.mu)
    
    if(Q.g > 0) {
      U.load    <- sqrt(abs(rnorm(Q.g, 0, g/2)) * diag(Q.g))
      z.load    <- matrix(rnorm(P * Q.g, 0, 1), nr=Q.g, ncol=P)
      v.load    <- t(U.load %*% z.load)
      mu.load   <- rnorm(Q.g, 0, Q.g/2)
      l.true    <- mu.load + v.load
    } else {
      l.true    <- matrix(, nr=P, nc=0)
    }
  
    shape.g     <- floor((g * max(Q.g, 2))/2)
    psi.true    <- 1/rgamma(n=P, shape=shape.g, rate=shape.g - 1)                            
    
  # Simulate data
    covmat      <- tcrossprod(l.true, l.true) + diag(psi.true)
    if(Q.g > 0) {
      U.om      <- chol(covmat)
    } else {
      U.om      <- sqrt(covmat)
    }
    z           <- matrix(rnorm(P * nn[g], 0, 1), nr=P, ncol=nn[g])
    var         <- U.om %*% z
    SimData     <- rbind(SimData, t(rep(mu.true, nn[g]) + var))
    names(mu.true)      <- vnames
    names(psi.true)     <- vnames
    dimnames(covmat)    <- list(vnames, vnames)
    dimnames(l.true)    <- list(vnames, if(Q.g > 0) paste0("Factor ", seq_len(Q.g)))
    true.mu[[g]]        <- mu.true
    true.l[[g]]         <- l.true
    true.psi[[g]]       <- psi.true
    true.cov[[g]]       <- covmat
  }
  
# Post-process data
  nn.seq        <- seq_along(nn)
  true.zlab     <- factor(rep(nn.seq, nn), labels=seq_along(nn.seq))
  permute       <- sample(seq_len(N), N, replace=FALSE)
  SimData       <- SimData[permute,]
  true.zlab     <- true.zlab[permute]
  dimnames(SimData)     <- list(paste0("Obs ", seq_len(N)), vnames)
  SimData       <- as.data.frame(SimData)
  attr(SimData,
       "Factors")       <- Q
  attr(SimData,
       "Groups")        <- G
  attr(SimData, 
       "Means")         <- true.mu
  attr(SimData, 
       "Loadings")      <- true.l
  attr(SimData, 
       "Uniquenesses")  <- true.psi
  attr(SimData, 
       "Labels")        <- true.zlab
  attr(SimData,
       "Covariance")    <- true.cov
  class(SimData)        <- c("data.frame", "IMIFA")
  return(SimData)
}