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
  if(!missing(nn)) {
    nn         <- as.integer(nn)
    if(any(length(nn) != G, 
           sum(nn)    !=N))              stop(paste0("nn must be an integer vector of length ", G, " which sums to ", N))
  } else {
    labs       <- factor(which(rmultinom(N, size=1, prob=rdirichlet(1, rep(0.5, G))) != 0, arr.ind=TRUE)[,1], levels=seq_len(G))
    nn         <- tabulate(labs, nbins=G)
  } 

  if(sum(nn == 0) > 0) {
    nn.ind     <- which(nn == 0)
    nn.max     <- which.max(nn)
    nn[nn.ind] <- nn[nn.ind] + 1
    nn[nn.max] <- nn[nn.max] - length(nn.ind)
  }                           
                                              
  vnames        <- paste0("Var ", seq_len(P))
  simdata       <- matrix(0, nr=0, nc=P)
  true.mu       <- setNames(vector("list", G), paste0("Group", seq_len(G)))
  true.l        <- true.mu
  true.psi      <- true.mu
  true.cov      <- true.mu
  
# Simulate true parameter values
  for(g in seq_len(G)) {
    Q.g         <- Q[g]
    mu.true     <- as.vector(rnorm(P, 0, max(Q.g, 1)) + sqrt(abs(rnorm(P, Q.g, g)) * diag(P)) %*% rnorm(P, 0, 1))
    
    if(Q.g > 0) {
      l.true    <- rnorm(Q.g, 0, Q.g/2) + t(sqrt(abs(rnorm(Q.g, 0, g/2)) * diag(Q.g)) %*% matrix(rnorm(P * Q.g, 0, 1), nr=Q.g, ncol=P))
    } else {
      l.true    <- matrix(, nr=P, nc=0)
    }
  
    shape.g     <- floor((g * max(Q.g, 2))/2)
    psi.true    <- 1/rgamma(n=P, shape=shape.g, rate=shape.g - 1)                            
    
  # Simulate data
    covmat      <- tcrossprod(l.true, l.true) + diag(psi.true)
    sigma       <- if(Q.g > 0) chol(covmat) else sqrt(covmat)
    simdata     <- rbind(simdata, t(rep(mu.true, nn[g]) + sigma %*% matrix(rnorm(P * nn[g], 0, 1), nr=P, ncol=nn[g])))
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
  true.zlab     <- factor(rep(nn.seq, nn), labels=nn.seq)
  permute       <- sample(seq_len(N), N, replace=FALSE)
  simdata       <- simdata[permute,]
  true.zlab     <- true.zlab[permute]
  dimnames(simdata)     <- list(paste0("Obs ", seq_len(N)), vnames)
  simdata       <- as.data.frame(simdata)
  attr(simdata,
       "Factors")       <- Q
  attr(simdata,
       "Groups")        <- G
  attr(simdata, 
       "Means")         <- true.mu
  attr(simdata, 
       "Loadings")      <- true.l
  attr(simdata, 
       "Uniquenesses")  <- true.psi
  attr(simdata, 
       "Labels")        <- true.zlab
  attr(simdata,
       "Covariance")    <- true.cov
  class(simdata)        <- c("data.frame", "IMIFA")
  return(simdata)
}