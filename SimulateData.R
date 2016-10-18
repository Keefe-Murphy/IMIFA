###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.IMIFA      <- function(N = 300, G = 3, P = 50, Q = rep(4, G), pis = rep(1/G, G), nn = NULL) {
  
  if(any(N  < 0, P  < 0, Q < 0, G <= 0)) stop("'N', 'P', and 'Q' must be strictly non-negative and 'G' must be strictly positive")
  if(any(Q >= P, Q >= N - 1))            stop(paste0("Cannot generate this many factors relative to N=", N, " and P=", P))
  if(length(Q) != G) {                   
    if(!missing(Q))  {
      if(length(Q) == 1) {
        Q      <- rep(Q, G)
      } else if(length(Q != G))          stop(paste0("'Q' must supplied for each of the G=", G, " groups"))  
    }
  }
  if(!missing(nn) && missing(pi.prop)) {
    nn         <- as.integer(nn)
    if(any(nn  == 0))                    stop("All 'nn' values be strictly positive; simulating empty groups not allowed")
    if(any(length(nn)  != G, 
           sum(nn)     != N))            stop(paste0("'nn' must be an integer vector of length G=", G, " which sums to N=", N))
  } else {
    if(any(length(pis) != G,
           sum(pis)    != 1))            stop(paste0("'pis' must be a vector of length G=", G, " which sums to ", 1))
    nn         <- rep(0,  G)
    while(any(nn < floor(N/(G * G))))  {
      labs     <- factor(which(rmultinom(N, size=1, prob=pis) != 0, arr.ind=TRUE)[,1], levels=seq_len(G))
      nn       <- tabulate(labs, nbins=G)  
    }
  } 
                                              
  vnames        <- paste0("Var ", seq_len(P))
  simdata       <- matrix(0, nr=0, nc=P)
  true.mu       <- setNames(vector("list", G), paste0("Group", seq_len(G)))
  true.l        <- true.mu
  true.psi      <- true.mu
  true.cov      <- true.mu
  
# Simulate true parameter values
  for(g in seq_len(G))   {
    Q.g         <- Q[g]
    mu.true     <- rnorm(P, rnorm(P, sample(-5:5, P, replace=TRUE), 1), 2)
    if(Q.g > 0) {
      l.true    <- matrix(rnorm(P * Q.g, 0, 1), nr=P, nc=Q.g)
    } else {
      l.true    <- matrix(, nr=P, nc=0)
    }
    psi.true    <- runif(P, 0, 0.5)
    
  # Simulate data
    covmat      <- tcrossprod(l.true, l.true) + diag(psi.true)
    if(!isSymmetric(covmat))             stop("Non-symmetric covariance matrix")
    sigma       <- if(Q.g > 0) chol(covmat) else sqrt(covmat)
    simdata     <- rbind(simdata, t(rep(mu.true, nn[g]) + sigma %*% matrix(rnorm(P * nn[g], 0, 1), nr=P, nc=nn[g])))
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