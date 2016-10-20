###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.IMIFA      <- function(N = 300, G = 3, P = 50, Q = rep(4, G), pis = rep(1/G, G), 
                           nn = NULL, method = c("conditional", "marginal")) {
  
  if(any(N  < 0, P  < 0, Q < 0, G <= 0)) stop("'N', 'P', and 'Q' must be strictly non-negative and 'G' must be strictly positive")
  if(any(N  < 2, N <= G))                stop("Must simulate more than one data-point and the number of groups cannot exceed N")
  if(any(Q >= P, Q >= N - 1))            stop(paste0("Cannot generate this many factors relative to N=", N, " and P=", P))
  if(length(Q) != G) {                   
    if(!missing(Q))  {
      if(length(Q) == 1) {
        Q      <- rep(Q, G)
      } else if(length(Q != G))          stop(paste0("'Q' must supplied for each of the G=", G, " groups"))  
    }
  }
  method       <- match.arg(method)
  Gseq         <- seq_len(G)
  Nseq         <- seq_len(N)
  Pseq         <- seq_len(P)
  nnames       <- paste0("Obs ", Nseq) 
  vnames       <- paste0("Var ", Pseq)
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
      labs     <- factor(which(rmultinom(N, size=1, prob=pis) != 0, arr.ind=TRUE)[,1], levels=Gseq)
      nn       <- tabulate(labs, nbins=G)  
    }
  } 
                                         
  simdata      <- matrix(0, nr=0, nc=P)
  nn.seq       <- seq_along(nn)
  true.mu      <- setNames(vector("list", G), paste0("Group", Gseq))
  true.l       <- true.mu
  true.psi     <- true.mu
  true.cov     <- true.mu
  
# Simulate true parameter values
  true.zlab    <- factor(rep(nn.seq, nn), labels=nn.seq)
  if(method == "conditional") {
    Q.max      <- max(Q)
    eta.true   <- matrix(rnorm(N * Q.max, 0, 1), nr=N, nc=Q.max) 
    dimnames(eta.true) <- list(nnames, paste0("Factor ", seq_len(Q.max)))
  }
  for(g in Gseq) {
    Q.g        <- Q[g]
    N.g        <- nn[g]
    mu.true    <- rnorm(P, rnorm(P, sample(-5:5, P, replace=TRUE), 1), 2)
    l.true     <- if(Q.g > 0) matrix(rnorm(P * Q.g, 0, 1), nr=P, nc=Q.g) else base::matrix(, nr=P, nc=0)
    psi.true   <- runif(P, 0, 0.5)
    
  # Simulate data
    covmat     <- diag(psi.true) + (if(method == "marginal") tcrossprod(l.true) else 0)
    if(!isSymmetric(covmat))             stop("Non-symmetric covariance matrix")
    sigma      <- if(any(Q.g > 0, method == "conditional")) chol(covmat) else sqrt(covmat)
    means      <- matrix(mu.true, nr=N.g, nc=P, byrow=TRUE) + (if(method == "conditional") tcrossprod(eta.true[true.zlab == g, seq_len(Q.g), drop=FALSE], l.true) else 0)
    simdata    <- rbind(simdata, means + matrix(rnorm(N.g * P), nr=N.g, nc=P) %*% sigma)
    names(mu.true)     <- vnames
    names(psi.true)    <- vnames
    dimnames(covmat)   <- list(vnames, vnames)
    dimnames(l.true)   <- list(vnames, if(Q.g > 0) paste0("Factor ", seq_len(Q.g)))
    true.mu[[g]]       <- mu.true
    true.l[[g]]        <- l.true
    true.psi[[g]]      <- psi.true
    true.cov[[g]]      <- covmat
  }
  
# Post-process data
  permute      <- sample(Nseq, N, replace=FALSE)
  simdata      <- simdata[permute,, drop=FALSE]
  true.zlab    <- true.zlab[permute]
  dimnames(simdata)    <- list(nnames, vnames)
  simdata      <- as.data.frame(simdata)
  attr(simdata,
       "Factors")      <- Q
  attr(simdata,
       "Groups")       <- G
  attr(simdata, 
       "Means")        <- true.mu
  if(method == "conditional") {
    eta.true   <- eta.true[permute,, drop=FALSE]
    attr(simdata,
         "Scores")     <- eta.true
  }
  attr(simdata, 
       "Loadings")     <- true.l
  attr(simdata, 
       "Uniquenesses") <- true.psi
  attr(simdata, 
       "Labels")       <- true.zlab
  attr(simdata,
       "Covariance")   <- true.cov
  class(simdata)       <- c("data.frame", "IMIFA")
  return(simdata)
}