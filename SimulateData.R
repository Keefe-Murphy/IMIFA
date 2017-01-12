######################################################
### Simulate Data (Single, Shrinkage & Group Case) ###
######################################################

sim_IMIFA      <- function(N = 300L, G = 3L, P = 50L, Q = rep(4L, G), pis = rep(1/G, G), 
                           nn = NULL, loc.diff = 1L, method = c("conditional", "marginal")) {
  
  N            <- as.integer(N)
  G            <- as.integer(G)
  P            <- as.integer(P)
  Q            <- as.integer(Q)
  if(any(N  < 0, P  < 0, Q < 0, G  <= 0)) stop("'N', 'P', and 'Q' must be strictly non-negative and 'G' must be strictly positive")
  if(any(length(N) != 1, length(loc.diff) != 1,
         length(G) != 1, length(P) != 1)) stop("'N', 'P', 'G', and 'loc.diff' must be of length 1")
  if(!is.numeric(loc.diff))               stop("'loc.diff' must be numeric")
  if(any(N  < 2, N <= G))                 stop("Must simulate more than one data-point and the number of groups cannot exceed N")
  if(any(Q >= P, Q >= N - 1))             stop(paste0("Cannot generate this many factors relative to N=", N, " and P=", P))
  if(length(Q) != G) {                   
    if(!missing(Q))  {
      if(length(Q) == 1) {
        Q      <- rep(Q, G)
      } else if(length(Q != G))           stop(paste0("'Q' must supplied for each of the G=", G, " groups"))  
    }
  }
  method       <- match.arg(method)
  Gseq         <- seq_len(G)
  Nseq         <- seq_len(N)
  Pseq         <- seq_len(P)
  nnames       <- paste0("Obs ", Nseq) 
  vnames       <- paste0("Var ", Pseq)
  if(!missing(nn) && missing(pis))     {
    nn         <- as.integer(nn)
    if(any(nn  == 0))                     stop("All 'nn' values must be strictly positive; simulating empty groups not allowed")
    if(any(length(nn)  != G, 
           sum(nn)     != N,
           !is.integer(nn)))              stop(paste0("'nn' must be an integer vector of length G=", G, " which sums to N=", N))
  } else {
    if(any(length(pis) != G,
           sum(pis)    != 1,
           !is.numeric(pis)))             stop(paste0("'pis' must be a numeric vector of length G=", G, " which sums to ", 1))
    nn         <- rep(0,  G)
    while(any(nn < floor(N/(G * G))))  {
      labs     <- .sim_z.p(N=N, prob.z=pis)
      nn       <- tabulate(labs, nbins=G)  
    }
  } 
                                         
  simdata      <- matrix(0, nr=0, nc=P)
  prior.mu     <- as.integer(scale(Gseq, center=TRUE, scale=FALSE))
  true.mu      <- stats::setNames(vector("list", G), paste0("Group", Gseq))
  true.l       <- true.mu
  true.psi     <- true.mu
  true.cov     <- true.mu
  
# Simulate true parameter values
  true.zlab    <- factor(rep(Gseq, nn), labels=Gseq)
  if(method    == "conditional") {
    Q.max      <- max(Q)
    eta.true   <- .sim_eta.p(Q=Q.max, N=N)
  }
  for(g in Gseq) {
    Q.g        <- Q[g]
    N.g        <- nn[g]
    mu.true    <- stats::setNames(.sim_mu.p(P=P, mu.zero=prior.mu[g] * loc.diff, sigma.mu=1), vnames)
    l.true     <- .sim_load.p(Q=Q.g, P=P, sigma.l=1)
    psi.true   <- stats::setNames(stats::runif(P, 0, 1), vnames)
    
  # Simulate data
    covmat     <- provideDimnames(diag(psi.true) + switch(method, marginal=tcrossprod(l.true), 0), base=list(vnames, vnames))
    if(!all(isSymmetric(covmat),
            is.double(covmat)))           stop("Invalid covariance matrix")
    if(!corpcor::is.positive.definite(covmat)) {
      covmat   <- corpcor::make.positive.definite(covmat)
    }
    sigma      <- if(any(Q.g > 0, method == "conditional")) .chol(covmat) else sqrt(covmat)
    means      <- matrix(mu.true, nr=N.g, nc=P, byrow=TRUE) + switch(method, conditional=tcrossprod(eta.true[true.zlab == g, seq_len(Q.g), drop=FALSE], l.true), 0)
    simdata    <- rbind(simdata, means + matrix(stats::rnorm(N.g * P), nr=N.g, nc=P) %*% sigma)
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
       "Means")        <- do.call(cbind, true.mu)
  if(method == "conditional") {
    eta.true   <- eta.true[permute,, drop=FALSE]
    dimnames(eta.true) <- list(nnames, if(Q.max > 0) paste0("Factor ", seq_len(Q.max)))
    attr(simdata,
         "Scores")     <- eta.true
  }
  attr(simdata, 
       "Loadings")     <- true.l
  attr(simdata,
       "Loc.Diff")     <- loc.diff
  attr(simdata, 
       "Uniquenesses") <- do.call(cbind, true.psi)
  attr(simdata, 
       "Labels")       <- true.zlab
  attr(simdata,
       "Covariance")   <- true.cov
  class(simdata)       <- c("data.frame", "IMIFA")
  return(simdata)
}