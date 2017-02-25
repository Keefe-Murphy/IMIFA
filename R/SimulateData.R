#' Simulating Data from a Mixture of Factor Analysers Structure
#'
#' Function to simulate data of any size and dimension from a mixture of (infinite) factor analysers structure.
#' @param N Desired overall number of observations in the simulated data set - a single integer.
#' @param G Desired number of clusters in the simulated data set - a single integer.
#' @param P Desired number of variables in the simulated dataset - a single integer.
#' @param Q Desired number of cluster-specific latent factors in the simulated data set. Can be specified either as a single integer if all clusters are to have the same number of factors, or a vector of length \code{G}.
#' @param pis Mixing proportions of the clusters in the dataset if \code{G} > 1. Must sum to 1. Defaults to \code{rep(1/G, G)}.
#' @param nn An alternative way to specify the size of each cluster, by giving the exact number of observations in each group explicitly. Must sum to \code{N}.
#' @param loc.diff A parameter to control the closeness of the clusters in terms of the difference in their location vectors. Defaults to 1.
#' @param method A switch indicating whether the mixture to be simulated from is the conditional distribution of the data given the latent variables (default), or simply the marginal distribution of the data.
#'
#' @return Invisibly returns a data.frame with \code{N} observations (rows) of \code{P} variables (columns). The true values of the parameters which generated these data are also stored.
#' @export
#' @importFrom corpcor "is.positive.definite" "make.positive.definite"
#'
#' @examples
#' # Simulate 100 observations from 3 balanced groups with cluster-specific numbers of latent factors
#' # sim_data <- sim_IMIFA_data(N=100, G=3, P=20, Q=c(2, 2, 5))
#' # names(attributes(sim_data))
#' # attr(sim_data, "Labels")
#' # tmp      <- mcmc_IMIFA(sim_data, method="MIFA", range.G=3, n.iters=5000)
#' @seealso The function \code{\link{mcmc_IMIFA}} for fitting an IMIFA related model to the simulated data set.
sim_IMIFA_data <- function(N = 300L, G = 3L, P = 50L, Q = rep(4L, G), pis = rep(1/G, G),
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
  if(any(Q  > .ledermann(N=N, P=P)))      stop(paste0("Cannot generate this many factors relative to N=", N, " and P=", P))
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
      labs     <- .sim_z_p(N=N, prob.z=pis)
      nn       <- tabulate(labs, nbins=G)
    }
  }

  simdata      <- matrix(0, nrow=0, ncol=P)
  prior.mu     <- as.integer(scale(Gseq, center=TRUE, scale=FALSE))
  true.mu      <- setNames(vector("list", G), paste0("Group", Gseq))
  true.l       <- true.mu
  true.psi     <- true.mu
  true.cov     <- true.mu

# Simulate true parameter values
  true.zlab    <- factor(rep(Gseq, nn), labels=Gseq)
  if(method    == "conditional") {
    Q.max      <- max(Q)
    eta.true   <- .sim_eta_p(Q=Q.max, N=N)
  }
  for(g in Gseq) {
    Q.g        <- Q[g]
    N.g        <- nn[g]
    mu.true    <- setNames(.sim_mu_p(P=P, mu.zero=prior.mu[g] * loc.diff, sig.mu.sqrt=1), vnames)
    l.true     <- .sim_load_p(Q=Q.g, P=P, sigma.l=1)
    psi.true   <- setNames(rgamma(P, 1, 1), vnames)

  # Simulate data
    covmat     <- provideDimnames(diag(psi.true) + switch(method, marginal=tcrossprod(l.true), 0), base=list(vnames))
    if(!all(isSymmetric(covmat),
            is.double(covmat)))           stop("Invalid covariance matrix")
    if(!is.positive.definite(covmat)) {
      covmat   <- make.positive.definite(covmat)
    }
    sigma      <- if(any(Q.g > 0, method == "conditional")) .chol(covmat) else sqrt(covmat)
    means      <- matrix(mu.true, nrow=N.g, ncol=P, byrow=TRUE) + switch(method, conditional=tcrossprod(eta.true[true.zlab == g, seq_len(Q.g), drop=FALSE], l.true), 0)
    simdata    <- rbind(simdata, means + matrix(rnorm(N.g * P), nrow=N.g, ncol=P) %*% sigma)
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
       "Proportions")  <- pis
  attr(simdata,
       "Covariance")   <- true.cov
  class(simdata)       <- c("data.frame")
    invisible(simdata)
}
