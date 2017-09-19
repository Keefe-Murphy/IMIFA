#' Simulate Data from a Mixture of Factor Analysers Structure
#'
#' Function to simulate data of any size and dimension from a mixture of (infinite) factor analysers structure.
#' @param N,G,P Desired overall number of observations, number of clusters, and number of variables in the simulated data set. All must be a single integer.
#' @param Q Desired number of cluster-specific latent factors in the simulated data set. Can be specified either as a single integer if all clusters are to have the same number of factors, or a vector of length \code{G}. Defaults to \code{floor(log(P))} in each cluster.
#' @param pis Mixing proportions of the clusters in the dataset if \code{G} > 1. Must sum to 1. Defaults to \code{rep(1/G, G)}.
#' @param mu True values of the mean parameters, either as a single value, a vector of length \code{G}, a vector of length \code{P}, or a \code{G * P} matrix. If \code{mu} is missing, \code{loc.diff} is invoked to simulate distinct means for each cluster.
#' @param psi True values of uniqueness parameters, either as a single value, a vector of length \code{G}, a vector of length \code{P}, or a \code{G * P} matrix. As such the user can specify uniquenesses as a diagonal or isotropic matrix, and further constrain uniquenesses across clusters if desired. If \code{psi} is missing, uniquenesses are simulated via \code{rgamma(P, 1, 1)} within each cluster.
#' @param nn An alternative way to specify the size of each cluster, by giving the exact number of observations in each cluster explicitly. Must sum to \code{N}.
#' @param loc.diff A parameter to control the closeness of the clusters in terms of the difference in their location vectors. Only relevant if \code{mu} is NOT supplied. Defaults to 1.
#' @param method A switch indicating whether the mixture to be simulated from is the conditional distribution of the data given the latent variables (default), or simply the marginal distribution of the data.
#'
#' @return Invisibly returns a \code{data.frame} with \code{N} observations (rows) of \code{P} variables (columns). The true values of the parameters which generated these data are also stored as attributes.
#' @export
#' @importFrom Rfast "is.symmetric"
#'
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @examples
#' # Simulate 100 observations from 3 balanced clusters with cluster-specific numbers of latent factors
#' # Specify isotropic uniquenesses within each cluster
#' sim_data <- sim_IMIFA_data(N=100, G=3, P=20, Q=c(2, 2, 5), psi=1:3)
#' names(attributes(sim_data))
#' labels   <- attr(sim_data, "Labels")
#'
#' # Visualise the data in two-dimensions
#' plot(cmdscale(dist(sim_data), k=2), col=labels)
#'
#' # Fit a MIFA model to this data
#' # tmp      <- mcmc_IMIFA(sim_data, method="MIFA", range.G=3, n.iters=5000)
#' @seealso \code{\link{mcmc_IMIFA}} for fitting an IMIFA related model to the simulated data set.
sim_IMIFA_data <- function(N = 300L, G = 3L, P = 50L, Q = rep(floor(log(P)), G), pis = rep(1/G, G), mu = NULL,
                           psi = NULL, nn = NULL, loc.diff = 1, method = c("conditional", "marginal")) {

  N            <- as.integer(N)
  G            <- as.integer(G)
  P            <- as.integer(P)
  Q            <- as.integer(Q)
  Q.warn       <- min(N - 1, Ledermann(P))
  if(any(N  < 0, P  < 0, Q < 0, G  <= 0)) stop("'N', 'P', and 'Q' must be strictly non-negative and 'G' must be strictly positive")
  if(any(length(N) != 1, length(loc.diff) != 1,
         length(G) != 1, length(P) != 1)) stop("'N', 'P', 'G', and 'loc.diff' must be of length 1")
  if(any(N  < 2, N <= G))                 stop("Must simulate more than one data-point and the number of clusters cannot exceed N")
  if(any(Q  > Q.warn))                    warning(paste0("Are you sure you want to generate this many factors relative to N=", N, " and P=", P, "?\nSuggested upper-bound = ", Q.warn), call.=FALSE)
  if(length(Q) != G) {
    if(!missing(Q))  {
      if(length(Q) == 1) {
        Q      <- rep(Q, G)
      } else if(length(Q != G))           stop(paste0("'Q' must supplied for each of the G=", G, " clusters"))
    }
  }
  if(!is.character(method))               stop("'method' must be a character vector of length 1")
  method       <- match.arg(method)
  Gseq         <- seq_len(G)
  Nseq         <- seq_len(N)
  Pseq         <- seq_len(P)
  nnames       <- paste0("Obs ", Nseq)
  vnames       <- paste0("Var ", Pseq)

  if(!missing(nn) && missing(pis))     {
    nn         <- as.integer(nn)
    if(!is.integer(nn)   ||
       any(length(nn)    != G,
           sum(nn)       != N))           stop(paste0("'nn' must be an integer vector of length G=", G, " which sums to N=", N))
    if(any(nn  <= 0))                     stop("All 'nn' values must be strictly positive; simulating empty clusters not allowed")
  } else {
    if(!is.numeric(pis)  ||
       any(length(pis)   != G,
           sum(pis)      != 1))           stop(paste0("'pis' must be a numeric vector of length G=", G, " which sums to ", 1))
    if(any(pis <= 0))                     stop("All 'pis' values must be strictly positive")
    nn         <- rep(0,  G)
    iter       <- 0
    while(any(nn   < floor(N/(G * G)),
              iter < 1000))            {
      labs     <- .sim_z_p(N=N, prob.z=pis)
      nn       <- tabulate(labs, nbins=G)
      iter     <- iter    + 1
    }
  }

  mu.miss      <- missing(mu)
  if(!mu.miss)  {
    musup      <- matrix(.len_check(as.matrix(mu),  switch0g = TRUE, method = ifelse(G > 1, "MFA", "FA"), P, G)[[1]], nrow=P, ncol=G, byrow=length(mu)  == G)
  } else {
    if(!is.numeric(loc.diff))             stop("'loc.diff' must be numeric")
    musup      <- as.integer(scale(Gseq, center=TRUE, scale=FALSE)) * loc.diff
  }
  psi.miss     <- missing(psi)
  if(!psi.miss) {
    psisup     <- matrix(.len_check(as.matrix(psi), switch0g = TRUE, method = ifelse(G > 1, "MFA", "FA"), P, G)[[1]], nrow=P, ncol=G, byrow=length(psi) == G)
  }
  simdata      <- base::matrix(0, nrow=0, ncol=P)
  true.mu      <- true.l   <-
  true.psi     <- true.cov <- stats::setNames(vector("list", G), paste0("Cluster", Gseq))
  sq_mat       <- if(P > 50)  function(x) diag(sqrt(diag(x))) else sqrt

# Simulate true parameter values
  true.zlab    <- factor(rep(Gseq, nn), labels=Gseq)
  if(method    == "conditional") {
    Q.max      <- max(Q)
    eta.true   <- .sim_eta_p(Q=Q.max, N=N)
  }
  for(g in Gseq) {
    Q.g        <- Q[g]
    N.g        <- nn[g]
    mu.true    <- stats::setNames(if(mu.miss) .sim_mu_p(P=P, mu.zero=musup[g], sig.mu.sqrt=1) else musup[,g], vnames)
    l.true     <- matrix(.sim_load_p(Q=Q.g, P=P, sigma.l=1), nrow=P, ncol=Q.g)
    psi.true   <- stats::setNames(if(psi.miss) stats::rgamma(P, 1, 1) else psisup[,g], vnames)

  # Simulate data
    covmat     <- provideDimnames({ if(P > 1) diag(psi.true) else as.matrix(psi.true) } + { if(Q.g > 0) switch(method, marginal=tcrossprod(l.true), 0) else 0}, base=list(vnames))
    if(!all(is.symmetric(covmat),
            is.double(covmat)))           stop("Invalid covariance matrix!")
    covmat     <- is.posi_def(covmat, make=TRUE)$X.new
    sigma      <- if(all(Q.g > 0, P > 1, method == "marginal")) .chol(covmat) else sq_mat(covmat)
    means      <- matrix(mu.true, nrow=N.g, ncol=P, byrow=TRUE) + switch(method, conditional=tcrossprod(eta.true[true.zlab == g, seq_len(Q.g), drop=FALSE], l.true), 0)
    simdata    <- rbind(simdata, means + matrnorm(N.g, P) %*% sigma)
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
       "Clusters")     <- G
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
       "Loc.Diff")     <- if(mu.miss) loc.diff
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
