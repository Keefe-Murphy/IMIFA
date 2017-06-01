#' Adaptive Gibbs Sampler for Nonparameteric Model-based Clustering using models from the IMIFA family
#'
#' Carries out Gibbs sampling for all models from the IMIFA family, facilitating model-based clustering with dimensionally reduced factor-analytic covariance structures, with automatic estimation of the number of clusters and cluster-specific factors as appropriate to the method employed. Factor analysis with one group (FA/IFA), finite mixtures (MFA/MIFA), overfitted mixtures (OMFA/OMIFA), infinite factor models which employ the multiplicative gamma process (MGP) shrinkage prior (IFA/MIFA/OMIFA/IMIFA), and infinite mixtures which employ Dirichlet Process Mixture Models (IMFA/IMIFA) are all provided. Creates a raw object of class 'IMIFA' from which the optimal/modal model can be extracted by \code{\link{get_IMIFA_results}}.
#'
#' @param dat A matrix or data frame such that rows correspond to observations (\code{N}) and columns correspond to variables (\code{P}). Non-numeric variables and rows with missing entries will be removed.
#' @param method An acronym for the type of model to fit where: \cr
#' \cr
#'  "\code{FA}" = Factor Analysis \cr
#'  "\code{IFA}" = Infinite Factor Analysis \cr
#'  "\code{MFA}" = Mixtures of Factor Analysers \cr
#'  "\code{MIFA}" = Mixtures of Infinite Factor Analysers \cr
#'  "\code{OMFA}" = Overfitted Mixtures of Factor Analysers \cr
#'  "\code{OMIFA}" = Overfitted Mixtures of Infinite Factor Analysers \cr
#'  "\code{IMFA}" = Infinite Mixtures of Factor Analysers \cr
#'  "\code{IMIFA}" = Infinite Mixtures of Infinite Factor Analysers \cr
#'  \cr
#'  The "\code{classify}" method is not yet implemented.
#' @param n.iters The number of iterations to run the Gibbs sampler for.
#' @param range.G Depending on the method employed, either the range of values for the number of clusters, or the conseratively high starting value for the number of clusters. Defaults to 1 for the "\code{FA}" and "\code{IFA}" methods. For the "\code{MFA}" and "\code{MIFA}" models this is to be given as a range of candidate models to explore. For the "\code{OMFA}", "\code{OMIFA}", "\code{IMFA}", and "\code{IMIFA}" models, this is the number of clusters with which the chain is to be initialised, in which case the default is \code{min(N - 1, max(25, ceiling(3 * log(N))))}. For the "\code{OMFA}", and "\code{OMIFA}" models this upper limit remains fixed for the entire length of the chain; \code{range.G} also doubles as the default \code{trunc.G} for the "\code{IMFA}" and "\code{IMIFA}" models. However, when \code{N < P}, or when this bound is close to or exceeds \code{N} for any of these overfitted/infinite mixture models, it is better to initialise at a value closer to the truth (i.e. \code{ceiling(log(N))} by default), though the upper bound remains the same - as a result the role of \code{range.G} when \code{N < P} is no longer to specify the upper bound (which can still be modified via \code{trunc.G}, at least for the "\code{IMFA}" and "\code{IMIFA}" methods) and the number of groups used for initialisation, but rather just the number of groups used for initialisation only. If \code{length(range.G) * length(range.Q)} is large, consider not storing unnecessary parameters, or breaking up the range of models to be explored into chunks, and sending each chunk to \code{\link{get_IMIFA_results}}.
#' @param range.Q Depending on the method employed, either the range of values for the number of latent factors, or, for methods ending in IFA the conservatively high starting value for the number of cluster-specific factors, in which case the default starting value is \code{floor(3 * log(P))}. For methods ending in IFA, different clusters can be modelled using different numbers of latent factors (incl. zero); for methods not ending in IFA it is possible to fit zero-factor models, corresponding to simple diagonal covariance structures. For instance, fitting the "\code{IMFA}" model with \code{range.Q=0} corresponds to a vanilla Dirichlet Process Mixture Model. If \code{length(range.G) * length(range.Q)} is large, consider not storing unnecessary parameters or breaking up the range of models to be explored into chunks, and sending each chunk to \code{\link{get_IMIFA_results}}.
#' @param burnin The number of burn-in iterations for the sampler. Defaults to \code{n.iters/5}. Note that chains can also be burned in later, using \code{\link{get_IMIFA_results}}.
#' @param thinning The thinning interval used in the simulation. Defaults to 2. No thinning corresponds to 1. Note that chains can also be thinned later, using \code{\link{get_IMIFA_results}}.
#' @param centering A logical value indicating whether mean centering should be applied to the data, defaulting to \code{TRUE}.
#' @param scaling The scaling to be applied - one of "\code{unit}", "\code{none}" or "\code{pareto}".
#' @param mu.zero The mean of the prior distribution for the mean parameter. Defaults to the sample mean of the data.
#' @param sigma.mu The covariance of the prior distribution for the mean parameter. Can be a scalar times the identity or a matrix of appropriate dimension. Defaults to the sample covariance matrix.
#' @param sigma.l The covariance of the prior distribution for the loadings. Defaults to 1. Only relevant for the finite factor methods.
#' @param alpha Depending on the method employed, either the hyperparameter of the Dirichlet prior for the cluster mixing proportions, or the Dirichlet process concentration parameter. Defaults to 0.5/range.G for the Overfitted methods - if supplied for "\code{OMFA}" and "\code{OMIFA}" methods, you are supplying the numerator of \code{alpha/range.G}, which should be less than half the dimension (per group!) of the free parameters of the smallest model considered in order to ensure superfluous clusters are emptied (for "\code{OMFA}", this corresponds to the smallest \code{range.Q}; for "\code{OMIFA}", this corresponds to a zero-factor model) [see: \code{\link{PGMM_dfree}} and Rousseau and Mengersen (2011)]. Defaults to 1 for the finite mixture models "\code{MFA}" and "\code{MIFA}". Defaults to \code{1 - discount} for the "\code{IMFA}" and "\code{IMIFA}" models if \code{learn.alpha=FALSE} or a simulation from the prior if \code{learn.alpha=TRUE}. Must be positive, unless \code{discount} is supplied for the "\code{IMFA}" or "\code{IMIFA}" methods.
#' @param psi.alpha The shape of the inverse gamma prior on the uniquenesses. Defaults to 2.5.
#' @param psi.beta The rate of the inverse gamma prior on the uniquenesses. Can be either a single parameter or a vector of variable specific rates.  If this is not supplied, \code{\link{psi_hyper}} is invoked to choose sensible values, depending on the value of \code{uni.prior}.
#' @param uni.type A switch indicating whether uniquenesses are to be "\code{unconstrained}" or "\code{isotropic}". Note that "\code{unconstrained}" here means variable-specific and group-specific, whereas "\code{isotropic}" here means isotropic but still group-specific. The "\code{isotropic}" constraint provides the link between factor analysis and the probabilistic principal component analysis model. Defaults to "\code{unconstrained}", but "\code{isotropic}" is recommended when \code{N < P}.
#' @param uni.prior A switch indicating whether uniquenesses rate hyperparameters are to be "\code{unconstrained}" or "\code{isotropic}". "\code{uni.prior}" must be "\code{isotropic}" if "\code{uni.type}" is "\code{isotropic}", but can take either value when "\code{uni.type}" is "\code{unconstrained}". Defaults to \code{uni.type} if that is supplied and \code{uni.prior} is not, otherwise defaults to "\code{unconstrained}", but "\code{isotropic}" is recommended when \code{N < P}. Only relevant when "\code{psi.beta}" is not supplied and \code{\link{psi_hyper}} is invoked.
#' @param z.init The method used to initialise the cluster labels. Defaults to \code{\link[mclust]{Mclust}}. Not relevant for the "\code{FA}" and "\code{"IFA"} methods.
#' @param z.list A user supplied list of cluster labels. Only relevant if \code{z.init == "z.list"}.
#' @param adapt A logical value indicating whether adaptation of the number of cluster-specific factors is to take place. Only relevant for methods ending in IFA, in which case the default is \code{TRUE}. Specifying \code{FALSE} and supplying \code{range.Q} provides a means to use the MGP prior in a finite factor context.
#' @param prop Proportion of elements within the neighbourhood \code{epsilon} of zero necessary to consider a loadings column redundant. Defaults to \code{floor(0.7 * P)/P}. Only relevant for methods ending in IFA.
#' @param epsilon Neighbourhood of zero within which a loadings entry is considered negligible according to \code{prop}. Defaults to 0.1. Only relevant for methods ending in IFA.
#' @param alpha.d1 Shape hyperparameter of the global shrinkage on the first column of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to 3. Only relevant for methods ending in IFA.
#' @param alpha.d2 Shape hyperparameter of the global shrinkage on subsequent columns of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to 6. Only relevant for methods ending in IFA.
#' @param beta.d1 Rate hyperparameter of the global shrinkage on the first column of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to 1. Only relevant for methods ending in IFA.
#' @param beta.d2 Rate hyperparameter of the global shrinkage on the first column of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to 1. Only relevant for methods ending in IFA.
#' @param nu Hyperparameter for the gamma prior on the local shrinkage parameters. Defaults to 2. Passed to \code{\link{MGP_check}} to ensure validity. Only relevant for methods ending in IFA.
#' @param nuplus1 Logical switch indicating whether the shape hyperparameter of the prior on the local shrinkage parameters is equal to \code{nu + 1}. If \code{FALSE}, it is simply equal to \code{nu}. Only relevant for methods ending in IFA.
#' @param adapt.at The iteration at which adaptation is to begin. Defaults to \code{burnin} for the "\code{IFA}" and "\code{MIFA}" methods, defaults to 0 for the "\code{OMIFA}" and "\code{IMIFA}". Cannot exceed \code{burnin}. Only relevant for methods ending in IFA.
#' @param b0 Intercept parameter for the exponentially decaying adaptation probability s.t. \code{p(iter) = 1/exp(b0 + b1 * (iter - adapt.at))}. Defaults to 0.1. Only relevant for methods ending in IFA.
#' @param b1 Slope parameter for the exponentially decaying adaptation probability s.t. \code{p(iter) = 1/exp(b0 + b1 * (iter - adapt.at))}. Defaults to 0.00005. Only relevant for methods ending in IFA.
#' @param trunc.G The maximum number of allowable and storable groups if the "\code{IMFA}" or "\code{IMIFA}" method is employed. Defaults to the same value as \code{range.G} (unless \code{N < P}, see \code{range.G} for details) and must be greater than or equal to this value. The number of active groups to be sampled at each iteration is adaptively truncated, with \code{trunc.G} as an upper limit for storage reasons. Note that large values of \code{trunc.G} may lead to memory capacity issues.
#' @param learn.alpha Logical indicating whether the Dirichlet process / Pitman concentration parameter is to be learned, or remain fixed for the duration of the chain. If being learned, a Ga(a, b) prior is assumed for \code{alpha}; updates take place via Gibbs sampling when \code{discount} is zero and via Metropolis-Hastings otherwise. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods, in which case the default is \code{TRUE}.
#' @param alpha.hyper A vector of length 2 giving hyperparameters for the Dirichlet process / Pitman-Yor concentration parameter \code{alpha}. If \code{isTRUE(learn.alpha)}, these are shape and rate parameter of a Gamma distribution. Defaults to Ga(2, 1). Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods, in which case the default is \code{TRUE}. The prior is shifted to have support on (-\code{discount}, \code{Inf}) when non-zero \code{discount} is supplied or \code{learn.d=TRUE}.
#' @param zeta Tuning parameter controlling the acceptance rate of the random-walk proposal for the Metropolis-Hastings steps when \code{learn.alpha=TRUE}. These steps are only invoked when either \code{discount} is non-zero or \code{learn.d=TRUE}, otherwise \code{alpha} is learned by Gibbs updates. Must be strictly positive. Defauts to 2.
#' @param ind.slice Logical indicitating whether the independent slice-efficient sampler is to be employed. If \code{FALSE} the dependent slice-efficient sampler is employed, whereby the slice sequence xi_1,...,xi_g is equal to the decreasingly ordered mixing proportions. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods. Defaults to \code{TRUE}.
#' @param rho Parameter controlling the rate of geometric decay for the independent slice-efficient sampler, s.t. xi = (1 - rho)rho^(g-1). Must lie in the interval (0, 1]. Higher values are associated with better mixing but longer run times. Defaults to 0.75, but 0.5 is an interesting special case which guarantees that the slice sequence xi_1,...,xi_g is equal to the \emph{expectation} of the decreasingly ordered mixing proportions. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods when \code{ind.slice} is \code{TRUE}.
#' @param IM.lab.sw Logial indicating whether the two forced label switching moves are to be implemented (defaults to \code{TRUE}) when running one of the infinite mixture models, with Dirichlet process or Pitman-Yor process priors. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param verbose Logical indicating whether to print output (e.g. run times) and a progress bar to the screen while the sampler runs. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param discount The discount parameter used when generalising the Dirichlet process to the Pitman-Yor process. Must lie in the interval [0, 1). If non-zero, \code{alpha} can be supplied greater than -discount. Defaults to 0. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param learn.d Logical indicating whether the \code{discount} parameter is to be updated via Metropolis-Hastings. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods, in which case the default is \code{FALSE}.
#' @param d.hyper Hyperparameters for the Beta(a,b) prior on the \code{discount} hyperparameter. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param kappa The prior distribution on the \code{discount} hyperparameter is assumed to be a mixture with point-mass at zero and a continuous Beta(a,b) distribution. \code{kappa} gives the weight of the point mass at zero. Must lie in the interval [0,1]. Defaults to 0.5. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param mu0g Logical indicating whether the \code{mu.zero} hyperparameter can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the "\code{MFA}" and "\code{MIFA}" methods when \code{z.list} is supplied.
#' @param psi0g Logical indicating whether the \code{psi.beta} hyperparameter(s) can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the "\code{MFA}" and "\code{MIFA}" methods when \code{z.list} is supplied.
#' @param delta0g Logical indicating whether the \code{alpha.d1}  and \code{alpha.d2} hyperparameters can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the "\code{MFA}" and "\code{MIFA}" methods when \code{z.list} is supplied.
#' @param mu.switch Logical indicating whether the means are to be stored (defaults to \code{TRUE}). May be useful not to store if memory is an issue. Warning: posterior inference won't be posssible.
#' @param score.switch Logical indicating whether the factor scores are to be stored. As the array containing each sampled scores matrix tends to be amongst the largest objects to be stored, this defaults to \code{FALSE} when \code{length(range.G) * length(range.Q) > 10}, otherwise the default is \code{TRUE}. May be useful not to store if memory is an issue - for the "\code{MIFA}", "\code{OMIFA}", and "\code{IMIFA}" methods, setting this switch to \code{FALSE} also offers a slight speed-up. Warning: posterior inference won't be posssible.
#' @param load.switch Logical indicating whether the factor loadings are to be stored (defaults to \code{TRUE}). May be useful not to store if memory is an issue. Warning: posterior inference won't be posssible.
#' @param psi.switch Logical indicating whether the uniquenesses are to be stored (defaults to \code{TRUE}). May be useful not to store if memory is an issue. Warning: posterior inference won't be posssible.
#' @param pi.switch Logical indicating whether the mixing proportions are to be stored (defaults to \code{TRUE}). May be useful not to store if memory is an issue. Warning: posterior inference won't be posssible.
#'
#' @return A list of lists of lists of class "IMIFA" to be passed to \code{\link{get_IMIFA_results}}. If the returned object is x, candidate models accesible via subsetting, where x is of the form x[[1:length(range.G)]][[1:length(range.Q)]]. However, these objects of class "IMIFA" should rarely if ever be manipulated by hand - use of the \code{\link{get_IMIFA_results}} function is \emph{strongly} advised. Dedicated \code{print} and \code{summary} functions exist for objects of class "\code{IMIFA}".
#' @export
#' @import stats
#' @importFrom utils "capture.output" "head" "setTxtProgressBar" "tail" "txtProgressBar"
#' @importFrom matrixStats "rowLogSumExps"
#' @importFrom Rfast "rowsums" "Order" "colVars" "rowmeans" "standardise" "sort_unique" "cora" "cova"
#' @importFrom e1071 "matchClasses"
#' @importFrom mvnfast "dmvn"
#' @importFrom slam "as.simple_sparse_array" "as.simple_triplet_matrix"
#' @importFrom mclust "Mclust" "mclustBIC"
#' @importFrom utils "memory.limit"
#'
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link{psi_hyper}}, \code{\link{MGP_check}}
#' @references
#' Murphy, K., Gormley, I. C. and Viroli, C. (2017) Infinite Mixtures of Infinite Factor Analysers: Nonparametric Model-Based Clustering via Latent Gaussian Models, \href{https://arxiv.org/abs/1701.07010}{arXiv:1701.07010}.
#'
#' Bhattacharya, A. and Dunson, D. B. (2011) Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291-306.
#'
#' Kalli, M., Griffin, J. E. and Walker, S. G. (2011) Slice sampling mixture models, \emph{Statistics and Computing}, 21(1): 93-105.
#'
#' Rousseau, J. and Mengersen, K. (2011) Asymptotic Behaviour of the posterior distribution in overfitted mixture models, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 73(5): 689-710.
#'
#' Tipping, M. E. and Bishop, C. M. (1999). Probabilistic principal component analysis, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 61(3): 611-622.
#'
#' @author Keefe Murphy
#'
#' @examples
#' # data(olive)
#' # data(coffee)
#'
#' # Fit an IMIFA model to the olive data. Accept all defaults.
#' # simIMIFA <- mcmc_IMIFA(olive, method="IMIFA")
#' # summary(simIMIFA)
#'
#' # Fit an IMIFA model assuming a Pitman-Yor prior, allowing the discount parameter to be learned.
#' # simPY    <- mcmc_IMIFA(olive, method="IMIFA", learn.d=TRUE)
#' # summary(simPY)
#'
#' # Fit a MFA model to the scaled olive data, with isotropic uniquenesses (i.e. MPPCA).
#' # Allow diagonal covariance as a special case where range.Q = 0. Accept all other defaults.
#' # simMFA   <- mcmc_IMIFA(olive, method="MFA", n.iters=10000, range.G=3:6,
#' #                        range.Q=0:3, centering=FALSE, uni.type="isotropic")
#'
#' # Fit a MIFA model to the centered & scaled coffee data, w/ cluster labels initialised by K-Means.
#' # Note that range.Q doesn't need to be specified. Allow IFA as a special case where range.G=1.
#' # simMIFA  <- mcmc_IMIFA(coffee, method="MIFA", n.iters=10000, range.G=1:3, z.init="kmeans")
#'
#' # Fit an IFA model to the centered and pareto scaled olive data.
#' # Note that range.G doesn't need to be specified. We can optionally supply a range.Q starting value.
#' # We can also enforce additional shrinkage using alpha.d1, alpha.d2, prop, and epsilon.
#' # simIFA   <- mcmc_IMIFA(olive, method="IFA", n.iters=10000, range.Q=4,
#' #                        alpha.d1=3.5, alpha.d2=7, prop=0.6, epsilon=0.12)
#'
#' # Fit an OMIFA model to the centered & scaled coffee data.
#' # Supply a sufficiently small alpha value. Try varying other hyperparameters.
#' # Accept the default value for the starting number of factors,
#' # but supply a value for the starting number of clusters.
#' # simOMIFA <- mcmc_IMIFA(coffee, method="OMIFA", range.G=10, psi.alpha=3, nu=3, alpha=0.8)
mcmc_IMIFA  <- function(dat = NULL, method = c("IMIFA", "IMFA", "OMIFA", "OMFA", "MIFA", "MFA", "IFA", "FA", "classify"), n.iters = 25000L, range.G = NULL, range.Q = NULL, burnin = n.iters/5,
                        thinning = 2L, centering = TRUE, scaling = c("unit", "pareto", "none"), mu.zero = NULL, sigma.mu = NULL, sigma.l = NULL, alpha = NULL, psi.alpha = NULL, psi.beta = NULL,
                        uni.type = c("unconstrained", "isotropic"), uni.prior = c("unconstrained", "isotropic"), z.init = c("mclust", "kmeans", "list", "priors"), z.list = NULL, adapt = TRUE,
                        prop = NULL, epsilon = NULL, alpha.d1 = NULL, alpha.d2 = NULL, beta.d1 = NULL, beta.d2 = NULL, nu = NULL, nuplus1 = TRUE, adapt.at = NULL, b0 = NULL, b1 = NULL,
                        trunc.G = NULL, learn.alpha = TRUE, alpha.hyper = NULL, zeta = NULL, ind.slice = TRUE, rho = NULL, IM.lab.sw = TRUE, verbose = interactive(), discount = NULL, learn.d = FALSE,
                        d.hyper = NULL, kappa = NULL, mu0g = FALSE, psi0g = FALSE, delta0g = FALSE, mu.switch = TRUE, score.switch = TRUE, load.switch = TRUE, psi.switch = TRUE, pi.switch = TRUE) {

  call      <- match.call()
  defopt    <- options()
  options(warn=1)
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
  if(!missing(method) && method == "classification") {
    method  <- "classify"
  }
  method    <- match.arg(method)
  scaling   <- match.arg(scaling)
  if(missing(dat))                  stop("Dataset must be supplied")
  dat.nam   <- gsub("[[:space:]]", "", deparse(substitute(dat)))
  nam.dat   <- gsub("\\[.*", "", dat.nam)
  pattern   <- c("(", ")")
  nam.x     <- gsub(".*\\[(.*)\\].*", "\\1)", dat.nam)
  if(!exists(nam.dat,
     envir=.GlobalEnv))             stop(paste0("Object ", match.call()$dat, " not found\n"))
  if(any(unlist(vapply(seq_along(pattern), function(p) grepl(pattern[p], nam.dat, fixed=TRUE), logical(1L))),
         !identical(dat.nam, nam.dat) && (any(grepl("[[:alpha:]]", gsub('c', '',  nam.x))) || grepl(":",
         nam.x,    fixed=TRUE))))   stop("Extremely inadvisable to supply 'dat' subsetted by any means other than row/column numbers or c() indexing: best to create new data object")
  zin.miss  <- missing(z.init)
  zli.miss  <- missing(z.list)
  if(!zli.miss) {
    z.nam   <- gsub("[[:space:]]", "", deparse(substitute(z.list)))
    nam.z   <- gsub("\\[.*", "", z.nam)
    nam.zx  <- gsub(".*\\[(.*)\\].*", "\\1)",   z.nam)
    if(!exists(nam.z,
               envir=.GlobalEnv))   stop(paste0("Object ", match.call()$z.list, " not found\n"))
    if(any(unlist(vapply(seq_along(pattern), function(p) grepl(pattern[p], nam.z, fixed=TRUE), logical(1L))),
           !identical(z.nam,   nam.z) && (any(grepl("[[:alpha:]]", gsub('c', '', nam.zx))) || grepl(":",
           nam.zx, fixed=TRUE))))   stop("Extremely inadvisable to supply 'z.list' subsetted by any means other than row/column numbers or c() indexing: best to create new object")
    if(!is.list(z.list))     z.list        <- lapply(list(z.list), as.factor)
    if(zin.miss &&
        z.init  != "list") { z.init        <- "list"
      if(verbose)                   message("'z.init' set to 'list' as 'z.list' was supplied")
    }
  }
  if(any(!is.logical(centering),
         length(centering) != 1))   stop("'centering' must be TRUE or FALSE")
  if(any(!is.logical(nuplus1),
         length(nuplus1)   != 1))   stop("'nuplus1' must be TRUE or FALSE")
  if(any(!is.logical(verbose),
         length(verbose)   != 1))   stop("'verbose' must be TRUE or FALSE")

# Remove non-numeric columns & apply centering & scaling if necessary
  burnin    <- as.integer(burnin)
  thinning  <- as.integer(thinning)
  n.iters   <- as.integer(n.iters)
  if(any(!is.integer(burnin),   burnin   < 0,
         length(burnin)   != 1))    stop("'burnin' must be a single integer")
  if(any(!is.integer(thinning), thinning < 1,
         length(thinning) != 1))    stop("'thinning' must be a single integer")
  if(any(!is.integer(n.iters),
         length(n.iters)  != 1))    stop("'n.iters' must be a single integer")
  n.iters   <- max(burnin + 1, as.integer(n.iters))
  iters     <- seq(from=burnin + 1, to=n.iters, by=thinning)
  iters     <- iters[iters > 0]
  raw.dat   <- as.data.frame(dat)
  num.check <- vapply(raw.dat, is.numeric, logical(1L))
  if(sum(num.check) != ncol(dat)) {
    if(verbose)                     message("Non-numeric columns removed")
    raw.dat <- raw.dat[num.check]
  }
  if(length(iters)  <= 1)           stop("Run a longer chain!")
  if(anyNA(raw.dat)) {
    if(verbose)                     message("Rows with missing values removed from data")
    raw.dat <- raw.dat[complete.cases(raw.dat),]
  }
  if(method != "classify") {
    scal    <- switch(scaling, none=FALSE, Rfast::colVars(as.matrix(raw.dat), std=TRUE))
    scal    <- switch(scaling, pareto=sqrt(scal), scal)
    dat     <- if(is.logical(scal)) standardise(as.matrix(raw.dat), center=centering, scale=scal) else scale(raw.dat, center=centering, scale=scal)
  } else   {
    dat     <- raw.dat
  }
  centered  <- switch(method, classify=all(round(colSums(dat)) == 0), any(centering, all(round(colSums(dat)) == 0)))
  N         <- as.integer(nrow(dat))
  P         <- as.integer(ncol(dat))
  lnN       <- log(N)
  NlP       <- N < P
  miss.uni  <- missing("uni.type")
  if(miss.uni) uni.type      <- "unconstrained"
  if(missing("uni.prior")) {
    uni.prior      <- uni.type
  }
  uni.type  <- match.arg(uni.type)
  uni.prior <- match.arg(uni.prior)
  if(all(uni.prior == "unconstrained",
         uni.type  == "isotropic")) stop("'uni.prior' can only be 'unconstrained' when 'uni.type' is 'unconstrained'")
  if(all(uni.prior == "unconstrained",
         uni.type  == "unconstrained",
         NlP, miss.uni, verbose))   message("Consider setting 'uni.type', or at least 'uni.prior', to 'isotropic' in N << P cases")

# Manage storage switches & warnings for other function inputs
  if(!missing(mu.switch)  && all(!mu.switch, ifelse(method == "classify",
     !centering, !centered)))       warning("Centering hasn't been applied - are you sure you want mu.switch=FALSE?", call.=FALSE)
  score.x   <- missing(score.switch)
  switches  <- c(mu.sw=mu.switch, s.sw=score.switch, l.sw=load.switch, psi.sw=psi.switch, pi.sw=pi.switch)
  if(any(length(switches) != 5,
         !is.logical(switches)))    stop("All logical parameter storage switches must be TRUE or FALSE")
  if(N < 2)                         stop("Must have more than one observation")
  G.x       <- missing(range.G)
  alpha.x   <- missing(learn.alpha)
  disc.x    <- missing(learn.d)
  if(any(!is.logical(learn.alpha),
         length(learn.alpha) != 1)) stop("'learn.alpha' must be TRUE or FALSE")
  if(any(!is.logical(learn.d),
         length(learn.d)     != 1)) stop("'learn.d' must be TRUE or FALSE")
  if(missing(d.hyper))       d.hyper       <- c(1L, 1L)
  if(length(d.hyper)         != 2)  stop("d.hyper' must be a vector of length 2")
  if(any(d.hyper   <= 0))           stop("'Discount Beta prior hyperparameters must be strictly positive")
  if(missing(kappa)) {
    kappa          <- 0.5
  }
  kappa            <- ifelse(all(!learn.d, discount == 0), 1, kappa)
  if(any(!is.numeric(kappa),
         length(kappa)       != 1)) stop("'kappa' must be a single number")
  if(kappa   <  0  || kappa   > 1)  stop("'kappa' must lie in the interval [0, 1]")
  if(kappa  ==  0) {
   if(all(!learn.d, discount == 0)) stop("'kappa' is zero and yet 'discount' is fixed at zero:\n either learn the discount parameter or specify a non-zero value")
  } else if(kappa  == 1) {
   if(any(learn.d,  discount != 0)) stop(paste0("'kappa' is exactly 1 and yet", ifelse(learn.d, " 'discount' is being learned ", if(discount != 0) " the discount is fixed at a non-zero value"), ":\n the discount should remain fixed at zero"))
  }
  discount         <- switch(method, IMFA=, IMIFA=ifelse(missing(discount), ifelse(learn.d, ifelse(kappa != 0 && runif(1) <= kappa, 0, rbeta(1, d.hyper[1], d.hyper[2])), 0), discount), 0)
  if(any(!is.numeric(discount),
         length(discount)  != 1))   stop("'discount' must be a single number")
  if(discount       < 0    ||
     discount      >= 1)            stop("'discount' must lie in the interval [0, 1)")
  if(!is.element(method, c("IMFA", "IMIFA"))) {
    if(learn.alpha) {
      learn.alpha  <- FALSE
      if(verbose   && !alpha.x)     message(paste0("'learn.alpha' forced to FALSE for the ", method, " method"))
    }
    if(learn.d)     {
      learn.d      <- FALSE
      if(verbose   && !disc.x)      message(paste0("'learn.d' must be FALSE for the ", method, " method"))
    }
  }
  if(!is.element(method, c("MFA", "MIFA")))      {
    if(length(range.G) > 1)         stop(paste0("Only one 'range.G' value can be specified for the ", method, " method"))
    if(all(!G.x, is.element(method, c("FA", "IFA"))) &&
       range.G  > 1)                warning(paste0("'range.G' must be 1 for the ", method, " method"), call.=FALSE)
    if(is.element(method, c("OMIFA", "OMFA", "IMFA", "IMIFA"))) {
      lnN2         <- ceiling(lnN)
      tmp.G        <- as.integer(min(N - 1, max(25, ceiling(3 * lnN))))
      if(G.x)   {
        range.G    <- G.init <- tmp.G
        if(NlP) {
          if(verbose)               message(paste0("Since N < P, the sampler will be initialised with a different default of ceiling(log(N)) = ", lnN2, " groups (unless 'range.G' is supplied)"))
          G.init   <- max(2, lnN2)
        }
      }
      if(all(!G.x, NlP))  {
        G.init     <- range.G
        range.G    <- tmp.G
      }
      if(range.G    < lnN2)         warning(paste0("'range.G' should be at least log(N) (=log(", N, "))", " for the ", method, " method"), call.=FALSE)
      if(is.element(method, c("IMFA", "IMIFA"))) {
        if(any(!is.logical(ind.slice),
           length(ind.slice) != 1)) stop("'ind.slice' must be TRUE or FALSE")
        if(any(!is.logical(IM.lab.sw),
           length(IM.lab.sw) != 1)) stop("'IM.lab.sw' must be TRUE or FALSE")
        if(missing(rho)) {
          rho      <- 0.75
        }
        if(all(length(rho) > 1,
           rho > 1 && rho <= 0))    stop("'rho' must be a single number in the interval (0, 1]")
        if(rho < 0.5)               warning("Are you sure 'rho' should be less than 0.5? This could adversely affect mixing", call.=FALSE)
        if(missing(alpha.hyper))    {
          alpha.hyper     <- if(learn.alpha) c(2L, 1L) else c(0L, 0L)
        }
        if(all(discount    > 0, !learn.d)) {
          alpha.hyper     <- unname(unlist(shift_GA(shape=alpha.hyper[1], rate=alpha.hyper[2], shift=-discount)))
        }
        if(missing(zeta)) {
          zeta     <- 2
        }
        if(any(!is.numeric(zeta),
               length(zeta) != 1,
               zeta < 0))           stop(paste0("'zeta' must be single strictly positive number"))
        if(all(length(alpha.hyper)  != 2,
           learn.alpha))            stop(paste0("'alpha.hyper' must be a vector of length 2, giving the shape and rate hyperparameters of the gamma prior for alpha when 'learn.alpha' is TRUE"))
        a.hyp1     <- alpha.hyper[1]
        a.hyp2     <- alpha.hyper[2]
        if(learn.alpha)   {
          if(a.hyp1   <= 0)         stop("The shape of the gamma prior for alpha must be strictly positive")
          if(a.hyp2   <= 0)         stop("The rate of the gamma prior for alpha must be strictly positive")
        }
        t.miss <- missing(trunc.G)
        if(t.miss)        {
          trunc.G  <- range.G
        }
        if(length(trunc.G) > 1)     stop("'trunc.G' must be a single number")
        if(all(ifelse(N > 50, trunc.G < 50,
           trunc.G  < N), !t.miss)) warning(paste0("'trunc.G' should only be less than min(N=", N, ", 50) for practical reasons in heavy computational/memory burden cases"), call.=FALSE)
        if(trunc.G  < range.G)      stop(paste0("'trunc.G' must be at least range.G=", range.G))
        if(trunc.G  > N)            stop(paste0("'trunc.G' cannot be greater than N=", N))
        if(trunc.G  > 50)           warning(paste0("'trunc.G' is large: this may lead to memory capacity issues"), call.=FALSE)
      }
    } else if(method == "classify") {
      if(!zin.miss &&
         z.init    != "list") {     stop("'z.init' must be set to 'list' for classification")
      } else z.init       <- "list"
      if(zli.miss)                  stop("Data labels must be supplied via 'z.list' for classification")
      levs         <- nlevels(unlist(z.list))
      if(length(z.list)    > 1)     stop("Only one set of labels can be supplied via 'z.list'")
      zlabels      <- unlist(z.list)
      if(length(zlabels)  != N)     stop(paste0("'z.list' must be a factor of length N=",  N))
      if(all(verbose, !G.x) && any(length(range.G > 1),
          range.G  != levs))   {    message("Forced 'range.G' equal to the number of levels in 'zlabels' for the 'classify' method")
      }
      range.G      <- levs
    } else {
      range.G      <- 1L
    }
    meth    <- method
  } else {
    alp3    <- 3L   * alpha
    if(G.x)                         stop("'range.G' must be specified")
    if(any(range.G  < 1))           stop("'range.G' must be strictly positive")
    if(any(range.G  > alp3 * lnN))  warning(paste0("'range.G' MUCH greater than log(N) (=log(", N, ")):\n Empty clusters are likely, consider running an overfitted or infinite mixture"), call.=FALSE)
    range.G <- G.init     <- sort_unique(range.G)
    meth    <- rep(method, length(range.G))
  }
  if(any(range.G >= N))             stop(paste0("'range.G' must be less than the number of observations N=", N))
  if(G.init[1]   == 1)     {
    if(is.element(method, c("IMIFA", "IMFA",
       "OMIFA", "OMFA")))  {        stop(paste0("'method' must be FA or IFA for a one group model under the ", method, " method"))
    } else {
      meth[1]    <- switch(method,  MFA=, FA="FA", MIFA=, IFA="IFA")
    }
    if(all(verbose, !is.element(method,
           c("FA", "IFA"))))        message(paste0("Forced use of ", meth[1], " method where 'range.G' is equal to 1"))
  }

# Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  cov.mat          <- if(P > 500) switch(scaling, unit=cora(as.matrix(dat)), cova(as.matrix(dat))) else switch(scaling, unit=cor(dat), cov(dat))
  datname          <- rownames(dat)
  if(any(length(unique(datname)) != N,
     is.null(datname)))      rownames(dat) <- seq_len(N)
  sigmu.miss       <- missing("sigma.mu")
  if(sigmu.miss)             sigma.mu      <- diag(cov.mat)
  if(scaling == "unit")      sigma.mu      <- sigma.mu[1]
  if(any(sigma.mu  <= 0, !is.numeric(sigma.mu),
     !is.element(length(sigma.mu),
     c(1, P))))                     stop(paste0("'sigma.mu' must be strictly positive, and of length 1 or P=", P))
  if(missing("psi.alpha"))   psi.alpha     <- 2.5
  if(any(psi.alpha <= 1, !is.numeric(psi.alpha),
     length(psi.alpha)   != 1))     stop("'psi.alpha' must be a single number strictly greater than 1 in order to bound uniquenesses away from zero")
  Q.miss    <- missing(range.Q)
  Q.min     <- min(ceiling(log(P)), ceiling(log(N)))
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(Q.miss)                      stop("'range.Q' must be specified")
    if(any(range.Q < 0))            stop(paste0("'range.Q' must be non-negative for the ", method, " method"))
    range.Q <- sort_unique(range.Q)
  } else {
    if(Q.miss)        range.Q    <- as.integer(min(ifelse(P > 500, 12 + floor(log(P)), floor(3 * log(P))), N - 1))
    if(any(!is.logical(adapt),
           length(adapt) != 1))     stop("'adapt' must be TRUE or FALSE")
    if(length(range.Q)    > 1)      stop(paste0("Only one starting value for 'range.Q' can be supplied for the ", method, " method"))
    if(range.Q    <= 0)             stop(paste0("'range.Q' must be strictly positive for the ", method, " method"))
    if(all(adapt, range.Q < Q.min)) stop(paste0("'range.Q' must be at least min(log(P), log(N)) for the ", method, " method when 'adapt' is TRUE"))
  }
  len.G     <- switch(method, classify=range.G, length(range.G))
  len.Q     <- length(range.Q)
  len.X     <- len.G * len.Q
  if(all(len.X > 10,
         suppressWarnings(memory.limit())  <= 16256,
         switches["s.sw"])) {
    if(!score.x)            {       warning(paste0("The large number of candidate models being explored (", len.X, ") could lead to memory issues\nConsider setting 'score.switch' to FALSE or breaking up the task into chunks and calling get_IMIFA_results() on each chunk"), call.=FALSE)
    } else                  {       warning(paste0("'score.switch' set to FALSE as too many candidate models are being explored (", len.X, ")\nPosterior inference on the scores will not be possible, though you can risk forcing storage by supplying score.switch=TRUE\nConsider breaking up the task into chunks and calling get_IMIFA_results() on each chunk"), call.=FALSE)
      switches["s.sw"]   <- FALSE
    }
  }
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(missing("sigma.l"))   sigma.l       <- 1L
    if(any(sigma.l <= 0, !is.numeric(sigma.l),
           length(sigma.l) != 1))   stop("'sigma.l' must be a single strictly positive number")
  } else {
    if(missing("nu"))        nu            <- 2L
    if(any(nu <= !nuplus1, !is.numeric(nu),
           length(nu) != 1))        stop(paste0("'nu' must be a single ", ifelse(nuplus1, "strictly positive number for the Ga(nu + 1, nu) parameterisation", "number strictly greater than 1 for the Ga(nu, nu) parameterisation")))
    if(missing("beta.d1"))   beta.d1       <- 1L
    if(missing("beta.d2"))   beta.d2       <- 1L
    if(any(!is.numeric(beta.d1),
           !is.numeric(beta.d1),
           length(beta.d1) != 1,
           length(beta.d2) != 1))   stop("'beta.d1' and 'beta.d2' must both be numberic and of length 1")
    if(missing("b0"))        b0            <- 0.1
    if(any(length(b0) != 1, !is.numeric(b0),
           b0  < 0))                stop("'b0' must be a non-negative scalar to ensure valid adaptation probability")
    if(missing("b1"))        b1            <- 0.00005
    if(any(length(b1) != 1, !is.numeric(b1),
           b1 <= 0))                stop("'b1' must be a single strictly positive scalar to ensure adaptation probability decreases")
    if(missing("prop"))      prop          <- floor(0.7 * P)/P
    if(missing("adapt.at"))  adapt.at      <- switch(method, IFA=, MIFA=burnin, 0)
    if(missing("epsilon"))   epsilon       <- ifelse(any(centered, centering), 0.1, 0.05)
    if(any(length(prop)    != 1, length(adapt.at) != 1,
           length(epsilon) != 1))   stop("'prop', 'adapt.at', and 'epsilon' must all be of length 1")
    if(any(!is.numeric(prop), !is.numeric(adapt.at),
           !is.numeric(epsilon)))   stop("'prop', 'adapt.at', and 'epsilon' must all be numeric")
    if(abs(prop - (1 - prop)) < 0)  stop("'prop' must be lie in the interval (0, 1)")
    if(adapt.at < 0 ||
       adapt.at > burnin)           stop("'adapt.at' must be lie in the interval [0, burnin]")
    if(epsilon <= 0 ||
       epsilon >= 1)                stop("'epsilon' must be lie in the interval (0, 1)")
  }
  Q.warn       <- min(N - 1, Ledermann(P))
  if(any(range.Q > Q.warn))   {
    if(all(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")),
       isTRUE(adapt)))        {     warning(paste0("Starting value for number of factors is greater than ", ifelse(any(range.Q > P), paste0("the number of variables (", P, ")"), paste0("the suggested Ledermann upper bound (", Q.warn, ")"))), call.=FALSE)
    } else if(any(is.element(method, c("FA", "MFA", "OMFA", "IMFA")),
              all(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")),
                  isTRUE(!adapt)))) warning(paste0("Number of factors is greater than ", ifelse(any(range.Q > P), paste0("the number of variables (", P, ")"), paste0("the suggested Ledermann upper bound (", Q.warn, ")"))), call.=FALSE)
  }
  if(any(all(method == "MFA",  any(range.G > 1)) && any(range.Q > 0),
         all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA",
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!switches["l.sw"],
           !switches["psi.sw"]))  {
                                    warning("Loadings & Uniquenesses not stored: will be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
    } else if(!switches["l.sw"])  { warning("Loadings not stored: may be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
    } else if(!switches["psi.sw"])  warning("Uniquenesses not stored: will be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
  }
  if(any(all(method == "MFA",  any(range.G > 1)),
         all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA",
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!switches["mu.sw"],
           !switches["psi.sw"]))  {
                                    warning("Means & Uniquenesses not stored: posterior mean estimates won't be available", call.=FALSE)
    } else if(!switches["mu.sw"]) { warning("Means not stored: posterior mean estimates won't be available", call.=FALSE)
    } else if(!switches["psi.sw"])  warning("Uniquenesses not stored: posterior mean estimates won't be available", call.=FALSE)
  }
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")) && all(range.Q == 0)) {
    if(all(switches[c("s.sw", "l.sw")]))  {
                                    warning("Scores & Loadings not stored as model has zero factors", call.=FALSE)
    } else if(switches["s.sw"])   { warning("Scores not stored as model has zero factors", call.=FALSE)
    } else if(switches["l.sw"])   { warning("Loadings not stored as model has zero factors", call.=FALSE)
    }
    switches[c("s.sw", "l.sw")]  <- FALSE
  } else {
    if(all(!switches[c("s.sw", "l.sw")])) {
                                    warning("Posterior Scores & Loadings won't be available as they're not being stored", call.=FALSE)
    } else if(!switches["s.sw"])  { warning("Posterior Scores won't be available as they're not being stored", call.=FALSE)
    } else if(!switches["l.sw"])  { warning("Posterior Loadings won't be available as they're not being stored", call.=FALSE)
    }
  }
  if(all(is.element(method, c("FA", "IFA")), verbose,
         !zin.miss ||
         !zli.miss))                message(paste0("z does not need to be initialised for the ", method, " method"))
  if(is.element(method, c("MFA", "MIFA", "classify"))) {
    if(any(!is.logical(mu0g),
           length(mu0g)    != 1))   stop("'mu0g' must be TRUE or FALSE")
    if(any(!is.logical(delta0g),
           length(delta0g) != 1))   stop("'delta0g' must be TRUE or FALSE")
    if(any(!is.logical(psi0g),
           length(psi0g)   != 1))   stop("'psi0g' must be TRUE or FALSE")
    if(all(method == "MFA",
           delta0g))                stop("'delta0g' cannot be TRUE for the 'MFA' method")
    if(method == "classify") mu0g          <- TRUE
  }
  dimension <- PGMM_dfree(Q=switch(method, FA=, classify=, MFA=, OMFA=, IMFA=range.Q,
               IFA=, MIFA=, OMIFA=, IMIFA=0), P=P, method=switch(uni.type, unconstrained="UUU", isotropic="UUC"))
  min.d2    <- min(dimension)/2
  min.d2G   <- min.d2 * G.init
  sw0gs     <- c(mu0g = mu0g, psi0g = psi0g, delta0g = delta0g)
  if(all(!is.element(method, c("MFA", "MIFA", "classify")),
         any(sw0gs)))               stop(paste0("'", names(which(sw0gs)), "' should be FALSE for the ", method, " method\n"))
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    if(missing("alpha"))   { alpha         <- switch(method, MFA=, MIFA=1, OMFA=, OMIFA=0.5/G.init, if(learn.alpha) max(1, rgamma(1, a.hyp1, a.hyp2)) - discount else 1 - discount)
    } else if(is.element(method,
      c("OMFA", "OMIFA")))   alpha         <- alpha/G.init
    if(length(alpha) != 1)          stop("'alpha' must be specified as a scalar to ensure an exchangeable prior")
    if(alpha <= -discount)          stop(paste0("'alpha' must be ",     ifelse(discount != 0, paste0("strictly greater than -discount (i.e. > ", - discount, ")"), "strictly positive")))
    if(all(is.element(method,  c("IMIFA",   "IMFA")),
       !learn.alpha))               warning(paste0("'alpha' fixed at ", ifelse(discount != 0, paste0("1 - 'discount' = "), ""), alpha, " as it's not being learned via Gibbs/Metropolis-Hastings updates"), call.=FALSE)
    if(all(is.element(method,  c("OMIFA",   "OMFA")),
       alpha >= min.d2))            warning(paste0("'alpha' over 'range.G' for the OMFA & OMIFA methods must be less than half the dimension (per group!)\n of the free parameters of the smallest model considered (= ", min.d2, "): consider suppling 'alpha' < ", min.d2G), call.=FALSE)
    if(any(all(is.element(method, c("MFA",  "MIFA")), alpha > 1),
           all(is.element(method, c("OMFA", "OMIFA")),
           alpha > 1/G.init)))      warning("Are you sure alpha should be greater than 1?", call.=FALSE)
                             z.init        <- match.arg(z.init)
    if(all(any(is.nan(rDirichlet(G=min(range.G), alpha = alpha))),
           is.element(method,  c("MFA",     "MIFA",
            "OMFA", "OMIFA"))))     stop("'alpha' is too small for simulated mixing proportions to be valid")
    if(all(is.element(method,  c("OMIFA",   "OMFA")), is.element(z.init,
       "priors")))                  stop(paste0("'z.init' cannot be set to 'priors' for the ", method, " method to ensure all groups are populated at the initialisation stage"))
    if(!zli.miss) {
      if(length(z.list)   != len.G)  {
                                    stop(paste0("'z.list' must be a list of length ", len.G))  }
                             list.levels   <- lapply(z.list, nlevels)
      if(!all(list.levels == G.init))             {
        if(!is.element(method, c("IMIFA",
                       "IMFA", "OMIFA", "OMFA"))) {
                                    stop(paste0("Each element of 'z.list' must have the same number of levels as 'range.G'"))
        } else                      stop(paste0("Only ", list.levels, " groups are populated according to z.list, but 'range.G' has been set to ", G.init, ":\n  Reset 'range.G' to this value to avoid redunandtly carrying around empty groups or supply a list with ", G.init, " levels"))
      }
      if(!all(lengths(z.list) == N)) {
                                    stop(paste0("Each element of 'z.list' must be a vector of length N=", N)) }
    }
    if(all(zli.miss, z.init   == "list"))         {
                                    stop(paste0("'z.list' must be supplied if 'z.init' is set to 'list'")) }
  }
  imifa     <- list(list())
  Gi        <- Qi  <- 1L
  gibbs.arg <- list(P = P, sigma.mu = sigma.mu, psi.alpha = psi.alpha, burnin = burnin, sw = switches,
                    thinning = thinning, iters = iters, verbose = verbose, uni.type = uni.type, uni.prior = uni.prior)
  if(is.element(method, c("IMIFA", "IMFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(trunc.G = trunc.G, rho = rho, ind.slice = ind.slice, learn.alpha = learn.alpha, learn.d = learn.d, kappa = kappa,
                                             zeta = zeta, IM.lab.sw = IM.lab.sw, a.hyper = alpha.hyper, discount = discount, d.hyper = d.hyper))
  }
  if(is.element(method, c("FA", "IFA", "MFA", "MIFA")))   {
    gibbs.arg      <- append(gibbs.arg, list(scaling = scaling))
  }
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(nu = nu, beta.d1 = beta.d1, beta.d2 = beta.d2, adapt.at = adapt.at, adapt = adapt,
                                             nuplus1 = nuplus1, b0 = b0, b1 = b1, prop = prop, epsilon = epsilon))
    temp.args      <- gibbs.arg
  } else {
    gibbs.arg      <- append(gibbs.arg, list(sigma.l = sigma.l))
  }

  init.start       <- proc.time()
  mu               <- list(colMeans(dat))
  beta.x           <- missing("psi.beta")
  mu0.x            <- missing("mu.zero")
  ad1.x            <- all(missing("alpha.d1"), is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA", "classify")))
  ad2.x            <- all(missing("alpha.d2"), is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA", "classify")))
  if(all(z.init != "list", any(sw0gs))) {
    if(delta0g)                     stop("'delta0g' can only be TRUE if z.init=list\n")
    if(all(!mu0.x, mu0g))           stop("'mu.zero' can only be supplied for each group if z.init=list")
    if(all(!beta.x, psi0g))         stop("'psi.beta' can only be supplied for each group if z.init=list")
  }
  if(beta.x) {
    psi.beta       <- temp.psi <- list(psi_hyper(shape=psi.alpha, covar=cov.mat, type=uni.prior))
  } else {
    psi.beta       <- .len_check(psi.beta, psi0g, method, P, G.init)
  }
  mu.zero          <- if(mu0.x) mu else .len_check(mu.zero, mu0g, method, P, G.init)
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    alpha.d1       <- if(ad1.x) list(3L) else .len_check(alpha.d1, delta0g, method, P, G.init, P.dim=FALSE)
    alpha.d2       <- if(ad2.x) list(6L) else .len_check(alpha.d2, delta0g, method, P, G.init, P.dim=FALSE)
    if(all(NlP, verbose,
           any(ad1.x, ad2.x)))      message("Consider applying more shrinkage with higher 'alpha.d1' and 'alpha.d2' hyperparameter values when N << P")
  }
  if(!is.element(method, c("FA", "IFA"))) {
    if(verbose)                     cat(paste0("Initialising...\n"))
    clust          <- list()
    pi.alpha       <- list()
    pi.prop        <- list()
    zi             <- list()
    for(g in seq_along(G.init)) {
      G            <- G.init[g]
      if(z.init    == "kmeans")      {
        k.res      <- try(suppressWarnings(kmeans(dat, centers=G, iter.max=20, nstart=100)), silent=TRUE)
        if(!inherits(k.res, "try-error"))  {
          zi[[g]]  <- as.integer(factor(k.res$cluster,        levels=seq_len(G)))
        } else                      stop("Cannot initialise cluster labels using kmeans. Try another z.init method")
      } else if(z.init  == "list")   {
        zi[[g]]    <- as.integer(z.list[[g]])
      } else if(z.init  == "mclust") {
        m.res      <- try(Mclust(dat, G, verbose=FALSE), silent=TRUE)
        if(!inherits(m.res, "try_error"))  {
          zi[[g]]  <- as.integer(factor(m.res$classification, levels=seq_len(G)))
        } else                      stop("Cannot initialise cluster labels using mclust. Try another z.init method")
      } else {
        zips       <- rep(1, N)
        if(!is.element(method, c("IMFA", "IMIFA"))) {
          while(all(length(unique(zips)) != G,
                any(prop.table(tabulate(zips, nbins=G)) < 1/G^2))) {
            pies   <- rDirichlet(alpha=alpha, G)
            zips   <- .sim_z_p(N=N, prob.z=pies)
          }
        } else {
          vies     <- .sim_vs_inf(alpha=alpha, N=N, nn=rep(N/G.init, G.init), discount=discount, len=G.init, lseq=seq_len(G.init))
          pies     <- .sim_pi_inf(vs=vies, len=G.init)
          zips     <- .sim_z_p(N=N, prob.z=pies)
        }
        zi[[g]]    <- as.integer(zips)
        rm(zips)
      }
      nngs         <- tabulate(zi[[g]], nbins=switch(method, IMFA=, IMIFA=trunc.G, G))
      pi.prop[[g]] <- prop.table(nngs)
      mu[[g]]      <- vapply(seq_len(G), function(gg) if(nngs[gg] > 0) colMeans(dat[zi[[g]] == gg,, drop=FALSE]) else rep(0, P), numeric(P))
      if(mu0.x)   {
        mu.zero[[g]]    <- if(mu0g) mu[[g]] else vapply(seq_len(G), function(gg) colMeans(dat), numeric(P))
      }
      if(beta.x)  {
        if(psi0g) {
          cov.gg   <- lapply(seq_len(G), function(gg, dat.gg = dat[zi[[g]] == gg,, drop=FALSE]) if(all(nngs[gg] > 1, P > nngs[g])) { if(P > 500) cova(as.matrix(dat.gg)) else cov(dat.gg) } else cov.mat)
          psi.beta[[g]] <- vapply(seq_len(G), function(gg) psi_hyper(shape=psi.alpha, covar=cov.gg[[gg]], type=uni.prior), numeric(P))
        } else {
          psi.beta[[g]] <- replicate(G, temp.psi[[1]])
        }
      }
      if(ad1.x)   {
        alpha.d1[[g]]   <- rep(unlist(alpha.d1), G)
      }
      if(ad2.x)   {
        alpha.d2[[g]]   <- rep(unlist(alpha.d2), G)
      }
      clust[[g]]   <- list(z = zi[[g]], pi.alpha = alpha, pi.prop = pi.prop[[g]])
      if(is.element(method, c("MFA", "MIFA"))) {
        sw0g.tmp   <- sw0gs
        if(all(g > 9, any(sw0gs))) {
          sw0g.tmp <- setNames(rep(FALSE, 4L), names(sw0gs))
                                    warning(paste0(names(which(sw0gs)), " set to FALSE where G > 9, as 'exact' label-switching is not possible in this case\n"), call.=FALSE)
        }
        clust[[g]] <- append(clust[[g]], list(l.switch = sw0g.tmp))
      }
      if(is.element(method, c("classify", "MIFA"))) {
        clust[[g]] <- append(clust[[g]], list(alpha.d1 = alpha.d1[[g]], alpha.d2 = alpha.d2[[g]]))
      }
    }
  }
  if(is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))) {
    mu.zero        <- if(all(lengths(mu.zero)  == 1)) list(mu.zero[[1]])  else list(mu.zero[[1]][,1])
    psi.beta       <- if(all(lengths(psi.beta) == 1)) list(psi.beta[[1]]) else list(psi.beta[[1]][,1])
    if(!is.element(method, c("OMFA", "IMFA"))) {
      alpha.d1     <- list(alpha.d1[[1]][1])
      alpha.d2     <- list(alpha.d2[[1]][1])
    }
  }
  if(all(round(vapply(mu.zero, sum, numeric(1L))) == 0)) {
    mu.zero        <- switch(method, classify=base::matrix(0, nr=1, nc=range.G), lapply(mu.zero, function(x) 0))
  }
  if(anyNA(unlist(psi.beta))) {
    psi.beta       <- lapply(psi.beta, function(x) replace(x, is.na(x), 0))
  }
  if(all(uni.type == "isotropic", unlist(lapply(psi.beta, function(x) {
    if(is.matrix(x)) any(apply(x, 2, function(y) {
      length(unique(round(y, nchar(y)))) }) != 1)  else  {
      length(unique(round(x,
      nchar(x)))) != 1 } } ))))     stop("'psi.beta' cannot be variable specific if 'uni.type' is 'isotropic'")
  if(any(unlist(psi.beta)   <= 0))  stop("'psi.beta' must be strictly positive")
  if(is.element(method, c("classify", "IFA", "MIFA", "IMIFA", "OMIFA"))) {
    if(!all(MGP_check(ad1=unlist(alpha.d1), ad2=unlist(alpha.d2), Q=unique(range.Q), nu=nu, bd1=beta.d1, bd2=beta.d2,
        plus1=nuplus1)[[1]]$valid)) stop("Invalid shrinkage hyperparameter values will not encourage loadings column removal.\n Try using the MGP_check() function in advance to ensure cumulative shrinkage property holds.")
    deltas         <- lapply(seq_along(G.init), function(g) list(alpha.d1 = alpha.d1[[g]], alpha.d2 = alpha.d2[[g]]))
  }
  init.time        <- proc.time() - init.start
  fac.time         <- 0
  G.range          <- switch(method, IMIFA=, IMFA=G.init, range.G)

  if(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA"))) {
    if(len.G == 1)  {
      start.time   <- proc.time()
      if(meth[Gi]  != "MIFA") {
        gibbs.arg  <- append(temp.args, deltas[[Gi]])
      }
      imifa[[Gi]][[Qi]]   <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = G.range, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "IFA") clust[[Gi]]), gibbs.arg))
      fac.time     <- fac.time + imifa[[Gi]][[Qi]]$time
    } else {
      start.time   <- proc.time()
      for(g in G.range)   {
        Gi         <- which(G.range == g)
        if(meth[Gi]  == "IFA") {
         gibbs.arg <- append(temp.args, deltas[[Gi]])
        }
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MIFA") clust[[Gi]]), gibbs.arg))
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && Gi  != len.G) cat(paste0("Model ", Gi, " of ", len.G, " complete"), "Initialising...", sep="\n")
      }
    }
  } else if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")))   {
    if(all(len.G == 1, len.Q == 1)) {
      start.time   <- proc.time()
      imifa[[Gi]][[Qi]]   <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = G.range, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] !=  "FA") clust[[Gi]]), gibbs.arg))
      fac.time     <- fac.time + imifa[[Gi]][[Qi]]$time
    } else if(len.G == 1) {
      start.time   <- proc.time()
      for(q in range.Q)   {
        Qi         <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = G.range, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] !=  "FA") clust[[Gi]]), gibbs.arg))
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && Qi  != len.Q) cat(paste0("Model ", Qi, " of ", len.Q, " complete"), "Initialising...", sep="\n")
      }
    } else if(len.Q == 1) {
      start.time   <- proc.time()
      for(g in G.range)   {
        Gi         <- which(G.range == g)
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && Gi  != len.G) cat(paste0("Model ", Gi, " of ", len.G, " complete"), "Initialising...", sep="\n")
      }
    } else {
      mi           <- 0
      start.time   <- proc.time()
      for(g in G.range)   {
        Gi         <- which(G.range == g)
        imifa[[Gi]]       <- list()
        for(q in range.Q) {
          Qi       <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        mi         <- mi + 1
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && mi  != len.X) cat(paste0("Model ", mi, " of ", len.X, " complete"), "Initialising...", sep="\n")
        }
      }
    }
  } else if(method == "classify") { stop("'classify' method not yet implemented")
    start.time     <- proc.time()
    if(centered)                    warning("Data supplied is globally centered, are you sure?", call.=FALSE)
    for(g in seq_len(range.G))  {
      tmp.dat      <- raw.dat[zlabels == levels(zlabels)[g],]
      scal         <- switch(scaling, none=FALSE, Rfast::colVars(as.matrix(tmp.dat), std=TRUE))
      scal         <- switch(scaling, pareto=sqrt(scal), scal)
      tmp.dat <- if(is.logical(scal)) standardise(as.matrix(tmp.dat), center=centering, scale=scal) else scale(tmp.dat, center=centering, scale=scal)
      if(sigmu.miss)   gibbs.arg$sigma.mu  <- if(nrow(tmp.dat) > 1) { if(P > 500) cova(as.matrix(tmp.dat)) else cov(tmp.dat) } else cov.mat
      imifa[[g]]          <- list()
      gibbs.arg    <- append(temp.args, lapply(deltas[[Gi]], "[[", g))
      imifa[[g]][[Qi]]    <- do.call(paste0(".gibbs_", "IFA"),
                                     args=append(list(data = tmp.dat, N = nrow(tmp.dat), mu = mu[[Gi]][,g], mu.zero = mu.zero[[Gi]][,g],
                                                      Q = range.Q, psi.beta = psi.beta[[Gi]][,g]), gibbs.arg))
      fac.time     <- fac.time + imifa[[g]][[Qi]]$time
      if(verbose   && g   != len.G) cat(paste0("Model ", g, " of ", len.G, " complete"), "Initialising...", sep="\n")
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/switch(method, MFA=len.X, MIFA=len.G, classify=range.G, len.Q)
  tot.time  <- tot.time    + init.time
  init.time <- init.time   + fac.time
  for(g in length(imifa)) {
   for(q in length(imifa[[g]])) {
     imifa[[g]][[q]]$time <- NULL
   }
  }

  imifa     <- switch(method, FA=, MFA=, OMFA=, IMFA={
     lapply(seq_along(imifa), function(x) setNames(imifa[[x]], paste0(range.Q, ifelse(range.Q == 1, "Factor", "Factors"))))
  }, lapply(seq_along(imifa), function(x) setNames(imifa[[x]], "IFA")))
  if(!is.element(method, c("FA", "IFA"))) {
    for(g in seq_along(G.init)) {
      attr(imifa[[g]],
           "Z.init")      <- factor(zi[[g]], levels=seq_len(G.init[g]))
    }
  }
  gnames    <- switch(method, classify=paste0("Group ", seq_len(range.G)), paste0(G.init, ifelse(G.init == 1, "Group", "Groups")))
  names(imifa)            <- gnames
  attr(imifa,
       "Alph.step")       <- learn.alpha
  attr(imifa, "Alpha")    <- if(!learn.alpha) alpha
  if(method == "classify") {
    attr(imifa,
         "Class.Props")   <- tabulate(z.list[[1]], range.G)/N
  }
  attr(imifa, "Call")     <- call
  attr(imifa, "Center")   <- any(centered, centering)
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa,
       "Disc.step")       <- learn.d
  attr(imifa, "Discount") <- if(!learn.d) discount
  attr(imifa, "Factors")  <- range.Q
  attr(imifa, "IM.labsw") <- all(is.element(method, c("IMFA", "IMIFA")), IM.lab.sw)
  attr(imifa,
       "Ind.Slice")       <- all(is.element(method, c("IMFA", "IMIFA")), ind.slice)
  attr(imifa, "Groups")   <- range.G
  if(!is.element(method, c("FA", "IFA"))) {
    attr(imifa, "Init.Z") <- z.init
    attr(imifa,
         "Label.Switch")  <- any(sw0gs)
  }
  method                  <- names(table(meth)[which.max(table(meth))])
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1, 1)),
                                    substr(method, 2, nchar(method)))
  attr(imifa, "Name")     <- dat.nam
  attr(imifa, "Nuplus1")  <- all(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")), nuplus1)
  attr(imifa, "Obs")      <- N
  attr(imifa, "Pitman")   <- all(is.element(method, c("IMFA", "IMIFA")), any(learn.d, discount > 0))
  attr(imifa, "Rho")      <- rho
  attr(imifa, "Scaling")  <- scal
  attr(attr(imifa,
  "Scaling"), "Method")   <- scaling
  attr(imifa, "Store")    <- length(iters)
  switches                <- c(switches, a.sw = learn.alpha, d.sw = learn.d)
  if(is.element(method, c("FA", "IFA"))) {
    switches["pi.sw"]     <- FALSE
  }
  attr(imifa, "Switch")   <- switches
  times                   <- lapply(list(Total = tot.time, Average = avg.time, Initialisation = init.time), round, 2)
  if(all(len.G  == 1,
         len.Q  == 1)) {
    times                 <- times[-2]
  }
  if(is.element(method, c("FA", "IFA", "classify"))) {
    times                 <- times[-length(times)]
  } else {
    attr(imifa, "G.init") <- G.init
  }
  class(times)            <- "listof"
  attr(imifa, "Time")     <- times
  attr(imifa, "Uni.Meth") <- c(Uni.Prior = uni.prior, Uni.Type = uni.type)
  attr(imifa, "Vars")     <- P
  if(verbose)                print(attr(imifa, "Time"))
  class(imifa)            <- "IMIFA"
    return(imifa)
}
