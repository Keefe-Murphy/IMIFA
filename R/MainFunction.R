#' Adaptive Gibbs Sampler for Nonparameteric Model-based Clustering using models from the IMIFA family
#'
#' Carries out Gibbs sampling for all models from the IMIFA family, facilitating model-based clustering with dimensionally reduced factor-analytic covariance structures, with automatic estimation of the number of clusters and cluster-specific factors as appropriate to the method employed. Factor analysis with one group (FA/IFA), finite mixtures (MFA/MIFA), overfitted mixtures (OMFA/OMIFA), infinite factor models which employ the multiplicative gamma process (MGP) shrinkage prior (IFA/MIFA/OMIFA/IMIFA), and infinite mixtures which employ Dirichlet/Pitman-Yor Process Mixture Models (IMFA/IMIFA) are all provided.
#'
#' @param dat A matrix or data frame such that rows correspond to observations (\code{N}) and columns correspond to variables (\code{P}). Non-numeric variables and rows with missing entries will be removed.
#' @param method An acronym for the type of model to fit where:
#' \describe{
#'  \item{"\code{FA}"}{Factor Analysis}
#'  \item{"\code{IFA}"}{Infinite Factor Analysis}
#'  \item{"\code{MFA}"}{Mixtures of Factor Analysers}
#'  \item{"\code{MIFA}"}{Mixtures of Infinite Factor Analysers}
#'  \item{"\code{OMFA}"}{Overfitted Mixtures of Factor Analysers}
#'  \item{"\code{OMIFA}"}{Overfitted Mixtures of Infinite Factor Analysers}
#'  \item{"\code{IMFA}"}{Infinite Mixtures of Factor Analysers}
#'  \item{"\code{IMIFA}"}{Infinite Mixtures of Infinite Factor Analysers}
#' }
#'  In principle, of course, one could overfit the "\code{MFA}" or "\code{MIFA}" models, but it is recommend to use the corresponding model options which begin with 'O' instead. Note that the "\code{classify}" method is not yet implemented.
#' @param n.iters The number of iterations to run the Gibbs sampler for.
#' @param range.G Depending on the method employed, either the range of values for the number of clusters, or the conseratively high starting value for the number of clusters. Defaults to 1 for the "\code{FA}" and "\code{IFA}" methods. For the "\code{MFA}" and "\code{MIFA}" models this is to be given as a range of candidate models to explore. For the "\code{OMFA}", "\code{OMIFA}", "\code{IMFA}", and "\code{IMIFA}" models, this is the number of clusters with which the chain is to be initialised, in which case the default is \code{min(N - 1, max(25, ceiling(3 * log(N))))}. For the "\code{OMFA}", and "\code{OMIFA}" models this upper limit remains fixed for the entire length of the chain; \code{range.G} also doubles as the default \code{trunc.G} for the "\code{IMFA}" and "\code{IMIFA}" models. However, when \code{N < P}, or when this bound is close to or exceeds \code{N} for any of these overfitted/infinite mixture models, it is better to initialise at a value closer to the truth (i.e. \code{ceiling(log(N))} by default), though the upper bound remains the same - as a result the role of \code{range.G} when \code{N < P} is no longer to specify the upper bound (which can still be modified via \code{trunc.G}, at least for the "\code{IMFA}" and "\code{IMIFA}" methods) and the number of clusters used for initialisation, but rather just the number of clusters used for initialisation only. If \code{length(range.G) * length(range.Q)} is large, consider not storing unnecessary parameters (via \code{\link{storeControl}}), or breaking up the range of models to be explored into chunks and sending each chunk to \code{\link{get_IMIFA_results}}.
#' @param range.Q Depending on the method employed, either the range of values for the number of latent factors, or, for methods ending in IFA the conservatively high starting value for the number of cluster-specific factors, in which case the default starting value is \code{floor(3 * log(P))}. For methods ending in IFA, different clusters can be modelled using different numbers of latent factors (incl. zero); for methods not ending in IFA it is possible to fit zero-factor models, corresponding to simple diagonal covariance structures. For instance, fitting the "\code{IMFA}" model with \code{range.Q=0} corresponds to a vanilla Dirichlet Process Mixture Model. If \code{length(range.G) * length(range.Q)} is large, consider not storing unnecessary parameters (via \code{\link{storeControl}}), or breaking up the range of models to be explored into chunks and sending each chunk to \code{\link{get_IMIFA_results}}.
#' @param burnin The number of burn-in iterations for the sampler. Defaults to \code{n.iters/5}. Note that chains can also be burned in later, using \code{\link{get_IMIFA_results}}.
#' @param thinning The thinning interval used in the simulation. Defaults to 2. No thinning corresponds to 1. Note that chains can also be thinned later, using \code{\link{get_IMIFA_results}}.
#' @param MGP A list of arguments pertaining to the multiplicative gamma process (MGP) shrinkage prior and adaptive Gibbs sampler (AGS). For use with the infinite factor models "\code{IFA}", "\code{MIFA}", "\code{OMIFA}", and "\code{IMIFA}" only. Defaults are set by a call to \code{\link{mgpControl}}, with further checking of validity by \code{\link{MGP_check}}
#' @param centering A logical value indicating whether mean centering should be applied to the data, defaulting to \code{TRUE}.
#' @param scaling The scaling to be applied - one of "\code{unit}", "\code{none}" or "\code{pareto}".
#' @param uni.type This argument specifies the type of constraint, if any, to be placed on the uniquenesses/idiosyncratic variances, i.e. whether a general diagonal matrix or isotropic diagonal matrix is to be assumed, and in turn whether these matrices are constrained to be equal across clusters. The default "\code{unconstrained}" corresponds to factor analysis (and mixtures thereof), whereas "\code{isotropic}" corresponds to probabilistic principal components analysers (and mixtures thereof). Constraints \emph{may} be particularly useful when \code{N < P}, though caution is advised when employing constraints for any of the infinite factor models, especially "\code{isotropic}" and "\code{single}", which may lead to overestimation of the number of clusters &/or factors if this specification is inappropriate. The four options correspond to the following 4 parsimonious Gaussian mixture models:
#' \describe{
#' \item{"\code{unconstrained}"}{\strong{UUU} - variable-specific and cluster-specific: \eqn{\Psi_g = \Psi_g}{Psi_g = Psi_g}.}
#' \item{"\code{isotropic}"}{\strong{UUC} - cluster-specific, equal across variables: \eqn{\Psi_g = \sigma_g^2 \mathcal{I}_p}{Psi_g = (sigma^2)_g I_p}.}
#' \item{"\code{constrained}"}{\strong{UCU} - variable-specific, equal across clusters: \eqn{\Psi_g = \Psi}{Psi_g = Psi}.}
#' \item{"\code{single}"}{\strong{UCC} - single value equal across clusters and variables: \eqn{\Psi_g = \sigma^2 \mathcal{I}_p}{Psi_g = sigma^2 I_p}.}
#' }
#'  The first letter \strong{U} here corresponds to constraints on loadings (not yet implemented), the second letter corresponds to constrained/unconstrained across clusters, and the third letter corresponds to the isotropic constraint. Of course, only the third letter is of relevance for the single-cluster "\code{FA}" and "\code{IFA}" models, such that "\code{unconstrained}" and "\code{constrained}" are equivalent for these models, and so too are "\code{isotropic}" and "\code{single}".
#' @param uni.prior A switch indicating whether uniquenesses rate hyperparameters are to be "\code{unconstrained}" or "\code{isotropic}", i.e. variable-specific or not. "\code{uni.prior}" must be "\code{isotropic}" if the last letter of "\code{uni.type}" is \strong{C}, but can take either value otherwise. Defaults to correspond to the last letter of \code{uni.type} if that is supplied and \code{uni.prior} is not, otherwise defaults to "\code{unconstrained}", but "\code{isotropic}" is recommended when \code{N < P}. Only relevant when "\code{psi.beta}" is not supplied and \code{\link{psi_hyper}} is therefore invoked.
#' @param alpha Depending on the method employed, either the hyperparameter of the Dirichlet prior for the cluster mixing proportions, or the Dirichlet process concentration parameter. Defaults to 0.5/range.G for the Overfitted methods - if supplied for "\code{OMFA}" and "\code{OMIFA}" methods, you are supplying the numerator of \code{alpha/range.G}, which should be less than half the dimension (per cluster!) of the free parameters of the smallest model considered in order to ensure superfluous clusters are emptied (for "\code{OMFA}", this corresponds to the smallest \code{range.Q}; for "\code{OMIFA}", this corresponds to a zero-factor model) [see: \code{\link{PGMM_dfree}} and Rousseau and Mengersen (2011)]. Defaults to 1 for the finite mixture models "\code{MFA}" and "\code{MIFA}". Defaults to \code{1 - discount} for the "\code{IMFA}" and "\code{IMIFA}" models if \code{learn.alpha=FALSE} or a simulation from the prior if \code{learn.alpha=TRUE}. Must be positive, unless \code{discount} is supplied for the "\code{IMFA}" or "\code{IMIFA}" methods.
#' @param storage A vector of named logical indicators governing storage of parameters of interest for all models in the IMIFA family. Defaults are set by a call to \code{\link{storeControl}}. It may be useful not to store certain parameters if memory is an issue.
#' @param psi.alpha The shape of the inverse gamma prior on the uniquenesses. Defaults to 2.5 if \code{uni.type} is one of "\code{unconstrained}" or "\code{constrained}", otherwise defaults to 3.5.
#' @param psi.beta The rate of the inverse gamma prior on the uniquenesses. Can be either a single parameter, a vector of variable specific rates, or a matrix of variable and cluster-specific rates. If this is not supplied, \code{\link{psi_hyper}} is invoked to choose sensible values, depending on the value of \code{uni.prior} and, for the "\code{MFA}" and "\code{MIFA}" models, the value of \code{psi0g}.
#' @param mu.zero The mean of the prior distribution for the mean parameter. Defaults to the sample mean of the data.
#' @param sigma.mu The covariance of the prior distribution for the mean parameter. Always assumed to be a diagonal matrix. Can be a scalar times the identity or a vector of appropriate dimension. If supplied as a matrix, only the diagonal elements will be extracted. Defaults to the diagonal entries of the sample covariance matrix.
#' @param sigma.l The covariance of the prior distribution for the loadings. Defaults to 1. Only relevant for the finite factor methods.
#' @param z.init The method used to initialise the cluster labels. Defaults to \code{\link[mclust]{Mclust}}. Other options include \code{kmeans}, hierarchical clustering via \code{\link[mclust]{hc}}, random initialisation via \code{priors}, and a user-supplied \code{list}. Not relevant for the "\code{FA}" and "\code{"IFA"} methods.
#' @param z.list A user supplied list of cluster labels. Only relevant if \code{z.init == "z.list"}.
#' @param trunc.G The maximum number of allowable and storable clusters if the "\code{IMFA}" or "\code{IMIFA}" method is employed. Defaults to the same value as \code{range.G} (unless \code{N < P}, see \code{range.G} for details) and must be greater than or equal to this value. The number of active clusters to be sampled at each iteration is adaptively truncated, with \code{trunc.G} as an upper limit for storage reasons. Note that large values of \code{trunc.G} may lead to memory capacity issues.
#' @param learn.alpha Logical indicating whether the Dirichlet process / Pitman concentration parameter is to be learned, or remain fixed for the duration of the chain. If being learned, a Ga(a, b) prior is assumed for \code{alpha}; updates take place via Gibbs sampling when \code{discount} is zero and via Metropolis-Hastings otherwise. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods, in which case the default is \code{TRUE}.
#' @param alpha.hyper A vector of length 2 giving hyperparameters for the Dirichlet process / Pitman-Yor concentration parameter \code{alpha}. If \code{isTRUE(learn.alpha)}, these are shape and rate parameter of a Gamma distribution. Defaults to Ga(2, 1). Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods, in which case the default is \code{TRUE}. The prior is shifted to have support on (-\code{discount}, \code{Inf}) when non-zero \code{discount} is supplied or \code{learn.d=TRUE}.
#' @param ind.slice Logical indicitating whether the independent slice-efficient sampler is to be employed. If \code{FALSE} the dependent slice-efficient sampler is employed, whereby the slice sequence xi_1,...,xi_g is equal to the decreasingly ordered mixing proportions. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods. Defaults to \code{TRUE}.
#' @param rho Parameter controlling the rate of geometric decay for the independent slice-efficient sampler, s.t. xi = (1 - rho)rho^(g-1). Must lie in the interval (0, 1]. Higher values are associated with better mixing but longer run times. Defaults to 0.75, but 0.5 is an interesting special case which guarantees that the slice sequence xi_1,...,xi_g is equal to the \emph{expectation} of the decreasingly ordered mixing proportions. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods when \code{ind.slice} is \code{TRUE}.
#' @param IM.lab.sw Logial indicating whether the two forced label switching moves are to be implemented (defaults to \code{TRUE}) when running one of the infinite mixture models, with Dirichlet process or Pitman-Yor process priors. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param discount The discount parameter used when generalising the Dirichlet process to the Pitman-Yor process. Must lie in the interval [0, 1). If non-zero, \code{alpha} can be supplied greater than -discount. Defaults to 0. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param learn.d Logical indicating whether the \code{discount} parameter is to be updated via Metropolis-Hastings. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods, in which case the default is \code{FALSE}.
#' @param d.hyper Hyperparameters for the Beta(a,b) prior on the \code{discount} hyperparameter. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods.
#' @param kappa The spike-and-slab prior distribution on the \code{discount} hyperparameter is assumed to be a mixture with point-mass at zero and a continuous Beta(a,b) distribution. \code{kappa} gives the weight of the point mass at zero (the 'spike'). Must lie in the interval [0,1]. Defaults to 0.5. Only relevant for the "\code{IMFA}" and "\code{IMIFA}" methods when \code{isTRUE(learn.d)}. A value of 0 ensures non-zero discount values (i.e. Pitman-Yor) at all times, and \emph{vice versa}.
#' @param zeta Tuning parameter controlling the acceptance rate of the random-walk proposal for the Metropolis-Hastings steps when \code{learn.alpha=TRUE}, where \code{2 * zeta} gives the full width of the uniform proposal distribution. These steps are only invoked when either \code{discount} is non-zero and fixed or \code{learn.d=TRUE}, otherwise \code{alpha} is learned by Gibbs updates. Must be strictly positive. Defauts to 2.
#' @param tune.zeta Used for tuning \code{zeta} & the width of the uniform proposal for \code{alpha} via diminishing Robbins-Monro type adaptation, when that parameter is learned via Metropolis-Hastings. Must be given in the form of a list with the following \emph{named} elements:
#' \describe{
#' \item{"\code{heat}"}{The initial adaptation intensity/step-size, such that larger values lead to larger updates. Must be strictly greater than zero. Defaults to 1 if not supplied but other elements of \code{tune.zeta} are.}
#' \item{"\code{lambda}"}{Iteration rescaling parameter which controls the speed at which adaptation diminishes, such that lower values cause the contribution of later iterations to diminish more slowly. Must lie in the interval (0.5, 1]. Defaults to 1 if not supplied but other elements of \code{tune.zeta} are.}
#' \item{"\code{target}"}{The target acceptance rate. Must lie in the interval [0, 1]. Defaults to 0.441, which is optimum for univariate targets, if not supplied but other elements of \code{tune.zeta} are.}
#' }
#'  \code{tune.zeta} is only relevant when \code{isTRUE(learn.alpha)} under the "\code{IMFA}" or "\code{IMIFA}" models, and either the \code{discount} remains fixed at a non-zero value, or when \code{isTRUE(learn.d)} and \code{kappa < 1}. Since Gibbs steps are invoked for updated \code{alpha} when \code{discount == 0}, adaption occurs according to a running count of the number of iterations with non-zero sampled \code{discount} values. If diminishing adaptation is invoked, the posterior mean \code{zeta} will be stored. Since caution is advised when employing adaptation, note that acceptance rates of between 10-50\% are generally considered adequate.
#' @param equal.pro Logical variable indicating whether or not the mixing mixing proportions are to be equal across clusters in the model (default = \code{FALSE}).
#' @param mu0g Logical indicating whether the \code{mu.zero} hyperparameter can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the "\code{MFA}" and "\code{MIFA}" methods when \code{z.list} is supplied.
#' @param psi0g Logical indicating whether the \code{psi.beta} hyperparameter(s) can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the "\code{MFA}" and "\code{MIFA}" methods when \code{z.list} is supplied, and only allowable when \code{uni.type} is one of \code{unconstrained} or \code{isotropic}.
#' @param verbose Logical indicating whether to print output (e.g. run times) and a progress bar to the screen while the sampler runs. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#'
#' @details Creates a raw object of class "\code{IMIFA}" from which the optimal/modal model can be extracted by \code{\link{get_IMIFA_results}}. Dedicated \code{print} and \code{summary} functions exist for objects of class "\code{IMIFA}".
#'
#' @note Further control over the specification of advanced function arguments can be obtained with recourse to the following functions:\cr
#' \itemize{
#' \item{\strong{\code{\link{mgpControl}}} - }{Supply arguments (with defaults) pertaining to the multiplicative gamma process (MGP) shrinkage prior and adaptive Gibbs sampler (AGS). For use with the infinite factor models "\code{IFA}", "\code{MIFA}", "\code{OMIFA}", and "\code{IMIFA}" only.}
#' \item{\strong{\code{\link{storeControl}}} - }{Supply logical indicators governing storage of parameters of interest for all models in the IMIFA family. It may be useful not to store certain parameters if memory is an issue (e.g. for large data sets or for a large number of MCMC iterations after burnin and thinning).}
#' }
#'
#' @return A list of lists of lists of class "\code{IMIFA}" to be passed to \code{\link{get_IMIFA_results}}. If the returned object is x, candidate models are accesible via subsetting, where x is of the form x[[1:length(range.G)]][[1:length(range.Q)]]. However, these objects of class "IMIFA" should rarely if ever be manipulated by hand - use of the \code{\link{get_IMIFA_results}} function is \emph{strongly} advised.
#' @export
#' @importFrom stats "acf" "complete.cases" "cor" "cov" "dbeta" "density" "dgamma" "dnorm" "factanal" "kmeans" "plogis" "pnorm" "qlogis" "quantile" "rbeta" "rexp" "rgamma" "rmultinom" "rnorm" "runif" "setNames"
#' @importFrom utils "capture.output" "setTxtProgressBar" "txtProgressBar" "memory.limit"
#' @importFrom matrixStats "rowLogSumExps"
#' @importFrom Rfast "Order" "colVars" "rowmeans" "standardise" "sort_unique" "cora" "cova" "Round" "groupcolVars"
#' @importFrom e1071 "matchClasses"
#' @importFrom mvnfast "dmvn"
#' @importFrom slam "as.simple_sparse_array" "as.simple_triplet_matrix"
#' @importFrom mclust "Mclust" "mclustBIC" "hc" "hclass"
#'
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link{psi_hyper}}, \code{\link{mgpControl}}, \code{\link{storeControl}}
#' @references
#' Murphy, K., Gormley, I. C. and Viroli, C. (2017) Infinite Mixtures of Infinite Factor Analysers: Nonparametric Model-Based Clustering via Latent Gaussian Models, <\href{https://arxiv.org/abs/1701.07010}{arXiv:1701.07010}>.
#'
#' Bhattacharya, A. and Dunson, D. B. (2011) Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291-306.
#'
#' Kalli, M., Griffin, J. E. and Walker, S. G. (2011) Slice sampling mixture models, \emph{Statistics and Computing}, 21(1): 93-105.
#'
#' McNicholas, P. D. and Murphy, T. B. (2008) Parsimonious Gaussian Mixture Models, \emph{Statistics and Computing}, 18(3): 285-296.
#'
#' Rousseau, J. and Mengersen, K. (2011) Asymptotic Behaviour of the posterior distribution in overfitted mixture models, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 73(5): 689-710.
#'
#' Tipping, M. E. and Bishop, C. M. (1999). Probabilistic principal component analysis, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 61(3): 611-622.
#'
#' @author Keefe Murphy - \href{keefe.murphy@ucd.ie}{<keefe.murphy@ucd.ie>}
#'
#' @examples
#' # data(olive)
#' # data(coffee)
#'
#' # Fit an IMIFA model to the olive data. Accept all defaults.
#' # simIMIFA <- mcmc_IMIFA(olive, method="IMIFA")
#' # summary(simIMIFA)
#'
#' # Fit an IMIFA model assuming a Pitman-Yor prior.
#' # Allow the alpha and discount parameter to be learned.
#' # Control the balance between the DP and PY priors using the kappa parameter.
#' # simPY    <- mcmc_IMIFA(olive, method="IMIFA", learn.d=TRUE, kappa=0.5)
#' # summary(simPY)
#'
#' # Fit a MFA model to the scaled olive data, with isotropic uniquenesses (i.e. MPPCA).
#' # Allow diagonal covariance as a special case where range.Q = 0.
#' # Don't store the scores. Accept all other defaults.
#' # simMFA   <- mcmc_IMIFA(olive, method="MFA", n.iters=10000, range.G=3:6,
#' #                        range.Q=0:3, storage=storeControl(store.switch=FALSE),
#' #                        centering=FALSE, uni.type="isotropic")
#'
#' # Fit a MIFA model to the centered & scaled coffee data, w/ cluster labels initialised by K-Means.
#' # Note that range.Q doesn't need to be specified. Allow IFA as a special case where range.G=1.
#' # simMIFA  <- mcmc_IMIFA(coffee, method="MIFA", n.iters=10000, range.G=1:3, z.init="kmeans")
#'
#' # Fit an IFA model to the centered and pareto scaled olive data.
#' # Note that range.G doesn't need to be specified. We can optionally supply a range.Q starting value.
#' # Enforce additional shrinkage using alpha.d1, alpha.d2, prop, and epsilon via mgpControl().
#' # simIFA   <- mcmc_IMIFA(olive, method="IFA", n.iters=10000, range.Q=4, scaling="pareto",
#' #                        MGP=mgpControl(alpha.d1=3.5, alpha.d2=7, prop=0.6, epsilon=0.12))
#'
#' # Fit an OMIFA model to the centered & scaled coffee data.
#' # Supply a sufficiently small alpha value. Try varying other hyperparameters.
#' # Accept the default value for the starting number of factors,
#' # but supply a value for the starting number of clusters.
#' # Try contraining uniquenesses to be common across both variables and clusters.
#' # simOMIFA <- mcmc_IMIFA(coffee, method="OMIFA", range.G=10, psi.alpha=3,
#' #                        MGP=mgpControl(nu=3), alpha=0.8, uni.type="single")
mcmc_IMIFA  <- function(dat = NULL, method = c("IMIFA", "IMFA", "OMIFA", "OMFA", "MIFA", "MFA", "IFA", "FA", "classify"), n.iters = 25000L, range.G = NULL, range.Q = NULL, burnin = n.iters/5, thinning = 2L, MGP = mgpControl(),
                        centering = TRUE, scaling = c("unit", "pareto", "none"), uni.type = c("unconstrained", "isotropic", "constrained", "single"), uni.prior = c("unconstrained", "isotropic"), alpha = NULL, storage = storeControl(),
                        psi.alpha = NULL, psi.beta = NULL, mu.zero = NULL, sigma.mu = NULL, sigma.l = NULL, z.init = c("mclust", "hc", "kmeans", "list", "priors"), z.list = NULL, trunc.G = NULL, learn.alpha = TRUE, alpha.hyper = NULL,
                        ind.slice = TRUE, rho = NULL, IM.lab.sw = TRUE, discount = NULL, learn.d = FALSE, d.hyper = NULL, kappa = NULL, zeta = NULL, tune.zeta = NULL, equal.pro = FALSE, mu0g = FALSE, psi0g = FALSE, verbose = interactive()) {

  call      <- match.call()
  defopt    <- options()
  options(warn=1)
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
  if(!missing(method) && method == "classification") {
    method  <- "classify"
  }
  if(!is.character(method))         stop("'method' must be a character vector of length 1")
  if(!is.character(scaling))        stop("'scaling' must be a character vector of length 1")
  method    <- match.arg(method)
  scaling   <- match.arg(scaling)
  if(missing(dat))                  stop("Dataset must be supplied")
  dat.nam   <- gsub("[[:space:]]", "", deparse(substitute(dat)))
  nam.dat   <- gsub("\\[.*", "", dat.nam)
  if(!exists(nam.dat,
     envir=.GlobalEnv))             stop(paste0("Object ", match.call()$dat, " not found\n"))
  if(any(!is.logical(centering),
         length(centering) != 1))   stop("'centering' must be a single logical indicator")
  if(any(!is.logical(verbose),
         length(verbose)   != 1))   stop("'verbose' must be a single logical indicator")

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
  if(length(iters)   < 10)          stop("Run a longer chain!")
  raw.dat   <- as.data.frame(dat)
  num.check <- vapply(raw.dat, is.numeric, logical(1L))
  if(anyNA(raw.dat)) {
    if(verbose)                     message("Rows with missing values removed from data")
    raw.dat <- raw.dat[complete.cases(raw.dat),, drop=FALSE]
  }
  if(sum(num.check) != ncol(raw.dat)) {
    if(verbose)                     message("Non-numeric columns removed from data")
    raw.dat <- raw.dat[,num.check, drop=FALSE]
  }
  if(any(dim(raw.dat) == 0))        stop("Empty data set after removal of ineligble rows/columns")
  if(method != "classify") {
    scal    <- switch(scaling, none=FALSE, Rfast::colVars(as.matrix(raw.dat), std=TRUE))
    scal    <- switch(scaling, pareto=sqrt(scal), scal)
    dat     <- if(is.logical(scal)) { if(any(centering, scal)) standardise(as.matrix(raw.dat), center=centering, scale=scal) else as.matrix(raw.dat) } else scale(raw.dat, center=centering, scale=scal)
  } else   {
    dat     <- as.matrix(raw.dat)
  }
  centered  <- switch(method, classify=all(Round(colSums(dat)) == 0), any(centering, all(Round(colSums(dat)) == 0)))
  if(!any(centered, centering) && all(verbose,
    scaling != "none"))             message("Are you sure you want to apply scaling without centering?")

  N         <- as.integer(nrow(dat))
  P         <- as.integer(ncol(dat))
  uni       <- P == 1
  lnN       <- log(N)
  NlP       <- N  < P
  miss.uni  <- missing(uni.type)
  miss.pri  <- missing(uni.prior)
  if(!is.character(uni.type))       stop("'uni.type' must be a character vector of length 1")
  if(!is.character(uni.prior))      stop("'uni.prior' must be a character vector of length 1")
  uni.type  <- ifelse(miss.uni, ifelse(uni, "isotropic", "unconstrained"), match.arg(uni.type))
  if(uni) {
    if(is.element(uni.type, c("unconstrained", "constrained"))) {
      uni.type     <- switch(uni.type, unconstrained=, isotropic="isotropic", constrained=, single="single")
      if(verbose)                   message(paste0("'uni.type' coerced to ", uni.type, " as the dataset is univariate"))
    }
  }
  uni.prior <- ifelse(miss.pri, switch(uni.type, constrained=, unconstrained="unconstrained", "isotropic"), match.arg(uni.prior))
  if(uni    && all(uni.prior  == "unconstrained", verbose)) {
    uni.prior     <- "isotropic";   message("'uni.prior' coerced to isotropic as the dataset is univariate")
  }
  if(all(is.element(uni.type, c("isotropic", "single")),
     uni.prior == "unconstrained")) stop("'uni.prior' can only be 'unconstrained' when 'uni.type' is 'unconstrained' or 'constrained'")
  if(all(is.element(uni.type, c("unconstrained", "constrained")), any(miss.uni, miss.pri),
         NlP, verbose))             message(paste0("Consider setting 'uni.type' to 'isotropic' or 'single'", ifelse(miss.pri && uni.prior == "unconstrained", ", or at least 'uni.prior' to 'isotropic', ", " "), "in N << P cases"))

# Manage storage switches & warnings for other function inputs
  store.x   <- attr(storage, "Missing")
  if(!store.x["mu.sw"] && all(!storage["mu.sw"], ifelse(method == "classify",
     !centering, !centered)))       warning("Centering hasn't been applied - are you sure you want mu.switch=FALSE?", call.=FALSE)
  if(N < 2)                         stop("Must have more than one observation")
  zin.miss  <- missing(z.init)
  zli.miss  <- missing(z.list)
  if(!is.character(z.init))         stop("'z.init' must be a character vector of length 1")
  z.init    <- match.arg(z.init)
  if(all(is.element(method, c("FA", "IFA")), verbose,
         !zin.miss || !zli.miss)) { message(paste0("z does not need to be initialised for the ", method, " method"))
  } else if(!zli.miss) {
    if(!is.list(z.list))     z.list        <- lapply(list(z.list), as.factor)
    if(zin.miss &&
       z.init   != "list") { z.init        <- "list"
       if(verbose)                   message("'z.init' set to 'list' as 'z.list' was supplied")
    }
  }

  G.x       <- missing(range.G)
  alpha.x   <- missing(learn.alpha)
  disc.x    <- missing(learn.d)
  if(any(!is.logical(learn.alpha),
         length(learn.alpha) != 1)) stop("'learn.alpha' must be a single logical indicator")
  if(any(!is.logical(learn.d),
         length(learn.d)     != 1)) stop("'learn.d' must be a single logical indicator")
  if(missing(d.hyper))       d.hyper       <- c(1L, 1L)
  if(length(d.hyper)         != 2)  stop("d.hyper' must be a vector of length 2")
  if(any(d.hyper   <= 0))           stop("'Discount Beta prior hyperparameters must be strictly positive")
  kappa            <- ifelse(all(!learn.d, discount == 0), 1, ifelse(missing(kappa), 0.5, kappa))
  if(any(!is.numeric(kappa),
         length(kappa)       != 1)) stop("'kappa' must be a single number")
  if(kappa   <  0  || kappa   > 1)  stop("'kappa' must lie in the interval [0, 1]")
  if(kappa  ==  0) {
   if(all(!learn.d, discount == 0)) stop("'kappa' is zero and yet 'discount' is fixed at zero:\neither learn the discount parameter or specify a non-zero value")
  } else if(kappa  == 1) {
   if(any(learn.d,  discount != 0)) stop(paste0("'kappa' is exactly 1 and yet", ifelse(learn.d, " 'discount' is being learned ", if(discount != 0) " the discount is fixed at a non-zero value"), ":\nthe discount should remain fixed at zero"))
  }
  discount         <- switch(method, IMFA=, IMIFA=ifelse(missing(discount), ifelse(learn.d, ifelse(kappa != 0 && runif(1) <= kappa, 0, rbeta(1, d.hyper[1], d.hyper[2])), 0), discount), 0)
  if(any(!is.numeric(discount),
         length(discount)    != 1)) stop("'discount' must be a single number")
  if(discount       < 0      ||
     discount      >= 1)            stop("'discount' must lie in the interval [0, 1)")
  if(!is.element(method, c("IMFA", "IMIFA"))) {
    if(learn.alpha) {
      learn.alpha  <- FALSE
      if(verbose   && !alpha.x)     message(paste0("'learn.alpha' forced to FALSE for the ", method, " method"))
    }
    if(learn.d)     {
      learn.d      <- FALSE
      if(verbose   && !disc.x)      message(paste0("'learn.d' forced to FALSE for the ",     method, " method"))
    }
  }

  if(!is.element(method, c("MFA", "MIFA")))      {
    if(length(range.G) > 1)         stop(paste0("Only one 'range.G' value can be specified for the ", method, " method"))
    if(all(!G.x, verbose, is.element(method, c("FA", "IFA"))) &&
       range.G  > 1)                message(paste0("'range.G' forced to 1 for the ", method, " method"))
    if(is.element(method, c("OMIFA", "OMFA", "IMFA", "IMIFA"))) {
      lnN2         <- ceiling(lnN)
      tmp.G        <- as.integer(min(N - 1, max(25, ceiling(3 * lnN))))
      if(G.x)   {
        range.G    <- G.init <- tmp.G
        if(NlP) {     if(verbose)   message(paste0("Since N << P, the sampler will be initialised with a different default of ceiling(log(N)) = ", lnN2, " clusters (unless 'range.G' is supplied)"))
          G.init   <- max(2, lnN2)
        }
      }
      if(!G.x)  {
        G.init     <- range.G
        if(NlP) {
          range.G  <- tmp.G
        }
      }
      if(range.G    < lnN2)         warning(paste0("'range.G' should be at least log(N) (=log(", N, "))", " for the ", method, " method"), call.=FALSE)
      if(is.element(method, c("IMFA", "IMIFA"))) {
        if(any(!is.logical(ind.slice),
           length(ind.slice) != 1)) stop("'ind.slice' must be a single logical indicator")
        x.IM.sw    <- missing(IM.lab.sw)
        if(any(!is.logical(IM.lab.sw),
           length(IM.lab.sw) != 1)) stop("'IM.lab.sw' must be a single logical indicator")
        if(all(!x.IM.sw, isTRUE(IM.lab.sw), !isTRUE(learn.alpha), verbose,
               !isTRUE(learn.d)))   message("May not be necessary to set 'IM.lab.sw' to TRUE when neither 'alpha' nor 'discount' are being learned")
        if(missing(rho)) {
          rho      <- 0.75
        }
        if(all(length(rho) > 1,
           rho > 1 || rho <= 0))    stop("'rho' must be a single number in the interval (0, 1]")
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
               length(zeta)  != 1,
               zeta < 0))           stop(paste0("'zeta' must be single strictly positive number"))
        if(all(length(alpha.hyper)  != 2,
           learn.alpha))            stop(paste0("'alpha.hyper' must be a vector of length 2, giving the shape and rate hyperparameters of the gamma prior for alpha when 'learn.alpha' is TRUE"))
        a.hyp1     <- alpha.hyper[1]
        a.hyp2     <- alpha.hyper[2]
        if(learn.alpha)   {
          if(a.hyp1   <= 0)         stop("The shape of the gamma prior for alpha must be strictly positive")
          if(a.hyp2   <= 0)         stop("The rate of the gamma prior for alpha must be strictly positive")
        }
        if((t.miss <- missing(trunc.G))) {
          trunc.G  <- range.G
        }
        if(length(trunc.G) > 1)     stop("'trunc.G' must be a single number")
        if(all(ifelse(N > 50, trunc.G < 50,
           trunc.G  < N), !t.miss)) message(paste0("Consider setting 'trunc.G' to min(N=", N, ", 50) unless practical reasons in heavy computational/memory burden cases prohibit it"))
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
      if(all(verbose, !G.x)  && any(length(range.G > 1),
          range.G  != levs))  {     message("Forced 'range.G' equal to the number of levels in 'zlabels' for the 'classify' method")
      }
     G.init <- range.G    <- levs
    } else {
     G.init <- range.G    <- 1L
     storage["pi.sw"]     <- FALSE
    }
    equal.pro      <- FALSE
    meth    <- method
  } else {
    alp3    <- 3L   * alpha
    if(length(equal.pro) > 1 ||
       !is.logical(equal.pro))      stop("'equal.pro' must be a single logical indicator")
    if(storage["pi.sw"]   && equal.pro) {
      if(all(!store.x["pi.sw"],
             verbose))              message("Forced non-storage of mixing proportions as 'equal.pro' is TRUE: only posterior mean estimates will be available")
      storage["pi.sw"]    <- FALSE
    }
    if(G.x)                         stop("'range.G' must be specified")
    if(any(range.G  < 1))           stop("'range.G' must be strictly positive")
    if(any(range.G  > alp3 * lnN))  warning(paste0("'range.G' MUCH greater than log(N) (=log(", N, ")):\nEmpty clusters are likely, consider running an overfitted or infinite mixture"), call.=FALSE)
    range.G <- G.init     <- sort_unique(range.G)
    meth    <- rep(method, length(range.G))
    if(any(!is.logical(mu0g),
           length(mu0g)   != 1))    stop("'mu0g' must be a single logical indicator")
    if(any(!is.logical(psi0g),
           length(psi0g)  != 1))    stop("'psi0g' must be a single logical indicator")
  }
  if(any(range.G >= N))             stop(paste0("'range.G' must be less than the number of observations N=", N))
  if(G.init[1]   == 1)     {
    if(is.element(method, c("IMIFA", "IMFA",
       "OMIFA", "OMFA")))  {        stop(paste0("'method' should be ", switch(method, IMFA=, OMFA="FA", OMIFA=, IMIFA="IFA"), " for a one group model under the ", method, " method"))
    } else {
      meth[1]    <- switch(method,  MFA=, FA="FA", MIFA=, IFA="IFA")
    }
    if(all(verbose, !is.element(method,
           c("FA", "IFA"))))        message(paste0("Forced use of ", meth[1], " method where 'range.G' is equal to 1"))
  }

# Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  datname          <- rownames(dat)
  if(any(length(unique(datname)) != N,
     is.null(datname)))      rownames(dat) <- seq_len(N)
  mu               <- list(as.matrix(colMeans(dat)))
  mu0.x            <- missing(mu.zero)
  mu.zero          <- if(mu0.x) mu  else     .len_check(mu.zero, mu0g, method, P, G.init)
  sigmu.miss       <- missing(sigma.mu)
  cov.mat          <- if(P > 500) switch(scaling, unit=cora(as.matrix(dat)), cova(as.matrix(dat))) else switch(scaling, unit=cor(dat), cov(dat))
  sigma.mu         <- if(sigmu.miss) diag(cov.mat) else if(is.matrix(sigma.mu)) diag(sigma.mu) else sigma.mu
  sigmu.len        <- length(unique(sigma.mu))
  if(sigmu.len     == 1)     sigma.mu      <- sigma.mu[1]
  if(any(sigma.mu  <= 0, !is.numeric(sigma.mu),
     !is.element(length(sigma.mu),
     c(1, P))))                     stop(paste0("'sigma.mu' must be strictly positive, and of length 1 or P=", P))
  if(missing(psi.alpha))     psi.alpha     <- 2.5 + is.element(uni.type, c("isotropic", "single"))
  if(any(psi.alpha <= 1, !is.numeric(psi.alpha),
     length(psi.alpha)     != 1))   stop("'psi.alpha' must be a single number strictly greater than 1 in order to bound uniquenesses away from zero")
  beta.x           <- missing(psi.beta)
  if(beta.x) {
    psi.beta       <- temp.psi   <- list(psi_hyper(shape=psi.alpha, covar=cov.mat, type=uni.prior))
  } else {
    psi.beta       <- lapply(.len_check(psi.beta, psi0g, method, P, G.init), matrix, nrow=P, ncol=G.init, byrow=length(psi.beta) == G.init)
  }
  Q.miss    <- missing(range.Q)
  Q.min     <- min(ceiling(log(P)), ceiling(log(N)))

  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(!missing(MGP))               message(paste0("'MGP' parameters not necessary for the ", method, " method"))
    if(missing(sigma.l))     sigma.l       <- 1L
    if(any(sigma.l <= 0, !is.numeric(sigma.l),
           length(sigma.l) != 1))   stop("'sigma.l' must be a single strictly positive number")
    if(Q.miss)                      stop("'range.Q' must be specified")
    if(any(range.Q  < 0))           stop(paste0("'range.Q' must be non-negative for the ", method, " method"))
    range.Q <- sort_unique(range.Q)
    delta0g <- FALSE
  } else {
    mgpmiss <- attr(MGP, "Missing")
    fQ0     <- uni || P    == 2
    MGP$prop   <- ifelse(fQ0 && mgpmiss$propx, 0, floor(MGP$prop * P)/P)
    adapt   <- MGP$adapt
    adaptat <- MGP$adaptat <- ifelse(mgpmiss$adaptatx, ifelse(fQ0, 0L, switch(method, IFA=, MIFA=burnin, 0L)), MGP$adaptat)
    nuplus1 <- MGP$nuplus1
    delta0g <- MGP$delta0g
    delta0x <- switch(method, MIFA=delta0g, FALSE)
    alpha.d1   <- .len_check(MGP$alpha.d1, delta0x, method, P, G.init, P.dim=FALSE)
    alpha.d2   <- .len_check(MGP$alpha.d2, delta0x, method, P, G.init, P.dim=FALSE)
    MGP     <- MGP[-c(1:3)]
    if(Q.miss)               range.Q       <- as.integer(ifelse(fQ0, 1, min(ifelse(P > 500, 12 + floor(log(P)), floor(3 * log(P))), N - 1)))
    if(length(range.Q)      > 1)    stop(paste0("Only one starting value for 'range.Q' can be supplied for the ", method, " method"))
    if(range.Q <= 0)                stop(paste0("'range.Q' must be strictly positive for the ", method, " method"))
    if(isTRUE(adapt)) {
     if(adaptat < 0        ||
        adaptat > burnin)           stop("'adapt.at' must be lie in the interval [0, burnin] if 'adapt' is TRUE")
     if(Q.min   > range.Q)          stop(paste0("'range.Q' must be at least min(log(P), log(N)) for the ", method, " method when 'adapt' is TRUE"))
    } else if(!fQ0 && is.element(method,
              c("OMIFA", "IMIFA"))) warning("'adapt=FALSE' is NOT recommended for the 'OMIFA' or 'IMIFA' methods", call.=FALSE)
  }

  len.G     <- switch(method, classify=range.G, length(range.G))
  len.Q     <- length(range.Q)
  len.X     <- len.G * len.Q
  if(all(len.X > 10,
         suppressWarnings(memory.limit())  <= 16256,
         storage["s.sw"]))  {
    if(!store.x["s.sw"])    {       warning(paste0("The large number of candidate models being explored (", len.X, ") could lead to memory issues\nConsider setting 'score.switch' to FALSE or breaking up the task into chunks and calling get_IMIFA_results() on each chunk"), call.=FALSE)
    } else                  {       warning(paste0("'score.switch' set to FALSE as too many candidate models are being explored (", len.X, ")\nPosterior inference on the scores will not be possible, though you can risk forcing storage by supplying score.switch=TRUE\nConsider breaking up the task into chunks and calling get_IMIFA_results() on each chunk"), call.=FALSE)
      storage["s.sw"]    <- FALSE
    }
  }
  Q.warn       <- min(N - 1, Ledermann(P))
  if(any(range.Q > Q.warn))   {
    if(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")) &&
       isTRUE(adapt))         {     warning(paste0("Starting value for number of factors is greater than ", ifelse(any(range.Q > P), paste0("the number of variables (", P, ")"), paste0("the suggested Ledermann upper bound (", Q.warn, ")"))), call.=FALSE)
    } else if(any(is.element(method, c("FA", "MFA", "OMFA", "IMFA")),
                  is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")) &&
                  isTRUE(!adapt)))  warning(paste0("Number of factors is greater than ", ifelse(any(range.Q > P), paste0("the number of variables (", P, ")"), paste0("the suggested Ledermann upper bound (", Q.warn, ")"))), call.=FALSE)
  }
  if(verbose   && any(all(method == "MFA",  any(range.G > 1)) && any(range.Q > 0),
                      all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA",
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!storage["l.sw"],
           !storage["psi.sw"]))   {
                                    message("Loadings & Uniquenesses not stored: will be unable to estimate covariance matrices and compute error metrics")
    } else if(!storage["l.sw"])   { message("Loadings not stored: will be unable to estimate covariance matrices and compute error metrics")
    } else if(!storage["psi.sw"])   message("Uniquenesses not stored: will be unable to estimate covariance matrices and compute error metrics")
  }
  if(all(storage["s.sw"], !storage["l.sw"],
     any(range.Q   != 0)))          message("Loadings not stored but scores are: Procrustes rotation of scores will not occur when passing results to get_IMIFA_results()")
  if(verbose   && is.element(method,
               c("FA", "IFA")))   {
    if(all(!storage["mu.sw"],
           !storage["psi.sw"]))   {
                                    message("Means & Uniquenesses not stored, but posterior mean estimates will still be available")
    } else if(!storage["mu.sw"])  { message("Means not stored, but posterior mean estimates will still be available")
    } else if(!storage["psi.sw"])   message("Uniquenesses not stored, but posterior mean estimates will still be available")
  }
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")) && any(range.Q == 0)) {
    if(all(storage[c("s.sw", "l.sw")]))   {
                                    message("Scores & Loadings not stored where 'range.Q=0' as model has zero factors", call)
    } else if(storage["s.sw"])    { message("Scores not stored where 'range.Q==0' as model has zero factors")
    } else if(storage["l.sw"])    { message("Loadings not stored where 'range.Q==0' as model has zero factors")
    }
    if(all(range.Q == 0))    storage[c("s.sw", "l.sw")] <- FALSE
  }

  sw0gs     <- c(mu0g = mu0g, psi0g = psi0g, delta0g = delta0g)
  if(is.element(method, c("MFA", "MIFA", "classify"))) {
    if(all(z.init != "list", any(sw0gs))) {
      if(isTRUE(delta0g))           stop("'delta0g' can only be TRUE if z.init='list' for the 'MIFA' method\n")
      if(all(!mu0.x,  mu0g))        stop("'mu.zero' can only be supplied for each cluster if z.init='list' for the 'MFA' & 'MIFA' methods")
      if(all(!beta.x, psi0g))       stop("'psi.beta' can only be supplied for each cluster if z.init='list' for the 'MFA' & 'MIFA' methods")
    }
    if(all(is.element(uni.type, c("constrained", "single")),
           isTRUE(psi0g)))  {       warning(paste0("'psi0g' forced to FALSE as uniquenesses are constrained across clusters (i.e. 'uni.type' = ", uni.type, ")"), call.=FALSE)
      psi0G <- FALSE
    }
    if(all(method == "MFA",
           isTRUE(delta0g)))        message("'delta0g' cannot be TRUE for the 'MFA' method")
    if(method == "classify") mu0g          <- TRUE
  } else if(any(sw0gs))             stop(paste0("'", names(which(sw0gs)), "' should be FALSE for the ", method, " method\n"))

  dimension <- PGMM_dfree(P=P, Q=switch(method, FA=, classify=, MFA=, OMFA=, IMFA=range.Q,
               IFA=, MIFA=, OMIFA=, IMIFA=0L), equal.pro=equal.pro, method=switch(uni.type, unconstrained="UUU", isotropic="UUC", constrained="UCU", single="UCC"))
  min.d2    <- min(dimension)/2
  min.d2G   <- min.d2 * G.init
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    if(missing(alpha))     { alpha         <- switch(method, MFA=, MIFA=1, OMFA=, OMIFA=0.5/G.init, ifelse(learn.alpha, max(1, rgamma(1, a.hyp1, a.hyp2)), 1))
    } else if(is.element(method,
      c("OMFA", "OMIFA")))   alpha         <- alpha/G.init
    if(length(alpha) != 1)          stop("'alpha' must be specified as a scalar to ensure an exchangeable prior")
    if(alpha <= -discount)          stop(paste0("'alpha' must be ",     ifelse(discount != 0, paste0("strictly greater than -discount (i.e. > ", - discount, ")"), "strictly positive")))
    if(all(is.element(method,  c("IMIFA",   "IMFA")),
       !learn.alpha))               warning(paste0("'alpha' fixed at ", Round(alpha, options()$digits), " as it's not being learned via Gibbs/Metropolis-Hastings updates"), call.=FALSE)
    if(all(is.element(method,  c("OMIFA",   "OMFA")),
       alpha >= min.d2))            warning(paste0("'alpha' over 'range.G' for the OMFA & OMIFA methods must be less than half the dimension (per cluster!)\nof the free parameters of the smallest model considered (= ", min.d2, "): consider suppling 'alpha' < ", min.d2G), call.=FALSE)
    if(any(all(is.element(method, c("MFA",  "MIFA")), alpha > 1),
           all(is.element(method, c("OMFA", "OMIFA")),
           alpha > 1/G.init)))      warning("Are you sure alpha should be greater than 1?", call.=FALSE)
    if(all(any(is.nan(rDirichlet(G=min(range.G), alpha = alpha))),
           is.element(method,  c("MFA",     "MIFA",
            "OMFA", "OMIFA"))))     stop("'alpha' is too small for simulated mixing proportions to be valid")
    if(all(is.element(method,  c("OMIFA",   "OMFA")), is.element(z.init,
       "priors")))                  stop(paste0("'z.init' cannot be set to 'priors' for the ", method, " method to ensure all clusters are populated at the initialisation stage"))
    if(!zli.miss) {
      if(length(z.list)   != len.G)  {
                                    stop(paste0("'z.list' must be a list of length ", len.G))  }
                             list.levels   <- lapply(z.list, nlevels)
      if(!all(list.levels == G.init))             {
        if(!is.element(method, c("IMIFA",
                       "IMFA", "OMIFA", "OMFA"))) {
                                    stop(paste0("Each element of 'z.list' must have the same number of levels as 'range.G'"))
        } else                      stop(paste0("Only ", list.levels, " clusters are populated according to z.list, but 'range.G' has been set to ", G.init, ":\nReset 'range.G' to this value to avoid redunandtly carrying around empty clusters or supply a list with ", G.init, " levels"))
      }
      if(!all(lengths(z.list) == N)) {
                                    stop(paste0("Each element of 'z.list' must be a vector of length N=", N)) }
    }
    if(all(zli.miss, z.init   == "list"))         {
                                    stop(paste0("'z.list' must be supplied if 'z.init' is set to 'list'")) }
  }

  imifa     <- list(list())
  Gi        <- Qi  <- 1L
  gibbs.arg <- list(P = P, sigma.mu = sigma.mu, psi.alpha = psi.alpha, burnin = burnin, sw = storage,
                    thinning = thinning, iters = iters, verbose = verbose, uni.type = uni.type, uni.prior = uni.prior)
  if(is.element(method, c("IMIFA", "IMFA"))) {
    gibbs.def      <- all(kappa     < 1,  learn.d, learn.alpha)
    gibbs.may      <- gibbs.def    &&     kappa > 0
    non.py         <- all(discount == 0, !learn.d, learn.alpha)
    def.py         <- all(discount != 0, !learn.d, learn.alpha)
    if(missing(tune.zeta)) {
      tune.zeta    <- list(heat = 0, lambda = NULL, target = NULL, do = FALSE)
    } else   {
      tz    <- tune.zeta
      if(!is.list(tz) || (!all(is.element(names(tz), c("heat", "lambda", "target"))) ||
         !all(lengths(tz) == 1)))   stop("'tune.zeta' must be a list with named elements 'heat', 'lambda' and 'target', all of length 1")
      if(is.null(tz$heat))   tz$heat   <- 1
      if(is.null(tz$lambda)) tz$lambda <- 1
      if(is.null(tz$target)) tz$target <- 0.441
      if(!all(vapply(tz,
         is.numeric, logical(1L)))) stop("Not all elements of 'tune.zeta' are numeric")
      tz$do <- any(gibbs.def, def.py)
      if(tz$heat    < 0)            stop("Invalid 'heat': must be >= 0")
      if(tz$target  < 0   ||
         tz$target  > 1)            stop("Invalid 'target': must lie in the interval [0, 1]")
      if(tz$lambda <= 0.5 ||
         tz$lambda  > 1)            stop("Invalid 'lambda': must lie in the interval (0.5, 1]")
      if(learn.alpha)  {
        if(tz$heat  > 0)   {
         if(isTRUE(gibbs.may))      warning("Are you sure you want to tune zeta?: Gibbs updates are possible as 'kappa' is between zero and 1", call.=FALSE)
         if(non.py && verbose)      message("Tuning zeta will have no effect as 'discount' is fixed at 0")
       } else if(tz$do)             warning("'heat' of 0 corresponds to no tuning: are you sure?", call.=FALSE)
      }
      tune.zeta    <- tz
    }
    gibbs.arg      <- append(gibbs.arg, list(trunc.G = trunc.G, rho = rho, ind.slice = ind.slice, learn.alpha = learn.alpha, learn.d = learn.d, kappa = kappa,
                                             zeta = zeta, IM.lab.sw = IM.lab.sw, a.hyper = alpha.hyper, discount = discount, d.hyper = d.hyper, tune.zeta = tune.zeta))
  }
  if(any(is.element(meth, c("FA",  "IFA")))) {
    gibbs.arg      <- append(gibbs.arg, list(scaling   = scaling))
  }
  if(is.element(method, c("MFA", "MIFA")))   {
    gibbs.arg      <- append(gibbs.arg, list(equal.pro = equal.pro))
  }
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    gibbs.arg      <- temp.args <- append(gibbs.arg, MGP)
  } else {
    gibbs.arg      <- append(gibbs.arg, list(sigma.l = sigma.l))
  }

  init.start       <- proc.time()
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
      } else if(z.init  == "hc")     {
        h.res      <- try(hc(dat, minclus=G), silent=TRUE)
        if(!inherits(h.res, "try-error"))  {
          zi[[g]]  <- as.integer(factor(hclass(h.res, G),     levels=seq_len(G)))
        } else                      stop("Cannot initialise cluster labels using hierarchical clustering. Try another z.init method")
      } else if(z.init  == "mclust") {
        m.res      <- try(Mclust(dat, G, verbose=FALSE), silent=TRUE)
        if(!inherits(m.res, "try-error"))  {
          zi[[g]]  <- as.integer(factor(m.res$classification, levels=seq_len(G)))
        } else                      stop("Cannot initialise cluster labels using mclust. Try another z.init method")
      } else {
        zips       <- rep(1, N)
        if(!is.element(method, c("IMFA", "IMIFA"))) {
          iter     <- 0
          while(all(length(unique(zips)) != G, iter < 100,
                any(prop.table(tabulate(zips, nbins=G)) < 1/G^2))) {
            pies   <- rDirichlet(alpha=alpha, G)
            zips   <- .sim_z_p(N=N, prob.z=pies)
            iter   <- iter + 1
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
      pi.prop[[g]] <- if(equal.pro) rep(1/G, G) else prop.table(nngs)
      mu[[g]]      <- vapply(seq_len(G), function(gg) if(nngs[gg] > 0) colMeans(dat[zi[[g]] == gg,, drop=FALSE]) else rep(0, P), numeric(P))
      mu[[g]]      <- if(uni)       t(mu[[g]])  else mu[[g]]
      if(mu0.x)   {
        mu.zero[[g]]    <- if(mu0g) mu[[g]]     else vapply(seq_len(G), function(gg) colMeans(dat), numeric(P))
      }
      mu.zero[[g]] <- if(uni && !mu0g) t(mu.zero[[g]]) else mu.zero[[g]]
      if(beta.x)  {
        if(psi0g) {
          cov.gg   <- lapply(seq_len(G), function(gg, dat.gg = dat[zi[[g]] == gg,, drop=FALSE]) if(all(nngs[gg] > 1, P <= nngs[g])) { if(P > 500) cova(as.matrix(dat.gg)) else cov(dat.gg) } else cov.mat)
          psi.beta[[g]] <- vapply(seq_len(G), function(gg) psi_hyper(shape=psi.alpha, covar=cov.gg[[gg]], type=uni.prior), numeric(P))
        } else {
          psi.beta[[g]] <- replicate(G, temp.psi[[1]])
        }
        psi.beta[[g]]   <- if(uni)  t(psi.beta[[g]])   else psi.beta[[g]]
      }
      if(is.element(method, c("MIFA", "OMIFA", "IMIFA")) && mgpmiss$ad1x) {
        alpha.d1[[g]]   <- rep(unlist(alpha.d1), G)
      }
      if(is.element(method, c("MIFA", "OMIFA", "IMIFA")) && mgpmiss$ad2x) {
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
    mu.zero        <- if(all(lengths(mu.zero)  == 1)) list(mu.zero[[1]])  else list(mu.zero[[1]][,1,  drop=FALSE])
    psi.beta       <- if(all(lengths(psi.beta) == 1)) list(psi.beta[[1]]) else list(psi.beta[[1]][,1, drop=FALSE])
    if(!is.element(method, c("OMFA", "IMFA")))   {
      alpha.d1     <- list(alpha.d1[[1]][1])
      alpha.d2     <- list(alpha.d2[[1]][1])
    }
  }
  if(all(Round(vapply(mu.zero, sum, numeric(1L))) == 0)) {
    mu.zero        <- switch(method, classify=base::matrix(0, nr=1, nc=range.G), lapply(mu.zero, function(x) 0))
  }
  if(mu0g && unlist(lapply(mu.zero,  function(x) {
    if(is.matrix(x)) any(apply(x, 1, function(y) {
      length(unique(Round(y, min(nchar(y))))) }) == 1)  else  {
      length(unique(Round(x,
      min(nchar(x))))) == 1 }})))   stop("'mu0g' must be FALSE if 'mu.zero' is not group-specific")
  if(anyNA(unlist(psi.beta))) {
    psi.beta       <- lapply(psi.beta, function(x) replace(x, is.na(x), 0))
  }
  if(all(is.element(uni.type, c("isotropic",   "single")),
         unlist(lapply(psi.beta,     function(x) {
    if(is.matrix(x)) any(apply(x, 2, function(y) {
      length(unique(Round(y, min(nchar(y))))) }) != 1)  else  {
      length(unique(Round(x,
      min(nchar(x))))) != 1 }}))))  stop("'psi.beta' cannot be variable-specific if 'uni.type' is 'isotropic' or 'single'")
  if(is.element(uni.type, c("constrained", "single"))) {
      if(unlist(lapply(psi.beta,     function(x) {
    if(is.matrix(x)) any(apply(x, 1, function(y) {
      length(unique(Round(y, min(nchar(y))))) }) != 1)  else  {
      length(unique(Round(x,
      min(nchar(x))))) != 1 }}))) { stop("'psi.beta' cannot be group-specific if 'uni.type' is 'constrained' or 'single'")
    } else if(isTRUE(psi0g) && is.element(method,
              c("MFA", "MIFA")))    stop("'psi0g' must be FALSE if 'psi.beta' is not group-specific")
  }
  if(any(unlist(psi.beta)   <= 0))  stop("'psi.beta' must be strictly positive")
  if(is.element(method, c("classify", "IFA", "MIFA", "IMIFA", "OMIFA"))) {
    if(!all(MGP_check(ad1=unlist(alpha.d1), ad2=unlist(alpha.d2), Q=unique(range.Q), nu=MGP$nu, bd1=MGP$beta.d1, bd2=MGP$beta.d2,
        plus1=nuplus1)[[1]]$valid)) stop("Invalid shrinkage hyperparameter values will not encourage loadings column removal.\nTry using the MGP_check() function in advance to ensure cumulative shrinkage property holds.")
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
      if(sigmu.miss)   gibbs.arg$sigma.mu  <- diag(if(nrow(tmp.dat) > 1) { if(P > 500) cova(as.matrix(tmp.dat)) else cov(tmp.dat) } else cov.mat)
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
  for(g  in length(imifa)) {
   for(q in length(imifa[[g]])) {
     imifa[[g]][[q]]$time <- NULL
   }
  }

  imifa     <- switch(method, FA=, MFA=, OMFA=, IMFA={
     lapply(seq_along(imifa), function(x) setNames(imifa[[x]], paste0(range.Q, ifelse(range.Q == 1, "Factor", "Factors"))))
  }, lapply(seq_along(imifa), function(x) setNames(imifa[[x]], "IFA")))
  gnames    <- switch(method, classify=paste0("Cluster ", seq_len(range.G)), paste0(G.init, ifelse(G.init == 1, "Cluster", "Clusters")))
  names(imifa)            <- gnames
  attr(imifa, "Adapt")    <- is.element(method, c("classify", "IFA", "MIFA", "OMIFA", "IMIFA")) && isTRUE(adapt)
  attr(imifa,
       "Alph.step")       <- learn.alpha
  attr(imifa, "Alpha")    <- if(!learn.alpha) alpha
  attr(imifa,
       "Class.Props")     <- if(method == "classify") tabulate(z.list[[1]], range.G)/N
  attr(imifa, "Call")     <- call
  attr(imifa, "Center")   <- any(centered, centering)
  attr(imifa, "Clusters") <- range.G
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa,
       "Disc.step")       <- learn.d
  attr(imifa, "Discount") <- if(!learn.d) discount
  attr(imifa, "Equal.Pi") <- equal.pro
  attr(imifa, "Factors")  <- range.Q
  attr(imifa, "G.init")   <- G.init
  attr(imifa, "IM.labsw") <- all(is.element(method, c("IMFA", "IMIFA")), IM.lab.sw)
  attr(imifa,
       "Ind.Slice")       <- all(is.element(method, c("IMFA", "IMIFA")), ind.slice)
  attr(imifa, "Init.Z")   <- if(!is.element(method, c("FA", "IFA")))     z.init
  attr(imifa,
       "Label.Switch")    <- if(!is.element(method, c("FA", "IFA")))     any(sw0gs)
  method                  <- names(table(meth)[which.max(table(meth))])
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1, 1)),
                                    substr(method, 2, nchar(method)))
  attr(imifa, "Name")     <- dat.nam
  attr(imifa, "Nuplus1")  <- is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")) && nuplus1
  attr(imifa, "Obs")      <- N
  attr(imifa, "Pitman")   <- all(is.element(method, c("IMFA", "IMIFA")), any(learn.d, discount > 0))
  attr(imifa, "Rho")      <- rho
  attr(imifa, "Scaling")  <- scal
  attr(attr(imifa,
  "Scaling"), "Method")   <- scaling
  attr(imifa, "Store")    <- length(iters)
  storage                 <- c(storage, a.sw = learn.alpha, d.sw = learn.d)
  attr(imifa, "Switch")   <- storage
  times                   <- lapply(list(Total = tot.time, Average = avg.time, Initialisation = init.time), round, 2)
  if(all(len.G  == 1,
         len.Q  == 1)) {
    times                 <- times[-2]
  }
  class(times)            <- "listof"
  attr(imifa, "Time")     <- if(is.element(method, c("FA", "IFA", "classify"))) times[-length(times)] else times
  attr(imifa, "TuneZeta") <- is.element(method, c("IMFA", "IMIFA")) && (tune.zeta$heat > 0 && tune.zeta$do)
  attr(imifa, "Uni.Meth") <- c(Uni.Prior = uni.prior, Uni.Type = uni.type)
  attr(imifa, "Vars")     <- P
  if(!is.element(method, c("FA", "IFA"))) {
    for(g in seq_along(G.init)) {
      attr(imifa[[g]],
           "Z.init")      <- factor(zi[[g]], levels=seq_len(G.init[g]))
    }
  }
  if(verbose)                print(attr(imifa, "Time"))
  class(imifa)            <- "IMIFA"
    return(imifa)
}
