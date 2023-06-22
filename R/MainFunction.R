#' Adaptive Gibbs Sampler for Nonparametric Model-based Clustering using models from the IMIFA family
#'
#' Carries out Gibbs sampling for all models from the IMIFA family, facilitating model-based clustering with dimensionally reduced factor-analytic covariance structures, with automatic estimation of the number of clusters and cluster-specific factors as appropriate to the method employed. Factor analysis with one group (FA/IFA), finite mixtures (MFA/MIFA), overfitted mixtures (OMFA/OMIFA), infinite factor models which employ the multiplicative gamma process (MGP) shrinkage prior (IFA/MIFA/OMIFA/IMIFA), and infinite mixtures which employ Pitman-Yor / Dirichlet Process Mixture Models (IMFA/IMIFA) are all provided.
#'
#' @param dat A matrix or data frame such that rows correspond to observations (\code{N}) and columns correspond to variables (\code{P}). Non-numeric variables will be discarded if they are explicitly coded as factors or ordinal factors; otherwise they will be treated as though they were continuous. Rows with missing entries will be also be automatically removed.
#' @param method An acronym for the type of model to fit where:
#' \describe{
#'  \item{\code{"FA"}}{Factor Analysis}
#'  \item{\code{"IFA"}}{Infinite Factor Analysis}
#'  \item{\code{"MFA"}}{Mixtures of Factor Analysers}
#'  \item{\code{"MIFA"}}{Mixtures of Infinite Factor Analysers}
#'  \item{\code{"OMFA"}}{Overfitted Mixtures of Factor Analysers}
#'  \item{\code{"OMIFA"}}{Overfitted Mixtures of Infinite Factor Analysers}
#'  \item{\code{"IMFA"}}{Infinite Mixtures of Factor Analysers}
#'  \item{\code{"IMIFA"}}{Infinite Mixtures of Infinite Factor Analysers}
#' }
#' In principle, of course, one could overfit the \code{"MFA"} or \code{"MIFA"} models, but it is recommend to use the corresponding model options which begin with `O' instead. Note that the \code{"classify"} method is not yet implemented.
#' @param range.G Depending on the method employed, either the range of values for the number of clusters, or the conservatively high starting value for the number of clusters. Defaults to (and must be!) \code{1} for the \code{"FA"} and \code{"IFA"} methods. For the \code{"MFA"} and \code{"MIFA"} models this is to be given as a range of candidate models to explore. For the \code{"OMFA"}, \code{"OMIFA"}, \code{"IMFA"}, and \code{"IMIFA"} models, this is the conservatively high number of clusters with which the chain is to be initialised (default = \code{max(25, ceiling(3 * log(N)))} for large N, or \code{min(N-1, ceiling(3 * log(N)))} for small N<=50).
#'
#' For the \code{"OMFA"}, and \code{"OMIFA"} models this upper limit remains fixed for the entire length of the chain; the upper limit for the for the \code{"IMFA"} and \code{"IMIFA"} models can be specified via \code{trunc.G} (see \code{\link{bnpControl}}), which must satisfy \code{range.G <= trunc.G < N}.
#'
#' If \code{length(range.G) * length(range.Q)} is large, consider not storing unnecessary parameters (via \code{\link{storeControl}}), or breaking up the range of models to be explored into chunks and sending each chunk to \code{\link{get_IMIFA_results}} separately.
#' @param range.Q Depending on the method employed, either the range of values for the number of latent factors or, for methods ending in IFA, the conservatively high starting value for the number of cluster-specific factors, in which case the default starting value is \code{round(3 * log(P))}.
#'
#' For methods ending in IFA, different clusters can be modelled using different numbers of latent factors (incl. zero); for methods not ending in IFA it is possible to fit zero-factor models, corresponding to simple diagonal covariance structures. For instance, fitting the \code{"IMFA"} model with \code{range.Q=0} corresponds to a vanilla Pitman-Yor / Dirichlet Process Mixture Model.
#'
#' If \code{length(range.G) * length(range.Q)} is large, consider not storing unnecessary parameters (via \code{\link{storeControl}}), or breaking up the range of models to be explored into chunks and sending each chunk to \code{\link{get_IMIFA_results}}.
#'
#' See \code{\link{Ledermann}} for bounds on \code{range.Q}; this is useful in both the finite factor and infinite factor settings, as one may wish to ensure the fixed number of factors, or upper limits on the number of factors, respectively, respects this bound to yield indentifiable solutions, particularly in low-dimensional settings.
#' @param MGP A list of arguments pertaining to the multiplicative gamma process (MGP) shrinkage prior and adaptive Gibbs sampler (AGS). For use with the infinite factor models \code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"} only. Defaults are set by a call to \code{\link{mgpControl}}, with further checking of validity by \code{\link{MGP_check}} (though arguments can also be supplied here directly).
#' @param BNP A list of arguments pertaining to the Bayesian Nonparametric Pitman-Yor / Dirichlet process priors, for use with the infinite mixture models \code{"IMFA"} and \code{"IMIFA"}, or select arguments related to the Dirichlet concentration parameter for the overfitted mixtures \code{"OMFA"} and \code{"OMIFA"}. Defaults are set by a call to \code{\link{bnpControl}} (though arguments can also be supplied here directly).
#' @param mixFA A list of arguments pertaining to \emph{all other} aspects of model fitting, e.g. MCMC settings, cluster initialisation, and hyperparameters common to every \code{method} in the \code{IMIFA} family. Defaults are set by a call to \code{\link{mixfaControl}} (though arguments can also be supplied here directly).
#' @param alpha Depending on the method employed, either the hyperparameter of the Dirichlet prior for the cluster mixing proportions, or the Pitman-Yor / Dirichlet process concentration parameter. Defaults to \code{1} for the finite mixture models \code{"MFA"} and \code{"MIFA"}, and must be a strictly positive scalar. Not relevant for the \code{"FA"} and \code{"IFA"} methods.
#' \describe{
#' \item{Under the \code{"IMFA"} and \code{"IMIFA"} models:}{\code{alpha} defaults to a simulation from the prior if \code{learn.alpha} is \code{TRUE}, otherwise \code{alpha} \emph{must} be specified. Must be positive, unless non-zero \code{discount} is supplied or \code{learn.d=TRUE} (the default), in which case it must be greater than \code{-discount}. Under certain conditions, \code{alpha} can remain fixed at \code{0} (see \code{\link{bnpControl}}). Additionally, when \code{discount} is negative, \code{alpha} must be a positive integer multiple of \code{abs(discount)} (default=\code{range.G * abs(discount)}).}
#' \item{Under the \code{"OMFA"} and \code{"OMIFA"} models:}{\code{alpha} defaults to a simulation from the prior if \code{learn.alpha} is \code{TRUE}, otherwise \code{alpha} defaults to \code{0.5/range.G}. If supplied, \code{alpha} must be positive, and you are supplying the numerator of \code{alpha/range.G}.
#'
#' If \code{alpha} remains fixed (i.e. \code{learn.alpha=FALSE}), \code{alpha} should be less than half the dimension (per cluster!) of the free parameters of the smallest model considered in order to ensure superfluous clusters are emptied (for \code{"OMFA"}, this corresponds to the smallest \code{range.Q}; for \code{"OMIFA"}, this corresponds to a zero-factor model) [see: \code{\link{PGMM_dfree}} and Rousseau and Mengersen (2011)].}
#' }
#' See \code{\link{bnpControl}} for further details of specifying \code{alpha} or specifying a prior for \code{alpha} under the \code{"IMFA"}, \code{"IMIFA"}, \code{"OMFA"}, or \code{"OMIFA"} methods.
#' @param storage A vector of named logical indicators governing storage of parameters of interest for all models in the IMIFA family. Defaults are set by a call to \code{\link{storeControl}}. It may be useful not to store certain parameters if memory is an issue.
#' @param ... An alternative means of passing control parameters directly via the named arguments of \code{\link{mixfaControl}}, \code{\link{mgpControl}}, \code{\link{bnpControl}}, and \code{\link{storeControl}}. Do not pass the output from calls to those functions here!
#' @param x,object Object of class \code{"IMIFA"}, for the \code{print.IMIFA} and \code{summary.IMIFA} functions, respectively.
#'
#' @details Creates a raw object of class \code{"IMIFA"} from which the optimal/modal model can be extracted by \code{\link{get_IMIFA_results}}. Dedicated \code{print} and \code{summary} functions exist for objects of class \code{"IMIFA"}.
#'
#' @note Further control over the specification of advanced function arguments can be obtained with recourse to the following functions:
#' \itemize{
#' \item{\strong{\code{\link{mgpControl}}} - }{Supply arguments (with defaults) pertaining to the multiplicative gamma process (MGP) shrinkage prior and adaptive Gibbs sampler (AGS). For use with the infinite factor models \code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"} only.}
#' \item{\strong{\code{\link{bnpControl}}} - }{Supply arguments (with defaults) pertaining to the Bayesian Nonparametric Pitman-Yor / Dirichlet process priors, for use with the infinite mixture models \code{"IMFA"} and \code{"IMIFA"}. Certain arguments related to the Dirichlet concentration parameter for the overfitted mixtures \code{"OMFA"} and \code{"OMIFA"} can be supplied in this manner also.}
#' \item{\strong{\code{\link{mixfaControl}}} - }{Supply arguments (with defaults) pertaining to \emph{all other} aspects of model fitting (e.g. MCMC settings, cluster initialisation, and hyperparameters common to every \code{method} in the \code{IMIFA} family.}
#' \item{\strong{\code{\link{storeControl}}} - }{Supply logical indicators governing storage of parameters of interest for all models in the IMIFA family. It may be useful not to store certain parameters if memory is an issue (e.g. for large data sets or for a large number of MCMC iterations after burnin and thinning).}
#' }
#' Note however that the named arguments of these functions can also be supplied directly. Parameter starting values are obtained by simulation from the relevant prior distribution specified in these control functions, though initial means and mixing proportions are computed empirically.
#' @return A list of lists of lists of class \code{"IMIFA"} to be passed to \code{\link{get_IMIFA_results}}. If the returned object is \code{x}, candidate models are accessible via subsetting, where \code{x} is of the following form:
#'
#' \code{x[[1:length(range.G)]][[1:length(range.Q)]]}.
#'
#' However, these objects of class "IMIFA" should rarely if ever be manipulated by hand - use of the \code{\link{get_IMIFA_results}} function is \emph{strongly} advised.
#' @keywords IMIFA main
#' @export
#' @importFrom matrixStats "colMeans2" "colSds" "colSums2" "colVars" "rowLogSumExps" "rowMeans2" "rowSums2"
#' @importFrom Rfast "matrnorm"
#' @importFrom mvnfast "dmvn"
#' @importFrom slam "as.simple_sparse_array" "as.simple_triplet_matrix"
#' @importFrom mclust "emControl" "Mclust" "mclustBIC" "mclustICL" "hc" "hclass" "hcE" "hcEEE" "hcEII" "hcV" "hcVII" "hcVVV"
#'
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link{mixfaControl}}, \code{\link{mgpControl}}, \code{\link{bnpControl}}, \code{\link{storeControl}}, \code{\link{Ledermann}}
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' Bhattacharya, A. and Dunson, D. B. (2011) Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291-306.
#'
#' Kalli, M., Griffin, J. E. and Walker, S. G. (2011) Slice sampling mixture models, \emph{Statistics and Computing}, 21(1): 93-105.
#'
#' Rousseau, J. and Mengersen, K. (2011) Asymptotic Behaviour of the posterior distribution in overfitted mixture models, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 73(5): 689-710.
#'
#' McNicholas, P. D. and Murphy, T. B. (2008) Parsimonious Gaussian mixture models, \emph{Statistics and Computing}, 18(3): 285-296.
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#'
#' @usage
#' mcmc_IMIFA(dat,
#'            method = c("IMIFA", "IMFA",
#'                       "OMIFA", "OMFA",
#'                       "MIFA", "MFA",
#'                       "IFA", "FA",
#'                       "classify"),
#'            range.G = NULL,
#'            range.Q = NULL,
#'            MGP = mgpControl(...),
#'            BNP = bnpControl(...),
#'            mixFA = mixfaControl(...),
#'            alpha = NULL,
#'            storage = storeControl(...),
#'            ...)
#' @examples
#' \donttest{# data(olive)
#' # data(coffee)
#'
#' # Fit an IMIFA model to the olive data. Accept all defaults.
#' # simIMIFA <- mcmc_IMIFA(olive, method="IMIFA")
#' # summary(simIMIFA)
#'
#' # Fit an IMIFA model assuming a Pitman-Yor prior.
#' # Control the balance between the DP and PY priors using the kappa parameter.
#' # simPY    <- mcmc_IMIFA(olive, method="IMIFA", kappa=0.75)
#' # summary(simPY)
#'
#' # Fit a MFA model to the scaled olive data, with isotropic uniquenesses (i.e. MPPCA).
#' # Allow diagonal covariance as a special case where range.Q = 0.
#' # Don't store the scores. Accept all other defaults.
#' # simMFA   <- mcmc_IMIFA(olive, method="MFA", n.iters=10000, range.G=3:6, range.Q=0:3,
#' #                        score.switch=FALSE, centering=FALSE, uni.type="isotropic")
#'
#' # Fit a MIFA model to the centered & scaled coffee data, w/ cluster labels initialised by K-Means.
#' # Note that range.Q doesn't need to be specified. Allow IFA as a special case where range.G=1.
#' # simMIFA  <- mcmc_IMIFA(coffee, method="MIFA", n.iters=10000, range.G=1:3, z.init="kmeans")
#'
#' # Fit an IFA model to the centered and pareto scaled olive data.
#' # Note that range.G doesn't need to be specified. We can optionally supply a range.Q starting value.
#' # Enforce additional shrinkage using alpha.d1, alpha.d2, prop, and eps (via mgpControl()).
#' # simIFA   <- mcmc_IMIFA(olive, method="IFA", n.iters=10000, range.Q=4, scaling="pareto",
#' #                        alpha.d1=2.5, alpha.d2=4, prop=0.6, eps=0.12)
#'
#' # Fit an OMIFA model to the centered & scaled coffee data.
#' # Supply a sufficiently small alpha value. Try varying other hyperparameters.
#' # Accept the default value for the starting number of factors,
#' # but supply a value for the starting number of clusters.
#' # Try constraining uniquenesses to be common across both variables and clusters.
#' # simOMIFA <- mcmc_IMIFA(coffee, method="OMIFA", range.G=10, psi.alpha=3,
#' #                        phi.hyper=c(2, 1), alpha=0.8, uni.type="single")}
mcmc_IMIFA  <- function(dat, method = c("IMIFA", "IMFA", "OMIFA", "OMFA", "MIFA", "MFA", "IFA", "FA", "classify"), range.G = NULL, range.Q = NULL,
                        MGP = mgpControl(...), BNP = bnpControl(...), mixFA = mixfaControl(...), alpha = NULL, storage = storeControl(...), ...) {

  call      <- match.call()
  defopt    <- options()
  options(warn=1)
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
  if(!missing(method) && method == "classification") {
    method  <- "classify"
  }
  if(!all(is.character(method)))    stop("'method' must be a character vector of length 1", call.=FALSE)
  method    <- match.arg(method)
  if(missing(dat))                  stop("Data set must be supplied", call.=FALSE)
  dat.nam   <- gsub("[[:space:]]", "", deparse(substitute(dat)))

# Remove non-numeric columns & apply centering & scaling if necessary
  mixfamiss <- attr(mixFA, "Missing")
  verbose   <- mixFA$verbose
  n.iters   <- mixFA$n.iters
  burnin    <- mixFA$burnin
  thinning  <- mixFA$thinning
  centering <- mixFA$centering
  scaling   <- mixFA$scaling
  iters     <- seq(from=burnin + 2L, to=n.iters + 1L, by=thinning)
  iters     <- iters[iters  > 0]
  if(length(iters)   < 10)          stop("Run a longer chain!", call.=FALSE)
  raw.dat   <- as.data.frame(dat)
  num.check <- vapply(raw.dat, is.numeric, logical(1L))
  if(anyNA(raw.dat)) {
    if(isTRUE(verbose))             message("Rows with missing values removed from data set\n")
    raw.dat <- raw.dat[stats::complete.cases(raw.dat),, drop=FALSE]
  }
  if(sum(num.check) != ncol(raw.dat)) {
    if(isTRUE(verbose))             message("Non-numeric columns removed from data set\n")
    raw.dat <- raw.dat[,num.check,    drop=FALSE]
  }
  raw.dat   <- data.matrix(raw.dat)
  glo.mean  <- colMeans2(raw.dat, refine=FALSE, useNames=FALSE)
  glo.scal  <- colSds(raw.dat,    refine=FALSE, useNames=FALSE, center=glo.mean)
  if(isTRUE(mixFA$drop0sd)) {
    sdx     <- glo.scal
    sd0ind  <- sdx  <= 0
    if(any(sd0ind))  { if(verbose)  message("Columns with standard deviation of zero removed from data set\n")
      raw.dat    <- raw.dat[,!sd0ind, drop=FALSE]
      sdx   <- sdx[!sd0ind]
    }
  }
  if(any(dim(raw.dat) == 0))        stop("Empty data set after removal of ineligble rows/columns", call.=FALSE)
  if(method != "classify")  {
    scal    <- switch(EXPR=scaling, none=FALSE, if(isTRUE(mixFA$drop0sd)) sdx else glo.scal)
    scal    <- switch(EXPR=scaling, pareto=sqrt(scal), scal)
    dat     <- if(centering) .scale2(raw.dat, center=glo.mean, scale=scal)    else .scale2(raw.dat, center=FALSE, scale=scal)
  }
  N         <- as.integer(nrow(dat))
  P         <- as.integer(ncol(dat))
  centered  <- switch(EXPR=method, classify=isTRUE(all.equal(colSums2(dat, useNames=FALSE), numeric(P))), (centering || isTRUE(all.equal(colSums2(dat, useNames=FALSE), numeric(P)))))
  if(!any(centered, centering) && all(verbose,
    scaling != "none"))             message("Are you sure you want to apply scaling without centering?\n")

  uni       <- P == 1
  lnN       <- log(N)
  NlP       <- N <= P
  miss.uni  <- mixfamiss$uni.type
  miss.pri  <- mixfamiss$uni.prior
  uni.type  <- ifelse(miss.uni, ifelse(uni, "isotropic", "unconstrained"), mixFA$uni.type)
  if(uni    && is.element(uni.type, c("unconstrained", "constrained"))) {
    uni.type     <- switch(EXPR=uni.type, unconstrained=, isotropic="isotropic", constrained=, single="single")
    if(isTRUE(verbose))             message(paste0("'uni.type' coerced to ", uni.type, " as the data set is univariate\n"))
  }
  uni.prior <- ifelse(miss.pri, switch(EXPR=uni.type, constrained=, unconstrained="unconstrained", "isotropic"), mixFA$uni.prior)
  if(uni    && uni.prior  == "unconstrained") {
    uni.prior    <- "isotropic"
    if(isTRUE(verbose))             message("'uni.prior' coerced to isotropic as the data set is univariate\n")
  }
  if(all(is.element(uni.type, c("isotropic", "single")),
     uni.prior == "unconstrained")) stop("'uni.prior' can only be 'unconstrained' when 'uni.type' is 'unconstrained' or 'constrained'", call.=FALSE)
  if(all(is.element(uni.type, c("unconstrained", "constrained")), any(miss.uni, miss.pri),
         NlP, verbose))             message(paste0("Consider setting 'uni.type' to 'isotropic' or 'single'", ifelse(miss.pri && uni.prior == "unconstrained", ", or at least 'uni.prior' to 'isotropic', ", " "), "in N <= P cases\n"))

# Manage storage switches & warnings for other function inputs
  if(!storage["u.sw"]  &&
    (!is.element(method, c("FA", "IFA")) ||
     !all(centering, centered)))    stop("'update.mu' can only be FALSE for the \"FA\" or \"IFA\" methods when centering is applied", call.=FALSE)
  if(!storage["u.sw"]  &&
     storage["mu.sw"])  {
    if(verbose)                     message("Forcing 'mu.switch' to FALSE as 'update.mu' is FALSE\n")
    storage["mu.sw"]   <- FALSE
  }
  store.x   <- attr(storage, "Missing")
  if(!store.x["mu.sw"] && all(!storage["mu.sw"], ifelse(method == "classify",
     !centering, !centered)))       warning("Centering hasn't been applied - are you sure you want mu.switch=FALSE?\n", call.=FALSE, immediate.=TRUE)
  if(N < 2)                         stop("Must have more than one observation", call.=FALSE)
  z.init    <- mixFA$z.init
  z.list    <- mixFA$z.list
  zin.miss  <- mixfamiss$z.init
  zli.miss  <- mixfamiss$z.list
  if(all(is.element(method, c("FA", "IFA")), verbose,
         !zin.miss || !zli.miss)) { message(paste0("z does not need to be initialised for the ", method, " method\n"))
  } else if(!zli.miss) {
    if(!inherits(z.list,
                 "list"))    z.list    <- lapply(list(z.list), as.factor)
    if(zin.miss &&
       z.init   != "list") { z.init    <- "list"
       if(isTRUE(verbose))          message("'z.init' set to 'list' as 'z.list' was supplied\n")
    }
  }
  dots      <- switch(EXPR = z.init,
                      mclust = mixFA$dots[c("modelName", "use", "modelNames", "criterion")],
                      hc = mixFA$dots[c("modelName", "use")],
                      kmeans = mixFA$dots[c("iter.max", "nstart")])
  if(is.element(z.init, c("hc", "mclust"))) {
    if(ifelse(z.init      ==    "mclust",
              is.null(dots[names(dots) != "modelNames"]$modelName),
              is.null(dots$modelName))) {
      dots$modelName      <- ifelse(NlP, "EII", "VVV")
    }
    if(is.element(dots$modelName,
       c("E", "V"))       && !uni)  stop("mclust's 'E' and 'V' hierarchical clustering models can only be employed for univariate data", call.=FALSE)
    dots$use              <- ifelse(is.null(dots$use), "VARS", dots$use)
  }
  if(z.init == "kmeans"   && is.null(dots$nstart))    {
    dots$nstart           <- 10L
  }
  dots      <- dots[!vapply(dots, is.null, logical(1L))]

  G.x       <- missing(range.G)
  bnpmiss   <- attr(BNP, "Missing")
  if(!is.element(method, c("MFA", "MIFA")))      {
    if(length(range.G) > 1)         stop(paste0("Only one 'range.G' value can be specified for the ", method, " method"), call.=FALSE)
    if(all(!G.x, verbose, is.element(method, c("FA", "IFA")))  &&
       range.G  > 1)                message(paste0("'range.G' forced to 1 for the ", method, " method\n"))
    if(is.element(method, c("OMIFA", "OMFA", "IMFA", "IMIFA"))) {
      lnN2         <- ceiling(lnN)
      if(G.x)   {
        G.init     <- range.G          <- as.integer(ifelse(N <= 50, min(N - 1L, ceiling(3 * lnN)), max(25L, ceiling(3 * lnN))))
      } else    {
        G.init     <- range.G
      }
      if(isTRUE(verbose))  {
        if(N <= 50 && G.x)          message("Consider initialising closer to the expected truth (~log(N)) when the sample size is small\n")
        if(range.G  < lnN2)         message(paste0("Suggestion:'range.G' should be at least log(N) (=log(", N, "))", " for the ", method, " method\n"))
      }
      if(G.init    >= N)            stop(paste0("'range.G' must less than N (", N, ")"), call.=FALSE)
      if(is.element(method, c("IMFA", "IMIFA")))  {
        if((tmiss  <- bnpmiss$trunc.G)) {
          trunc.G  <- BNP$trunc.G      <- max(min(50L, N - 1L), range.G)
        } else trunc.G    <- BNP$trunc.G
        if(all(verbose, ifelse(N <= 50, trunc.G   < N - 1L,
           trunc.G  < 50), !tmiss)) message(paste0("Consider setting 'trunc.G' to min(N-1=", N - 1L, ", 50) unless practical reasons in heavy computational/memory burden cases prohibit it\n"))
        if(trunc.G  < range.G)      stop(paste0("'trunc.G' must be at least range.G=", range.G), call.=FALSE)
        if(trunc.G >= N)            stop(paste0("'trunc.G' cannot be greater than N-1 (", N - 1L, ")"), call.=FALSE)
        if(trunc.G  > 50)           warning(paste0("'trunc.G' is large: this may lead to memory capacity issues\n"), call.=FALSE, immediate.=TRUE)
      }
    } else if(method      == "classify")   {
      if(!zin.miss &&
         z.init    != "list") {     stop("'z.init' must be set to 'list' for classification", call.=FALSE)
      } else z.init       <- "list"
      if(zli.miss)                  stop("Data labels must be supplied via 'z.list' for classification", call.=FALSE)
      levs         <- nlevels(unlist(z.list))
      if(length(z.list)    > 1)     stop("Only one set of labels can be supplied via 'z.list'", call.=FALSE)
      zlabels      <- unlist(z.list)
      if(length(zlabels)  != N)     stop(paste0("'z.list' must be a factor of length N=",  N), call.=FALSE)
      if(all(verbose, !G.x)    && any(length(range.G > 1),
          range.G  != levs))  {     message("Forced 'range.G' equal to the number of levels in 'zlabels' for the 'classify' method\n")
      }
     G.init <- range.G    <- levs
    } else {
     G.init <- range.G    <- 1L
     storage["pi.sw"]     <- FALSE
    }
    meth    <- method
    equal.pro      <- FALSE
  } else {
    equal.pro      <- mixFA$equal.pro
    if(storage["pi.sw"]   && equal.pro) {
      if(all(!store.x["pi.sw"],
             verbose))              message("Forced non-storage of mixing proportions as 'equal.pro' is TRUE\n")
      storage["pi.sw"]    <- FALSE
    }
    if(G.x)                         stop("'range.G' must be specified",                   call.=FALSE)
    if(any(floor(range.G) != range.G)  ||
       any(range.G  < 1))           stop("'range.G' must be a strictly positive integer", call.=FALSE)
    range.G <- G.init     <- sort(unique(range.G))
    meth    <- rep(method, length(range.G))
  }
  if(any(range.G >= N))             stop(paste0("'range.G' must be less than the number of observations N=", N), call.=FALSE)
  if(G.init[1L]  == 1)     {
    if(is.element(method, c("IMIFA", "IMFA",
       "OMIFA", "OMFA")))  {        stop(paste0("'method' should be ", switch(EXPR=method, IMFA=, OMFA="FA", OMIFA=, IMIFA="IFA"), " for a one group model under the ", method, " method"), call.=FALSE)
    } else {
      meth[1]    <- switch(EXPR=method, MFA=, FA="FA", MIFA=, IFA="IFA")
    }
    if(all(verbose, !is.element(method,
           c("FA", "IFA"))))        message(paste0("Forced use of ", meth[1L], " method where 'range.G' is equal to 1\n"))
  }

  if(!is.element(method, c("IMFA", "IMIFA", "OMFA", "OMIFA"))) {
    if(verbose   &&
      (!missing(BNP)         ||
        any(!unlist(bnpmiss))))      message(paste0("'bnpControl()' parameters not necessary for the ", method, " method\n"))
  } else          {
    learn.a      <- BNP$learn.alpha
    tune.zeta    <- BNP$tune.zeta
    BNP$zeta     <- ifelse(bnpmiss$zeta, switch(EXPR=method, IMFA=, IMIFA=2L, 0.75), BNP$zeta)
    tune.zeta$do <- switch(EXPR=method, OMFA=, OMIFA=tune.zeta$do, tune.zeta$do && attr(tune.zeta, "IM.Need"))
    if(tune.zeta$do)   {
      if(tune.zeta$stop.zeta <=
         tune.zeta$start.zeta)      stop(paste0("'stop.zeta' must be greater than 'start.zeta' (=", tune.zeta$start.zeta, ")"), call.=FALSE)
      if(!(tune.zeta$do      <-
         tune.zeta$start.zeta <
         n.iters)     &&
         isTRUE(verbose))     {     message("Diminishing adaptation to tune zeta not invoked as 'start.zeta' is not less than 'n.iters'\n")
      } else if(!attr(tune.zeta, "stopx") && isTRUE(verbose) &&
                tune.zeta$stop.zeta       >=
                n.iters)            message("'stop.zeta' not invoked as it is not less than 'n.iters'\n")
    }
    if(is.element(method, c("IMFA", "IMIFA"))) {
      discount   <- BNP$discount
      learn.d    <- BNP$learn.d
      IM.lab.sw  <- BNP$IM.lab.sw
      if(diff(BNP$a.hyper)    < 0)  warning("The rate hyperparameter for the prior on alpha should be >= to the shape hyperparameter, in order to encourage clustering\n", call.=FALSE, immediate.=TRUE)
      if(all(discount         > 0,
             !learn.d, learn.a)) {
        BNP$a.hyper          <- unname(shift_GA(shape=BNP$a.hyper[1L], rate=BNP$a.hyper[2L], shift=-discount))
      }
      if(all(!bnpmiss$IM.lab.sw, IM.lab.sw, !learn.a,
             verbose, !learn.d))    message("May not be necessary to set 'IM.lab.sw' to TRUE when neither 'alpha' nor 'discount' are being learned\n")
      if(all(BNP$exchange,
             IM.lab.sw))         {
        if(verbose)                 message("'IM.lab.sw' forced to FALSE as 'exchange' is TRUE\n")
        IM.lab.sw  <-
        BNP$IM.lab.sw        <- FALSE
      }
      if(all(any(BNP$exchange, BNP$thresh),
             BNP$ind.slice))     {
        if(verbose)                 message("'ind.slice' forced to FALSE as 'exchange' &/or 'thresh' is TRUE\n")
        BNP$ind.slice        <- FALSE
      }
    } else BNP$a.hyper[2L]   <- BNP$a.hyper[2L] * G.init
  }

# Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  datname          <- rownames(dat)
  if(any(length(unique(datname)) != N,
     is.null(datname)))  rownames(dat) <- seq_len(N)
  cmeans           <- colMeans2(dat, refine=FALSE, useNames=FALSE)
  mu               <- list(as.matrix(cmeans))
  mu0.x            <- mixfamiss$mu.zero
  mu0g             <- mixFA$mu0g
  mu.zero          <- if(mu0.x) mu  else  .len_check(mixFA$mu.zero, mu0g, method, P, G.init)
  sigma.mu         <- mixFA$sigma.mu
  sigmu.miss       <- mixfamiss$sigma.mu
  psi.alpha        <- mixFA$psi.alpha
  psi0g            <- mixFA$psi0g
  beta.x           <- mixfamiss$psi.beta
  if(beta.x && psi.alpha <= 1)      stop("'psi.alpha' must be strictly greater than 1 when invoking the default for 'psi.beta' in order to bound uniquenesses away from zero",    call.=FALSE)
  if(psi.alpha     <= 1)            warning("'psi.alpha' is not strictly greater than 1; uniquenesses may not be sufficiently bounded away from zero, algorithm may terminate\n", call.=FALSE)
  obsnames         <- rownames(dat)
  varnames         <- colnames(dat)
  covmat           <- switch(EXPR=scaling, unit=stats::cor(dat), stats::cov(dat))
  dimnames(covmat) <- list(varnames, varnames)
  if(anyNA(covmat))                 warning(paste0("Covariance matrix cannot be estimated: ", ifelse(beta.x || (sigmu.miss && scaling != "unit"), "deriving mean/uniqueness hyperparameters may not be possible, neither will certain posterior predictive checks\n", "certain posterior predictive checks will not be possible\n")), call.=FALSE)
  sigma.mu         <- if(sigmu.miss && scaling == "unit") 1L else if(sigmu.miss) diag(covmat) else if(is.matrix(sigma.mu)) diag(sigma.mu) else sigma.mu
  sigma.mu         <- sigma.mu/mixFA$prec.mu
  if(anyNA(sigma.mu))               stop(ifelse(sigmu.miss, "Not possible to derive default 'sigma.mu': this argument now must be supplied", "NA in 'sigma.mu'"), call.=FALSE)
  sigmu.len        <- length(unique(sigma.mu))
  if(sigmu.len     == 1) sigma.mu      <- sigma.mu[1L]
  if(any(sigma.mu  <= 0, !is.numeric(sigma.mu),
     !is.element(length(sigma.mu),
     c(1, P))))                     stop(paste0("'sigma.mu' must be strictly positive, and of length 1 or P=", P))
  if(beta.x) {
    psi.beta       <- temp.psi   <- tryCatch(list(psi_hyper(shape=psi.alpha, dat=dat, type=uni.prior, covar=covmat, ...)), error=function(e) message(paste0(e)))
    if(is.null(psi.beta))           stop("Not possible to derive default 'psi.beta': try supplying this argument directly", call.=FALSE)
  } else {
    psi.beta       <- mixFA$psi.beta
    psi.beta       <- lapply(.len_check(psi.beta, psi0g, method, P, G.init), matrix, nrow=P, ncol=G.init, byrow=length(psi.beta) == G.init)
  }
  if(any(unlist(psi.beta) < 1E-03,
         psi.alpha  < 1E-03))       warning("Excessively small values for the uniquenesses hyperparameters may lead to critical numerical issues & should thus be avoided\n", call.=FALSE, immediate.=TRUE)
  Q.miss    <- missing(range.Q)
  Q.min     <- min(ceiling(log(P)), ceiling(log(N)))

  mgpmiss   <- attr(MGP, "Missing")
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
   if(verbose      &&
      !missing(MGP)        ||
      any(!unlist(mgpmiss)))        message(paste0("'mgpControl()' parameters not necessary for the ", method, " method\n"))
   if(Q.miss)                       stop("'range.Q' must be specified", call.=FALSE)
   if(any(floor(range.Q)   != range.Q) ||
      any(range.Q   < 0))           stop(paste0("'range.Q' must be a non-negative integer for the ", method, " method"), call.=FALSE)
   range.Q  <- sort(unique(range.Q))
   delta0g  <- FALSE
  } else {
   fQ0      <- uni || P    == 2
   MGP$prop <- ifelse(fQ0  && mgpmiss$propx, 0.5, floor(MGP$prop * P)/P)
   adapt    <- MGP$adapt
   truncate <- MGP$truncated
   delta0g  <- MGP$delta0g && method   == "MIFA"
   MGP$nu1  <- nu1         <- MGP$phi.hyper[1L]
   MGP$nu2  <- nu2         <- MGP$phi.hyper[2L]
   MGP$rho1 <- rho1        <- if(MGP$cluster.shrink && method != "IFA") MGP$sigma.hyper[1L]
   MGP$rho2 <- rho2        <- if(MGP$cluster.shrink && method != "IFA") MGP$sigma.hyper[2L]
   alpha.d1 <- .len_check(MGP$alpha.d1, delta0g, method, P, G.init, P.dim=FALSE)
   alpha.d2 <- .len_check(MGP$alpha.d2, delta0g, method, P, G.init, P.dim=FALSE)
   MGP      <- MGP[-seq_len(5L)]
   start.AGS       <-  MGP$start.AGS   <- ifelse(mgpmiss$startAGSx, pmin(burnin, ifelse(fQ0, 2L, switch(EXPR=method, IFA=, MIFA=burnin, 2L))), MGP$start.AGS)
   if(Q.miss)                range.Q   <- as.integer(ifelse(fQ0, 1L, min(ifelse(P > 500, 12L + round(log(P)), round(3 * log(P))), N - 1L, P - 1L)))
   if(length(range.Q)       > 1)    stop(paste0("Only one starting value for 'range.Q' can be supplied for the ", method, " method"), call.=FALSE)
   if(range.Q      <= 0)            stop(paste0("'range.Q' must be strictly positive for the ", method, " method"), call.=FALSE)
   if(isTRUE(adapt))  {
     if(start.AGS     > burnin)     stop("'start.AGS' must be <= 'burnin' if 'adapt' is TRUE", call.=FALSE)
     if(MGP$stop.AGS <= start.AGS)  stop(paste0("'stop.AGS' must be greater than 'start.AGS' (=", start.AGS, ")"),  call.=FALSE)
     if(!mgpmiss$stopAGSx  && verbose  &&
        MGP$stop.AGS >= n.iters)    message("'stop.AGS' not invoked as it is not less than 'n.iters'\n")
     if(Q.min   > range.Q)          stop(paste0("'range.Q' must be at least min(log(P), log(N)) for the ", method, " method when 'adapt' is TRUE"), call.=FALSE)
    } else if(!fQ0 && is.element(method,
              c("OMIFA", "IMIFA"))) warning("'adapt=FALSE' is NOT recommended for the 'OMIFA' or 'IMIFA' methods\n", call.=FALSE, immediate.=TRUE)
  }

  len.G     <- switch(EXPR=method, classify=range.G, length(range.G))
  len.Q     <- length(range.Q)
  len.X     <- len.G * len.Q
  if(all(len.X > 10,
         suppressWarnings(utils::memory.limit()) <= 16256,
         storage["s.sw"]))  {
    if(!store.x["s.sw"])    {       warning(paste0("The large number of candidate models being explored (", len.X, ") could lead to memory issues\nConsider setting 'score.switch' to FALSE or breaking up the task into chunks and calling get_IMIFA_results() on each chunk\n"), call.=FALSE, immediate.=TRUE)
    } else                  {       warning(paste0("'score.switch' set to FALSE as too many candidate models are being explored (", len.X, ")\nPosterior inference on the scores will not be possible, though you can risk forcing storage by supplying score.switch=TRUE\nConsider breaking up the task into chunks and calling get_IMIFA_results() on each chunk\n"), call.=FALSE, immediate.=TRUE)
      storage["s.sw"]    <- FALSE
    }
  }
  Q.warn       <- min(N - 1L, Ledermann(P, isotropic=is.element(uni.type, c("isotropic", "single"))))
  if(any(range.Q > Q.warn))   {
    if(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")) &&
       isTRUE(adapt))         {     warning(paste0("Starting value for number of factors is greater than ", ifelse(any(range.Q > P), paste0("the number of variables (", P, ")"), paste0("the suggested Ledermann upper bound (", Q.warn, ")\n"))), call.=FALSE, immediate.=TRUE)
    } else if(any(is.element(method, c("FA", "MFA", "OMFA", "IMFA")),
                  is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")) &&
                  isFALSE(adapt)))  warning(paste0("Number of factors is greater than ", ifelse(any(range.Q > P), paste0("the number of variables (", P, ")"), paste0("the suggested Ledermann upper bound (", Q.warn, ")\n"))), call.=FALSE, immediate.=TRUE)
  }
  if(verbose   && !all(storage))  {
    if(any(!storage))             {
      if(all(storage["s.sw"], !storage["l.sw"],
         any(range.Q   != 0)))      message("Loadings not stored but scores are: Procrustes rotation of scores will not occur when passing results to get_IMIFA_results()\n")
      sX       <- c("mu.sw", "l.sw", "psi.sw")
      if(any(!storage[sX]))       {
        if(!is.element(method,
           c("FA", "IFA"))      ||
          (storage[6L] &&
           any(!storage[sX[-1L]]))) message(paste0("Non-storage of parameters means posterior predictive checking error metrics almost surely will not be available", ifelse(is.element(max(which(!storage[sX])), c(1L, 4L)), paste0(":\nposterior mean parameter estimates of the ", ifelse(!storage["mu.sw"], ifelse(storage["psi.sw"], "means", "means and uniquenesses"), "uniquenesses"), " will still be available\n"), "\n")))
      } else if(any(all(method  == "MFA",  any(range.G > 1)) && any(range.Q > 0),
                    all(method  == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))) &&
          (!equal.pro  &&
           !storage["pi.sw"]))      message("Non-storage of mixing proportions parameters means posterior predictive checking error metrics will almost surely not all be available\n")
      if(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")) && isFALSE(adapt) &&
         !storage["l.sw"])          message("Non-storage of loadings means post-hoc adaptation of the number of active factors will not be possible\n")
    }
  } else if(verbose && is.element(method, c("FA", "MFA", "OMFA", "IMFA")) && any(range.Q == 0)) {
    if(all(storage[c("s.sw", "l.sw")]))   {
                                    message("Scores & Loadings not stored where 'range.Q=0' as model has zero factors\n")
    } else if(storage["s.sw"])    { message("Scores not stored where 'range.Q==0' as model has zero factor\n")
    } else if(storage["l.sw"])    { message("Loadings not stored where 'range.Q==0' as model has zero factors\n")
    }
    if(all(range.Q  == 0))   storage[c("s.sw", "l.sw")] <- FALSE
  }

  sw0gs     <- c(mu0g = mu0g, psi0g = psi0g, delta0g = delta0g)
  if(is.element(method, c("MFA", "MIFA", "classify"))) {
    if(all(z.init != "list", any(sw0gs))) {
      if(isTRUE(delta0g))           stop("'delta0g' can only be TRUE if z.init='list' for the 'MIFA' method\n", call.=FALSE)
      if(all(!mu0.x,  mu0g))        stop("'mu.zero' can only be supplied for each cluster if z.init='list' for the 'MFA' & 'MIFA' methods", call.=FALSE)
      if(all(!beta.x, psi0g))       stop("'psi.beta' can only be supplied for each cluster if z.init='list' for the 'MFA' & 'MIFA' methods", call.=FALSE)
    }
    if(all(is.element(uni.type, c("constrained", "single")),
           isTRUE(psi0g)))  {       warning(paste0("'psi0g' forced to FALSE as uniquenesses are constrained across clusters (i.e. 'uni.type' = ", uni.type, ")\n"), call.=FALSE)
      psi0g <- FALSE
    }
    if(method == "classify") mu0g      <- TRUE
  } else if(any(sw0gs))             stop(paste0("'", names(which(sw0gs)), "' should be FALSE for the ", method, " method\n"), call.=FALSE)

  if(!is.element(method, c("FA", "IFA", "classify"))) {
    if(missing(alpha))     { alpha     <- switch(EXPR=method, MFA=, MIFA=1L, OMFA=, OMIFA=ifelse(learn.a, min(1/G.init, stats::rgamma(1, BNP$a.hyper[1L], BNP$a.hyper[2L])), 0.5/G.init), ifelse(learn.a && discount >= 0, max(1L, stats::rgamma(1, BNP$a.hyper[1L], BNP$a.hyper[2L])), ifelse(learn.a, range.G * abs(discount), 1L)))
      if(is.element(method, c("IMIFA", "IMFA"))      &&
         !learn.a)                  stop("'alpha' must be specified if it is to remain fixed when 'learn.alpha' is FALSE, as it's not being learned via Gibbs/Metropolis-Hastings updates", call.=FALSE)
    } else if(is.element(method,
      c("IMFA", "IMIFA"))) {
      if(discount < 0     &&
        (alpha   <= 0     ||
       !.IntMult(alpha, discount))) stop("'alpha' must be a positive integer multiple of 'abs(discount)' when 'discount' is negative", call.=FALSE)
    } else if(is.element(method,
      c("OMFA", "OMIFA")))   alpha     <- alpha/G.init
    if(length(alpha) != 1)          stop("'alpha' must be specified as a scalar to ensure an exchangeable prior", call.=FALSE)
    if(is.element(method,   c("IMIFA", "IMFA")))      {
      if(kappa0      <- alpha    <= 0)  {
        discount     <- BNP$discount   <- ifelse(ifelse(learn.d, discount == 0, bnpmiss$discount), pmin(pmax(stats::rbeta(1, BNP$d.hyper[1L], BNP$d.hyper[2L]), .Machine$double.eps - alpha), 1 - .Machine$double.eps), discount)
        if(kappa0    <- kappa0   && !learn.a)         {
         if(!learn.d &&
            discount == 0)          stop("Set 'learn.d'=TRUE or fix a non-zero 'discount' value if fixing 'alpha' at <= 0", call.=FALSE)
         if(learn.d  && bnpmiss$kappa)  {
           BNP$kappa <- 0L
         } else if(BNP$kappa     != 0  &&
                   learn.d)          stop("Set 'kappa'=0 if fixing 'alpha' at <= 0 and 'learn.d'=TRUE", call.=FALSE)
        }
      }
      if(alpha       <= -discount)  stop(paste0("'alpha' must be ",     ifelse(discount != 0, paste0("strictly greater than -discount (i.e. > ", - discount, ")"), "strictly positive")), call.=FALSE)
    }
    if(is.element(method, c("OMIFA",   "OMFA")) && !learn.a) {
      min.d2         <- 0.5 * PGMM_dfree(P=P, Q=switch(EXPR=method, OMFA=min(range.Q), OMIFA=0L), equal.pro=equal.pro,
                                         method=switch(EXPR=uni.type, unconstrained="UUU", isotropic="UUC", constrained="UCU", single="UCC"))
      if(alpha       >= min.d2)     warning(paste0("'alpha' over 'range.G' for the OMFA & OMIFA methods when 'learn.alpha=FALSE', should be less than half the dimension (per cluster!)\nof the free parameters of the smallest model considered (= ", min.d2, "): consider suppling 'alpha' < ", min.d2 * G.init, "\n"), call.=FALSE, immediate.=TRUE)
    }
    if(any(all(is.element(method, c("MFA",  "MIFA")),  alpha > 1),
           all(is.element(method, c("OMFA", "OMIFA")) && !learn.a,
           alpha > 1/G.init)))      warning("Are you sure alpha should be greater than 1?\n", call.=FALSE, immediate.=TRUE)
    if(all(is.element(method,     c("MFA", "MIFA")),
           range.G    > lnN,
           alpha      < 1))         warning("Empty clusters are likely with small 'alpha' and large 'range.G', consider running an overfitted or infinite mixture\n", call.=FALSE, immediate.=TRUE)
    if(!zli.miss)     {
      if(length(z.list)   != len.G)  {
                                    stop(paste0("'z.list' must be a list of length ", len.G), call.=FALSE)  }
      list.levels    <- lapply(z.list, nlevels)
      if(!all(list.levels == G.init))             {
        if(!is.element(method, c("IMIFA",
                       "IMFA", "OMIFA", "OMFA"))) {
                                    stop(paste0("Each element of 'z.list' must have the same number of levels as 'range.G'"), call.=FALSE)
        } else                      stop(paste0("Only ", list.levels, " clusters are populated according to z.list, but 'range.G' has been set to ", G.init, ":\nReset 'range.G' to this value to avoid redunandtly carrying around empty clusters or supply a list with ", G.init, " levels"), call.=FALSE)
      }
      if(!all(lengths(z.list) == N)) {
                                    stop(paste0("Each element of 'z.list' must be a vector of length N=", N), call.=FALSE) }
    }
    if(all(zli.miss, z.init   == "list"))         {
                                    stop(paste0("'z.list' must be supplied if 'z.init' is set to 'list'"),    call.=FALSE) }
  }

  imifa     <- list(list())
  Gi        <- Qi  <- 1L
  gibbs.arg <- list(P = P, sigma.mu = sigma.mu, psi.alpha = psi.alpha, burnin = burnin, sw = storage, col.mean = cmeans,
                    thinning = thinning, iters = iters, verbose = verbose, uni.type = uni.type, uni.prior = uni.prior)
  if(is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))) {
    gibbs.arg      <- append(gibbs.arg, BNP)
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
    gibbs.arg      <- append(gibbs.arg, list(sigma.l   = mixFA$sigma.l))
  }

  init.start       <- proc.time()
  if(!is.element(method,  c("FA",    "IFA"))   &&
     !all(G.init == 1))         {
    if(isTRUE(verbose))             message(paste0("Initialising...\n"))
    clust          <- list()
    pi.alpha       <- list()
    pi.prop        <- list()
    zi             <- list()
    if(is.element(z.init, c("hc", "mclust")))   {
      hcG          <- G.init[G.init > 1]
      Zhc          <- tryCatch(hc(data=dat, modelName=dots$modelName, minclus=min(hcG), use=dots$use), error=function(e) stop(paste0("Hierarchical clustering initialisation failed", ifelse(z.init == "hc", ". ", " (when initialising using 'mclust'). "), "Try another z.init method", ifelse(length(dots) > 1, " or supply different 'hc' arguments via '...' to mixfaControl", "")), call.=FALSE))
      if(z.init    == "mclust") {
        mcarg      <- list(data=dat, G=hcG, verbose=FALSE, control=emControl(equalPro=equal.pro), initialization=list(hcPairs=Zhc), unlist(mixFA$dots["modelNames"]))
        mcl        <- suppressWarnings(if(identical(dots$criterion, "icl")) do.call(mclustICL, mcarg) else do.call(mclustBIC, mcarg))
        class(mcl) <- "mclustBIC"
        if(!is.null(dots$criterion) && verbose &&
           !is.element(dots$criterion,
           c("bic", "icl")))        message("Using 'bic' to determine optimal mclust model to initialise with\n")
      } else        {
        hc1        <- any(G.init   == 1)
        hcZ        <- hclass(hcPairs=Zhc, G=hcG)
      }
    }
    for(g in seq_along(G.init)) {
      G            <- G.init[g]
      Gseq         <- seq_len(G)
      if(G == 1)    {
        zi[[g]]    <- rep(1L, N)
      } else if(z.init   == "mclust") {
        m.res      <- try(suppressWarnings(Mclust(data=dat, G=G, x=mcl)), silent=TRUE)
        if(!inherits(m.res, "try-error")  &&
           !is.null(m.res))     {
          zi[[g]]  <- as.integer(factor(m.res$classification, levels=Gseq))
        } else                      stop(paste0("Cannot initialise cluster labels using mclust. Try another z.init method", ifelse(length(dots) > 1, " or supply different 'mclust' arguments via '...' to mixfaControl", "")), call.=FALSE)
      } else if(z.init   == "hc")     {
        zi[[g]]    <- as.integer(factor(hcZ[,g - hc1],  levels=Gseq))
      } else if(z.init   == "kmeans") {
        k.res      <- try(suppressWarnings(stats::kmeans(dat, centers=G, dots)), silent=TRUE)
        if(!inherits(k.res, "try-error"))  {
          zi[[g]]  <- as.integer(factor(k.res$cluster, levels=Gseq))
        } else                      stop(paste0("Cannot initialise cluster labels using kmeans. Try another z.init method", ifelse(length(dots) > 1, " or supply different 'kmeans' arguments via '...' to mixfaControl", "")), call.=FALSE)
      } else if(z.init   == "list")   {
        zi[[g]]    <- as.integer(z.list[[g]])
      } else {
        zips       <- rep(1L, N)
        iter       <- 0L
        n0         <- TRUE
        switch(EXPR=method,
               IMFA=, IMIFA=          {
          while(n0 && iter < 100)     {
            G.sim  <- G.init - 1L
            vies   <- .sim_vs_inf(alpha=alpha, discount=discount, len=G.sim, lseq=seq_len(G.sim), N=N, nn=rep(N/G.init, G.sim))
            pies   <- .sim_pi_inf(vs=c(vies, 1L), len=G.init)
            zips   <- .sim_z_p(N=N, prob.z=pies)
            iter   <- iter + 1L
            n0     <- any(tabulate(zips, nbins=G.init) <= 0)
          }
        },          {
          a.tmp    <- switch(EXPR=method, OMFA=, OMIFA=alpha * G.init, alpha)
          while(n0 && iter < 100)     {
            pies   <- rDirichlet(alpha=a.tmp, G.init)
            if(any(is.nan(pies)))   stop("'alpha' is too small to initialise cluster labels from the prior", call.=FALSE)
            zips   <- .sim_z_p(N=N, prob.z=pies)
            iter   <- iter + 1L
            n0     <- any(tabulate(zips, nbins=G.init) <= 0)
          }
        })
        if(isTRUE(n0))              stop("Empty clusters after simulating labels from the priors; try another 'z.init' method", call.=FALSE)
        zi[[g]]    <- as.integer(zips)
        rm(zips)
      }
      nngs         <- tabulate(zi[[g]], nbins=switch(EXPR=method, IMFA=, IMIFA=trunc.G, G))
      if(verbose   && zli.miss)     message(paste0("G=", G, " - initial cluster sizes: ", paste(nngs[Gseq], collapse=", "), "\n"))
      pi.prop[[g]] <- if(equal.pro) rep(1/G, G) else prop.table(nngs)
      mu[[g]]      <- vapply(Gseq, function(gg) if(nngs[gg] > 0) colMeans2(dat[zi[[g]] == gg,, drop=FALSE], refine=FALSE, useNames=FALSE) else integer(P), numeric(P))
      mu[[g]]      <- if(uni)       t(mu[[g]])  else mu[[g]]
      if(mu0.x)   {
        mu.zero[[g]]    <- if(mu0g) mu[[g]]     else replicate(G, cmeans, simplify="array")
      }
      mu.zero[[g]] <- if(uni && !mu0g) t(mu.zero[[g]]) else mu.zero[[g]]
      if(beta.x)  {
        if(psi0g) {
          dat.gg   <- lapply(Gseq, function(gg) dat[zi[[g]] == gg,, drop=FALSE])
          cov.gg   <- lapply(Gseq, function(gg) if(all(nngs[gg] > 1, P <= nngs[g])) stats::cov(dat.gg[[gg]]) else covmat)
          psi.beta[[g]] <- try(vapply(Gseq, function(gg) psi_hyper(shape=psi.alpha, dat=dat.gg[[gg]], type=uni.prior, covar=cov.gg[[gg]], ...), numeric(P)), silent=TRUE)
        }
        if(!psi0g  || inherits(psi.beta[[g]], "try-error")) {
          psi.beta[[g]] <- replicate(G, temp.psi[[1L]])
        }
        if(any(psi.beta[[g]]  <
               1E-03))              stop("Excessively small values for the uniquenesses hyperparameters will lead to critical numerical issues & should thus be avoided", call.=FALSE)
        psi.beta[[g]]   <- if(uni)  t(psi.beta[[g]])   else psi.beta[[g]]
      }
      clust[[g]]   <- list(z = zi[[g]], pi.alpha = alpha, pi.prop = pi.prop[[g]])
      if(is.element(method, c("MFA", "MIFA"))) {
        clust[[g]] <- append(clust[[g]], list(l.switch = sw0gs))
      }
      if(is.element(method, c("classify", "MIFA"))) {
        clust[[g]] <- append(clust[[g]], list(alpha.d1 = alpha.d1[[g]], alpha.d2 = alpha.d2[[g]]))
      }
    }
  }

  if(is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))) {
    mu.zero        <- if(all(lengths(mu.zero)  == 1L)) list(mu.zero[[1L]])  else list(mu.zero[[1L]][,1L,  drop=FALSE])
    psi.beta       <- if(all(lengths(psi.beta) == 1L)) list(psi.beta[[1L]]) else list(psi.beta[[1L]][,1L, drop=FALSE])
    if(!is.element(method, c("OMFA", "IMFA")))  {
      alpha.d1     <- list(alpha.d1[[1L]][1L])
      alpha.d2     <- list(alpha.d2[[1L]][1L])
    }
  }
  if(isTRUE(all.equal(vapply(mu.zero, sum, numeric(1L)),
                      numeric(length(G.init)))))   {
    mu.zero        <- switch(EXPR=method, classify=base::matrix(0L, nrow=1L, ncol=range.G), lapply(mu.zero, function(x) 0L))
  }
  if(mu0g && unlist(lapply(mu.zero,   function(x)  {
    if(is.matrix(x)) any(apply(x, 1L, function(y)  {
     length(unique(round(y, min(.ndeci(y))))) })  == 1) else  {
     length(unique(round(x,
     min(.ndeci(x))))) == 1 }})))   stop("'mu0g' must be FALSE if 'mu.zero' is not group-specific", call.=FALSE)
  if(anyNA(unlist(psi.beta))) {
    psi.beta       <- lapply(psi.beta, function(x) replace(x, is.na(x), 0L))
  }
  if(all(is.element(uni.type, c("isotropic",   "single")),
         unlist(lapply(psi.beta,      function(x)  {
    if(is.matrix(x)) any(apply(x, 2L, function(y)  {
     length(unique(round(y, min(.ndeci(y))))) })  != 1) else  {
     length(unique(round(x,
     min(.ndeci(x))))) != 1 }}))))  stop("'psi.beta' cannot be variable-specific if 'uni.type' is 'isotropic' or 'single'", call.=FALSE)
  if(is.element(uni.type, c("constrained", "single")))        {
      if(unlist(lapply(psi.beta,      function(x)  {
    if(is.matrix(x)) any(apply(x, 1L, function(y)  {
     length(unique(round(y, min(.ndeci(y))))) })  != 1) else  {
     !is.element(method, c("FA", "IFA"))          &&
     length(unique(round(x,
     min(.ndeci(x))))) != 1 }}))) { stop("'psi.beta' cannot be group-specific if 'uni.type' is 'constrained' or 'single'",  call.=FALSE)
    } else if(isTRUE(psi0g) && is.element(method,
              c("MFA", "MIFA")))    stop("'psi0g' must be FALSE if 'psi.beta' is not group-specific", call.=FALSE)
  }
  if(any(unlist(psi.beta)   <= 0))  stop("'psi.beta' must be strictly positive", call.=FALSE)
  if(is.element(method, c("classify", "IFA", "MIFA", "IMIFA", "OMIFA"))) {
    ad1uu          <- unique(unlist(alpha.d1))
    ad2uu          <- unique(unlist(alpha.d2))
    check.mgp      <- suppressWarnings(MGP_check(ad1=ad1uu, ad2=ad2uu, Q=unique(range.Q), phi.shape=nu1, phi.rate=nu2, sigma.shape=rho1, sigma.rate=rho2, truncated=truncate, bd1=MGP$beta.d1, bd2=MGP$beta.d2))
    if(!all(check.mgp$valid))       stop("Invalid shrinkage hyperparameter values WILL NOT encourage loadings column removal.\nTry using the MGP_check() function in advance to ensure the cumulative shrinkage property holds.", call.=FALSE)
    if(any(attr(check.mgp, "Warning"))) {
      if(any(ad2uu <= ad1uu))       warning("Column shrinkage hyperparameter values MAY NOT encourage loadings column removal.\n'alpha.d2' should be moderately large relative to 'alpha.d1'\n", call.=FALSE, immediate.=TRUE)
      if(any(nu2    > nu1 - 1))     warning("Expectation of local shrinkage hyperprior is not less than 1\n", call.=FALSE, immediate.=TRUE)
    }
    deltas         <- lapply(seq_along(G.init), function(g) list(alpha.d1 = alpha.d1[[g]], alpha.d2 = alpha.d2[[g]]))
  }
  init.time        <- proc.time() - init.start
  fac.time         <- 0
  G.range          <- switch(EXPR=method, IMIFA=, IMFA=G.init, range.G)

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
        if(verbose && Gi  != len.G) message(paste0("Model ", Gi, " of ", len.G, " complete"), "\nInitialising...", sep="\n")
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
        if(verbose && Qi  != len.Q) message(paste0("Model ", Qi, " of ", len.Q, " complete"), "\nInitialising...", sep="\n")
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
        if(verbose && Gi  != len.G) message(paste0("Model ", Gi, " of ", len.G, " complete"), "\nInitialising...", sep="\n")
      }
    } else {
      mi           <- 0L
      start.time   <- proc.time()
      for(g in G.range)   {
        Gi         <- which(G.range == g)
        imifa[[Gi]]       <- list()
        for(q in range.Q) {
          Qi       <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0(".gibbs_", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        mi         <- mi + 1L
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && mi  != len.X) message(paste0("Model ", mi, " of ", len.X, " complete"), "\nInitialising...", sep="\n")
        }
      }
    }
  } else if(method == "classify") { stop("'classify' method not yet implemented", call.=FALSE)
    start.time     <- proc.time()
    if(centered)                    warning("Data supplied is globally centered, are you sure?\n", call.=FALSE)
    for(g in seq_len(range.G))  {
      tmp.dat      <- raw.dat[zlabels == levels(zlabels)[g],]
      scal         <- switch(EXPR=scaling, none=FALSE, colSds(tmp.dat, refine=FALSE, useNames=FALSE))
      scal         <- switch(EXPR=scaling, pareto=sqrt(scal), scal)
      tmp.dat      <- .scale2(as.matrix(tmp.dat), center=centering, scale=scal)
      if(sigmu.miss) {
       gibbs.arg$sigma.mu <- diag(if(nrow(tmp.dat) > 1) stats::cov(tmp.dat) else covmat)
      }
      imifa[[g]]          <- list()
      gibbs.arg    <- append(temp.args, lapply(deltas[[Gi]], "[[", g))
      imifa[[g]][[Qi]]    <- do.call(paste0(".gibbs_", "IFA"),
                                     args=append(list(data = tmp.dat, N = nrow(tmp.dat), mu = mu[[Gi]][,g], mu.zero = mu.zero[[Gi]][,g],
                                                      Q = range.Q, psi.beta = psi.beta[[Gi]][,g]), gibbs.arg))
      fac.time     <- fac.time + imifa[[g]][[Qi]]$time
      if(verbose   && g   != len.G) message(paste0("Model ", g, " of ", len.G, " complete"), "\nInitialising...", sep="\n")
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/switch(EXPR=method, MFA=len.X, MIFA=len.G, classify=range.G, len.Q)
  tot.time  <- tot.time    + init.time
  init.time <- init.time   + fac.time
  for(g  in length(imifa)) {
   for(q in length(imifa[[g]])) {
     imifa[[g]][[q]]$time <- NULL
   }
  }

  imifa     <- switch(EXPR=method, FA=, MFA=, OMFA=, IMFA={
     lapply(seq_along(imifa), function(x) stats::setNames(imifa[[x]], paste0(range.Q, ifelse(range.Q == 1, "Factor", "Factors"))))
  }, lapply(seq_along(imifa), function(x) stats::setNames(imifa[[x]], "IFA")))
  gnames    <- switch(EXPR=method, classify=paste0("Cluster ", seq_len(range.G)), paste0(G.init, ifelse(G.init == 1, "Cluster", "Clusters")))
  names(imifa)            <- gnames
  attr(imifa, "Adapt")    <- is.element(method, c("classify", "IFA", "MIFA", "OMIFA", "IMIFA")) && isTRUE(adapt)
  attr(imifa,
       "Alph.step")       <- is.element(method, c("IMFA", "IMIFA", "OMFA", "OMIFA")) && learn.a
  attr(imifa, "Alpha")    <- if(!attr(imifa, "Alph.step")) alpha
  attr(imifa, "C.Shrink") <- is.element(method, c("MIFA", "OMIFA", "IMIFA"))         && MGP$cluster.shrink
  attr(imifa,
       "Class.Props")     <- if(method == "classify") tabulate(z.list[[1L]], range.G)/N
  attr(imifa, "Call")     <- call
  attr(imifa, "Center")   <- any(centered, centering)
  attr(imifa, "Cov.Emp")  <- covmat
  attr(imifa, "Clusters") <- range.G
  attr(imifa, "Dataset")  <- dat
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa,
       "Disc.step")       <- is.element(method, c("IMFA", "IMIFA"))    && learn.d
  attr(imifa, "Discount") <- if(is.element(method, c("IMFA", "IMIFA")) && !learn.d) discount
  attr(imifa, "Equal.Pi") <- equal.pro
  attr(imifa, "Exchange") <- BNP$exchange
  attr(imifa, "Factors")  <- range.Q
  attr(imifa, "ForceQg")  <- MGP$forceQg && is.element(method, c("MIFA", "OMIFA", "IMIFA"))
  attr(imifa, "G.init")   <- G.init
  attr(imifa, "G.Mean")   <- if(attr(imifa, "Center")) glo.mean
  attr(imifa, "G.Scale")  <- if(scaling != "none")     switch(EXPR=scaling, pareto=sqrt(glo.scal), glo.scal)
  attr(imifa, "IM.labsw") <- is.element(method, c("IMFA", "IMIFA"))    && IM.lab.sw
  attr(imifa,
       "Ind.Slice")       <- is.element(method, c("IMFA", "IMIFA"))    && BNP$ind.slice
  attr(imifa, "Init.Z")   <- if(!is.element(method, c("FA", "IFA")))      z.init
  attr(imifa, "Kappa0")   <- is.element(method, c("IMFA", "IMIFA"))    && kappa0
  attr(imifa,
       "Label.Switch")    <- if(!is.element(method, c("FA", "IFA")))      any(sw0gs)
  method                  <- names(table(meth)[which.max(table(meth))])
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1L, 1L)), substr(method, 2L, nchar(method)))
  attr(imifa, "Name")     <- dat.nam
  attr(imifa, "Obs")      <- N
  attr(imifa, "Obsnames") <- obsnames
  attr(imifa, "Pitman")   <- is.element(method,    c("IMFA", "IMIFA")) && any(learn.d, discount > 0)
  attr(imifa, "Rho")      <- if(is.element(method, c("IMFA", "IMIFA")))   BNP$rho
  attr(imifa, "Scale")    <- scal
  attr(imifa, "Scaling")  <- scaling
  attr(imifa, "Sd0.drop") <- if(isTRUE(mixFA$drop0sd) && any(sd0ind)) sd0ind
  attr(imifa, "Store")    <- length(iters)
  attr(imifa, "Switch")   <- c(storage, a.sw = attr(imifa, "Alph.step"), d.sw = attr(imifa, "Disc.step"))
  times                   <- lapply(list(Total = tot.time, Average = avg.time, Initialisation = init.time), round, 2L)
  if(all(len.G  == 1,
         len.Q  == 1)) {
    times                 <- times[-2L]
  }
  class(times)            <- "listof"
  attr(imifa, "Thresh")   <- BNP$thresh
  attr(imifa, "Time")     <- if(is.element(method, c("FA", "IFA", "classify"))) times[-length(times)] else times
  attr(imifa, "Truncate") <- is.element(method, c("classify", "IFA", "MIFA", "OMIFA", "IMIFA")) && isTRUE(truncate)
  attr(imifa, "TuneZeta") <- is.element(method, c("IMFA", "IMIFA", "OMFA", "OMIFA")) && tune.zeta$do
  attr(imifa, "Uni.Meth") <- c(Uni.Prior = uni.prior, Uni.Type = uni.type)
  attr(imifa, "Varnames") <- varnames
  attr(imifa, "Vars")     <- P
  if(!is.element(method, c("FA", "IFA"))) {
    for(g in seq_along(G.init)) {
      attr(imifa[[g]],
           "Z.init")      <- factor(zi[[g]], levels=seq_len(G.init[g]))
    }
  }
  if(isTRUE(verbose))        cat("\n"); print(attr(imifa, "Time"))
  class(imifa)            <- "IMIFA"
    return(imifa)
}
