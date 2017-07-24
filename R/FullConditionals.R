###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals

  # Means
    .sim_mu      <- function(N, P, mu.sigma, psi.inv, sum.data, sum.eta, lmat, mu.zero) {
      mu.omega   <- 1/(mu.sigma + N * psi.inv)
        mu.omega  * (psi.inv * (sum.data - lmat %*% sum.eta) + mu.sigma * mu.zero) + sqrt(mu.omega) * rnorm(P)
    }

  # Scores
    .sim_score   <- function(N, Q, lmat, psi.inv, c.data, Q1) {
      load.psi   <- lmat * psi.inv
      u.eta      <- diag(Q) + crossprod(load.psi, lmat)
      u.eta      <- if(Q1) sqrt(u.eta) else .chol(u.eta)
      mu.eta     <- c.data %*% (load.psi %*% if(Q1) 1/(u.eta * u.eta) else chol2inv(u.eta))
        mu.eta    + t(backsolve(u.eta, matrix(rnorm(length(mu.eta)), nrow=Q, ncol=N)))
    }

  # Loadings
    .sim_load    <- function(l.sigma, Q, c.data, eta, psi.inv, EtE, Q1)  {
      u.load     <- l.sigma + psi.inv * EtE
      u.load     <- if(Q1) sqrt(u.load) else .chol(u.load)
        psi.inv   * (if(Q1) 1/(u.load * u.load) else chol2inv(u.load)) %*% crossprod(eta, c.data) + backsolve(u.load, rnorm(Q))
    }

    .sim_load_s  <- function(Q, c.data, eta, phi, tau, psi.inv, EtE, Q1) {
      u.load     <- diag(phi * tau, Q) + psi.inv * EtE
      u.load     <- if(Q1) sqrt(u.load) else .chol(u.load)
        psi.inv   * (if(Q1) 1/(u.load  * u.load) else chol2inv(u.load)) %*% crossprod(eta, c.data) + backsolve(u.load, rnorm(Q))
    }

  # Uniquenesses
    .sim_psi_uu  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat, Q0) {
      S.mat      <- c.data  - if(Q0) tcrossprod(eta, lmat) else 0
        rgamma(P, shape=N/2 + psi.alpha, rate=colSums(S.mat * S.mat)/2 + psi.beta)
    }

    .sim_psi_uc  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat, Q0) {
      S.mat      <- c.data  - if(Q0) tcrossprod(eta, lmat) else 0
        rep(rgamma(1, shape=(N * P)/2 + psi.alpha, rate=sum(S.mat * S.mat)/2 + psi.beta), P)
    }

    .sim_psi_cu  <- function(u.shape, psi.beta, S.mat, V) {
      rgamma(V, shape=u.shape, rate=colSums(do.call(rbind, S.mat))/2 + psi.beta)
    }

    .sim_psi_cc  <- function(u.shape, psi.beta, S.mat, V = 1) {
      rgamma(V, shape=u.shape, rate=sum(unlist(S.mat))/2 + psi.beta)
    }

    .sim_psi_u1  <- function(u.shape, psi.beta, S.mat, V) {
      rgamma(V, shape=u.shape, rate=colSums(S.mat * S.mat)/2 + psi.beta)
    }

    .sim_psi_c1  <- function(u.shape, psi.beta, S.mat, V = 1) {
      rgamma(V, shape=u.shape, rate=sum(S.mat * S.mat)/2 + psi.beta)
    }

  # Local Shrinkage
    .sim_phi     <- function(Q, P, nu, tau, load.2, plus1) {
        base::matrix(rgamma(P * Q, shape=0.5 + nu + plus1, rate=(nu + sweep(load.2, 2, tau, FUN="*"))/2), nrow=P, ncol=Q)
    }

  # Global Shrinkage
    .sim_delta1  <- function(Q, P, alpha.d1, delta.1, beta.d1, tau, sum.term) {
        rgamma(1, shape=alpha.d1 + P * Q/2, rate=beta.d1 + 0.5/delta.1 * tau %*% sum.term)
    }

    .sim_deltak  <- function(Q, P, k, alpha.d2, beta.d2, delta.k, tau.kq, sum.term.kq) {
        rgamma(1, shape=alpha.d2 + P/2 * (Q - k + 1), rate=beta.d2 + 0.5/delta.k * tau.kq %*% sum.term.kq)
    }

  # Mixing Proportions
#' Simulate Mixing Proportions from a Dirichlet Distribution
#'
#' Generates samples from the Dirichlet distrubution with parameter \code{alpha} efficiently by simulating Gamma(\code{alpha}, 1) random variables and normalising them.
#' @param G The number of clusters for which weights need to be sampled.
#' @param alpha The Dirichlet hyperparameter, either of length 1 or \code{G}. When the length of \code{alpha} is 1, this amounts to assuming an exchangeable prior. Be warned that this will be recycled if necessary.
#' @param nn A vector giving the number of observations in each of G clusters so that Dirichlet posteriors rather than priors can be sampled from. This defaults to 0, i.e. simulation from the prior. Be warned that this will be recycled if necessary.
#'
#' @return A Dirichlet vector of \code{G} weights which sum to 1.
#'
#' @note Though the function is available for standalone use, note that no checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}.
#'
#' @references Devroye, L. (1986) \emph{Non-Uniform Random Variate Generation}, Springer-Verlag, New York, p. 594.
#' @export
#'
#' @examples
#' (prior     <- rDirichlet(G=5, alpha=1))
#' (posterior <- rDirichlet(G=5, alpha=1, nn=c(20, 41, 32, 8, 12)))
    rDirichlet   <- function(G, alpha, nn = 0) {
      tmp        <- rgamma(G, shape=alpha + nn, rate=1)
        tmp/sum(tmp)
    }

    .sim_vs_inf  <- function(alpha, nn = 0, N = sum(nn), discount, len, lseq = NULL) {
        if(discount == 0) rbeta(len, 1 + nn, alpha + N - cumsum(nn)) else rbeta(len, 1 - discount + nn, alpha + lseq * discount + N - cumsum(nn))
    }

    .sim_pi_inf  <- function(vs, len, init = 0) {
        vs * cumprod(1 - c(init, vs[-len]))
    }

  # Cluster Labels
#' Simulate Cluster Labels from Unnormalised Log-Probabilities using the Gumbel-Max Trick
#'
#' Samples cluster labels for N observations from G clusters efficiently using log-probabilities and the so-called Gumbel-Max trick, without requiring that the log-probabilities be normalised; thus redunant computation can be avoided. Computation takes place on the log scale for stability/underflow reasons (to ensure negligible probabilities won't cause computational difficulties); in any case, many functions for calculating multivariate normal densities already output on the log scale.
#' @param probs An N x G matrix of unnormalised probabilities on the log scale, where N is he number of observations that require labels to be sampled and G is the number of active clusters s.t. sampled labels can take values in \code{1:G}.
#' @param slice A logical indicating whether or not the indicator correction for slice sampling has been applied to \code{probs}. Defaults to \code{FALSE} but is \code{TRUE} for the "\code{IMIFA}" and "\code{IMFA}" methods under \code{\link{mcmc_IMIFA}}. Details of this correction are given in Murphy et. al. (2017). When set to \code{TRUE}, this results in a speed-improvement when \code{probs} contains non-finite values (e.g. \code{-Inf}, corresponding to zero on the probability scale).
#' @return A vector of N sampled cluster labels, with the largest label no greater than G.
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link[matrixStats]{rowLogSumExps}}
#'
#' @note Though the function is available for standalone use, note that no checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}.\cr
#' If the normalising constant is required for another reason, e.g. to compute the log-likelihood, it can be calculated by summing the output obtained by calling \code{\link[matrixStats]{rowLogSumExps}} on \code{probs}.
#'
#' @references Murphy, K., Gormley, I. C. and Viroli, C. (2017) Infinite Mixtures of Infinite Factor Analysers: Nonparametric Model-Based Clustering via Latent Gaussian Models, <\href{https://arxiv.org/abs/1701.07010}{arXiv:1701.07010}>.
#'
#' Yellot, J. I. Jr. (1977) The relationship between Luce's choice axiom, Thurstone's theory of comparative judgment, and the double exponential distribution, \emph{Journal of Mathematical Psychology}, 15: 109-144.
#' @export
#'
#' @author Keefe Murphy
#'
#' @examples
#' # Create weights for 3 components
#'   G         <- 3
#'   weights   <- seq_len(G)
#'
#' # Call gumbel_max() repeatedly to obtain samples of the labels, zs
#'   iters     <- 10000
#'   zs        <- vapply(seq_len(iters), function(i)
#'                gumbel_max(probs=log(weights)), numeric(1L))
#'
#' # Compare answer to the normalised weights
#'   tabulate(zs, nbins=G)/iters
#'   (normalised <- as.numeric(weights/sum(weights)))
#'
#' # Simulate a matrix of dirichlet weights & the associated vector of N labels
#'   N       <- 400
#'   G       <- 8
#'   sizes   <- seq(from=85, to=15, by=-10)
#'   weights <- matrix(rDirichlet(N * G, alpha=1, nn=sizes), byrow=TRUE, nrow=N, ncol=G)
#'   (zs     <- gumbel_max(probs=log(weights)))
    gumbel_max   <- function(probs, slice = FALSE) {
     if(is.vector(probs)) probs <- t(probs)
     if(isTRUE(slice)) {
      fps        <- is.finite(probs)
      probs[fps] <- probs[fps] - log(rexp(sum(fps)))
     } else   {
      probs      <- probs - log(rexp(length(probs)))
     }
      max.col(probs)
    }

  # Alpha
    .sim_alpha_g <- function(alpha, shape, rate, G, N) {
      shape2     <- shape  + G - 1
      rate2      <- rate   - log(rbeta(1, alpha + 1, N))
      weight     <- shape2/(shape2 + N * rate2)
        weight    * rgamma(1, shape=shape2 + 1, rate=rate2) + (1 - weight) * rgamma(1, shape=shape2, rate=rate2)
    }

    .log_palpha  <- function(alpha, discount, alpha.shape, alpha.rate, N, G) {
      l.prior    <- dgamma(alpha +  discount, shape=alpha.shape, rate=alpha.rate, log=TRUE)
        lgamma(alpha + 1)  - lgamma(alpha + N) + sum(log(alpha + discount  * seq_len(G - 1))) + l.prior
    }

    .sim_alpha_m <- function(alpha, discount, alpha.shape, alpha.rate, N, G, zeta) {
      inter      <- c(max( - discount, alpha   - zeta),  alpha + zeta)
      propa      <- runif(1, inter[1], inter[2])
      cprob      <- .log_palpha(alpha, discount, alpha.shape,    alpha.rate, N, G)
      pprob      <- .log_palpha(propa, discount, alpha.shape,    alpha.rate, N, G)
      propinter  <- c(max( - discount, propa   - zeta),  propa + zeta)
      logpr      <- pprob  - cprob   - log(diff(propinter))    + log(diff(inter))
      acpt       <- logpr >= 0  ||   - rexp(1) < logpr
        return(list(alpha  = ifelse(acpt, propa, alpha), rate  = acpt, l.prob = logpr))
    }

  # Adaptively Tune Zeta
    .tune_zeta   <- function(zeta, time, l.rate, heat = 1, target = 0.441, lambda = 1) {
      exp(heat/time^lambda * (exp(min(0, l.rate)) - target)) * zeta
    }

  # Discount
    .log_pdisc   <- function(discount, alpha, disc.shape1, disc.shape2, N, G, kappa, unif, nn) {
      l.prior    <- ifelse(discount == 0, log(kappa),  log1p(- kappa) + ifelse(unif, 0, dbeta(discount, shape1=disc.shape1, shape2=disc.shape2, log=TRUE)))
        sum(log(alpha + discount * seq_len(G - 1)))  + sum(lgamma(nn  - discount) -  lgamma(1  - discount)) + l.prior
    }

    .sim_disc_mh <- function(discount, alpha, disc.shape1, disc.shape2, N, G, kappa, unif, nn) {
      propd      <- ifelse(alpha > 0,  ifelse(kappa != 0 && runif(1) <= kappa, 0, ifelse(unif,
                           runif(1),   rbeta(1,  disc.shape1, disc.shape2))), runif(1, max(0,  - alpha), 1))
      if(identical(discount, propd)) {
          return(list(disc = discount, rate    = 0L))
      } else {
        cprob    <- .log_pdisc(discount,  alpha, disc.shape1, disc.shape2, N, G,  kappa, unif, nn)
        pprob    <- .log_pdisc(propd,     alpha, disc.shape1, disc.shape2, N, G,  kappa, unif, nn)
        logpr    <- pprob  - cprob
        acpt     <- logpr >= 0  ||   - rexp(1) < logpr
          return(list(disc = ifelse(acpt, propd, discount),   rate   = acpt))
      }
    }

# Priors
  # Means
    .sim_mu_p    <- function(P, mu.zero, sig.mu.sqrt) {
      sig.mu.sqrt * rnorm(P) + mu.zero
    }

  # Scores
    .sim_eta_p   <- function(Q, N) {
        matrix(rnorm(N * Q), nrow=N, ncol=Q)
    }

  # Loadings
    .sim_load_p  <- function(Q, P, sigma.l) {
        sqrt(sigma.l) * rnorm(P * Q)
    }

    .sim_load_ps <- function(Q, sigma.l, phi, tau) {
        sqrt(1/(phi * tau)) * rnorm(Q)
    }

  # Uniquenesses
    .sim_psi_ipu <- function(P, psi.alpha, psi.beta) {
        rgamma(n=P, shape=psi.alpha, rate=psi.beta)
    }

    .sim_psi_ipc <- function(P, psi.alpha, psi.beta) {
        rep(rgamma(1, shape=psi.alpha, rate=psi.beta), P)
    }

  # Local Shrinkage
    .sim_phi_p   <- function(Q, P, nu, plus1) {
        base::matrix(rgamma(n=P * Q, shape=nu + plus1, rate=nu), nrow=P, ncol=Q)
    }

  # Global Shrinkage
    .sim_delta_p <- function(Q = 2L, alpha, beta) {
        rgamma(n=Q - 1, shape=alpha, rate=beta)
    }

  # Cluster Labels
    .sim_z_p     <- function(N, prob.z) {
        which(rmultinom(N, size=1, prob=prob.z) != 0, arr.ind=TRUE)[,1]
    }

# Other Functions
  # Uniqueness Hyperparameters
#' Find sensible inverse gamma hyperparameters for variance/uniqueness parameters
#'
#' Takes a shape hyperparameter and covariance matrix, and finds data-driven rate hyperparameters in such a way that Heywood problems are avoided for factor analysis or probabilistic principal components analysis (and mixtures thereof).
#' @param shape A positive shape hyperparameter.
#' @param covar A square, positive-semidefinite covariance matrix.
#' @param type A switch indicating whether a single rate (\code{isotropic}) or variable-specific rates (\code{unconstrained}) are to be derived. The isotropic constraint provides the link between factor analysis and the probabilistic principal components analysis model. Uniquenesses are only allowed to be variable specific under the factor analysis model.
#'
#' @details Rates are allowed to be variable-specific or a single value under the factor analysis model, but \emph{must} be a single value for the PPCA model. Used internally by \code{\link{mcmc_IMIFA}} when its argument \code{psi_beta} is not supplied.
#'
#' @return Either a single rate hyperparameter or \code{ncol(covar)} variable specific hyperparameters.
#' @export
#'
#' @seealso \code{\link{mcmc_IMIFA}}
#' @references Fruwirth-Schnatter, S. and Lopes, H. F. (2010). Parsimonious Bayesian factor analysis when the number of factors is unknown, \emph{Technical Report}. The University of Chicago Booth School of Business.
#'
#' Tipping, M. E. and Bishop, C. M. (1999). Probabilistic principal component analysis, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 61(3): 611-622.
#'
#' @author Keefe Murphy
#'
#' @examples
#' data(olive)
#' olive2 <- olive[,-(1:2)]
#' (rates <- psi_hyper(shape=2.5, covar=cov(olive2), type="isotropic"))
#'
#' olive_scaled <- scale(olive2, center=TRUE, scale=TRUE)
#' (rate  <- psi_hyper(shape=3, covar=cov(olive_scaled), type="unconstrained"))
    psi_hyper   <- function(shape, covar, type=c("unconstrained", "isotropic")) {
      if(!all(is.posi_def(covar, semi=TRUE),
              is.symmetric(covar),
              is.double(covar)))           stop("Invalid covariance matrix supplied")
      if(any(!is.numeric(shape),
             length(shape) != 1))          stop("'shape' must be a single digit")
      inv.cov   <- try(base::solve(covar), silent=TRUE)
      if(inherits(inv.cov, "try-error"))   {
        covsvd  <- svd(covar)
        posi    <- covsvd$d > max(sqrt(.Machine$double.eps) * covsvd$d[1L], 0)
        inv.cov <- if(all(posi)) covsvd$v %*% (t(covsvd$u)/covsvd$d) else if(!any(posi))
                   array(0, dim(covar)[2L:1L]) else covsvd$v[,posi, drop=FALSE] %*% (t(covsvd$u[,posi, drop=FALSE])/covsvd$d[posi])
      }
        unname((shape - 1)/switch(match.arg(type), unconstrained=diag(inv.cov),
                                  isotropic=rep(exp(mean(log(diag(inv.cov)))), ncol(covar))))
    }

  # Alpha/Discount Shifted Gamma Hyperparameters
#' Moment Matching Parameters of Shifted Gamma Distributions
#'
#' This function takes shape and rate parameters of a Gamma distribution and modifies them to achieve the same expected value and variance when the left extent of the support of the distribution is shifted up or down.
#' @param shape Shape parameter a of a Gamma(a, b) distribution. Must be strictly positive.
#' @param rate Rate parameter b of a Gamma(a, b) distribution. Must be strictly positive.
#' @param shift Modifier, such that the Gamma distribution has support on (\code{shift}, \eqn{\infty}). Can be positive or negative, though typically negative and small.
#' @param param Switch controlling whether the supplied \code{rate} parameter is indeed a rate, or actually a scale parameter. Also governs whether the output is given in terms of rate or scale. Defaults to "\code{rate}".
#'
#' @return A list of length 2, containing the modified shape and rate parameters, respectively.
#' @export
#'
#' @author Keefe Murphy
#'
#' @examples
#' # Shift a Ga(shape=4, rate=2) distribution to the left by 1;
#' # achieving the same expected value of 2 and variance of 1.
#' shift_GA(4, 2, -1)
    shift_GA    <- function(shape, rate, shift = 0L, param = c("rate", "scale")) {
      if(length(shape) > 1 ||
        !is.numeric(shape) || shape <= 0) stop("Argument 'shape' must be a single strictly positive number")
      if(length(rate)  > 1 ||
        !is.numeric(rate)  || rate  <= 0) stop("Argument 'rate' must be a single strictly positive number")
      if(length(shift) > 1 ||
        !is.numeric(shift))               stop("Argument 'shift' must be a single number")
      param     <- match.arg(param)
      rate      <- switch(param, rate=rate, 1/rate)
      exp       <- shape/rate
      if(shift  >= exp)                   warning("This expected value is not achievable with the supplied 'shift'", call.=FALSE)
      var       <- exp/rate
      exp       <- pmax(var * rate   - shift, 0)
      rate      <- exp/var
      shape     <- rate * exp
        return(list(shape   = shape, rate = switch(param, rate=rate, 1/rate)))
    }

  # Check Shrinkage Hyperparemeters
#' Check the validity of Multiplicative Gamma Process (MGP) hyperparameters
#'
#' Checks the hyperparameters for the multiplicative gamma process (MGP) shrinkage prior in order to ensure that the property of cumulative shrinkage holds.
#' @param ad1 Shape hyperparameter for \eqn{\delta_1}{delta_1}.
#' @param ad2 Shape hyperparameter for \eqn{\delta_2}{delta_2}.
#' @param Q Number of latent factors. Defaults to 3, which is enough to check if the cumulative shrinkage property holds. Supply \code{Q} if the actual \emph{a priori} expected shrinkage factors are of interest.
#' @param nu Hyperparameter for the local shrinkage parameters. Defaults to 2. Not necessary for checking if the cumulative shrinkage property holds, but worth supplying if the actual \emph{a priori} expected shrinkage factors are of interest.
#' @param bd1 Rate hyperparameter for \eqn{\delta_1}{delta_1}. Defaults to 1.
#' @param bd2 Rate hyperparameter for \eqn{\delta_2}{delta_2}. Defaults to 1.
#' @param plus1 Logical indicator for whether the Gamma prior on the local shrinkage parameters is of the form Ga(\code{nu + 1, nu}), the default, or Ga(\code{nu, nu}).
#' @param inverse Logical indicator for whether the cumulative shrinkage property is assessed against the induced Inverse Gamma prior, the default, or in terms of the Gamma prior (which is incorrect). This is always \code{TRUE} when used inside \code{\link{mcmc_IMIFA}}: the \code{FALSE} option exists only for demonstration purposes.
#'
#' @details This is called inside \code{\link{mcmc_IMIFA}} for the "\code{IFA}", "\code{MIFA}", "\code{OMIFA}" and "\code{IMIFA}" methods. The arguments \code{ad1, ad2, nu, bd1} and \code{bd2} are vectorised.
#'
#' @return A list of length 2 containing the following objects:
#' \describe{
#'   \item{expectation}{The vector of actual expected shrinkage factors, \emph{a priori}.}
#'   \item{valid}{A logical indicating whether the cumulative shrinkage property holds.}
#' }
#' @export
#' @seealso \code{\link{mcmc_IMIFA}}
#' @references
#' Bhattacharya, A. and Dunson, D. B. (2011). Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291-306.
#'
#' Durante, D. (2017). A note on the multiplicative gamma process, \emph{Statistics & Probability Letters}, 122: 198-204.
#'
#' @author Keefe Murphy
#'
#' @examples
#' # Check if expected shrinkage under the MGP increases with the column index (WRONG approach!).
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, nu=2, inverse=FALSE)[[1]]$valid
#'
#' # Check if the induced IG prior on the MGP global shrinkage parameters
#' # is stochastically increasing, thereby inducing cumulative shrinkage (CORRECT approach!).
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, nu=2, inverse=TRUE)[[1]]$valid
#'
#' # Check again with a parameterisation that IS valid and examine the expected shrinkage values.
#' (shrink <- MGP_check(ad1=1.5, ad2=2.8, Q=10, nu=2, inverse=TRUE)[[1]])
    MGP_check   <- Vectorize(function(ad1, ad2, Q = 3, nu = 2, bd1 = 1L, bd2 = 1L, plus1 = TRUE, inverse = TRUE) {
      if(any(!is.logical(plus1),
             length(plus1)    != 1))       stop("'plus1' must be TRUE or FALSE")
      if(any(!is.logical(inverse),
             length(inverse)  != 1))       stop("'inverse' must be TRUE or FALSE")
      if(missing(ad1) || missing(ad2))     stop("Shrinkage shape hyperparameters 'ad1' and 'ad2' must be supplied")
      if(missing(nu))                      stop("Local shrinkage parameter 'nu' must be supplied")
      if(missing(Q))                       stop("Number of latent factors 'Q' must be supplied")
      if(any(nu <= !plus1,
             !is.numeric(nu)))             stop(paste0("'nu' must be a single ", ifelse(plus1,
                                                "strictly positive number for the Ga(nu + 1, nu) parameterisation",
                                                "number strictly greater than 1 for the Ga(nu, nu) parameterisation")))
      if(any(c(ad1, ad2)  < 1))            stop("All shrinkage shape hyperparameter values must be at least 1")
      if(any(c(bd1, bd2) <= 0))            stop("All shrinkage rate hyperparameter values must be strictly positive")
      rate      <- nu
      shape     <- ifelse(plus1, rate   + 1, rate)
      if(inverse) {
        ad1     <- ifelse(ad1 == 1, ad1 + .Machine$double.eps, ad1)
        ad2     <- ifelse(ad2 == 1, ad2 + .Machine$double.eps, ad2)
        exp.seq <- rate/(shape - 1) * bd1/(ad1 - 1) * (bd2/(ad2 - 1))^(seq_len(Q) - 1)
        check   <- is.unsorted(exp.seq)
      } else {
        exp.seq <- shape/rate * ad1/bd1 * (ad2/bd2)^(seq_len(Q) - 1)
        check   <- !is.unsorted(exp.seq)
      }
        return(list(expectation = exp.seq, valid = ifelse(Q < 2, TRUE, check)))
    }, vectorize.args = c("ad1", "ad2", "nu", "bd1", "bd2"), SIMPLIFY = FALSE)

  # Number of 'free' parameters
#' Estimate the Number of Free Parameters in Finite Factor Analytic Mixture Models (PGMM)
#'
#' Estimates the dimension of the 'free' parameters in fully finite factor analytic mixture models, otherwise known as Parsimonious Gaussian Mixture Models (PGMM), typically necessary for the penalty term of various model selection criteria.
#' @param Q The number of latent factors (which can be 0, corresponding to a model with diagonal covariance). This argument is vectorised.
#' @param P The number of variables.
#' @param G The number of clusters. This defaults to 1.
#' @param method By default, calculation assumes the \code{UUU} model with unconstrained loadings and unconstrained isotropic uniquesses. The other seven models detailed in McNicholas and Murphy (2008) are also given. The first letter denotes whether loadings are constrained/unconstrained across clusters; the second letter denotes the same for the uniquenesses; the final letter denotes whether uniquenesses are in turn constrained to be isotropic.
#' @param equal.pro Logical variable indicating whether or not the mixing mixing proportions are equal across clusters in the model (default = \code{FALSE}).
#'
#' @details This function is used to calculate the penalty terms for the \code{aic.mcmc} and \code{bic.mcmc} model selection criteria implemented in \code{\link{get_IMIFA_results}} for \emph{finite} factor models (though \code{\link{mcmc_IMIFA}} currently only implements the \code{UUU}, \code{UUC}, \code{UCU}, and \code{UCC} covariance structures). The function is vectorized with respect to the argument \code{Q}.

#' @return A vector of length \code{length(Q)}.
#'
#' @note Though the function is available for standalone use, note that no checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}.
#'
#' @export
#' @references McNicholas, P. D. and Murphy, T. B. (2008) Parsimonious Gaussian Mixture Models, \emph{Statistics and Computing}, 18(3): 285-296.
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link{mcmc_IMIFA}}
#'
#' @author Keefe Murphy
#'
#' @examples
#' (UUU <- PGMM_dfree(Q=4:5, P=50, G=3, method="UUU"))
#' (CCC <- PGMM_dfree(Q=4:5, P=50, G=3, method="CCC", equal.pro=TRUE))
    PGMM_dfree   <- Vectorize(function(Q, P, G = 1, method = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC"), equal.pro = FALSE) {
      if(length(equal.pro) > 1 ||
         !is.logical(equal.pro))           stop("'equal.pro' must be a single logical indicator")
      meth       <- unlist(strsplit(match.arg(method), ""))
      lambda     <- P * Q - 0.5 * Q * (Q - 1)
      lambda     <- switch(meth[1], C=lambda, U=G  * lambda)
      psi        <- switch(meth[2], C=1,      U=G)
      psi        <- switch(meth[3], C=1,      U=P) * psi
        as.integer(ifelse(equal.pro, 0, G - 1) + G * P + lambda + psi) },  vectorize.args = "Q")

  # Label Switching
    .lab_switch <- function(z.new, z.old) {
      tab       <- table(z.new, z.old, dnn=NULL)
      tab.tmp   <- tab[rowsums(tab) != 0,colSums(tab) != 0, drop=FALSE]
      Gs        <- seq_len(max(unique(as.numeric(z.new))))
      nc        <- ncol(tab.tmp)
      nr        <- nrow(tab.tmp)
      ng        <- tabulate(z.new, length(Gs))
      if(nc > nr) {
        tmp.mat <- matrix(rep(0, nc), nrow=nc - nr, ncol=nc)
        rownames(tmp.mat) <- setdiff(as.numeric(colnames(tab.tmp)), as.numeric(rownames(tab.tmp)))[seq_len(nc - nr)]
        tab.tmp <- rbind(tab.tmp, tmp.mat)
      } else if(nr > nc) {
        tmp.mat <- matrix(rep(0, nr), nrow=nr, ncol=nr - nc)
        colnames(tmp.mat) <- setdiff(as.numeric(rownames(tab.tmp)), as.numeric(colnames(tab.tmp)))[seq_len(nr - nc)]
        tab.tmp <- cbind(tab.tmp, tmp.mat)
      }
      if(nr == 1) {
        z.perm  <- setNames(as.numeric(colnames(tab.tmp)), as.numeric(colnames(tab.tmp)))
      } else if(nc == 1) {
        z.perm  <- setNames(as.numeric(colnames(tab.tmp)), as.numeric(colnames(tab.tmp)))
      } else {
        z.perm  <- tryCatch(suppressWarnings(matchClasses(tab.tmp, method="exact",  verbose=FALSE)),
          error=function(e) suppressWarnings(matchClasses(tab.tmp, method="greedy", verbose=FALSE)))
        z.perm  <- setNames(as.numeric(z.perm), names(z.perm))
      }
      if(length(Gs) > length(z.perm)) {
        z.perm  <- c(z.perm, setNames(setdiff(Gs, z.perm), setdiff(Gs, names(z.perm))))
      }
      z.names   <- as.numeric(names(z.perm))
      z.perm    <- z.perm[Order(z.names)]
      z.sw      <- factor(z.new, labels=z.perm[which(ng > 0)])
        return(list(z = as.integer(levels(z.sw))[z.sw], z.perm = z.perm))
    }

  # Similarity matrix and 'average' clustering
#' Summarise MCMC samples of clustering labels with a similarity matrix and find the 'average' clustering
#'
#' This function takes a Monte Carlo sample of cluster labels, computes an average similarity matrixand returns the clustering with minimum squared distance to this average.
#' @param zs A matrix containing samples of clustering labels where the columns correspond to the number of observations (N) and the rows correspond to the number of iterations (M).
#'
#' @details This function takes a Monte Carlo sample of cluster labels, converts them to adjacency matrices, and computes a similarity matrix as an average of the adjacency matrices. The dimension of the similarity matrix is invariant to label switching and the number of clusters in each sample, desirable features when summarising partitions of Bayesian nonparametric models such as IMIFA. As a summary of the posterior clustering, the clustering with minimum squared distance to this 'average' clustering is reported.\cr
#'
#' A heatmap of \code{z.sim} may provide a useful visualisation, if appropriately ordered. The user is also invited to perform hierarchical clustering using \code{\link[stats]{hclust}} after first converting this similarity matrix to a distance matrix - "complete" linkage is recommended.
#'
#' @return A list containing three elements:
#' \describe{
#' \item{z.avg}{The 'average' clustering, with minimum squared distance to \code{z.sim}.}
#' \item{z.sim}{The N x N similary matrix, in a sparse format (see \code{\link[slam]{simple_triplet_matrix}}).}
#' \item{dist.z}{A vector of length M recording the distances between each clustering and the 'average' clustering.}
#' }
#' @export
#' @importFrom mcclust "cltoSim" "comp.psm"
#'
#' @note This is liable to take quite some time to run, especially if the number of observations &/or number of iterations is large. Depending on how distinct the clusters are, \code{z.sim} may be stored better in a non-sparse format. This function can optionally be called inside \code{\link{get_IMIFA_results}}.
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link[slam]{simple_triplet_matrix}}, \code{\link[stats]{hclust}}, \code{\link[mcclust]{comp.psm}}, \code{\link[mcclust]{cltoSim}}
#'
#' @author Keefe Murphy
#'
#' @examples
#' # Run a IMIFA model and extract the sampled cluster labels
#' # data(olive)
#' # sim    <- mcmc_IMIFA(olive, method="IMIFA", n.iters=5000)
#' # zs     <- sim[[1]][[1]]$z.store
#'
#' # Get the similarity matrix and visualise it
#' # zsimil <- Zsimilarity(zs)
#' # z.sim  <- as.matrix(zsimil$z.sim)
#' # z.sim2 <- replace(z.sim, z.sim == 0, NA)
#' # image(z.sim2, col=heat.colors(30)[30:1]); box(lwd=2)
#'
#' # Extract the clustering with minimum squared distance to this
#' # 'average' and evaluate its performance against the true labels
#' # z.avg  <- zsimil$z.avg
#' # table(z.avg, olive$area)
#'
#' # Perform hierarchical clustering on the distance matrix
#' # Hcl    <- hclust(as.dist(1 - z.sim), method="complete")
#' # plot(Hcl)
#' # hier.z <- cutree(Hcl, k=3)
#' # table(hier.z, olive$area)
    Zsimilarity <- function(zs) {
      if(!is.matrix(zs))                   stop("'zs' must be a matrix with rows corresponding to the number of observations and columns corresponding to the number of iterations")
      zsim      <- comp.psm(zs)
      dist.z    <- vapply(seq_len(nrow(zs)), function(i, x=cltoSim(zs[i,]) - zsim) tryCatch(suppressWarnings(sum(x * x)), error=function(e) Inf), numeric(1L))
      Z.avg     <- zs[which.min(dist.z),]
      attr(Z.avg, "Distance")  <-  min(dist.z)
        return(list(z.avg = Z.avg, z.sim = as.simple_triplet_matrix(zsim), dist.z = dist.z))
    }

  # Move 1
    .lab_move1  <- function(nn.ind, pi.prop, nn) {
      sw        <- sample(nn.ind, 2L)
      log.pis   <- log(pi.prop[sw])
      nns       <- nn[sw]
      a.prob    <- (nns[1] - nns[2]) * (log.pis[1]     - log.pis[2])
        return(list(rate1  = a.prob >= 0 || - rexp(1)  < a.prob, sw = sw))
    }

  # Move 2
    .lab_move2  <- function(G, Vs, nn) {
      sw        <- sample(G, 1L, prob=c(rep(1, G - 2), 0.5, 0.5))
      sw        <- if(is.element(sw, c(G, G - 1))) c(G - 1, G) else c(sw, sw + 1)
      nns       <- nn[sw]
      log.vs    <- log1p( - Vs[sw])
      a.prob    <- nns[1] * log.vs[2]       - nns[2]   * log.vs[1]
      a.prob[is.nan(a.prob)]       <-       - Inf
        return(list(rate2 = a.prob >= 0  || - rexp(1)  < a.prob, sw = sw))
    }

  # Positive-(Semi)Definite Checker
    #' Check Postive-(Semi)definiteness of a matrix
    #'
    #' Tests whether all eigenvalues of a symmetric matrix are positive (or strictly non-negative) to check for positive-definiteness and positive-semidefiniteness, respectively. If the supplied matrix doesn't satisfy the test, the nearest matrix which does can optionally be returned.
    #' @param x A matrix, assumed to be real and symmetric.
    #' @param tol Tolerance for singular values and for absolute eigenvalues - only those with values larger than tol are considered non-zero (default: tol = \code{max(dim(x))*max(E)*.Machine$double.eps}, where \code{E} is the vector of absolute eigenvalues).
    #' @param semi Logical switch to test for positive-semidefiniteness when \code{TRUE} or positive-definiteness when \code{FALSE} (the default).
    #' @param make Logical switch to return the nearest matrix which satisifies the test - if the test has been passed, this is of course just \code{x} itself, otherwise the nearest positive-(semi)definite matrix. Note that for reasons due to finite precision arithmetic, finding the nearest positive-definite and nearest positive-semidefinite matrices are effectively equivalent tasks.
    #'
    #' @return If \code{isTRUE(make)}, a list with two components:
    #' \describe{
    #' \item{check}{A logical value indicating whether the matrix satisfies the test.}
    #' \item{X.new}{The nearest matrix which satisfies the test (which may just be the input matrix itself.)}
    #' }
    #' Otherwise, only the logical value indicating whether the matrix satisfies the test is returned.
    #'
    #' @export
    #'
    #' @examples
    #' x    <- cov(matrix(rnorm(100), nrow=10, ncol=10))
    #' is.posi_def(x)
    #' is.posi_def(x, semi=TRUE)
    #'
    #' Xnew <- is.posi_def(x, semi=FALSE, make=TRUE)$X.new
    #' identical(x, Xnew)
    #' identical(x, is.posi_def(x, semi=TRUE, make=TRUE)$X.new)
    is.posi_def <- function(x, tol = NULL, semi = FALSE, make = FALSE)  {
      if(!is.matrix(x)     &&
        nrow(x) != ncol(x))                stop("argument x is not a square matrix")
      if(!is.symmetric(x))                 stop("argument x is not a symmetric matrix")
      if(!is.double(x))                    stop("argument x is not a numeric matrix")
      if(!is.logical(semi) ||
         length(semi) > 1)                 stop("argument semi is not a single logical indicator")
      if(!is.logical(make) ||
         length(make) > 1)                 stop("argument make is not a single logical indicator")
      d         <- nrow(x)
      eigs      <- eigen(x, symmetric = TRUE)
      eval      <- eigs$values
      abseigs   <- abs(eval)
      tol       <- if(missing(tol)) max(abseigs) * d * .Machine$double.eps else tol
      if(length(tol)  > 1  ||
         !is.numeric(tol))                 stop("argument tol is not a single number")
      test      <- replace(eval, abseigs < tol, 0)
      check     <- !any(if(isTRUE(semi)) test < 0 else test <= 0)
      if(isTRUE(make))  {
        evec    <- eigs$vectors
        return(list(check = check, X.new = if(all(check)) x else x + evec %*% tcrossprod(diag(pmax(ifelse(isTRUE(semi), 0, .Machine$double.eps), 2 * tol - eval), d), evec)))
      } else check
    }

  # Ledermann Bound
#' Ledermann Bound
#'
#' Returns the maximum possible number of latent factors in a factor analysis model for data of dimension \code{P}. This Ledermann bound is given by the largest integer smaller than or equal to the solution \eqn{k}{k} of \eqn{(M - k)^2 \geq M + k}{(M - k)^2 >= M + k}.
#' @param P Integer number of variables in data set.
#'
#' @return The Ledermann bound, a non-negative integer.
#' @export
#'
#' @examples
#' Ledermann(25)
    Ledermann   <- function(P) {
      P         <- as.integer(P)
      if(length(P)   > 1  || P <= 0)       stop('argument P is a not a single positive integer')
      R         <- P + 0.5 * (1 - sqrt(8 * P  + 1))
        as.integer(floor(ifelse(1e-10 > abs(R - round(R)), round(R), R)))
    }

  # Procrustes Transformation
#' Procrustes Transformation
#'
#' This function performs a Procrustes transformation on a matrix \code{X} to minimize the squared distance between \code{X} and another comparable matrix \code{Xstar}.
#' @param X The matrix to be transformed.
#' @param Xstar The target matrix.
#' @param translate Logical value indicating whether \code{X} should be translated (defaults to \code{FALSE}).
#' @param dilate Logical value indicating whether \code{X} should be dilated (defaults to \code{FALSE}).
#' @param sumsq Logical value indicating whether the sum of squared differences between \code{X} and \code{Xstar} should be calculated and returned.
#'
#' @details{
#'    \code{R}, \code{tt}, and \code{d} are chosen so that:
#'
#'    \deqn{d \times \mathbf{X} \mathbf{R} + 1\hspace*{-3pt}1 \underline{t}^\top \approx X^\star}{d X R + 1 t' approximately Xstar}
#'
#'    \code{X.new} is given by:
#'
#'    \deqn{X_{\textrm{new}} = d \times \mathbf{X} \mathbf{R} + 1\hspace*{-3pt}1 \underline{t}^\top}{X.new = d X R + 1 t'}
#'}
#'
#' @return A list containing:
#' \describe{
#' \item{X.new}{The matrix that is the Procrustes transformed version of \code{X}.}
#' \item{R}{The rotation matrix.}
#' \item{t}{The translation vector (if \code{isTRUE(translate)}).}
#' \item{d}{The scaling factor (is \code{isTRUE(dilate)}).}
#' \item{ss}{The sum of squared differences (if \code{isTRUE(sumsq)}).}
#' }
#' @export
#'
#' @references Borg, I. and Groenen, P. J. F. (1997) \emph{Modern Multidimensional Scaling}. Springer-Verlag, New York, pp. 340-342.
#'
#' @examples
#' # Match two matrices, allowing translation and dilation
#' mat1     <- diag(rnorm(10))
#' mat2     <- 0.05 * matrix(rnorm(100), 10, 10) + mat1
#' proc     <- Procrustes(X=mat1, Xstar=mat2, translate=TRUE, dilate=TRUE, sumsq=TRUE)
#'
#' # Extract the transformed matrix, rotation matrix, translation vector and scaling factor
#' mat_new  <- proc$X.new
#' mat_rot  <- proc$R
#' mat_t    <- proc$t
#' mat_d    <- proc$d
#'
#' # Compare the sum of squared differences to a Procestean transformation with rotation only
#' mat_ss   <- proc$ss
#' mat_ss2  <- Procrustes(X=mat1, Xstar=mat2, sumsq=TRUE)$ss
    Procrustes  <- function(X, Xstar, translate = FALSE, dilate = FALSE, sumsq = FALSE) {
      if((N <- nrow(X)) != nrow(Xstar))    stop("X and Xstar do not have the same number of rows")
      if((P <- ncol(X)) != ncol(Xstar))    stop("X and Xstar do not have the same number of columns")
      J         <- if(translate) diag(N) - matrix(1/N, N, N)                         else diag(N)
      C         <- crossprod(Xstar, J) %*% X
      svdX      <- svd(C)
      R         <- tcrossprod(svdX$v, svdX$u)
      d         <- if(dilate)    sum(diag(C %*% R))/sum(diag(crossprod(X, J) %*% X)) else 1
      tt        <- if(translate) crossprod(Xstar - d * X %*% R, matrix(1, N, 1))/N   else 0
      X.new     <- d * X %*% R + if(translate) matrix(tt, N, P, byrow = TRUE)        else tt
        return(c(list(X.new = X.new), list(R = R), if(translate) list(t = tt),
                 if(dilate) list(d = d), if(sumsq) list(ss = sum((X - X.new)^2))))
    }

  # Length Checker
    .len_check  <- function(obj0g, switch0g, method, P, range.G, P.dim = TRUE) {
      V         <- ifelse(P.dim, P, 1L)
      rGseq     <- seq_along(range.G)
      obj.name  <- deparse(substitute(obj0g))
      sw.name   <- deparse(substitute(switch0g))
      if(!is.list(obj0g))        obj0g  <- list(obj0g)
      if(length(obj0g) != length(range.G))    {
        if(!P.dim)             {
          obj0g <- replicate(length(range.G), obj0g)
        } else                             stop(paste0(obj.name, " must be a list of length ", length(range.G)))
      }
      len       <- lengths(obj0g)
      if(is.element(method, c("FA", "IFA")))  {
        if(any(!is.element(len, c(1, V)))) stop(paste0(obj.name, " must be list of length 1 containing a scalar", ifelse(P.dim, paste0(" or a vector of length P=", V), ""), " for a 1-group model"))
      } else {
        if(any(is.element(len,
           c(1, range.G, V)))) {
          if(all(len == range.G)) obj0g <- if(switch0g) lapply(rGseq, function(g) matrix(obj0g[[g]], nrow=1))  else stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"))
          if(all(len == V))       obj0g <- if(V == 1)   lapply(rGseq, function(g) rep(obj0g[[g]], range.G[g])) else lapply(rGseq, function(g) matrix(obj0g[[g]], nrow=V, ncol=range.G[g]))
        } else if(!all(vapply(rGseq, function(g) is.matrix(obj0g[[g]]) && any(identical(dim(obj0g[[g]]), c(1, range.G[g])), identical(dim(obj0g[[g]]), c(V, range.G[g]))), logical(1L)))) {
                                           stop(paste0(ifelse(length(range.G) > 1, "Each element of ", ""), obj.name, " must be either of length 1, ", ifelse(P.dim, paste0("P=", V, ", or it's corresponding range.G, or a matrix with P rows and it's corresponding range.G columns"), paste0("or G=", range.G))))
        } else if(all(vapply(obj0g, is.matrix, logical(1L)), !switch0g) && any(vapply(rGseq, function(g) any(dim(obj0g[[g]]) == range.G[g]), logical(1L)))) {
                                           stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"))
        }
      }
      if(all(length(unique(unlist(obj0g))) > 1,
             !switch0g, !P.dim))           stop(paste0(obj.name, " must be a scalar if ", sw.name, " is FALSE"))
        obj0g
    }

  # Moments of Dirichlet / Pitman-Yor Processes
#' 1st Moment of the Dirichlet / Pitman-Yor processes
#'
#' Calculates the expected number of clusters under a Dirichlet process or Pitman-Yor process prior for a sample of size \code{N} at given values of the concentration parameter \code{alpha} and optionally also the \code{discount} parameter. Useful for soliciting sensible priors (or fixed values) for \code{alpha} or \code{discount} under the "\code{IMFA}" and "\code{IMIFA}" methods for \code{\link{mcmc_IMIFA}}.
#' @param N The sample size.
#' @param alpha The concentration parameter. Must be specified and must be strictly greater than \code{-discount}.
#' @param discount The discount parameter for the Pitman-Yor process. Must lie in the interval [0, 1). Defaults to 0 (i.e. the Dirichlet process).
#'
#' @details All arguments are vectorised. Users can also consult \code{\link{G_variance}} and \code{\link{G_priorDensity}} in order to solicit sensible priors.
#'
#' @return The expected number of clusters under the specified prior conditions.
#' @export
#'
#' @note Requires use of the \code{Rmpfr} and \code{gmp} libraries for non-zero \code{discount} values.
#'
#' @seealso \code{\link{G_variance}}, \code{\link{G_priorDensity}}, \code{\link[Rmpfr]{Rmpfr}}
#'
#' @author Keefe Murphy
#'
#' @examples
#' G_expected(N=50, alpha=19.23356)
#'
#' # require("Rmpfr")
#' # G_expected(N=50, alpha=c(19.23356, 12.21619, 1), discount=c(0, 0.25, 0.7300045))
    G_expected  <- Vectorize(function(N, alpha, discount = 0L) {
      if(!all(is.numeric(N), is.numeric(discount),
         is.numeric(alpha)))               stop("All inputs must be numeric")
      if(discount  < 0  || discount >= 1)  stop("'discount' must lie in the interval [0,1)")
      if(alpha  <= - discount)             stop("'alpha' must be strictly greater than -discount")
      if(alpha  == 0)                      stop("'alpha' equal to zero not yet implemented")
      if(suppressMessages(requireNamespace("Rmpfr", quietly=TRUE))) {
        mpfrind <- TRUE
        on.exit(.detach_pkg("Rmpfr"))
        on.exit(.detach_pkg("gmp"), add=TRUE)
        alpha   <- Rmpfr::mpfr(alpha, precBits=256)
      } else if(discount != 0)             stop("'Rmpfr' package not installed")
      if(discount == 0) {
        exp     <- alpha * (digamma(alpha + N) - digamma(alpha))
        if(mpfrind)     {
          gmp::asNumeric(exp)
        } else {
          exp
        }
      } else {
        adx     <- alpha/discount
          gmp::asNumeric(adx * Rmpfr::pochMpfr(alpha + discount, N)/Rmpfr::pochMpfr(alpha, N) - adx)
      }
    })

#' 2nd Moment of Dirichlet / Pitman-Yor processes
#'
#' Calculates the variance in the number of clusters under a Dirichlet process or Pitman-Yor process prior for a sample of size \code{N} at given values of the concentration parameter \code{alpha} and optionally also the \code{discount} parameter. Useful for soliciting sensible priors (or fixed values) for \code{alpha} or \code{discount} under the "\code{IMFA}" and "\code{IMIFA}" methods for \code{\link{mcmc_IMIFA}}.
#' @inheritParams G_expected
#' @details All arguments are vectorised. Users can also consult \code{\link{G_expected}} or \code{\link{G_priorDensity}} in order to solicit sensible priors.
#'
#' @return The variance of the number of clusters under the specified prior conditions.
#' @export
#'
#' @note Requires use of the \code{Rmpfr} and \code{gmp} libraries for non-zero \code{discount} values.
#' @seealso \code{\link{G_expected}}, \code{\link{G_priorDensity}}, \code{\link[Rmpfr]{Rmpfr}}
#'
#' @author Keefe Murphy
#'
#' @examples
#' G_variance(N=50, alpha=19.23356)
#'
#' # require("Rmpfr")
#' # G_variance(N=50, alpha=c(19.23356, 12.21619, 1), discount=c(0, 0.25, 0.7300045))
    G_variance  <- Vectorize(function(N, alpha, discount = 0L) {
      if(!all(is.numeric(N), is.numeric(discount),
         is.numeric(alpha)))               stop("All inputs must be numeric")
      if(discount  < 0  || discount >= 1)  stop("'discount' must lie in the interval [0,1)")
      if(alpha  <= - discount)             stop("'alpha' must be strictly greater than -discount")
      if(alpha  == 0)                      stop("'alpha' equal to zero not yet implemented")
      if(suppressMessages(requireNamespace("Rmpfr", quietly=TRUE))) {
        mpfrind <- TRUE
        on.exit(.detach_pkg(Rmpfr))
        on.exit(.detach_pkg(gmp), add=TRUE)
        alpha   <- Rmpfr::mpfr(alpha, precBits=256)
      } else if(discount != 0)             stop("'Rmpfr' package not installed")
      alpha2    <- alpha * alpha
      if(discount == 0) {
        var     <- alpha  * (digamma(alpha + N) - digamma(alpha))
        if(mpfrind)     {
          alpha <- gmp::asNumeric(alpha)
          gmp::asNumeric(var + alpha2 * (trigamma(alpha + N) - trigamma(alpha)))
        } else {
          var + alpha2 * (trigamma(alpha + N) - trigamma(alpha))
        }
      } else {
        sum.ad  <- alpha + discount
        poch.a  <- Rmpfr::pochMpfr(alpha, N)
        poch.ad <- Rmpfr::pochMpfr(sum.ad, N)
        subterm <- alpha/discount * poch.ad/poch.a
          gmp::asNumeric((alpha * sum.ad)/(discount * discount) * Rmpfr::pochMpfr(sum.ad + discount, N)/poch.a - subterm - subterm * subterm)
      }
    })

  # Print functions
#' @method print IMIFA
#' @export
    print.IMIFA <- function(x, ...) {
      meth      <- attr(x, "Method")
      name      <- attr(x, "Name")
      fac       <- attr(x, "Factors")
      grp       <- attr(x, "Clusters")
      Qmsg      <- Gmsg <- msg   <- NULL
      for(i in seq_along(fac[-length(fac)])) {
        Qmsg    <- c(Qmsg, (paste0(fac[i], ifelse(i + 1 < length(fac), ", ", " "))))
      }
      for(i in seq_along(grp[-length(grp)])) {
        Gmsg    <- c(Gmsg, (paste0(grp[i], ifelse(i + 1 < length(grp), ", ", " "))))
      }
      Qmsg      <- if(length(fac) > 1) paste(c(Qmsg, paste0("and ", fac[length(fac)])), sep="", collapse="") else fac
      Gmsg      <- if(length(grp) > 1) paste(c(Gmsg, paste0("and ", grp[length(grp)])), sep="", collapse="") else grp
      Qmsg      <- paste0(" with ", Qmsg, " factor", ifelse(length(fac) == 1, "", "s"))
      Gmsg      <- paste0(" with ", Gmsg, " group",  ifelse(length(grp) == 1, "", "s"))
      if(is.element(meth, c("FA", "OMFA", "IMFA"))) {
        msg     <- Qmsg
      } else {
        msg     <- switch(meth, MFA=paste0(Gmsg, " and", Qmsg), MIFA=Gmsg)
      }
        cat(paste0(meth, " simulations for '", name, "' dataset", msg, " to be passed to get_IMIFA_results(...)\n"))
    }
#' @method summary IMIFA
#' @export
    summary.IMIFA        <- function(object, ...) {
      meth      <- attr(object, "Method")
      name      <- attr(object, "Name")
      call      <- attr(object, "Call")
      fac       <- attr(object, "Factors")
      grp       <- attr(object, "Clusters")
      Qmsg      <- Gmsg <- msg   <- NULL
      for(i in seq_along(fac[-length(fac)])) {
        Qmsg    <- c(Qmsg, (paste0(fac[i], ifelse(i + 1 < length(fac), ", ", " "))))
      }
      for(i in seq_along(grp[-length(grp)])) {
        Gmsg    <- c(Gmsg, (paste0(grp[i], ifelse(i + 1 < length(grp), ", ", " "))))
      }
      Qmsg      <- if(length(fac) > 1) paste(c(Qmsg, paste0("and ", fac[length(fac)])), sep="", collapse="") else fac
      Gmsg      <- if(length(grp) > 1) paste(c(Gmsg, paste0("and ", grp[length(grp)])), sep="", collapse="") else grp
      Qmsg      <- paste0(" with ", Qmsg, " factor", ifelse(length(fac) == 1, "", "s"))
      Gmsg      <- paste0(" with ", Gmsg, " group",  ifelse(length(grp) == 1, "", "s"))
      if(is.element(meth, c("FA", "OMFA", "IMFA"))) {
        msg     <- Qmsg
      } else {
        msg     <- switch(meth, MFA=paste0(Gmsg, " and", Qmsg), MIFA=Gmsg)
      }
        cat("Call:\t"); print(call); cat("\n")
        cat(paste0(meth, " simulations for '", name, "' dataset", msg, " to be passed to get_IMIFA_results(...)\n"))
    }

#' @method print Results_IMIFA
#' @export
    print.Results_IMIFA  <- function(x, ...) {
      method    <- attr(x, "Method")
      G         <- x$GQ.results$G
      Q         <- x$GQ.results$Q
      if(is.element(method, c("FA", "IFA")))  {
        msg     <- paste0("The chosen ", method, " model has ", Q, " factor", ifelse(Q == 1, "", "s"))
      } else if(is.element(method, c("MFA", "OMFA", "IMFA"))) {
        msg     <- paste0("The chosen ", method, " model has ", G, " group",  ifelse(G == 1, " with ", "s, each with "), unique(Q), " factor", ifelse(unique(Q) == 1, "", "s"))
      } else {
        Q.msg   <- NULL
        for(i in seq_along(Q[-length(Q)])) {
          Q.msg <- c(Q.msg, (paste0(Q[i], ifelse(i + 1 < length(Q), ", ", " "))))
        }
        Q.msg   <- if(length(Q) > 1) paste(c(Q.msg, paste0("and ", Q[length(Q)])), sep="", collapse="") else Q
        msg     <- paste0("The chosen ", method, " model has ", G, " group",  ifelse(G == 1, " with ", "s, with "), Q.msg, " factor", ifelse(G == 1 && Q == 1, "", paste0("s", ifelse(G == 1, "", " respectively"))), sep="")
      }
        cat(paste0(msg, ": this Results_IMIFA object can be passed to plot(...)\n"))
    }

#' @method summary Results_IMIFA
#' @export
    summary.Results_IMIFA <- function(object, ...) {
      criterion <- unlist(strsplit(toupper(attr(object$GQ.results, "Criterion")), "[.]"))
      criterion <- ifelse(length(criterion) > 1, ifelse(criterion[1] != "LOG", paste0(criterion[1], ".", tolower(criterion[2])), "LogIntegratedLikelihood"), criterion)
      crit.mat  <- object$GQ.results[[paste0(criterion, "s")]]
      call      <- attr(object, "Call")
      msg       <- NULL
      if(any(dim(crit.mat) > 1)) {
        msg     <- paste0(", and ", ifelse(substr(criterion, 1, 1) == "A", "an ", "a "),  criterion, " of ", round(max(crit.mat), 2), "\n")
      }
        cat("Call:\t"); print(call); cat("\n")
        cat(paste0(capture.output(print.Results_IMIFA(object)), msg))
    }

  # Control functions
#' Set storage values for use with the IMIFA family of models
#'
#' Supplies a list of values for logical switches indicating whether parameters of interest (means, scores, loadings, uniquenesses, and mixing proportions) should be stored when running models from the IMIFA family via \code{\link{mcmc_IMIFA}}.
#' @param load.switch Logical indicating whether the factor loadings are to be stored (defaults to \code{TRUE}).
#' @param mu.switch Logical indicating whether the means are to be stored (defaults to \code{TRUE}).
#' @param pi.switch Logical indicating whether the mixing proportions are to be stored (defaults to \code{TRUE}).
#' @param psi.switch Logical indicating whether the uniquenesses are to be stored (defaults to \code{TRUE}).
#' @param score.switch Logical indicating whether the factor scores are to be stored. As the array containing each sampled scores matrix tends to be amongst the largest objects to be stored, this defaults to \code{FALSE} inside \code{\link{mcmc_IMIFA}} when \code{length(range.G) * length(range.Q) > 10}, otherwise the default is \code{TRUE}. For the "\code{MIFA}", "\code{OMIFA}", and "\code{IMIFA}" methods, setting this switch to \code{FALSE} also offers a slight speed-up.
#'
#' @details \code{\link{storeControl}} is provided for assigning values for IMIFA models within \code{\link{mcmc_IMIFA}}. It may be useful not to store certain parameters if memory is an issue. Warning: posterior inference and plotting won't be posssible for parameters not stored. In particular, when loadings and uniquenesses are not stored, it will not be possible to estimate covariance matrices and compute error metrics.
#'
#' @note Further warning messages may appear when \code{\link{mcmc_IMIFA}} is called depending on the particularities of the data set and the IMIFA method employed etc. as additional checks occur.

#' @return A named list in which the names are the names of the storage switches and the values are logicals indicating whether that parameter is to be stored. The list also contains as an attribute a logical for each switch indicating whether it was actually supplied (\code{TRUE}) or the default was accepted (\code{FALSE}).
#' @export
#'
#' @seealso \code{\link{mcmc_IMIFA}}
#' @author Keefe Murphy
#' @examples
#' storeControl(score.switch=FALSE)
   storeControl <- function(load.switch = TRUE, mu.switch = TRUE, pi.switch = TRUE, psi.switch = TRUE, score.switch = TRUE) {
      switches  <- c(l.sw=load.switch, mu.sw=mu.switch, pi.sw=pi.switch, psi.sw=psi.switch, s.sw=score.switch)
      attr(switches, "Missing") <- c(l.sw=missing(load.switch), mu.sw=missing(mu.switch), pi.sw=missing(pi.switch), psi.sw=missing(psi.switch), s.sw=missing(score.switch))
      if(any(length(switches)   != 5,
             !is.logical(switches)))       stop("All logical parameter storage switches must be TRUE or FALSE")
        switches
   }

  # Other Hidden Functions
    .chol       <- function(x) tryCatch(chol(x), error=function(e) {
      d         <- nrow(x)
      eigs      <- eigen(x, symmetric = TRUE)
      eval      <- eigs$values
      evec      <- eigs$vectors
        return(chol(x + evec %*% tcrossprod(diag(pmax(0, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec)))
      }
    )

    .detach_pkg <- function(pkg, character.only = FALSE) {
      searches  <- paste("package", if(!character.only) deparse(substitute(pkg)) else pkg, sep=":")
      while(searches %in% search()) {
        detach(searches, unload=TRUE, character.only=TRUE)
      }
    }

    .empty_mat <- function(nr) {
      base::matrix(0L, nrow=nr, ncol=0)
    }

    .ent_exit  <- function(opts = options()) {
      ent      <- readline("Hit <Return> to see next plot or type 'EXIT'/hit <Esc> to exit: ")
      options(show.error.messages=FALSE)
      on.exit(suppressWarnings(options(opts)), add=TRUE)
        if(ent  %in% c("exit", "EXIT"))    stop()
    }

    .logdensity     <- function(x,  left = 0) { # export?
      d        <- tryCatch(density(x, bw = "SJ"),  error  = function(e) density(x))
      h        <- d$bw
      w        <- 1/pnorm(left,  mean = x, sd = h, lower.tail = FALSE)
        return(suppressWarnings(density(x, bw = h, kernel = "gaussian", weights = w/length(x))))
    }

    .logitdensity   <- function(x)  { # export?
      y         <- qlogis(x[x  > 0  &   x < 1])
      g         <- tryCatch(density(y, bw = "SJ"), error  = function(e) density(y))
      xgrid     <- plogis(g$x)
      g$y       <- g$y/(xgrid  * (1 - xgrid))
      g$x       <- xgrid
        return(g)
    }

    .power2     <- function(x) x * x

    .which0     <- function(x) which(x == 0)
    #
