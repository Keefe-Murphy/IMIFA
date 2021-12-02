###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals

  # Means
    .sim_mu      <- function(N, P, mu.sigma, psi.inv, sum.data, sum.eta, lmat, mu.prior) {
      mu.omega   <- sqrt(mu.sigma + N * psi.inv)
        (psi.inv/mu.omega * (sum.data - lmat %*% sum.eta) + mu.prior + stats::rnorm(P))/mu.omega
    }

  # Scores
    .sim_score   <- function(N, Q, lmat, psi.inv, c.data, Q1) {
      load.psi   <- lmat * psi.inv
      if(Q1)      {
        u.eta2   <- sum(load.psi * lmat) + 1
        u.eta    <- sqrt(u.eta2)
        mu.eta   <- (c.data %*% load.psi)/u.eta2
          mu.eta  + stats::rnorm(Q)/u.eta
      } else      {
        u.eta    <- diag(Q) + crossprod(load.psi, lmat)
        u.eta    <- .chol(u.eta)
        mu.eta   <- c.data %*% (load.psi %*% chol2inv(u.eta))
          mu.eta  + t(backsolve(u.eta, matrnorm(Q, N), k=Q))
      }
    }

  # Loadings
    .sim_load    <- function(l.sigma, Q, c.data, eta, psi.inv, EtE, Q1)  {
      u.load     <- l.sigma  + psi.inv * EtE
      if(Q1)      {
        u.load   <- sqrt(u.load)
          (psi.inv/u.load * crossprod(eta, c.data) + stats::rnorm(Q))/u.load
      } else      {
        u.load   <- .chol(u.load)
          backsolve(u.load, backsolve(u.load, psi.inv * crossprod(eta, c.data), transpose=TRUE, k=Q) + stats::rnorm(Q), k=Q)
      }
    }

    .sim_load_s  <- function(Q, c.data, eta, phi, tau, psi.inv, EtE, Q1, sigma = 1L) {
      u.load     <- diag(phi * tau * sigma, Q) + psi.inv * EtE
      if(Q1)      {
        u.load   <- sqrt(u.load)
          (psi.inv/u.load * crossprod(eta, c.data) + stats::rnorm(Q))/u.load
      } else      {
        u.load   <- .chol(u.load)
          backsolve(u.load, backsolve(u.load, psi.inv * crossprod(eta, c.data), transpose=TRUE, k=Q) + stats::rnorm(Q), k=Q)
      }
    }

  # Uniquenesses
  #' @importFrom matrixStats "colSums2"
    .sim_psi_uu  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat, Q0) {
      S.mat      <- c.data  - if(Q0) tcrossprod(eta, lmat) else 0L
        stats::rgamma(P, shape=N/2L + psi.alpha, rate=colSums2(S.mat^2)/2 + psi.beta)
    }

    .sim_psi_uc  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat, Q0) {
      S.mat      <- c.data  - if(Q0) tcrossprod(eta, lmat) else 0L
        rep(stats::rgamma(1, shape=(N * P)/2L + psi.alpha, rate=sum(S.mat^2)/2 + psi.beta), P)
    }

  #' @importFrom matrixStats "colSums2"
    .sim_psi_cu  <- function(u.shape, psi.beta, S.mat, V) {
        stats::rgamma(V, shape=u.shape, rate=colSums2(do.call(rbind, S.mat))/2 + psi.beta)
    }

    .sim_psi_cc  <- function(u.shape, psi.beta, S.mat, V = 1L) {
        stats::rgamma(V, shape=u.shape, rate=sum(unlist(S.mat))/2 + psi.beta)
    }

  #' @importFrom matrixStats "colSums2"
    .sim_psi_u1  <- function(u.shape, psi.beta, S.mat, V) {
        stats::rgamma(V, shape=u.shape, rate=colSums2(S.mat^2)/2  + psi.beta)
    }

    .sim_psi_c1  <- function(u.shape, psi.beta, S.mat, V = 1L) {
        stats::rgamma(V, shape=u.shape, rate=sum(S.mat^2)/2 + psi.beta)
    }

  # Local Shrinkage
    .sim_phi     <- function(Q, P, nu1.5, nu2, tau, load.2, sigma = 1L) {
        matrix(stats::rgamma(P * Q, shape=nu1.5, rate=nu2 + (sigma * sweep(load.2, 2L, tau, FUN="*", check.margin=FALSE))/2), nrow=P, ncol=Q)
    }

  # Column Shrinkage
    .sim_delta1  <- function(Q, P.5, alpha.d1, delta.1, beta.d1, tau, sum.term, sigma = 1L) {
        stats::rgamma(1, shape=alpha.d1 + P.5 * Q, rate=beta.d1 + (sigma * 0.5)/delta.1 * tau %*% sum.term)
    }

    .sim_deltak  <- function(Q, P.5, k, alpha.d2, beta.d2, delta.k, tau.kq, sum.term.kq, sigma = 1L) {
        stats::rgamma(1, shape=alpha.d2 + P.5 * (Q - k + 1L), rate=beta.d2 + (sigma * 0.5)/delta.k * tau.kq %*% sum.term.kq)
    }

    .sim_deltaKT <- function(Q, P.5, k, alpha.d2, beta.d2, delta.k, tau.kq, sum.term.kq, sigma = 1L) {
        rltrgamma(1,     shape=alpha.d2 + P.5 * (Q - k + 1L), rate=beta.d2 + (sigma * 0.5)/delta.k * tau.kq %*% sum.term.kq)
    }

  # Cluster Shrinkage
    .sim_sigma   <- function(G, P.5, Qs, rho1, rho2, sum.terms, tau) {
        stats::rgamma(G, shape=rho1 + P.5 * Qs, rate=rho2 + mapply("%*%", sum.terms, tau)/2)
    }

  # Mixing Proportions
#' Simulate Mixing Proportions from a Dirichlet Distribution
#'
#' Generates samples from the Dirichlet distribution with parameter \code{alpha} efficiently by simulating Gamma(\code{alpha}, 1) random variables and normalising them.
#' @param G The number of clusters for which weights need to be sampled.
#' @param alpha The Dirichlet hyperparameter, either of length 1 or \code{G}. When the length of \code{alpha} is 1, this amounts to assuming an exchangeable prior, which doesn't favour one component over another. Be warned that this will be recycled if necessary. Larger values have the effect of making the returned samples more equal.
#' @param nn A vector giving the number of observations in each of G clusters so that Dirichlet posteriors rather than priors can be sampled from. This defaults to 0, i.e. simulation from the prior. Must be non-negative. Be warned that this will be recycled if necessary.
#'
#' @return A Dirichlet vector of \code{G} weights which sum to 1.
#'
#' @note Though the function is available for standalone use, note that few checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}. In particular, \code{alpha} and \code{nn} may be invisibly recycled.
#'
#' While small values of \code{alpha} have the effect of increasingly concentrating the mass onto fewer components, note that this function may return \code{NaN} for excessively small values of \code{alpha}, when \code{nn=0}; see the details of \code{rgamma} for small \code{shape} values.
#'
#' @references Devroye, L. (1986) \emph{Non-Uniform Random Variate Generation}, Springer-Verlag, New York, p. 594.
#' @keywords utility
#' @export
#' @usage
#' rDirichlet(G,
#'            alpha,
#'            nn = 0L)
#' @examples
#' (prior      <- rDirichlet(G=5, alpha=1))
#' (posterior  <- rDirichlet(G=5, alpha=1, nn=c(20, 41, 32, 8, 12)))
#' (asymmetric <- rDirichlet(G=5, alpha=c(3,4,5,1,2), nn=c(20, 41, 32, 8, 12)))
    rDirichlet   <- function(G, alpha, nn = 0L) {
      if(!all(0   < alpha))                stop("'alpha' must be strictly positive",  call.=FALSE)
      if(!all(0  <= nn))                   stop("'nn' must be strictly non-negative", call.=FALSE)
      shape      <- alpha + nn
      tmp        <- if(all(shape == 1)) stats::rexp(G, rate=1L) else stats::rgamma(G, shape=shape, rate=1L)
        tmp/sum(tmp)
    }

    .sim_vs_inf  <- function(alpha, nn = 0L, N = sum(nn), discount = 0, len, lseq = NULL) {
        if(discount == 0) stats::rbeta(len, 1L + nn, alpha + N    - cumsum(nn))                             else
        if(discount  > 0) stats::rbeta(len, 1  - discount  + nn, alpha + lseq * discount  + N - cumsum(nn)) else {
          stats::rbeta(len, 1 - discount + nn, pmax(alpha  + lseq * discount  + N - cumsum(nn), 0L))
         #stats::rbeta(len, alpha + nn, (if(len == 1) 1 else len - lseq) * alpha + N - cumsum(nn))
        }
    }

    .sim_pi_inf  <- function(vs, len, prev.prod = 0) {
        vs * cumprod(1 - c(prev.prod, vs[-len]))
    }

  # Cluster Labels
#' Simulate Cluster Labels from Unnormalised Log-Probabilities using the Gumbel-Max Trick
#'
#' Samples cluster labels for N observations from G clusters efficiently using log-probabilities and the so-called Gumbel-Max trick, without requiring that the log-probabilities be normalised; thus redundant computation can be avoided.
#' @param probs An N x G matrix of unnormalised probabilities on the log scale, where N is he number of observations that require labels to be sampled and G is the number of active clusters s.t. sampled labels can take values in \code{1:G}. Typically \code{N > G}.
#' @param slice A logical indicating whether or not the indicator correction for slice sampling has been applied to \code{probs}. Defaults to \code{FALSE} but is \code{TRUE} for the \code{"IMIFA"} and \code{"IMFA"} methods under \code{\link{mcmc_IMIFA}}. Details of this correction are given in Murphy et. al. (2020). When set to \code{TRUE}, this results in a speed-improvement when \code{probs} contains non-finite values (e.g. \code{-Inf}, corresponding to zero on the probability scale).
#' @return A vector of N sampled cluster labels, with the largest label no greater than G.
#' @importFrom Rfast "rowMaxs"
#' @keywords utility
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link[matrixStats]{rowLogSumExps}}
#'
#' @details Computation takes place on the log scale for stability/underflow reasons (to ensure negligible probabilities won't cause computational difficulties); in any case, many functions for calculating multivariate normal densities already output on the log scale.
#'
#' @note Though the function is available for standalone use, note that no checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}.
#'
#' If the normalising constant is required for another reason, e.g. to compute the log-likelihood, it can be calculated by summing the output obtained by calling \code{\link[matrixStats]{rowLogSumExps}} on \code{probs}.
#'
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' Yellott, J. I. Jr. (1977) The relationship between Luce's choice axiom, Thurstone's theory of comparative judgment, and the double exponential distribution, \emph{Journal of Mathematical Psychology}, 15(2): 109-144.
#' @export
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' gumbel_max(probs,
#'            slice = FALSE)
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
#' # Simulate a matrix of Dirichlet weights & the associated vector of N labels
#'   N       <- 400
#'   G       <- 8
#'   sizes   <- seq(from=85, to=15, by=-10)
#'   weights <- matrix(rDirichlet(N * G, alpha=1, nn=sizes), byrow=TRUE, nrow=N, ncol=G)
#'   (zs     <- gumbel_max(probs=log(weights)))
    gumbel_max   <- function(probs, slice = FALSE) {
     if(anyNA(probs))                      stop("Missing values not allowed in 'probs'", call.=FALSE)
     if(length(slice) > 1 ||
        !is.logical(slice))                stop("'slice' must be a single logical indicator", call.=FALSE)
     if(is.vector(probs)) probs <- t(probs)
     if(isTRUE(slice)) {
      fps        <- is.finite(probs)
      probs[fps] <- probs[fps] - log(stats::rexp(sum(fps)))
     } else   {
      probs      <- probs - log(stats::rexp(length(probs)))
     }
        Rfast::rowMaxs(probs)  # i.e. max.col(probs)
    }

  # Alpha
    .sim_alpha_g <- function(alpha, shape, rate, G, N) {
      shape2     <- shape  + G  - 1L
      rate2      <- rate   - log(stats::rbeta(1, alpha + 1, N))
      weights    <- shape2 / (shape2  + N * rate2)
        sum(c(weights,   1 - weights) * c(stats::rgamma(2L, shape=c(shape2 + 1, shape2), rate=rate2)))
    }

    .log_palpha  <- function(alpha, discount, alpha.shape, alpha.rate, N, G) {
      l.prior    <- stats::dgamma(alpha   + discount, shape=alpha.shape, rate=alpha.rate, log=TRUE)
        lgamma(alpha  + 1) - lgamma(alpha + N) + sum(log(alpha + discount    * seq_len(G - 1L))) + l.prior
    }

    .sim_alpha_m <- function(alpha, discount, alpha.shape, alpha.rate, N, G, zeta) {
      propinter  <- c(max( - discount, alpha   - zeta), alpha  + zeta)
      propa      <- stats::runif(1,    propinter[1L], propinter[2L])
      cprob      <- .log_palpha(alpha, discount,  alpha.shape,   alpha.rate, N, G)
      pprob      <- .log_palpha(propa, discount,  alpha.shape,   alpha.rate, N, G)
      currinter  <- c(max( - discount, propa   - zeta), propa  + zeta)
      logpr      <- pprob  - cprob    - log(diff(currinter))   + log(diff(propinter))
      acpt       <- logpr >= 0 ||     - stats::rexp(1L) <= logpr
        return(list(alpha  = ifelse(acpt, propa, alpha), rate  = acpt, l.prob = logpr))
    }

    .log_Oalpha  <- function(x, G, N, nn, shape, rate)  {
      G          <- G * x
        lgamma(G) - lgamma(N + G)         +
        sum(lgamma(nn + x)   - lgamma(x)) +
        stats::dgamma(x, shape=shape, rate=rate, log=TRUE)
    }

    .sim_alpha_o <- function(alpha, zeta = 1, G, N, nn, shape, rate) {
      propa      <- exp(log(alpha)    + stats::rnorm(1, 0, zeta))
      logpr      <- .log_Oalpha(propa, G, N, nn, shape, rate) - .log_Oalpha(alpha, G, N, nn, shape, rate) + log(propa) - log(alpha)
      acpt       <- logpr >= 0 ||     - stats::rexp(1L) <= logpr
        return(list(alpha  = ifelse(acpt, propa, alpha), rate = acpt, l.prob = logpr))
    }

  # Adaptively Tune Zeta
    .tune_zeta   <- function(zeta, time, l.rate, heat = 1L, target = 0.441, lambda = 1L) {
        exp(heat/time^lambda * (exp(min(0, l.rate)) - target)) * zeta
    }

  # Discount
    .log_pdslab  <- function(discount, disc.shape1, disc.shape2, G, unif, nn)            {
      l.prior    <- ifelse(unif, 0, stats::dbeta(discount, shape1=disc.shape1, shape2=disc.shape2, log=TRUE))
        sum(log(discount   * seq_len(G - 1L))) +  sum(lgamma(nn - discount)) - lgamma(1  - discount) * G     + l.prior
    }

    .log_pdspike <- function(discount, disc.shape1, disc.shape2, G, unif, nn, alpha, kappa)     {
      l.prior    <- ifelse(discount == 0, log(kappa),     log1p(- kappa) + ifelse(unif, 0, stats::dbeta(discount, shape1=disc.shape1, shape2=disc.shape2, log=TRUE)))
        sum(log(alpha + discount  * seq_len(G  - 1L)))  + sum(lgamma(nn  - discount) - lgamma(1 - discount)) + l.prior
    }

    .sim_d_slab  <- function(discount, disc.shape1, disc.shape2, G, unif, nn, ...)   {
      propd      <- ifelse(unif, stats::runif(1), stats::rbeta(1, disc.shape1, disc.shape2))
      propd      <- ifelse(propd == 1, propd   - .Machine$double.eps,      propd)
      cprob      <- .log_pdslab(discount,  disc.shape1, disc.shape2, G, unif, nn)
      pprob      <- .log_pdslab(propd,     disc.shape1, disc.shape2, G, unif, nn)
      logpr      <- pprob  - cprob
      acpt       <- logpr >= 0   ||  - stats::rexp(1L) <= logpr
        return(list(disc = ifelse(acpt, propd, discount),   rate = acpt))
    }

    .sim_d_spike <- function(discount, disc.shape1, disc.shape2, G, unif, nn, alpha, kappa)     {
      propd      <- ifelse(alpha  > 0,  ifelse(kappa   != 0 && stats::runif(1) <= kappa, 0L, ifelse(unif,
                           stats::runif(1), stats::rbeta(1, disc.shape1, disc.shape2))), stats::runif(1, max(0, - alpha), 1))
      propd      <- ifelse(propd == 1, propd - .Machine$double.eps, propd)
      if(identical(discount, propd)) {
          return(list(disc = discount, rate = 0L))
      } else {
        cprob    <- .log_pdspike(discount, disc.shape1, disc.shape2, G, unif, nn, alpha, kappa)
        pprob    <- .log_pdspike(propd,    disc.shape1, disc.shape2, G, unif, nn, alpha, kappa)
        logpr    <- pprob  - cprob
        acpt     <- logpr >= 0   ||  - stats::rexp(1L) <= logpr
          return(list(disc = ifelse(acpt, propd, discount), rate = acpt))
      }
    }

# Priors
  # Means
    .sim_mu_p    <- function(P, mu.zero, sig.mu.sqrt) {
        stats::rnorm(P, mu.zero, sig.mu.sqrt)
    }

  # Scores
    #' @importFrom Rfast "matrnorm"
    .sim_eta_p   <- function(N, Q) {
        matrnorm(N, Q)
    }

  # Loadings
    .sim_load_p  <- function(Q, P, sig.l.sqrt) {
        stats::rnorm(P * Q, sd=sig.l.sqrt)
    }

    .sim_load_ps <- function(Q, phi, tau, sigma = 1L) {
        stats::rnorm(Q)/sqrt(phi * tau * sigma)
    }

  # Uniquenesses
    .sim_psi_ipu <- function(P, psi.alpha, psi.beta) {
        .rgamma0(n=P, shape=psi.alpha, rate=psi.beta)
    }

    .sim_psi_ipc <- function(P, psi.alpha, psi.beta) {
        rep(.rgamma0(1, shape=psi.alpha, rate=psi.beta), P)
    }

  # Local Shrinkage
    .sim_phi_p   <- function(Q, P, nu1, nu2) {
        base::matrix(.rgamma0(n=P * Q, shape=nu1, rate=nu2), nrow=P, ncol=Q)
    }

  # Column Shrinkage
    .sim_delta_p <- function(Q = 2L, alpha, beta) {
        stats::rgamma(n=pmax(0L, Q - 1L), shape=alpha, rate=beta)
    }

    .sim_deltaPT <- function(Q = 2L, alpha, beta) {
      rltrgamma(n=pmax(0L, Q - 1L),       shape=alpha, rate=beta)
    }

  # Cluster Shrinkage
    .sim_sigma_p <- function(G = 1L, rho1, rho2)  {
        .rgamma0(n=G, shape=rho1, rate=rho2)
    }

  # Cluster Labels
    .sim_z_p     <- function(N, prob.z) {
        sample(seq_along(prob.z), size=N, replace=TRUE, prob=prob.z)
    }

# Other Functions
  # Uniqueness Hyperparameters
#' Find sensible inverse gamma hyperparameters for variance/uniqueness parameters
#'
#' Takes an inverse-Gamma shape hyperparameter, and an inverse covariance matrix (or estimate thereof), and finds data-driven scale hyperparameters in such a way that Heywood problems are avoided for factor analysis or probabilistic principal components analysis (and mixtures thereof).
#' @param shape A positive shape hyperparameter.
#' @param dat The data matrix for which the inverse covariance matrix is to be estimated. If data are to be centered &/or scaled within \code{\link{mcmc_IMIFA}}, then \code{dat} must also be standardised in the same way.
#' @param type A switch indicating whether a single scale (\code{isotropic}) or variable-specific scales (\code{unconstrained}) are to be derived. Both options are allowed under models in \code{\link{mcmc_IMIFA}} with "constrained" or "unconstrained" uniquenesses, but only a single scale can be specified for models with "isotropic" or "single" uniquenesses.
#' @param beta0 See Note below. Must be a strictly positive numeric scalar. Defaults to \code{3}, but is only invoked when explicitly supplied or when the number of observations \code{N} is not greater than the number of variables \code{P} (or when inverting the sample covariance matrix somehow otherwise fails). In some cases, \code{beta0} should not be so small as to render the estimate of the inverse unstable.
#' @param ... Catches unused arguments. Advanced users can also supply the sample covariance matrix of \code{dat}, if available, through \code{...} via the argument \code{covar}.
#'
#' @details Constraining uniquenesses to be isotropic provides the link between factor analysis and the probabilistic PCA model. When used in conjunction with \code{\link{mcmc_IMIFA}} with "isotropic" or "single" uniquenesses, \code{type} must be \code{isotropic}, but for "unconstrained" or "constrained" uniquenesses, it's possible to specify either a single scale (\code{type="isotropic"}) or variable-specific scales (\code{type="unconstrained"}).
#'
#' Used internally by \code{\link{mcmc_IMIFA}} when its argument \code{psi_beta} is not supplied.
#'
#' @note When \code{N > P}, where \code{N} is the number of observations and \code{P} is the number of variables, the inverse of the sample covariance matrix is used by default.
#'
#' When \code{N <= P}, the inverse either does not exist or the estimate thereof is highly unstable. Thus, an estimate of the form \eqn{\left(\beta_0 + \frac{N}{2}\right)\left(\beta_0\mathcal{I}_p + 0.5\sum_{i=1}^N x_i x_i^\top\right)^{-1}}{(beta0 + N/2) * solve(diag(beta0, P) + 0.5 * crossprod(dat))} is used instead.
#'
#' For unstandardised data, the estimate is instead constructed using a standardised version of the data, and the resulting inverse \emph{correlation} matrix estimate is scaled appropriately by the diagonal entries of the sample covariance matrix of the original data.
#'
#' This estimate can also be used in \code{N > P} cases by explicitly supplying \code{beta0}. It will also be used if inverting the sample covariance matrix fails in \code{N > P} cases.
#'
#' The optional argument \code{beta0} can be supplied to \code{\link{mcmc_IMIFA}} via the control function \code{\link{mixfaControl}}.
#'
#' @return Either a single scale hyperparameter or \code{ncol(dat)} variable-specific scale hyperparameters.
#' @keywords utility
#' @importFrom Rfast "is.symmetric"
#' @export
#'
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link{mixfaControl}}
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' Fruwirth-Schnatter, S. and Lopes, H. F. (2010). Parsimonious Bayesian factor analysis when the number of factors is unknown, \emph{Technical Report}. The University of Chicago Booth School of Business.
#'
#' Fruwirth-Schnatter, S. and Lopes, H. F. (2018). Sparse Bayesian factor analysis when the number of factors is unknown, \emph{to appear}. <\href{https://arxiv.org/abs/1804.04231}{arXiv:1804.04231}>.
#'
#' Tipping, M. E. and Bishop, C. M. (1999). Probabilistic principal component analysis, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 61(3): 611-622.
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' psi_hyper(shape,
#'           dat,
#'           type = c("unconstrained", "isotropic"),
#'           beta0 = 3,
#'           ...)
#' @examples
#' data(olive)
#' olive2  <- olive[,-(1:2)]
#' shape   <- 2.5
#' (scale1 <- psi_hyper(shape=shape, dat=olive2))
#'
#' # Try again with scaled data
#' olive_S <- scale(olive2, center=TRUE, scale=TRUE)
#'
#' # Use the inverse of the sample covariance matrix
#' (scale2 <- psi_hyper(shape=shape, dat=olive_S))
#'
#' # Use the estimated inverse covariance matrix
#' (scale3 <- psi_hyper(shape=shape, dat=olive_S, beta0=3))
#'
#' # In the normalised example, the mean uniquenesses (given by scale/(shape - 1)),
#' # can be interpreted as the prior proportion of the variance that is idiosyncratic
#' (prop1   <- scale1/(shape - 1))
#' (prop2   <- scale2/(shape - 1))
    psi_hyper    <- function(shape, dat, type=c("unconstrained", "isotropic"), beta0 = 3, ...) {
      if(!all(is.character(type)))         stop("'type' must be a character vector of length 1",  call.=FALSE)
      if(any(!is.numeric(shape),
             length(shape) != 1))          stop("'shape' must be a single digit",                 call.=FALSE)
      if(shape   <= 1)                     stop("Scale parameters not defined when 'shape' <= 1", call.=FALSE)
      dat        <- as.matrix(dat)
      N          <- nrow(dat)
      P          <- ncol(dat)
      if(is.null(covar     <- list(...)[["covar"]]))         {
        covar    <- stats::cov(dat)
      } else if(!all(is.symmetric(covar), is.double(covar)) ||
        !is.posi_def(covar, semi=TRUE))    stop("Invalid covariance matrix supplied", call.=FALSE)

      if((noest  <- N  > P && missing(beta0)))               {
        inv.cov  <- try(chol2inv(chol(covar)), silent=TRUE)
        if((noest     <- !inherits(inv.cov, "try-error")))   {
          return(unname((shape - 1)/switch(EXPR=match.arg(type), unconstrained=diag(inv.cov),
                                                isotropic=rep(max(diag(inv.cov)), P))))
        }
      }

      if(!noest)  {
        if(!is.numeric(beta0)   ||
           length(beta0)   != 1 ||
           beta0 <= 0)                     stop("'beta0' must be a strictly positive numeric scalar", call.=FALSE)
        if(beta0  < P &&
           beta0 != floor(beta0))          message("Suggestion:'beta0' should be an integer when less than P\n")
        if(!all((cvar <- diag(covar)) == 1)) {
          dat    <- .scale2(dat, TRUE, TRUE)
        } else cvar   <- 1L
        inv.cov  <- (beta0  + N/2L) * diag(tryCatch(chol2inv(chol(diag(beta0, P) + 0.5 * crossprod(dat))), error=function(e) {
                                           stop("Cannot estimate inverse covariance matrix: try increasing 'beta0'", call.=FALSE) }))
        inv.cov  <- 1/cvar  * inv.cov
          return(unname((shape - 1)/switch(EXPR=match.arg(type), unconstrained=inv.cov,
                                                isotropic=rep(max(inv.cov), P))))
      }
    }

  # Alpha/Discount Shifted Gamma Hyperparameters
#' Moment Matching Parameters of Shifted Gamma Distributions
#'
#' This function takes shape and rate parameters of a Gamma distribution and modifies them to achieve the same expected value and variance when the left extent of the support of the distribution is shifted up or down.
#' @param shape,rate Shape and rate parameters a and b, respectively, of a Gamma(a, b) distribution. Both must be strictly positive.
#' @param shift Modifier, such that the Gamma distribution has support on (\code{shift}, \eqn{\infty}). Can be positive or negative, though typically negative and small.
#' @param param Switch controlling whether the supplied \code{rate} parameter is indeed a rate, or actually a scale parameter. Also governs whether the output is given in terms of rate or scale. Defaults to \code{"rate"}.
#'
#' @note This function is invoked within \code{\link{mcmc_IMIFA}} when \code{discount} is \emph{fixed} at a non-zero value and \code{learn.alpha=TRUE}.
#' @return A named vector of length 2, containing the modified shape and rate parameters, respectively.
#' @keywords utility
#' @export
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' shift_GA(shape,
#'          rate,
#'          shift = 0,
#'          param = c("rate", "scale"))
#' @examples
#' # Shift a Ga(shape=4, rate=2) distribution to the left by 1;
#' # achieving the same expected value of 2 and variance of 1.
#' shift_GA(4, 2, -1)
    shift_GA     <- function(shape, rate, shift = 0, param = c("rate", "scale")) {
      if(length(shape) > 1 ||
        !is.numeric(shape) || shape <= 0)  stop("Argument 'shape' must be a single strictly positive number", call.=FALSE)
      if(length(rate)  > 1 ||
        !is.numeric(rate)  || rate  <= 0)  stop("Argument 'rate' must be a single strictly positive number", call.=FALSE)
      if(length(shift) > 1 ||
        !is.numeric(shift))                stop("Argument 'shift' must be a single number", call.=FALSE)
      if(!all(is.character(param)))        stop("'param' must be a character vector of length 1", call.=FALSE)
      param      <- match.arg(param)

      rate       <- switch(EXPR=param, rate=rate, 1/rate)
      exp        <- shape/rate
      if(shift   >= exp)                   warning("This expected value is not achievable with the supplied 'shift'\n", call.=FALSE, immediate.=TRUE)
      var        <- exp/rate
      exp        <- max(exp - shift, 0L)
      rate       <- exp/var
      shape      <- rate    * exp
        return(c(shape = shape, rate = switch(EXPR=param, rate=rate, 1/rate)))
    }

  # Check Shrinkage Hyperparameters
#' Check the validity of Multiplicative Gamma Process (MGP) hyperparameters
#'
#' Checks the hyperparameters for the multiplicative gamma process (MGP) shrinkage prior in order to ensure that the property of cumulative shrinkage (in expectation) holds, i.e. checks whether growing mass is assigned to small neighbourhoods of zero as the column index increases.
#' @param ad1,ad2 Shape hyperparameters for \eqn{\delta_1}{delta_1} and \eqn{\delta_k \forall k\ge 2}{delta_k for all k >= 2}, respectively.
#' @param Q Number of latent factors. Defaults to \code{3}, which is enough to check if the cumulative shrinkage property holds. Supply \code{Q} if the actual \emph{a priori} expected shrinkage factors are of interest.
#' @param phi.shape,phi.rate The shape and rate hyperparameters for the gamma prior on the local shrinkage parameters. Not necessary for checking if the cumulative shrinkage property holds, but worth supplying \emph{both} if the actual \emph{a priori} expected shrinkage factors are of interest. The default value(s) depends on the value of \code{inverse}, but are chosen in such a way that the local shrinkage has no effect on the expectation unless both are supplied. Cannot be incorporated into the expectation if \code{phi.shape < 1} and \code{isTRUE(inverse)}.
#' @param sigma.shape,sigma.rate The shape and rate hyperparameters for the gamma prior on the cluster shrinkage parameters. Not necessary for checking if the cumulative shrinkage property holds, but worth supplying \emph{both} if the actual \emph{a priori} expected shrinkage factors are of interest. The default value(s) depends on the value of \code{inverse}, but are chosen in such a way that the cluster shrinkage has no effect on the expectation unless both are supplied. Cannot be incorporated into the expectation if \code{sigma.shape < 1} and \code{isTRUE(inverse)}.
#' @param bd1,bd2 Rate hyperparameters for \eqn{\delta_1}{delta_1} and \eqn{\delta_k \forall k\ge 2}{delta_k for all k >= 2}, respectively. Both default to \code{1}.
#' @param truncated A logical value indicating whether the version of the MGP prior based on left-truncated gamma distributions is invoked (see \code{\link{ltrgamma}} and the Zhang et al. reference below). Defaults to \code{FALSE}. Note that, when \code{TRUE}, only the shrinkage parameters for the first loadings column are affected and the conditions needed to pass this check are much less strict. Moreover, more desirable shrinkage properties are easily obtained.
#' @param inverse Logical indicator for whether the cumulative shrinkage property is assessed against the induced Inverse Gamma prior, the default, or in terms of the Gamma prior (which is incorrect). This is always \code{TRUE} when used inside \code{\link{mcmc_IMIFA}}: the \code{FALSE} option exists only for demonstration purposes.
#'
#' @details This is called inside \code{\link{mcmc_IMIFA}} for the \code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"} and \code{"IMIFA"} methods. This function is vectorised with respect to the arguments \code{ad1}, \code{ad2}, \code{phi.shape}, \code{phi.rate}, \code{sigma.shape}, \code{sigma.rate}, \code{bd1} and \code{bd2}.
#'
#' @return A list of length 2 containing the following objects:
#' \itemize{
#'   \item{\strong{expectation} - }{The vector (or list of vectors) of actual expected \emph{a priori} shrinkage factors.}
#'   \item{\strong{valid} - }{A logical (or vector of logicals) indicating whether the cumulative shrinkage property holds (in expectation).}
#' }
#' @export
#' @note It is \emph{recommended} that \code{ad2} be moderately large relative to \code{ad1}, even if \code{valid} can sometimes be \code{TRUE} when this is not the case (e.g. when \code{truncated=TRUE}). Similarly, satisfying this condition is no guarantee that \code{valid} will be \code{TRUE}, unless \code{truncated=TRUE}. Therefore, a warning is returned if \code{ad1 <= ad2}, regardless of the value taken by \code{valid}, when \code{truncated=FALSE} (the default).
#' @keywords control
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link{ltrgamma}}
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' Durante, D. (2017). A note on the multiplicative gamma process, \emph{Statistics & Probability Letters}, 122: 198-204.
#'
#' Bhattacharya, A. and Dunson, D. B. (2011). Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291-306.
#'
#' Zhang, X., Dunson, D. B., and Carin, L. (2011) Tree-structured infinite sparse factor model. In Getoor, L. and Scheffer, T. (Eds.), \emph{Proceedings of the 28th International Conference on Machine Learning}, ICML'11, Madison, WI, USA, pp. 785-792. Omnipress.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' MGP_check(ad1,
#'           ad2,
#'           Q = 3L,
#'           phi.shape = NULL,
#'           phi.rate = NULL,
#'           sigma.shape = NULL,
#'           sigma.rate = NULL,
#'           bd1 = 1,
#'           bd2 = 1,
#'           truncated = FALSE,
#'           inverse = TRUE)
#' @examples
#' # Check if expected shrinkage under the MGP increases with the column index (WRONG approach!).
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, phi.shape=3, inverse=FALSE)$valid   #TRUE
#'
#' # Check if the induced IG prior on the MGP column shrinkage parameters
#' # is stochastically increasing, thereby inducing cumulative shrinkage (CORRECT approach!).
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, phi.shape=3, inverse=TRUE)$valid    #FALSE
#'
#' # Check again with a parameterisation that IS valid and examine the expected shrinkage values
#' (shrink <- MGP_check(ad1=1.5, ad2=2.8, Q=10, phi.shape=2, phi.rate=0.5, inverse=TRUE))
#'
#' # Check previously invalid parameterisation again using truncated version of the MGP prior
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, phi.shape=3, truncated=FALSE)$valid # TRUE
    MGP_check    <- function(ad1, ad2, Q = 3L, phi.shape = NULL, phi.rate = NULL, sigma.shape = NULL, sigma.rate = NULL, bd1 = 1, bd2 = 1, truncated = FALSE, inverse = TRUE) {
      if(length(truncated)     != 1  ||
         !is.logical(truncated))           stop("'truncated' must be a single logical indicator", call.=FALSE)
      if(length(inverse) != 1  ||
         !is.logical(inverse))             stop("'inverse' must be a single logical indicator",   call.=FALSE)
      phi.shape  <- if(is.null(phi.shape))   1L + inverse        else phi.shape
      phi.rate   <- if(is.null(phi.rate))    phi.shape - inverse else phi.rate
      sig.shape  <- if(is.null(sigma.shape)) 1L + inverse        else sigma.shape
      sig.rate   <- if(is.null(sigma.rate))  sig.shape - inverse else sigma.rate

      if(length(Q)       != 1  ||
         !is.numeric(Q)  || floor(Q) != Q) stop("'Q' must be a single integer value", call.=FALSE)
      if(missing(ad1) || missing(ad2))     stop("Column shrinkage shape hyperparameters 'ad1' and 'ad2' must be supplied",         call.=FALSE)
      max.len    <- max(length(ad1), length(ad2), length(Q), length(phi.shape), length(phi.rate), length(sig.shape), length(sig.rate), length(bd1), length(bd2))
      if(!is.element(length(ad1),
                     c(1, max.len)))       stop(paste0("'ad1' must be of length 1 or ", max.len, " for proper recycling"),         call.=FALSE)
      if(!is.element(length(ad2),
                     c(1, max.len)))       stop(paste0("'ad2' must be of length 1 or ", max.len, " for proper recycling"),         call.=FALSE)
      if(!is.element(length(bd1),
                     c(1, max.len)))       stop(paste0("'bd1' must be of length 1 or ", max.len, " for proper recycling"),         call.=FALSE)
      if(!is.element(length(bd2),
                     c(1, max.len)))       stop(paste0("'bd2' must be of length 1 or ", max.len, " for proper recycling"),         call.=FALSE)
      if(!is.element(length(phi.shape),
                     c(1, max.len)))       stop(paste0("'phi.shape' must be of length 1 or ", max.len, " for proper recycling"),   call.=FALSE)
      if(!is.element(length(phi.rate),
                     c(1, max.len)))       stop(paste0("'phi.rate' must be of length 1 or ", max.len, " for proper recycling"),    call.=FALSE)
      if(!is.element(length(sig.shape),
                     c(1, max.len)))       stop(paste0("'sigma.shape' must be of length 1 or ", max.len, " for proper recycling"), call.=FALSE)
      if(!is.element(length(sig.rate),
                     c(1, max.len)))       stop(paste0("'sigma.rate' must be of length 1 or ", max.len, " for proper recycling"),  call.=FALSE)

      if(inverse &
         any(phi.shape   <= 1) -> PX)      warning("Can't incorporate local shrinkage parameters into the expectation when 'phi.shape' is not strictly greater than 1\n",     call.=FALSE, immediate.=TRUE)
      if(inverse &
         any(sig.shape   <= 1) -> SX)      warning("Can't incorporate cluster shrinkage parameters into the expectation when 'sigma.shape' is not strictly greater than 1\n", call.=FALSE, immediate.=TRUE)
      if(any(phi.rate    <= 0))            stop("All local shrinkage rate hyperparameter values must be strictly positive",        call.=FALSE)
      if(any(sig.rate    <= 0))            stop("All cluster shrinkage rate hyperparameter values must be strictly positive",      call.=FALSE)
      if(any(c(ad1, ad2)  < 1))            stop("All column shrinkage shape hyperparameter values must be at least 1",             call.=FALSE)
      if(any(c(bd1, bd2) <= 0))            stop("All column shrinkage rate hyperparameter values must be strictly positive",       call.=FALSE)
      if(any(WX  <- ad1  >= ad2 &
             !truncated))                  warning("'ad2' should be moderately large relative to 'ad1' to encourage loadings column removal\n", call.=FALSE, immediate.=TRUE)

      Qseq       <- seq_len(Q)  - 1L
      ML         <- seq_len(max.len)
      nu1        <- if(PX) 1L   + inverse  else phi.shape
      nu2        <- if(PX) 1L              else phi.rate
      rho1       <- if(SX) 1L   + inverse  else sig.shape
      rho2       <- if(SX) 1L              else sig.rate
      if(isTRUE(inverse)) {
        if(any(XW        <-
           phi.rate       >
           phi.shape      - 1L))           warning("A priori expectation of the induced inverse-gamma local shrinkage hyperprior is not <=1\n", call.=FALSE, immediate.=TRUE)
        ad1      <- ifelse(ad1 == 1, ad1 + .Machine$double.eps, ad1)
        ad2      <- ifelse(ad2 == 1, ad2 + .Machine$double.eps, ad2)
        exp.Q1   <- nu2/(nu1    - 1)     * bd1/(ad1  - 1) * rho2/(rho1 - 1)
        exp.Qk   <- ifelse(truncated, exp_ltrgamma(ad2, bd2, inverse=TRUE),  bd2/(ad2 - 1))
        exp.Q1   <- if(length(exp.Q1)    < length(exp.Qk)) rep(exp.Q1, max.len) else exp.Q1
        exp.Qk   <- if(length(exp.Qk)    < length(exp.Q1)) rep(exp.Qk, max.len) else exp.Qk
        exp.seq  <- lapply(ML, function(i) exp.Q1[i] * exp.Qk[i]^Qseq)
        check    <- vapply(exp.seq,  is.unsorted, logical(1L))
      } else {
        exp.Q1   <- nu1/nu2     * ad1/bd1            * rho1/rho2
        exp.Qk   <- ifelse(truncated, exp_ltrgamma(ad2, bd2, inverse=FALSE), ad2/bd2)
        exp.Q1   <- if(length(exp.Q1)    < length(exp.Qk)) rep(exp.Q1, max.len) else exp.Q1
        exp.Qk   <- if(length(exp.Qk)    < length(exp.Q1)) rep(exp.Qk, max.len) else exp.Qk
        exp.seq  <- lapply(ML, function(i) exp.Q1[i] * exp.Qk[i]^Qseq)
        check    <- !vapply(exp.seq, is.unsorted, logical(1L))
      }
      exp.seq    <- if(length(exp.seq) == 1) exp.seq[[1L]] else exp.seq
      res        <- list(expectation = exp.seq, valid = Q < 2 | check)
      attr(res, "Warning")    <- WX     | (inverse  && XW)
        return(res)
    }

  # Number of 'free' parameters
#' Estimate the Number of Free Parameters in Finite Factor Analytic Mixture Models (PGMM)
#'
#' Estimates the dimension of the 'free' parameters in fully finite factor analytic mixture models, otherwise known as Parsimonious Gaussian Mixture Models (PGMM), typically necessary for the penalty term of various model selection criteria.
#' @param Q The number of latent factors (which can be 0, corresponding to a model with diagonal covariance). This argument is vectorised.
#' @param P The number of variables. Must be a single strictly positive integer.
#' @param G The number of clusters. This defaults to 1. Must be a single strictly positive integer.
#' @param method By default, calculation assumes the \code{UUU} model with unconstrained loadings and unconstrained diagonal uniquesseses (i.e. the factor analysis model). The seven other models detailed in McNicholas and Murphy (2008) are given too (of which currently the first four are accommodated within \code{\link{mcmc_IMIFA}}). The first letter denotes whether loadings are constrained/unconstrained across clusters; the second letter denotes the same for the uniquenesses; the final letter denotes whether uniquenesses are in turn constrained to be isotropic. Finally, the 4 extra 4-letter models from the EPGMM family (McNicholas and Murphy, 2010), are also included.
#' @param equal.pro Logical variable indicating whether or not the mixing mixing proportions are equal across clusters in the model (default = \code{FALSE}).
#'
#' @return A vector of length \code{length(Q)} giving the total number of parameters, including means and mixing proportions, and not only covariance parameters. Set \code{equal.pro} to \code{FALSE} and subtract \code{G * P} from the result to determine the number of covariance parameters only.
#' @keywords utility
#'
#' @note This function is used to calculate the penalty terms for the \code{aic.mcmc} and \code{bic.mcmc} model selection criteria implemented in \code{\link{get_IMIFA_results}} for \emph{finite} factor models (though \code{\link{mcmc_IMIFA}} currently only implements the \code{UUU}, \code{UUC}, \code{UCU}, and \code{UCC} covariance structures). The function is vectorised with respect to the argument \code{Q}.
#'
#' Though the function is available for standalone use, note that no checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}.
#' @export
#' @references McNicholas, P. D. and Murphy, T. B. (2008) Parsimonious Gaussian mixture models, \emph{Statistics and Computing}, 18(3): 285-296.
#'
#' McNicholas, P. D. and Murphy, T. B. (2010) Model-Based clustering of microarray expression data via latent Gaussian mixture models, \emph{Bioinformatics}, 26(21): 2705-2712.
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link{mcmc_IMIFA}}
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' PGMM_dfree(Q,
#'            P,
#'            G = 1L,
#'            method = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC",
#'                       "CCU", "CCC", "CCUU", "UCUU", "CUCU", "UUCU"),
#'            equal.pro = FALSE)
#' @examples
#' (UUU <- PGMM_dfree(Q=0:5, P=50, G=3, method="UUU"))
#' (CCC <- PGMM_dfree(Q=0:5, P=50, G=3, method="CCC", equal.pro=TRUE))
    PGMM_dfree   <- function(Q, P, G = 1L, method = c("UUU", "UUC", "UCU", "UCC", "CUU", "CUC", "CCU", "CCC", "CCUU", "UCUU", "CUCU", "UUCU"), equal.pro = FALSE) {
      if(any(Q    < 0, floor(Q)  != Q))    stop("'Q' must consist of strictly non-negative integers", call.=FALSE)
      if(any(P   <= 0, floor(P)  != P,
        length(P) > 1))                    stop("'P' must be a single strictly positive integer", call.=FALSE)
      if(any(G   <= 0, floor(G)  != G,
        length(G) > 1))                    stop("'G' must be a single strictly positive integer", call.=FALSE)
      if(!all(is.character(method)))       stop("'method' must be a single character string", call.=FALSE)
      if(length(equal.pro)  > 1  ||
         !is.logical(equal.pro))           stop("'equal.pro' must be a single logical indicator", call.=FALSE)
      meth       <- unlist(strsplit(match.arg(method), ""))
      lambda     <- P * Q   - 0.5 * Q * (Q  - 1L)
      lambda     <- switch(EXPR=meth[1L], C = lambda,  U = G    * lambda)
      if(length(meth) < 4)  {
        psi      <- switch(EXPR=meth[2L], C = 1L,      U = G)
        psi      <- switch(EXPR=meth[3L], C = 1L,      U = P)   * psi
      } else      {
        epgmm    <- paste(meth[-1L], sep="", collapse="")
        psi      <- switch(EXPR=epgmm, CUU  = G   + P - 1L, UCU = 1L + G * (P - 1L))
      }
        as.integer(ifelse(equal.pro, 0L, G  - 1L) + G * P  + lambda  + psi)
    }

  # Label Switching
  #' @importFrom matrixStats "colSums2" "rowSums2"
    .lab_switch  <- function(z.new, z.old) {
      tab        <- table(z.new, z.old, dnn=NULL)
      tab.tmp    <- tab[rowSums2(tab) != 0, colSums2(tab) != 0, drop=FALSE]
      nc         <- ncol(tab.tmp)
      nr         <- nrow(tab.tmp)
      ng         <- table(z.new)
      Gs         <- .as_numchar(names(ng))
      if(nc > nr) {
        tmp.mat  <- matrix(0L, nrow=nc - nr, ncol=nc)
        rownames(tmp.mat) <- setdiff(.as_numchar(colnames(tab.tmp)), .as_numchar(rownames(tab.tmp)))[seq_len(nc - nr)]
        tab.tmp  <- rbind(tab.tmp, tmp.mat)
      } else if(nr > nc)   {
        tmp.mat  <- matrix(0L, nrow=nr, ncol=nr - nc)
        colnames(tmp.mat) <- setdiff(.as_numchar(rownames(tab.tmp)), .as_numchar(colnames(tab.tmp)))[seq_len(nr - nc)]
        tab.tmp  <- cbind(tab.tmp, tmp.mat)
      }

      if(length(tab.tmp)  == 1)       {
        z.perm   <- stats::setNames(.as_numchar(colnames(tab.tmp)), .as_numchar(rownames(tab.tmp)))
      } else if(nr == 1)   {
        z.perm   <- stats::setNames(.as_numchar(colnames(tab.tmp)), .as_numchar(colnames(tab.tmp)))
      } else if(nc == 1)   {
        z.perm   <- stats::setNames(.as_numchar(rownames(tab.tmp)), .as_numchar(rownames(tab.tmp)))
      } else {
        z.perm   <- tryCatch(suppressWarnings(.match_classes(tab.tmp, method="exact",  verbose=FALSE)),
          error=function(e) suppressWarnings(.match_classes(tab.tmp,  method="greedy", verbose=FALSE)))
        z.perm   <- stats::setNames(.as_numchar(z.perm), names(z.perm))
      }
      if(length(Gs) > length(z.perm)) {
        z.perm   <- c(z.perm, stats::setNames(setdiff(Gs, z.perm), setdiff(Gs, names(z.perm))))
      }
      z.names    <- .as_numchar(names(z.perm))
      z.perm     <- z.perm[order(z.names)]
      z.sw       <- factor(z.new, labels=z.perm[which(ng > 0)])
        return(list(z = as.integer(levels(z.sw))[z.sw], z.perm = z.perm))
    }

  # Similarity matrix and 'average' clustering
#' Summarise MCMC samples of clustering labels with a similarity matrix and find the 'average' clustering
#'
#' This function takes a Monte Carlo sample of cluster labels, computes an average similarity matrix and returns the clustering with minimum mean squared error to this average. The \code{\link[mcclust]{mcclust}} package \strong{must} be loaded.
#' @param zs A matrix containing samples of clustering labels where the columns correspond to the number of observations (N) and the rows correspond to the number of iterations (M).
#'
#' @details This function takes a Monte Carlo sample of cluster labels, converts them to adjacency matrices, and computes a similarity matrix as an average of the adjacency matrices. The dimension of the similarity matrix is invariant to label switching and the number of clusters in each sample, desirable features when summarising partitions of Bayesian nonparametric models such as IMIFA. As a summary of the posterior clustering, the clustering with minimum mean squared error to this 'average' clustering is reported.
#'
#' A heatmap of \code{z.sim} may provide a useful visualisation, if appropriately ordered. The user is also invited to perform hierarchical clustering using \code{\link[stats]{hclust}} after first converting this similarity matrix to a distance matrix - "complete" linkage is recommended. Alternatively, \code{\link[mclust]{hc}} could be used.
#' @return A list containing three elements:
#' \item{z.avg}{The 'average' clustering, with minimum squared distance to \code{z.sim}.}
#' \item{z.sim}{The N x N similarity matrix, in a sparse format (see \code{\link[slam]{simple_triplet_matrix}}).}
#' \item{MSE.z}{A vector of length M recording the MSEs between each clustering and the 'average' clustering.}
#' @export
#' @keywords utility
#' @references Carmona, C., Nieto-barajas, L. and Canale, A. (2018) Model based approach for household clustering with mixed scale variables. \emph{Advances in Data Analysis and Classification}, 13(2): 559-583.
#' @importFrom slam "as.simple_triplet_matrix"
#'
#' @note The \code{\link[mcclust]{mcclust}} package \strong{must} be loaded.
#'
#' This is liable to take quite some time to run, especially if the number of observations &/or number of iterations is large. Depending on how distinct the clusters are, \code{z.sim} may be stored better in a non-sparse format. This function can optionally be called inside \code{\link{get_IMIFA_results}}.
#' @seealso \code{\link{get_IMIFA_results}}, \code{\link[slam]{simple_triplet_matrix}}, \code{\link[stats]{hclust}}, \code{\link[mclust]{hc}}, \code{\link[mcclust]{comp.psm}}, \code{\link[mcclust]{cltoSim}}
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
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
#' # z.col  <- mat2cols(z.sim, cols=heat.colors(30, rev=TRUE))
#' # z.col[z.sim == 0] <- NA
#' # plot_cols(z.col, na.col=par()$bg); box(lwd=2)
#'
#' # Extract the clustering with minimum squared distance to this
#' # 'average' and evaluate its performance against the true labels
#' # table(zsimil$z.avg, olive$area)
#'
#' # Perform hierarchical clustering on the distance matrix
#' # Hcl    <- hclust(as.dist(1 - z.sim), method="complete")
#' # plot(Hcl)
#' # table(cutree(Hcl, k=3), olive$area)
    Zsimilarity  <- function(zs) {
      has.pkg    <- suppressMessages(requireNamespace("mcclust", quietly=TRUE))
      if(!has.pkg)                         stop("'mcclust' package not installed", call.=FALSE)
      if(!is.matrix(zs))                   stop("'zs' must be a matrix with rows corresponding to the number of observations and columns corresponding to the number of iterations", call.=FALSE)
      if(anyNA(zs))                        stop("Missing values are not allowed in 'zs'", call.=FALSE)
      zsim       <- mcclust::comp.psm(zs)
      mse.z      <- vapply(seq_len(nrow(zs)), function(i, x=mcclust::cltoSim(zs[i,]) - zsim) tryCatch(suppressWarnings(mean(x * x)), error=function(e) Inf), numeric(1L))
      Z.avg      <- zs[which.min(mse.z),]
      attr(Z.avg, "MSE")  <- min(mse.z)
        return(list(z.sim  = as.simple_triplet_matrix(zsim), z.avg = Z.avg, MSE.z = mse.z))
    }

#' Posterior Confusion Matrix
#'
#' For a (\code{N * G}) matrix of posterior cluster membership probabilities, this function creates a (\code{G * G}) posterior confusion matrix, whose hk-th entry gives the average probability that observations with maximum posterior allocation h will be assigned to cluster k.
#' @param z A (\code{N * G}) matrix of posterior cluster membership probabilities whose (ig)-th entry gives the posterior probability that observation i belongs to cluster g. Entries must be valid probabilities in the interval [0,1]; missing values are not allowed.
#'
#' Otherwise, a list of such matrices can be supplied, where each matrix in the list has the same dimensions.
#' @param scale A logical indicator whether the PCM should be rescaled by its row sums. When \code{TRUE} (the default), the benchmark matrix for comparison is the identity matrix of order \code{G}, corresponding to a situation with no uncertainty in the clustering. When \code{FALSE}, the row sums give the number of observations in each cluster.
#' @return A (\code{G * G}) posterior confusion matrix, whose hk-th entry gives the average probability that observations with maximum posterior allocation h will be assigned to cluster k. When \code{scale=TRUE}, the benchmark matrix for comparison is the identity matrix of order \code{G}, corresponding to a situation with no uncertainty in the clustering.
#' @export
#' @keywords utility
#' @importFrom matrixStats "rowSums2"
#' @importFrom Rfast "rowOrder" "rowSort"
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{get_IMIFA_results}
#' @references Ranciati, S., Vinciotti, V. and Wit, E., (2017) Identifying overlapping terrorist cells from the Noordin Top actor-event network, \emph{to appear}. <\href{https://arxiv.org/abs/1710.10319v1}{arXiv:1710.10319v1}>.
#'
#' @examples
#' # data(olive)
#' # sim  <- mcmc_IMIFA(olive, n.iters=1000)
#' # res  <- get_IMIFA_results(sim)
#' # (PCM <- post_conf_mat(res$Clust$post.prob))
#'
#' # par(mar=c(5.1, 4.1, 4.1, 3.1))
#' # PCM  <- replace(PCM, PCM == 0, NA)
#' # plot_cols(mat2cols(PCM, col=heat.colors(30, rev=TRUE), na.col=par()$bg)); box(lwd=2)
#' # heat_legend(PCM, cols=heat.colors(30, rev=TRUE))
#' # par(mar=c(5.1, 4.1, 4.1, 2.1))
   post_conf_mat <- function(z, scale = TRUE) {
      if(inherits(z, "list"))     {
        if(!all(vapply(z, is.matrix,
                        logical(1L))))     stop("Elements of the list 'z' must be matrices", call.=FALSE)
        if(any(vapply(z, anyNA,
                        logical(1L))))     stop("Missing values are not allowed in 'z'", call.=FALSE)
        if(any(vapply(z, function(x)
           any(x  < 0) ||
           any(x  > 1), logical(1L))))     stop("Values in 'z' must be valid probabilities in the interval [0,1]", call.=FALSE)
        nit      <- length(z)
        Ns       <- vapply(z, nrow, numeric(1L))
        Gs       <- vapply(z, ncol, numeric(1L))
        N        <- Ns[1L]
        G        <- Gs[1L]
        if(any(Ns      != N)     ||
           any(Gs      != G)     ||
           N      < G)                     stop("All matrices in the list 'z' must have the same dimensions, with more columns than rows", call.=FALSE)
        tX       <- lapply(z, rowSort,  descending=TRUE)
        rX       <- lapply(z, rowOrder, descending=TRUE)
      } else      {
        if(!is.matrix(z))                  stop("'z' must be a matrix", call.=FALSE)
        if(anyNA(z))                       stop("Missing values are not allowed in 'z'", call.=FALSE)
        if(any(z  < 0) ||
           any(z  > 1))                    stop("Values in 'z' must be valid probabilities in the interval [0,1]", call.=FALSE)
        nit      <- 1L
        N        <- nrow(z)
        G        <- ncol(z)
        if(N      < G)                     stop("'z' must have more rows than columns",  call.=FALSE)
        tX       <- list(rowSort(z,  descending=TRUE))
        rX       <- list(rowOrder(z, descending=TRUE))
      }
      if(length(scale) > 1       ||
         !is.logical(scale))               stop("'scale' must be a single logical indicator", call.=FALSE)
      PCM        <- matrix(0, nrow=G, ncol=G)
      for(n      in seq_len(nit)) {
        for(k    in seq_len(G))   {
          for(i  in seq_len(N))   {
            PCM[rX[[n]][i,1L], rX[[n]][i,k]] <- PCM[rX[[n]][i,1L], rX[[n]][i,k]] + tX[[n]][i,k]
           }
        }
      }
      PCM        <- PCM/nit
        if(scale)   PCM/rowSums2(PCM) else PCM
    }

  # Move 1
    .lab_move1   <- function(nn.ind, pi.prop, nn) {
      sw         <- sample(nn.ind, 2L)
      log.pis    <- log(pi.prop[sw])
      nns        <- nn[sw]
      a.prob     <- (nns[1L] - nns[2L]) * (log.pis[1L] - log.pis[2L])
        return(list(rate1    = a.prob  >= 0 || - stats::rexp(1L) <= a.prob, sw = sw))
    }

  # Move 2
    .lab_move2   <- function(G, Vs, nn)    {
      sw         <- sample.int(G, 1L, prob=c(rep(1, G - 2L), 0.5, 0.5))
      sw         <- if(is.element(sw, c(G, G - 1L))) c(G - 1L, G) else c(sw, sw + 1L)
      nns        <- nn[sw]
      if(nns[1L] == 0)      {
        return(list(rate2   = TRUE,
               sw = sw))
      } else      {
        log.vs   <- log1p(  - Vs[sw])
        a.prob   <- nns[2L] * log.vs[1L]
        a.prob   <- nns[1L] * log.vs[2L]       - ifelse(is.nan(a.prob), 0L, a.prob)
        return(list(rate2   = a.prob >= 0   || - stats::rexp(1L) <= a.prob, sw = sw))
      }
    }

# Positive-(Semi)Definite Checker
#' Check Positive-(Semi)definiteness of a matrix
#'
#' Tests whether all eigenvalues of a symmetric matrix are positive (or strictly non-negative) to check for positive-definiteness and positive-semidefiniteness, respectively. If the supplied matrix doesn't satisfy the test, the nearest matrix which does can optionally be returned.
#' @param x A matrix, assumed to be real and symmetric.
#' @param tol Tolerance for singular values and for absolute eigenvalues - only those with values larger than tol are considered non-zero.
#'
#' (default: tol = \code{max(dim(x))*max(E)*.Machine$double.eps}, where \code{E} is the vector of absolute eigenvalues).
#' @param semi Logical switch to test for positive-semidefiniteness when \code{TRUE} or positive-definiteness when \code{FALSE} (the default).
#' @param make Logical switch to return the nearest matrix which satisfies the test - if the test has been passed, this is of course just \code{x} itself, otherwise the nearest positive-(semi)definite matrix. Note that for reasons due to finite precision arithmetic, finding the nearest positive-definite and nearest positive-semidefinite matrices are effectively equivalent tasks.
#'
#' @return If \code{isTRUE(make)}, a list with two components:
#' \item{\code{check}}{A logical value indicating whether the matrix satisfies the test.}
#' \item{\code{X.new}}{The nearest matrix which satisfies the test (which may just be the input matrix itself.)}
#' Otherwise, only the logical value indicating whether the matrix satisfies the test is returned.
#' @keywords utility
#' @export
#' @usage
#' is.posi_def(x,
#'             tol = NULL,
#'             semi = FALSE,
#'             make = FALSE)
#' @examples
#' x    <- cov(matrix(rnorm(100), nrow=10, ncol=10))
#' is.posi_def(x)                                           #FALSE
#' is.posi_def(x, semi=TRUE)                                #TRUE
#'
#' Xnew <- is.posi_def(x, semi=FALSE, make=TRUE)$X.new
#' identical(x, Xnew)                                       #FALSE
#' identical(x, is.posi_def(x, semi=TRUE, make=TRUE)$X.new) #TRUE
    is.posi_def  <- function(x, tol = NULL, semi = FALSE, make = FALSE)  {
      if(!is.matrix(x)     &&
        nrow(x)  != ncol(x))               stop("argument 'x' is not a square matrix",    call.=FALSE)
      if(anyNA(x))                         stop("argument 'x' contains missing values",   call.=FALSE)
      if(!is.symmetric(x))                 stop("argument 'x' is not a symmetric matrix", call.=FALSE)
      if(!is.numeric(x))                   stop("argument 'x' is not a numeric matrix",   call.=FALSE)
      if(!is.logical(semi) ||
         length(semi) > 1)                 stop("argument 'semi' is not a single logical indicator", call.=FALSE)
      if(!is.logical(make) ||
         length(make) > 1)                 stop("argument 'make' is not a single logical indicator", call.=FALSE)

      d          <- nrow(x)
      eigs       <- eigen(x, symmetric = TRUE)
      eval       <- eigs$values
      abseigs    <- abs(eval)
      tol        <- if(missing(tol)) max(abseigs) * d * .Machine$double.eps else tol
      if(length(tol)  > 1  ||
         !is.numeric(tol))                 stop("argument 'tol' is not a single number", call.=FALSE)
      test       <- replace(eval, abseigs < tol, 0L)
      check      <- all(if(isTRUE(semi)) test >= 0 else test > 0)
      if(isTRUE(make))  {
        evec     <- eigs$vectors
        return(list(check = check, X.new = if(check) x else x + evec %*% tcrossprod(diag(pmax.int(2L * tol - eval, ifelse(isTRUE(semi), 0L, .Machine$double.eps)), d), evec)))
      } else check
    }

  # Ledermann Bound
#' Ledermann Bound
#'
#' Returns the maximum number of latent factors in a factor analysis model for data of dimension \code{P} which actually achieves dimension reduction in terms of the number of covariance parameters. This Ledermann bound is given by the largest integer smaller than or equal to the solution \eqn{k}{k} of \eqn{(M - k)^2 \geq M + k}{(M - k)^2 >= M + k}.
#' @param P Integer number of variables in data set. This argument is vectorised.
#' @param isotropic Logical indicating whether uniquenesses are constrained to be isotropic, in which case the bound is simply \eqn{P-1}{P-1}. Defaults to \code{FALSE}.
#'
#' @return The Ledermann bound, a non-negative integer, or a vector of \code{length(P)} such bounds.
#' @keywords utility
#' @usage
#' Ledermann(P,
#'           isotropic = FALSE)
#' @export
#'
#' @examples
#' Ledermann(c(25, 50, 100))
#'
#' data(olive)
#' Ledermann(ncol(olive[,-c(1,2)]))
    Ledermann    <- function(P, isotropic = FALSE) { # heteroscedastic factors
      if(!is.numeric(P)   ||
         any(P   <= 0, floor(P) != P))      stop("'P' must be a strictly positive integer", call.=FALSE)
      if(length(isotropic) > 1  ||
         !is.logical(isotropic))            stop("'isotropic' must be a single logical indicator", call.=FALSE)
      if(isTRUE(isotropic))      {
          as.integer(P - 1L)
      } else      {
        R        <- P + 0.5 * (1 - sqrt(8L * P  + 1L))
          as.integer(floor(ifelse(1e-10 > abs(R - round(R)), round(R), R)))
      }
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
#'    \deqn{d \times \mathbf{X} \mathbf{R} + 1\hspace*{-3.5pt}1 \underline{t}^\top \approx X^\star}{d X R + 1 t' approximately Xstar}
#'
#'    \code{X.new} is given by:
#'
#'    \deqn{X_{\textrm{new}} = d \times \mathbf{X} \mathbf{R} + 1\hspace*{-3.5pt}1 \underline{t}^\top}{X.new = d X R + 1 t'}
#' }
#'
#' @note \code{X} is padded out with columns containing \code{0} if it has fewer columns than \code{Xstar}.
#' @return A list containing:
#' \item{X.new}{The matrix that is the Procrustes transformed version of \code{X}.}
#' \item{R}{The rotation matrix.}
#' \item{t}{The translation vector (if \code{isTRUE(translate)}).}
#' \item{d}{The scaling factor (is \code{isTRUE(dilate)}).}
#' \item{ss}{The sum of squared differences (if \code{isTRUE(sumsq)}).}
#' @keywords utility
#' @export
#'
#' @references Borg, I. and Groenen, P. J. F. (1997) \emph{Modern Multidimensional Scaling}. Springer-Verlag, New York, pp. 340-342.
#' @usage
#' Procrustes(X,
#'            Xstar,
#'            translate = FALSE,
#'            dilate = FALSE,
#'            sumsq = FALSE)
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
#' # Compare the sum of squared differences to a Procrustean transformation with rotation only
#' mat_ss   <- proc$ss
#' mat_ss2  <- Procrustes(X=mat1, Xstar=mat2, sumsq=TRUE)$ss
    Procrustes   <- function(X, Xstar, translate = FALSE, dilate = FALSE, sumsq = FALSE) {
      if(any(!is.logical(translate),
             length(translate) != 1))      stop("'translate' must be a single logical indicator", call.=FALSE)
      if(any(!is.logical(dilate),
             length(dilate)    != 1))      stop("'dilate' must be a single logical indicator",    call.=FALSE)
      if(any(!is.logical(sumsq),
             length(sumsq)     != 1))      stop("'sumsq' must be a single logical indicator",     call.=FALSE)
      if((N <- nrow(X)) != nrow(Xstar))    stop("X and Xstar do not have the same number of rows",       call.=FALSE)
      if((P <- ncol(X)) != (P2 <- ncol(Xstar)))  {
        if(P < P2)       {                 warning("X padded out to same number of columns as Xstar\n",  call.=FALSE, immediate.=TRUE)
          X <- cbind(X, matrix(0L, nrow=N, ncol=P2 - P))
        } else                             stop("X cannot have more columns than Xstar",          call.=FALSE)
      }
      if(P2 == 0)                          stop("Xstar must contain at least one column",         call.=FALSE)
      if(anyNA(Xstar)   || anyNA(X))       stop("X and Xstar are not allowed to contain missing values", call.=FALSE)
      J          <- if(translate) diag(N) - 1/N                                      else diag(N)
      C          <- if(translate) crossprod(Xstar, J) %*% X                          else crossprod(Xstar, X)
      if(!all(P  == 1,
              P2 == 1))  {
        svdX     <- svd(C)
        R        <- tcrossprod(svdX$v, svdX$u)
      } else R   <- 1L
      d          <- if(dilate)    sum(C * R)/sum(crossprod(J, X) * X)                else 1L
      tt         <- if(translate) crossprod(Xstar - d * X %*% R, matrix(1L, N, 1))/N else 0L
      X.new      <- d * X %*% R + if(translate) matrix(tt, N, P2, byrow = TRUE)      else tt
        return(c(list(X.new = X.new), list(R = R), if(translate) list(t = tt),
                 if(dilate) list(d = d), if(sumsq) list(ss = sum((X[,seq_len(P2), drop=FALSE] - X.new)^2L))))
    }

  # Length Checker
    .len_check   <- function(obj0g, switch0g, method, P, range.G, P.dim = TRUE) {
      V          <- ifelse(P.dim, P, 1L)
      rGseq      <- seq_along(range.G)
      obj.name   <- deparse(substitute(obj0g))
      obj.name   <- ifelse(grepl("$", obj.name, fixed=TRUE), vapply(strsplit(obj.name, "$", fixed=TRUE), "[[", character(1L), 2L), obj.name)
      sw.name    <- deparse(substitute(switch0g))
      if(!inherits(obj0g,
                   "list"))       obj0g <- list(obj0g)
      if(length(obj0g) != length(range.G))    {
        if(!P.dim)             {
          obj0g  <- replicate(length(range.G), obj0g)
        } else                             stop(paste0(obj.name, " must be a list of length ", length(range.G)), call.=FALSE)
      }
      len        <- lengths(obj0g)

      if(is.element(method, c("FA", "IFA"))  ||
         all(range.G == 1))    {
        if(any(!is.element(len, c(1, V)))) stop(paste0(obj.name, " must be a scalar", ifelse(P.dim, paste0(" or a vector of length P=", V), ""), " for a 1-group model"), call.=FALSE)
      } else {
        if(any(is.element(len,
           c(1, range.G, V)))) {
          if(length(unique(len)) > 1    &&
             !switch0g)                    stop(paste0(sw.name, " must be TRUE if the dimension/length of ", obj.name, " varies"), call.=FALSE)
          if(all(len == range.G)) obj0g <- if(switch0g) lapply(rGseq, function(g) matrix(obj0g[[g]], nrow=1L))  else stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"), call.=FALSE)
          if(all(len == V))       obj0g <- if(V == 1)   lapply(rGseq, function(g) rep(obj0g[[g]], range.G[g]))  else lapply(rGseq, function(g) matrix(obj0g[[g]], nrow=V, ncol=range.G[g]))
        } else if(!all(vapply(rGseq, function(g, dimG=as.integer(dim(obj0g[[g]]))) is.matrix(obj0g[[g]]) && any(identical(dimG, as.integer(c(1, range.G[g]))), identical(dimG, as.integer(c(V, range.G[g])))), logical(1L)))) {
                                           stop(paste0(ifelse(length(rGseq) > 1, "Each element of ", ""), obj.name, " must be of length 1", ifelse(P.dim, paste0(", P=", V, ifelse(length(rGseq) == 1, paste0(", or G=",  range.G),  ", or its corresponding range.G"),
                                                ", or a matrix with P rows and ", ifelse(length(rGseq) == 1, "G columns", "its corresponding range.G number of columns")),  ifelse(switch0g,  ifelse(length(rGseq) == 1, paste0( " or G=",  range.G),  " or its corresponding range.G"), ""))), call.=FALSE)
        } else if(all(vapply(obj0g, is.matrix, logical(1L)), !switch0g) && any(vapply(rGseq, function(g) any(dim(obj0g[[g]]) == range.G[g]), logical(1L)))) {
                                           stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"), call.=FALSE)
        }
      }
      if(all(length(unique(unlist(obj0g))) > 1,
             !switch0g, !P.dim))           stop(paste0(obj.name, " must be a scalar if ", sw.name, " is FALSE"), call.=FALSE)
        obj0g
    }

  # Moments of Pitman-Yor / Dirichlet Processes
#' 1st & 2nd Moments of the Pitman-Yor / Dirichlet Processes
#'
#' Calculate the \emph{a priori} expected number of clusters (\code{G_expected}) or the variance of the number of clusters (\code{G_variance}) under a PYP or DP prior for a sample of size \code{N} at given values of the concentration parameter \code{alpha} and optionally also the Pitman-Yor \code{discount} parameter. Useful for soliciting sensible priors (or fixed values) for \code{alpha} or \code{discount} under the \code{"IMFA"} and \code{"IMIFA"} methods for \code{\link{mcmc_IMIFA}}. Additionally, for a given sample size \code{N} and given expected number of clusters \code{EG}, \code{G_calibrate} elicits a value for the concentration parameter \code{alpha} \strong{or} the \code{discount} parameter.
#' @param N The sample size.
#' @param alpha The concentration parameter. Must be specified (though not for \code{G_calibrate}) and must be strictly greater than \code{-discount}. The case \code{alpha=0} is accommodated. When \code{discount} is negative \code{alpha} must be a positive integer multiple of \code{abs(discount)}. See \strong{Details} for behaviour for \code{G_calibrate}.
#' @param discount The discount parameter for the Pitman-Yor process. Must be less than 1, but typically lies in the interval [0, 1). Defaults to 0 (i.e. the Dirichlet process). When \code{discount} is negative \code{alpha} must be a positive integer multiple of \code{abs(discount)}. See \strong{Details} for behaviour for \code{G_calibrate}.
#' @param MPFR Logical indicating whether the high-precision libraries \code{\link[Rmpfr]{Rmpfr}} and \code{gmp} are invoked, at the expense of run-time. Defaults to \code{TRUE} and \strong{must} be \code{TRUE} for \code{G_expected} when \code{alpha=0} or \code{G_variance} when \code{discount} is non-zero. For \code{G_calibrate}, it is \emph{strongly recommended} to use \code{MPFR=TRUE} when \code{discount} is non-zero and strictly necessary when \code{alpha=0} is supplied. See \strong{\code{Note}}.
#' @param EG The prior expected number of clusters. Must exceed \code{1} and be less than \code{N}.
#' @param ... Additional arguments passed to \code{\link[stats]{uniroot}}, e.g. \code{maxiter}.
#'
#' @details All arguments are vectorised. Users can also consult \code{\link{G_priorDensity}} in order to solicit sensible priors.
#'
#' For \code{G_calibrate}, \strong{only one} of \code{alpha} or \code{discount} can be supplied, and the function elicits a value for the opposing parameter which achieves the desired expected number of clusters \code{EG} for the given sample size \code{N}. By default, a value for \code{alpha} subject to \code{discount=0} (i.e. the Dirichlet process) is elicited. Note that \code{alpha} may not be a positive integer multiple of \code{discount} as it should be if \code{discount} is negative. See \strong{Examples} below.
#' @return The expected number of clusters under the specified prior conditions (\code{G_expected}), or the variance of the number of clusters (\code{G_variance}), or the concentration parameter \code{alpha} \strong{or} \code{discount} parameter achieving a particular expected number of clusters (\code{G_calibrate}).
#' @keywords utility
#' @export
#' @name G_moments
#' @rdname G_moments
#'
#' @note \code{G_variance} requires use of the \code{\link[Rmpfr]{Rmpfr}} and \code{gmp} libraries for non-zero \code{discount} values. \code{G_expected} requires these libraries only for the \code{alpha=0} case. These libraries are \emph{strongly recommended} (but they are not required) for \code{G_calbirate} when \code{discount} is non-zero, but they are required when \code{alpha=0} is supplied. Despite the high precision arithmetic used, the functions can still be unstable for large \code{N} and/or extreme values of \code{alpha} and/or \code{discount}. See the argument \code{MPFR}.
#'
#' @seealso \code{\link{G_priorDensity}}, \code{\link[Rmpfr]{Rmpfr}}, \code{\link[stats]{uniroot}}
#' @references De Blasi, P., Favaro, S., Lijoi, A., Mena, R. H., Prunster, I., and Ruggiero, M. (2015) Are Gibbs-type priors the most natural generalization of the Dirichlet process?, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 37(2): 212-229.
#'
#' Yamato, H. and Shibuya, M. (2000) Moments of some statistics of Pitman sampling formula, \emph{Bulletin of Informatics and Cybernetics}, 32(1): 1-10.
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' G_expected(N,
#'            alpha,
#'            discount = 0,
#'            MPFR = TRUE)
#' @examples
#' G_expected(N=50, alpha=19.23356, MPFR=FALSE)
#' G_variance(N=50, alpha=19.23356, MPFR=FALSE)
#'
#' G_expected(N=50, alpha=c(19.23356, 12.21619, 1),
#'            discount=c(0, 0.25, 0.7300045), MPFR=FALSE)
#' # require("Rmpfr")
#' # G_variance(N=50, alpha=c(19.23356, 12.21619, 1),
#' #            discount=c(0, 0.25, 0.7300045), MPFR=c(FALSE, TRUE, TRUE))
#'
#' # Examine the growth rate of the DP
#' DP   <- sapply(c(1, 5, 10), function(i) G_expected(1:200, alpha=i, MPFR=FALSE))
#' matplot(DP, type="l", xlab="N", ylab="G")
#'
#' # Examine the growth rate of the PYP
#' # PY <- sapply(c(0.25, 0.5, 0.75), function(i) G_expected(1:200, alpha=1, discount=i))
#' # matplot(PY, type="l", xlab="N", ylab="G")
#'
#' # Other special cases of the PYP are also facilitated
#' # G_expected(N=50, alpha=c(27.1401, 0), discount=c(-27.1401/100, 0.8054448))
#' # G_variance(N=50, alpha=c(27.1401, 0), discount=c(-27.1401/100, 0.8054448))
#'
#' # Elicit values for alpha under a DP prior
#' G_calibrate(N=50, EG=25)
#'
#' # Elicit values for alpha under a PYP prior
#' # require("Rmpfr")
#' # G_calibrate(N=50, EG=25, discount=c(-27.1401/100, 0.25, 0.7300045))
#'
#' # Elicit values for discount under a PYP prior
#' # G_calibrate(N=50, EG=25, alpha=c(12.21619, 1, 0), maxiter=2000)
    G_expected   <- Vectorize(function(N, alpha, discount = 0, MPFR = TRUE) {
      if(!all(is.numeric(N), is.numeric(discount),
         is.numeric(alpha)))               stop("All inputs must be numeric",     call.=FALSE)
      if(discount >= 1)                    stop("'discount' must be less than 1", call.=FALSE)
      if(discount > 0    &&
         alpha   <= - discount)            stop("'alpha' must be strictly greater than -discount",     call.=FALSE)
      if(discount < 0    &&
        (alpha   <= 0    ||
        !.IntMult(alpha,    discount)))    stop("'alpha' must be a positive integer multiple of 'abs(discount)' when 'discount' is negative", call.=FALSE)
      if(alpha   == 0    && discount <= 0) stop("'discount' must be strictly positive when 'alpha'=0", call.=FALSE)
      if(alpha   == 0    && !isTRUE(MPFR)) stop("'MPFR' must be TRUE when 'alpha' == 0", call.=FALSE)
      igmp       <- isNamespaceLoaded("Rmpfr")
      if(mpfrind <- (isTRUE(MPFR)    &&
                     suppressMessages(requireNamespace("Rmpfr", quietly=TRUE)) &&
                     .version_above("gmp", "0.5-4")))       {
        if(isFALSE(igmp)) {
          on.exit(.detach_pkg("Rmpfr"))
          on.exit(.detach_pkg("gmp"), add=TRUE)
        }
        alpha    <- Rmpfr::mpfr(alpha, precBits=256)
      }   else if(isTRUE(MPFR))            warning("'Rmpfr' package not installed\n", call.=FALSE, immediate.=TRUE)

      if(alpha    == 0)   {
        if(mpfrind)       {
          return(gmp::asNumeric(Rmpfr::pochMpfr(discount + 1, N - 1L)/gamma(Rmpfr::mpfr(N, precBits=256))))
        } else                             stop("'Rmpfr' must be installed when 'alpha'=0", call.=FALSE)
      }
      if(discount == 0)   {
        exp      <- alpha * (digamma(alpha + N) - digamma(alpha))
       #exp      <- sum(alpha/(alpha  +  0L:(N  - 1L)))
        if(mpfrind)       {
          gmp::asNumeric(exp)
        } else    {
          exp
        }
      } else      {
        adx      <- alpha/discount
        if(mpfrind)       {
          gmp::asNumeric(adx  * Rmpfr::pochMpfr(alpha + discount, N)/Rmpfr::pochMpfr(alpha, N) - adx)
        } else    {
          adx * (prod(discount/(alpha + 0L:(N - 1L))  + 1L) - 1L)
        }
      }
    })

#' @keywords utility
#' @rdname G_moments
#' @usage
#' G_variance(N,
#'            alpha,
#'            discount = 0,
#'            MPFR = TRUE)
#' @export
    G_variance   <- Vectorize(function(N, alpha, discount = 0, MPFR = TRUE) {
      if(!all(is.numeric(N), is.numeric(discount),
         is.numeric(alpha)))               stop("All inputs must be numeric",     call.=FALSE)
      if(discount >= 1)                    stop("'discount' must be less than 1", call.=FALSE)
      if(discount > 0    &&
         alpha   <= - discount)            stop("'alpha' must be strictly greater than -discount",     call.=FALSE)
      if(discount < 0    &&
        (alpha   <= 0    ||
        !.IntMult(alpha,    discount)))    stop("'alpha' must be a positive integer multiple of 'abs(discount)' when 'discount' is negative", call.=FALSE)
      if(alpha   == 0    && discount <= 0) stop("'discount' must be strictly positive when 'alpha'=0", call.=FALSE)
      if(discount != 0   && !isTRUE(MPFR)) stop("'MPFR' must be TRUE when 'discount' is non-zero",     call.=FALSE)
      igmp       <- isNamespaceLoaded("Rmpfr")
      if(mpfrind <- (isTRUE(MPFR)    &&
                     suppressMessages(requireNamespace("Rmpfr", quietly=TRUE)) &&
                     .version_above("gmp", "0.5-4"))) {
        if(isFALSE(igmp)) {
          on.exit(.detach_pkg("Rmpfr"))
          on.exit(.detach_pkg("gmp"), add=TRUE)
        }
        alpha    <- Rmpfr::mpfr(alpha,   precBits=256)
      }   else    {
        if(discount      != 0)        {    stop("'Rmpfr' package not installed",      call.=FALSE)
        } else if(isTRUE(MPFR))            warning("'Rmpfr' package not installed\n", call.=FALSE, immediate.=TRUE)
      }

      alpha2     <- alpha * alpha
      if(discount == 0)   {
        var      <- alpha * (digamma(alpha  + N) - digamma(alpha))
        if(mpfrind)       {
          alpha  <- gmp::asNumeric(alpha)
          gmp::asNumeric(var + alpha2 * (trigamma(alpha + N) - trigamma(alpha)))
        } else {
          var  + alpha2   * (trigamma(alpha + N) - trigamma(alpha))
        }
      } else if(alpha    == 0) {
        poch.1   <- gamma(Rmpfr::mpfr(N, precBits=256))
        subterm  <- Rmpfr::pochMpfr(discount + 1, N - 1L)/poch.1
          gmp::asNumeric(Rmpfr::pochMpfr(2 * discount, N)/(discount * poch.1) - subterm  * (1L + subterm))
      } else   {
        sum.ad   <- alpha + discount
        poch.a   <- Rmpfr::pochMpfr(alpha, N)
        poch.ad  <- Rmpfr::pochMpfr(sum.ad, N)
        subterm  <- alpha/discount * poch.ad/poch.a
          gmp::asNumeric((alpha * sum.ad)/(discount * discount) * Rmpfr::pochMpfr(sum.ad + discount, N)/poch.a - subterm * (1L + subterm))
      }
    })

#' @keywords utility
#' @rdname G_moments
#' @usage
#' G_calibrate(N,
#'             EG,
#'             alpha = NULL,
#'             discount = 0,
#'             MPFR = TRUE,
#'             ...)
#' @export
    G_calibrate  <- Vectorize(function(N, EG, alpha = NULL, discount = 0, MPFR = TRUE, ...) {
      if(!all(is.numeric(N), is.numeric(discount), is.null(alpha) || is.numeric(alpha),
         is.numeric(EG)))                  stop("All inputs must be numeric",     call.=FALSE)
      if(discount >= 1)                    stop("'discount' must be less than 1", call.=FALSE)
      if(EG       <= 1)                    stop("'EG' must be greater than 1",    call.=FALSE)
      if(EG       >= N)                    stop("'EG' must be less than 'N'",     call.=FALSE)
      igmp       <- isNamespaceLoaded("Rmpfr")
      if(mpfrind <- (isTRUE(MPFR)    &&
                     suppressMessages(requireNamespace("Rmpfr", quietly=TRUE)) &&
                     .version_above("gmp", "0.5-4"))) {
        if(isFALSE(igmp))   {
          on.exit(.detach_pkg("Rmpfr"))
          on.exit(.detach_pkg("gmp"), add=TRUE)
        }
      } else if(isTRUE(MPFR))              warning("'Rmpfr' package not installed\n", call.=FALSE, immediate.=TRUE)
      if(isTRUE(mpfrind)   &&
        (!is.null(alpha)   || discount != 0))          {
        if(!is.null(alpha) && alpha    == 0)           {
          p1     <- gamma(Rmpfr::mpfr(N, precBits=256))
          RFA    <- function(N, discount, p1)          {
            gmp::asNumeric(Rmpfr::pochMpfr(discount + 1, N - 1L)/p1)
          }
        } else    {
          RFA    <- function(N, alpha, discount)       {
            x    <- gmp::asNumeric((Rmpfr::pochMpfr(alpha  + discount, N)/Rmpfr::pochMpfr(alpha, N)  - 1)    * alpha/discount)
            ifelse(x == 0, N, ifelse(is.infinite(x)   && x < 0, 1, x))
          }
        }
      } else RFA <- function(N, alpha, discount) alpha/discount * (prod(1 + discount/(alpha  + 0L:(N - 1L))) - 1)

      if(is.null(alpha))    {
        if(discount  == 0)  {
          inter  <- c(.Machine$double.eps,    .Machine$double.xmax)
          X      <- try(suppressWarnings(stats::uniroot(function(x) sum(x/(x + 0L:(N - 1L))) - EG, interval=inter, ...)), silent=TRUE)
        } else    {
          inter  <- if(isTRUE(mpfrind)) c(-discount + .Machine$double.eps, .Machine$double.xmax) else c(-discount + 0.000001,  100000)
          X      <- try(suppressWarnings(stats::uniroot(function(x) RFA(N, x, discount)      - EG, interval=inter, ...)), silent=TRUE)
        }
        if(inherits(X, "try-error")) {     warning(paste0("uniroot failed to elicit a discount value", ifelse(isFALSE(MPFR), ": consider setting MPFR=TRUE\n","\n")), call.=FALSE, immediate.=TRUE)
          Y      <- stats::setNames(NA,         "alpha")
        } else Y <- stats::setNames(X$root,     "alpha")
      }   else if(missing(discount) ||
                  discount == 0)     {
        if(alpha == 0)               {
          if(!isTRUE(MPFR))                stop("'MPFR' must be TRUE when 'alpha' == 0",          call.=FALSE)
          inter  <- c(.Machine$double.eps,   1 - .Machine$double.eps)
          X      <- try(suppressWarnings(stats::uniroot(function(x) RFA(N, x,    p1)         - EG, interval=inter, ...)), silent=TRUE)
        } else    {
          inter  <- c(-.Machine$double.xmax, 1 - .Machine$double.eps)
          X      <- try(suppressWarnings(stats::uniroot(function(x) RFA(N, alpha, x)         - EG, interval=inter, ...)), silent=TRUE)
        }
        if(inherits(X, "try-error")) {     warning(paste0("uniroot failed to elicit a discount value", ifelse(isFALSE(MPFR), ": consider setting MPFR=TRUE\n","\n")), call.=FALSE, immediate.=TRUE)
          Y      <- stats::setNames(NA,      "discount")
        } else Y <- stats::setNames(X$root,  "discount")
      }   else                             stop("'alpha' and 'discount' cannot both be supplied", call.=FALSE)
      dots       <- list(...)
      maxiter    <- ifelse(length(dots) > 0 && any(names(dots) %in% "maxiter"), dots$maxiter, 1000)
      if(!inherits(X, "try-error")  &&
         X$iter  == maxiter)               warning(paste0("uniroot failed to converge in ", maxiter, " iterations\n"), call.=FALSE)
        return(Y)
    },  vectorize.args = c("N", "EG", "alpha", "discount", "MPFR"))

  # Print functions
#' @method print IMIFA
#' @rdname mcmc_IMIFA
#' @usage
#' \method{print}{IMIFA}(x,
#'       ...)
#' @include MainFunction.R
#' @export
    print.IMIFA  <- function(x, ...) {
      meth       <- attr(x, "Method")
      name       <- attr(x, "Name")
      fac        <- attr(x, "Factors")
      grp        <- attr(x, "Clusters")
      Qmsg       <- Gmsg <- msg   <- NULL
      for(i in seq_along(fac[-length(fac)])) {
        Qmsg     <- c(Qmsg, (paste0(fac[i], ifelse(i + 1 < length(fac), ", ", " "))))
      }
      for(i in seq_along(grp[-length(grp)])) {
        Gmsg     <- c(Gmsg, (paste0(grp[i], ifelse(i + 1 < length(grp), ", ", " "))))
      }

      Qmsg       <- if(length(fac) > 1) paste(c(Qmsg, paste0("and ", fac[length(fac)])), sep="", collapse="") else fac
      Gmsg       <- if(length(grp) > 1) paste(c(Gmsg, paste0("and ", grp[length(grp)])), sep="", collapse="") else grp
      Qmsg       <- paste0(" with ", Qmsg, " factor", ifelse(length(fac) > 1 || fac > 1, "s", ""))
      Gmsg       <- paste0(" with ", Gmsg, " group",  ifelse(length(grp) > 1 || grp > 1, "s", ""))
      if(is.element(meth, c("FA", "OMFA", "IMFA"))) {
        msg      <- Qmsg
      } else {
        msg      <- switch(EXPR=meth, MFA=paste0(Gmsg, " and", Qmsg), MIFA=Gmsg)
      }
      cat(paste0(meth, " simulations for '", name, "' data set", msg, " to be passed to get_IMIFA_results(...)\n"))
        invisible()
    }

#' @method summary IMIFA
#' @rdname mcmc_IMIFA
#' @usage
#' \method{summary}{IMIFA}(object,
#'         ...)
#' @include MainFunction.R
#' @export
    summary.IMIFA        <- function(object, ...) {
      meth       <- attr(object, "Method")
      name       <- attr(object, "Name")
      call       <- attr(object, "Call")
      fac        <- attr(object, "Factors")
      grp        <- attr(object, "Clusters")
      Qmsg       <- Gmsg <- msg   <- NULL
      for(i in seq_along(fac[-length(fac)])) {
        Qmsg     <- c(Qmsg, (paste0(fac[i], ifelse(i + 1 < length(fac), ", ", " "))))
      }
      for(i in seq_along(grp[-length(grp)])) {
        Gmsg     <- c(Gmsg, (paste0(grp[i], ifelse(i + 1 < length(grp), ", ", " "))))
      }

      Qmsg       <- if(length(fac) > 1) paste(c(Qmsg, paste0("and ", fac[length(fac)])), sep="", collapse="") else fac
      Gmsg       <- if(length(grp) > 1) paste(c(Gmsg, paste0("and ", grp[length(grp)])), sep="", collapse="") else grp
      Qmsg       <- paste0(" with ", Qmsg, " factor", ifelse(length(fac) > 1 || fac > 1, "s", ""))
      Gmsg       <- paste0(" with ", Gmsg, " group",  ifelse(length(grp) > 1 || grp > 1, "s", ""))
      if(is.element(meth, c("FA", "OMFA", "IMFA"))) {
        msg      <- Qmsg
      } else {
        msg      <- switch(EXPR=meth, MFA=paste0(Gmsg, " and", Qmsg), MIFA=Gmsg)
      }
      summ       <- list(call = call, details = paste0(meth, " simulations for '", name, "' data set", msg, " to be passed to get_IMIFA_results(...)\n"))
      class(summ)        <- "summary_IMIFA"
        summ
    }

#' @method print summary_IMIFA
#' @export
  print.summary_IMIFA    <- function(x, ...) {
    cat("Call:\t")
    print(x$call)
    cat(paste0("\n", x$details))
  }

#' @method print Results_IMIFA
#' @rdname get_IMIFA_results
#' @usage
#' \method{print}{Results_IMIFA}(x,
#'       ...)
#' @export
    print.Results_IMIFA  <- function(x, ...) {
      method     <- attr(x, "Method")
      adapt      <- attr(x, "Adapt") || !is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA"))
      choice     <- ifelse(attr(x, "Choice"), paste0(" (", switch(EXPR=method, MFA="both ", ""), "chosen by ", toupper(attr(x, "Criterion")), ")"), "")
      G          <- x$GQ.results$G
      Q          <- x$GQ.results$Q
      switch(EXPR=method, FA=, IFA={
        msg      <- paste0("The chosen ", method, " model has ", Q, " factor", ifelse(Q == 1, "", "s"), choice, ifelse(adapt, "", " (no adaptation took place)"))
      }, MFA=, OMFA=, IMFA= {
        msg      <- paste0("The chosen ", method, " model has ", G, " group",  ifelse(G == 1, " with ", "s, each with "), unique(Q), " factor", ifelse(unique(Q) == 1, "", "s"), choice)
      }, {
        Q.msg    <- NULL
        for(i in seq_along(Q[-length(Q)])) {
          Q.msg  <- c(Q.msg, (paste0(Q[i], ifelse(i + 1 < length(Q), ", ", " "))))
        }
        Q.msg    <- if(!adapt) paste0(unique(Q)) else { if(length(Q) > 1) paste(c(Q.msg, paste0("and ", Q[length(Q)])), sep="", collapse="") else Q }
        msg      <- paste0("The chosen ", method, " model has ", G, " group", ifelse(G  == 1, "", "s"), choice, ifelse(G == 1, " with ", paste0(ifelse(adapt, "", " each"), " with ")), Q.msg, " factor", ifelse(all(Q == 1), "", paste0("s", ifelse(G == 1, "", " respectively"))), ifelse(adapt, "", " (no adaptation took place)"), sep="")
      })
      cat(paste0(msg, ":\nthis Results_IMIFA object can be passed to plot(...)\n"))
      if(!isTRUE(attr(x,  "Nowarn.G"))) {  cat("\n");  message(attr(x, "Nowarn.G"))
      }
      if(!isTRUE(attr(x,  "Nowarn.Q"))) {
        if(isTRUE(attr(x, "Nowarn.G"))) {  cat("\n")}; message(attr(x, "Nowarn.Q"))
      }
        invisible()
    }

#' @method summary Results_IMIFA
#' @rdname get_IMIFA_results
#' @usage
#' \method{summary}{Results_IMIFA}(object,
#'         MAP = TRUE,
#'         ...)
#' @export
    summary.Results_IMIFA <- function(object, MAP = TRUE, ...) {
      if(any(!is.logical(MAP),
             length(MAP) != 1))            stop("'MAP' must be a single logical indicator", call.=FALSE)
      criterion  <- unlist(strsplit(toupper(attr(object$GQ.results, "Criterion")), "[.]"))
      criterion  <- ifelse(length(criterion) > 1, ifelse(criterion[1L] != "LOG", paste0(criterion[1L], ".", tolower(criterion[2L])), "LogIntegratedLikelihood"), criterion)
      crit.mat   <- object$GQ.results[[paste0(criterion, "s")]]
      call       <- attr(object, "Call")
      msg        <- NULL
      if(any(dim(crit.mat) > 1)) {
        msg      <- paste0(", and ", ifelse(substr(criterion, 1L, 1L) == "A", "an ", "a "),  criterion, " of ", round(max(crit.mat), 2L), "\n")
      }
      summ       <- list(call = call, details = paste0(utils::capture.output(print(object)), msg),
                         MAP=object$Clust$MAP, showMAP = MAP)
      class(summ)        <- "summary_ResultsIMIFA"
        summ
    }

#' @method print summary_ResultsIMIFA
#' @export
    print.summary_ResultsIMIFA  <- function(x, ...) {
      cat("Call:\t")
      print(x$call)
      cat(paste0("\n", x$details))
      if(isTRUE(x$showMAP))      {
        cat("\n\nClustering table :")
        print(table(x$MAP), row.names = FALSE)
      }
    }

  # Control functions
#' Control settings for the IMIFA family of factor analytic mixtures
#'
#' Supplies a list of arguments for use in \code{\link{mcmc_IMIFA}} pertaining to \emph{ALL} methods in the \code{IMIFA} family: e.g. MCMC settings, cluster initialisation, generic hyperparameters for factor-analytic mixtures, etc.
#' @param n.iters The number of iterations to run the sampler for. Defaults to 25000.
#' @param burnin The number of burn-in iterations for the sampler. Defaults to \code{n.iters/5}. Note that chains can also be burned in later, using \code{\link{get_IMIFA_results}}.
#' @param thinning The thinning interval used in the simulation. Defaults to 2. No thinning corresponds to 1. Note that chains can also be thinned later, using \code{\link{get_IMIFA_results}}.
#' @param centering A logical value indicating whether mean centering should be applied to the data, defaulting to \code{TRUE}.
#' @param scaling The scaling to be applied - one of \code{"unit"}, \code{"none"} or \code{"pareto"}. Defaults to \code{"unit"}.
#' @param uni.type This argument specifies the type of constraint, if any, to be placed on the uniquenesses/idiosyncratic variances, i.e. whether a general diagonal matrix or isotropic diagonal matrix is to be assumed, and in turn whether these matrices are constrained to be equal across clusters. The default \code{"unconstrained"} corresponds to factor analysis (and mixtures thereof), whereas \code{"isotropic"} corresponds to probabilistic principal components analysers (and mixtures thereof).
#'
#' Constraints \emph{may} be particularly useful when \code{N <= P}, though caution is advised when employing constraints for any of the infinite factor models, especially \code{"isotropic"} and \code{"single"}, which may lead to overestimation of the number of clusters &/or factors if this specification is inappropriate. The four options correspond to the following 4 parsimonious Gaussian mixture models:
#' \describe{
#' \item{\code{"unconstrained"}}{(\strong{UUU}) - variable-specific and cluster-specific: \eqn{\Psi_g = \Psi_g}{Psi_g = Psi_g}.}
#' \item{\code{"isotropic"}}{(\strong{UUC}) - cluster-specific, equal across variables: \eqn{\Psi_g = \psi\mathcal{I}_p}{Psi_g = (sigma^2)_g I_p}.}
#' \item{\code{"constrained"}}{(\strong{UCU}) - variable-specific, equal across clusters: \eqn{\Psi_g = \Psi}{Psi_g = Psi}.}
#' \item{\code{"single"}}{(\strong{UCC}) - single value equal across clusters and variables: \eqn{\Psi_g = \psi\mathcal{I}_p}{Psi_g = sigma^2 I_p}.}
#' }
#' The first letter \strong{U} here corresponds to constraints on loadings (not yet implemented), the second letter corresponds to uniquenesses constrained/unconstrained across clusters, and the third letter corresponds to the isotropic constraint on the uniquenesses. Of course, only the third letter is of relevance for the single-cluster \code{"FA"} and \code{"IFA"} models, such that \code{"unconstrained"} and \code{"constrained"} are equivalent for these models, and so too are \code{"isotropic"} and \code{"single"}.
#' @param psi.alpha The shape of the inverse gamma prior on the uniquenesses. Defaults to 2.5. Must be greater than 1 if \code{psi.beta} is \emph{not} supplied. Otherwise be warned that values less than or equal to 1 may not bound uniquenesses sufficiently far away from 0, and the algorithm may therefore terminate. Also, excessively small values may lead to critical numerical issues and should thus be avoided.
#' @param psi.beta The scale of the inverse gamma prior on the uniquenesses. Can be either a single parameter, a vector of variable specific scales, or (if \code{psi0g} is \code{TRUE}) a matrix of variable and cluster-specific scales. If this is not supplied, \code{\link{psi_hyper}} is invoked to choose sensible values, depending on the value of \code{uni.prior} and the data size and dimension, for the \code{"MFA"} and \code{"MIFA"} models only, the value of \code{psi0g}. Excessively small values may lead to critical numerical issues and should thus be avoided.
#'
#' Note that optional arguments to \code{psi_hyper} can be supplied via the \code{...} construct here.
#' @param mu.zero The mean of the prior distribution for the mean parameter. Either a scalar of a vector of appropriate dimension. Defaults to the sample mean of the data.
#' @param sigma.mu The covariance of the prior distribution for the cluster mean parameters. Always assumed to be a diagonal matrix, and set to the identity matrix by default. Can also be a scalar by which the identity is multiplied, a vector of appropriate dimension; if supplied as a matrix, only the diagonal elements will be extracted. Specifying \code{sigma.mu=NULL} will use the diagonal entries of the sample covariance matrix: for unit-scaled data this is simply the identity again. See \code{prec.mu} for further control over the hypercovariance in the prior for the means.
#' @param prec.mu A scalar controlling the degree of flatness of the prior for the cluster means by scaling \code{sigma.mu} (i.e. multiplying every element of \code{sigma.mu} by \code{1/prec.mu}). Lower values lead to a more diffuse prior. Defaults to \code{0.01}, such that the prior is relatively non-informative by default. Of course, \code{prec.mu=1} nullifies any effect of this argument. The user can supply a scaled \code{sigma.mu} directly, but this argument is especially useful when specifying \code{sigma.mu=NULL}, such that the diagonal entries of the sample covariance matrix are used.
#' @param sigma.l A scalar controlling the diagonal covariance of the prior distribution for the loadings. Defaults to \code{1}, i.e. the identity; otherwise a diagonal matrix with non-zero entries all equal to \code{sigma.l} Only relevant for the finite factor methods.
#' @param z.init The method used to initialise the cluster labels. Defaults to model-based agglomerative hierarchical clustering via \code{"\link[mclust]{hc}"}. Other options include \code{"\link[stats]{kmeans}"} (with 10 random starts, by default), \code{\link[mclust]{Mclust}} via \code{"mclust"}, random initialisation via \code{"priors"}, and a user-supplied \code{"list"} (\code{z.list}). Not relevant for the \code{"FA"} and \code{"IFA"} methods. Arguments for the relevant functions can be passed via the \code{...} construct. For \code{"\link[mclust]{hc}"}, \code{VVV} is used by default, unless the data is high-dimensional, in which case the default is \code{EII}. The option \code{"priors"} may lead to empty components at initialisation, which will return an error.
#'
#' In any case, unless \code{z.list} is explicitly supplied, or \code{verbose} is \code{FALSE}, the initial cluster sizes will be printed to the screen to alert users to potentially bad initialisiations (e.g. heavily imbalanced initial cluster sizes).
#' @param z.list A user supplied list of cluster labels. Only relevant if \code{z.init == "z.list"}.
#' @param equal.pro Logical variable indicating whether or not the mixing mixing proportions are to be equal across clusters in the model (default = \code{FALSE}). Only relevant for the \code{"MFA"} and \code{"MIFA"} methods.
#' @param uni.prior A switch indicating whether uniquenesses scale hyperparameters are to be \code{"unconstrained"} or \code{"isotropic"}, i.e. variable-specific or not. \code{"uni.prior"} must be \code{"isotropic"} if the last letter of \code{uni.type} is \strong{C}, but can take either value otherwise. Defaults to correspond to the last letter of \code{uni.type} if that is supplied and \code{uni.prior} is not, otherwise defaults to \code{"unconstrained"} (though \code{"isotropic"} is recommended when \code{N <= P}). Only relevant when \code{psi.beta} is not supplied and \code{\link{psi_hyper}} is therefore invoked (with optional arguments passable via the \code{...} construct).
#' @param mu0g Logical indicating whether the \code{mu.zero} hyperparameter can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the \code{"MFA"} and \code{"MIFA"} methods when \code{z.list} is supplied.
#' @param psi0g Logical indicating whether the \code{psi.beta} hyperparameter(s) can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the \code{"MFA"} and \code{"MIFA"} methods when \code{z.list} is supplied, and only allowable when \code{uni.type} is one of \code{unconstrained} or \code{isotropic}.
#' @param drop0sd Logical indicating whether to drop variables with no standard deviation (defaults to \code{TRUE}). This is \emph{strongly} recommended, especially a) when \code{psi.beta} is not supplied &/or \code{sigma.mu=NULL}, and either/both are therefore estimated using the empirical covariance matrix, &/or b) if some form of posterior predictive checking is subsequently desired when calling \code{\link{get_IMIFA_results}}.
#' @param verbose Logical indicating whether to print output (e.g. run times) and a progress bar to the screen while the sampler runs. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param ... Also catches unused arguments. A number of optional arguments can be also supplied here:
#' \itemize{
#' \item{The additional \code{\link{psi_hyper}} argument \code{beta0}, especially when \code{N <= P}.}
#' \item{Additional arguments to be passed to \code{\link[mclust]{hc}} (\code{modelName} & \code{use} only), to \code{\link[mclust]{Mclust}} (\code{modelNames}, and the arguments for \code{\link[mclust]{hc}} with which \code{\link[mclust]{Mclust}} is itself initialised - \code{modelName} & \code{use}), or to \code{\link[stats]{kmeans}} (\code{iter.max} and \code{nstart} only), depending on the value of \code{z.init}.}
#' \item{Additionally, when \code{z.init="mclust"}, \code{criterion} can be passed here (can be \code{"icl"} or \code{"bic"}, the default) to control how the optimum \code{\link[mclust]{Mclust}} model to initialise with is determined.}
#' }
#' @return A named list in which the names are the names of the arguments and the values are the values of the arguments.
#' @export
#' @note Users should be careful to note that data are mean-centered (\code{centering=TRUE}) and unit-scaled (\code{scaling="unit"}) by default when supplying other parameters among the list above, especially those related in any way to \code{psi.hyper}, or to the other control functions \code{\link{mgpControl}} and \code{\link{bnpControl}}.
#' @keywords control
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link{psi_hyper}}, \code{\link[mclust]{hc}}, \code{\link[stats]{kmeans}}, \code{\link[mclust]{Mclust}}, \code{\link{mgpControl}}, \code{\link{bnpControl}}, \code{\link{storeControl}}
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' McNicholas, P. D. and Murphy, T. B. (2008) Parsimonious Gaussian mixture models, \emph{Statistics and Computing}, 18(3): 285-296.
#'
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' mixfaControl(n.iters = 25000L,
#'              burnin = n.iters/5L,
#'              thinning = 2L,
#'              centering = TRUE,
#'              scaling = c("unit", "pareto", "none"),
#'              uni.type = c("unconstrained", "isotropic",
#'                           "constrained", "single"),
#'              psi.alpha = 2.5,
#'              psi.beta = NULL,
#'              mu.zero = NULL,
#'              sigma.mu = 1L,
#'              prec.mu = 0.01,
#'              sigma.l = 1L,
#'              z.init = c("hc", "kmeans", "list", "mclust", "priors"),
#'              z.list = NULL,
#'              equal.pro = FALSE,
#'              uni.prior = c("unconstrained", "isotropic"),
#'              mu0g = FALSE,
#'              psi0g = FALSE,
#'              drop0sd = TRUE,
#'              verbose = interactive(),
#'              ...)
#' @examples
#' mfctrl <- mixfaControl(n.iters=200, prec.mu=1E-03, sigma.mu=NULL,
#'                        beta0=1, uni.type="constrained")
#'
#' # data(olive)
#' # sim  <- mcmc_IMIFA(olive, "IMIFA", mixFA=mfctrl)
#'
#' # Alternatively specify these arguments directly
#' # sim  <- mcmc_IMIFA(olive, "IMIFA", n.iters=200, prec.mu=1E-03,
#' #                    sigma.mu=NULL, beta0=1, uni.type="constrained")
  mixfaControl   <- function(n.iters = 25000L, burnin = n.iters/5L, thinning = 2L, centering = TRUE, scaling = c("unit", "pareto", "none"),
                             uni.type = c("unconstrained", "isotropic", "constrained", "single"), psi.alpha = 2.5, psi.beta = NULL, mu.zero = NULL,
                             sigma.mu = 1L, prec.mu = 0.01, sigma.l = 1L, z.init = c("hc", "kmeans", "list", "mclust", "priors"), z.list = NULL, equal.pro = FALSE,
                             uni.prior = c("unconstrained", "isotropic"), mu0g = FALSE, psi0g = FALSE, drop0sd = TRUE, verbose = interactive(), ...) {
    miss.args    <- list(uni.type = missing(uni.type), psi.beta = missing(psi.beta), mu.zero = missing(mu.zero),
                         sigma.mu = is.null(sigma.mu), z.init = missing(z.init), z.list = missing(z.list), uni.prior = missing(uni.prior))
    burnin       <- as.integer(burnin)
    n.iters      <- max(burnin + 1L, as.integer(n.iters))
    thinning     <- as.integer(thinning)
    if(any(!is.integer(n.iters),
       length(n.iters)     != 1))          stop("'n.iters' must be a single integer", call.=FALSE)
    if(any(!is.integer(burnin),    burnin     < 0,
       length(burnin)      != 1))          stop("'burnin' must be a single non-negative integer", call.=FALSE)
    if(any(!is.integer(thinning),  thinning   < 1,
       length(thinning)    != 1))          stop("'thinning' must be a single strictly positive integer", call.=FALSE)
    if(any(!is.logical(centering),
       length(centering)   != 1))          stop("'centering' must be a single logical indicator", call.=FALSE)
    if(!all(is.character(scaling)))        stop("'scaling' must be a character vector of length 1", call.=FALSE)
    scaling      <- match.arg(scaling)
    if(any(!is.numeric(prec.mu),   prec.mu   <= 0,
       length(prec.mu)     != 1))          stop("'prec.mu' must be a single strictly positive number", call.=FALSE)
    if(!all(is.character(uni.type)))       stop("'uni.type' must be a character vector of length 1", call.=FALSE)
    uni.type     <- match.arg(uni.type)
    if(any(!is.numeric(psi.alpha), psi.alpha <= 0,
       length(psi.alpha)   != 1))          stop("'psi.alpha' must be a single strictly positive number", call.=FALSE)
    if(any(!is.numeric(sigma.l),   sigma.l   <= 0,
       length(sigma.l)     != 1))          stop("'sigma.l' must be a single strictly positive number", call.=FALSE)
    if(!all(is.character(z.init)))         stop("'z.init' must be a character vector of length 1", call.=FALSE)
    z.init       <- match.arg(z.init)
    if(length(equal.pro)    > 1 ||
       !is.logical(equal.pro))             stop("'equal.pro' must be a single logical indicator", call.=FALSE)
    if(!all(is.character(uni.prior)))      stop("'uni.prior' must be a character vector of length 1", call.=FALSE)
    uni.prior    <- match.arg(uni.prior)
    if(any(!is.logical(mu0g),
       length(mu0g)        != 1))          stop("'mu0g' must be a single logical indicator", call.=FALSE)
    if(any(!is.logical(psi0g),
       length(psi0g)       != 1))          stop("'psi0g' must be a single logical indicator", call.=FALSE)
    if(any(!is.logical(drop0sd),
       length(drop0sd)     != 1))          stop("'drop0sd' must be a single logical indicator", call.=FALSE)
    if(any(!is.logical(verbose),
       length(verbose)     != 1))          stop("'verbose' must be a single logical indicator", call.=FALSE)
    mixfa        <- list(n.iters = n.iters, burnin = burnin, thinning = thinning, centering = centering, scaling = scaling, uni.type = uni.type, psi.alpha = psi.alpha,
                         psi.beta = psi.beta, mu.zero = mu.zero, sigma.mu = sigma.mu, prec.mu = prec.mu, sigma.l = sigma.l, z.init = z.init, z.list = z.list,
                         equal.pro = equal.pro, uni.prior = uni.prior, mu0g = mu0g, psi0g = psi0g, drop0sd = drop0sd, verbose = verbose)
    dots         <- list(...)
    dots         <- dots[unique(names(dots))]
    mixfa        <- c(mixfa, list(dots = dots[!(names(dots) %in% names(mixfa))]))
    attr(mixfa, "Missing") <- miss.args
      return(mixfa)
  }

#' Control settings for the Bayesian Nonparametric priors for infinite mixture models (or shrinkage priors for overfitted mixtures)
#'
#' Supplies a list of arguments for use in \code{\link{mcmc_IMIFA}} pertaining to the use of the Bayesian Nonparametric Pitman-Yor / Dirichlet process priors with the infinite mixture models \code{"IMFA"} and \code{"IMIFA"}. Certain arguments related to the Dirichlet concentration parameter for the overfitted mixtures \code{"OMFA"} and \code{"OMIFA"} can be supplied in this manner also.
#' @param learn.alpha
#' \describe{
#' \item{For the \code{"IMFA"} and \code{"IMIFA"} methods:}{A logical indicating whether the Pitman-Yor / Dirichlet process concentration parameter is to be learned (defaults to \code{TRUE}), or remain fixed for the duration of the chain. If being learned, a Ga(a, b) prior is assumed for \code{alpha}; updates take place via Gibbs sampling when \code{discount} is zero and via Metropolis-Hastings when \code{discount > 0}. If not being learned, \code{alpha} \emph{must} be supplied.
#'
#' In the special case of \code{discount < 0}, \code{alpha} must be supplied as a positive integer multiple of \code{abs(discount)}; in this instance, \code{learn.alpha} is forced to \code{TRUE} and \code{alpha} is updated with the changing number of components as the positive integer.}
#' \item{For the \code{"OMFA"} and \code{"OMIFA"} methods:}{A logical indicating whether the Dirichlet concentration parameter is to be learned (defaults to \code{TRUE}) or remain fixed for the duration of the chain. If being learned, a Ga(a, b * G) is assumed for \code{alpha}, where G is the number of mixture components \code{range.G}, and updates take place via Metropolis-Hastings. If not being learned \code{alpha} \emph{must} be supplied.}
#' }
#' @param alpha.hyper
#' \describe{
#' \item{For the \code{"IMFA"} and \code{"IMIFA"} methods:}{A vector of length 2 giving hyperparameters for the prior on the Pitman-Yor / Dirichlet process concentration parameter \code{alpha}. If \code{isTRUE(learn.alpha)}, these are shape and rate parameters of a Gamma distribution. Defaults to Ga(\code{2}, \code{4}). Choosing a larger rate is particularly important, as it encourages clustering. The prior is shifted to have support on (\code{-discount}, \code{Inf}) when non-zero \code{discount} is supplied and remains fixed (i.e. \code{learn.d=FALSE}) or when \code{learn.d=TRUE}.}
#' \item{For the \code{"OMFA"} and \code{"OMIFA"} methods:}{A vector of length 2 giving hyperparameters a and b for the prior on the Dirichlet concentration parameter \code{alpha}. If \code{isTRUE(learn.alpha)}, these are shape and rate parameters of a Gamma distribution. Defaults to Ga(2, 4). Note that the supplied rate will be multiplied by \code{range.G}, to encourage clustering, such that the form of the prior is Ga(a, b * G).}
#' }
#' @param discount The discount parameter used when generalising the Dirichlet process to the Pitman-Yor process. Defaults to 0, but typically must lie in the interval [0, 1). If greater than zero, \code{alpha} can be supplied greater than \code{-discount}. By default, Metropolis-Hastings steps are invoked for updating this parameter via \code{learn.d}.
#'
#' The special case of \code{discount < 0} is allowed, in which case \code{learn.d=FALSE} is forced and \code{alpha} must be supplied as a positive integer multiple of \code{abs(discount)}. Fixing \code{discount > 0.5} is discouraged (see \code{learn.alpha}).
#' @param learn.d Logical indicating whether the \code{discount} parameter is to be updated via Metropolis-Hastings (defaults to \code{TRUE}, unless \code{discount} is supplied as a negative value).
#' @param d.hyper Hyperparameters for the Beta(a,b) prior on the \code{discount} parameter. Defaults to Beta(1,1), i.e. Uniform(0,1).
#' @param ind.slice Logical indicating whether the independent slice-efficient sampler is to be employed (defaults to \code{TRUE}). If \code{FALSE} the dependent slice-efficient sampler is employed, whereby the slice sequence \eqn{\xi_1,\ldots,\xi_g}{xi_1,...,xi_g} is equal to the decreasingly ordered mixing proportions.
#' @param rho Parameter controlling the rate of geometric decay for the independent slice-efficient sampler, s.t. \eqn{\xi=(1-\rho)\rho^{g-1}}{xi = (1 - rho)rho^(g-1)}. Must lie in the interval [0, 1). Higher values are associated with better mixing but longer run times. Defaults to 0.75, but 0.5 is an interesting special case which guarantees that the slice sequence \eqn{\xi_1,\ldots,\xi_g}{xi_1,...,xi_g} is equal to the \emph{expectation} of the decreasingly ordered mixing proportions. Only relevant when \code{ind.slice} is \code{TRUE}.
#' @param trunc.G The maximum number of allowable and storable clusters under the \code{"IMIFA"} and \code{"IMFA"} models. The number of active clusters to be sampled at each iteration is adaptively truncated, with \code{trunc.G} as an upper limit for storage reasons. Defaults to \code{max(min(N-1, 50), range.G))} and must satisfy \code{range.G <= trunc.G < N}. Note that large values of \code{trunc.G} may lead to memory capacity issues.
#' @param kappa The spike-and-slab prior distribution on the \code{discount} hyperparameter is assumed to be a mixture with point-mass at zero and a continuous Beta(a,b) distribution. \code{kappa} gives the weight of the point mass at zero (the 'spike'). Must lie in the interval [0,1]. Defaults to 0.5. Only relevant when \code{isTRUE(learn.d)}. A value of 0 ensures non-zero discount values (i.e. Pitman-Yor) at all times, and \emph{vice versa}. Note that \code{kappa} will default to exactly 0 if \code{alpha<=0} and \code{learn.alpha=FALSE}.
#' @param IM.lab.sw Logical indicating whether the two forced label switching moves are to be implemented (defaults to \code{TRUE}) when running one of the infinite mixture models.
#' @param zeta
#' \describe{
#' \item{For the \code{"IMFA"} and \code{"IMIFA"} methods:}{Tuning parameter controlling the acceptance rate of the random-walk proposal for the Metropolis-Hastings steps when \code{learn.alpha=TRUE}, where \code{2 * zeta} gives the full width of the uniform proposal distribution. These steps are only invoked when either \code{discount} is non-zero and fixed or \code{learn.d=TRUE}, otherwise \code{alpha} is learned by Gibbs updates. Must be strictly positive (if invoked). Defaults to \code{2}.}
#' \item{For the \code{"OMFA"} and \code{"OMIFA"} methods:}{Tuning parameter controlling the standard deviation of the log-normal proposal for the Metropolis-Hastings steps when \code{learn.alpha=TRUE}. Must be strictly positive (if invoked). Defaults to \code{0.75}.}
#' }
#' @param tune.zeta A list with the following named arguments, used for tuning \code{zeta} (which is either the width of the uniform proposal for the \code{"IMFA"} or \code{"IMIFA"} methods or the standard deviation of the log-normal proposal for the \code{"OMFA"} or \code{"OMIFA"} methods) for \code{alpha}, via diminishing Robbins-Monro type adaptation, when the \code{alpha} parameter is learned via Metropolis-Hastings steps:
#' \describe{
#' \item{\code{heat}}{The initial adaptation intensity/step-size, such that larger values lead to larger updates. Must be strictly greater than zero. Defaults to 1 if not supplied but other elements of \code{tune.zeta} are.}
#' \item{\code{lambda}}{Iteration rescaling parameter which controls the speed at which adaptation diminishes, such that lower values cause the contribution of later iterations to diminish more slowly. Must lie in the interval (0.5, 1]. Defaults to 1 if not supplied but other elements of \code{tune.zeta} are.}
#' \item{\code{target}}{The target acceptance rate. Must lie in the interval [0, 1]. Defaults to 0.441, which is optimum for univariate targets, if not supplied but other elements of \code{tune.zeta} are.}
#' \item{\code{start.zeta}}{The iteration at which diminishing adaptation begins. Defaults to \code{100}.}
#' \item{\code{stop.zeta}}{The iteration at which diminishing adaptation is to stop completely. Defaults to \code{Inf}, such that diminishing adaptation is never explicitly made to stop. Must be greater than \code{start.zeta}.}
#' }
#' At least one \code{tune.zeta} argument must be supplied for diminishing adaptation to be invoked. \code{tune.zeta} arguments are only relevant when \code{learn.alpha} is \code{TRUE} (and, for the \code{"IMFA"} and \code{"IMIFA"} methods, when either of the following is also true: the \code{discount} remains fixed at a non-zero value, or when \code{learn.d} is \code{TRUE} and \code{kappa < 1}). Since Gibbs steps are invoked for updating \code{alpha} when \code{discount == 0} under the \code{"IMFA"} or \code{"IMIFA"} methods, adaption occurs according to a running count of the number of iterations with non-zero sampled \code{discount} values for those methods.
#'
#' If diminishing adaptation is invoked, the posterior mean \code{zeta} will be stored. Since caution is advised when employing adaptation, note that acceptance rates of between 10-50\% are generally considered adequate.
#' @param ... Catches unused arguments.
#'
#' @details The crucial concentration parameter \code{alpha} is documented within the main \code{\link{mcmc_IMIFA}} function, and is relevant to all of the \code{"IMIFA"}, \code{"IMFA"}, \code{"OMIFA"}, and \code{"OMFA"} methods.
#'
#' All arguments here are relevant to the \code{"IMFA"} and \code{"IMIFA"} methods, but the following are also related to the \code{"OMFA"} and \code{"OMIFA"} methods, and may behave differently in those instances: \code{learn.alpha}, \code{alpha.hyper}, \code{zeta}, and \code{tune.zeta}.
#'
#' @return A named list in which the names are the names of the arguments related to the BNP prior(s) and the values are the values supplied to the arguments.
#' @note Certain supplied arguments will be subject to further checks within \code{\link{mcmc_IMIFA}}. \code{\link{G_priorDensity}} and \code{\link{G_moments}} can help with soliciting sensible DP/PYP priors.
#'
#' Under the \code{"IMFA"} and \code{"IMIFA"} methods, a Pitman-Yor process prior is specified by default. A Dirichlet process prior can be easily invoked when the \code{discount} is fixed at \code{0} and \code{learn.d=FALSE}. The normalized stable process can also be specified as a prior distribution, as a special case of the Pitman-Yor process, when \code{alpha} remains fixed at \code{0} and \code{learn.alpha=FALSE} (provided the \code{discount} is fixed at a strictly positive value or \code{learn.d=TRUE}). The special case of the Pitman-Yor process with negative \code{discount} is also allowed as an experimental feature for which caution is advised, though \code{learn.d} and \code{learn.alpha} are forced to \code{FALSE} and \code{TRUE}, respectively, in this instance.
#' @keywords control
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' Kalli, M., Griffin, J. E. and Walker, S. G. (2011) Slice sampling mixture models, \emph{Statistics and Computing}, 21(1): 93-105.
#' @export
#'
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link{G_priorDensity}}, \code{\link{G_moments}}, \code{\link{mixfaControl}}, \code{\link{mgpControl}}, \code{\link{storeControl}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' bnpControl(learn.alpha = TRUE,
#'            alpha.hyper = c(2L, 4L),
#'            discount = 0,
#'            learn.d = TRUE,
#'            d.hyper = c(1L, 1L),
#'            ind.slice = TRUE,
#'            rho = 0.75,
#'            trunc.G = NULL,
#'            kappa = 0.5,
#'            IM.lab.sw = TRUE,
#'            zeta = NULL,
#'            tune.zeta = list(...),
#'            ...)
#' @examples
#' bnpctrl <- bnpControl(learn.d=FALSE, ind.slice=FALSE, alpha.hyper=c(3, 3))
#'
#' # data(olive)
#' # sim   <- mcmc_IMIFA(olive, "IMIFA", n.iters=5000, BNP=bnpctrl)
#'
#' # Alternatively specify these arguments directly
#' # sim   <- mcmc_IMIFA(olive, "IMIFA", n.iters=5000, learn.d=FALSE,
#' #                     ind.slice=FALSE, alpha.hyper=c(3, 3))
  bnpControl     <- function(learn.alpha = TRUE, alpha.hyper = c(2L, 4L), discount = 0, learn.d = TRUE, d.hyper = c(1L, 1L),
                             ind.slice = TRUE, rho = 0.75, trunc.G = NULL, kappa = 0.5, IM.lab.sw = TRUE, zeta = NULL, tune.zeta = list(...), ...) {
    miss.args    <- list(discount = missing(discount), IM.lab.sw = missing(IM.lab.sw), kappa = missing(kappa),
                         trunc.G = missing(trunc.G), zeta = missing(zeta), learn.d = missing(learn.d), learn.alpha = missing(learn.alpha))
    if(any(!is.logical(learn.alpha),
           length(learn.alpha)    != 1))   stop("'learn.alpha' must be a single logical indicator", call.=FALSE)
    if(all(length(alpha.hyper)    != 2,
           learn.alpha))                   stop(paste0("'alpha.hyper' must be a vector of length 2, giving the shape and rate hyperparameters of the gamma prior for alpha when 'learn.alpha' is TRUE"), call.=FALSE)
    if(learn.alpha       &&
       any(alpha.hyper   <= 0))            stop("The shape and rate of the gamma prior for alpha must both be strictly positive", call.=FALSE)
    if(any(!is.logical(learn.d),
           length(learn.d)        != 1))   stop("'learn.d' must be a single logical indicator", call.=FALSE)
    if(all(length(d.hyper)        != 2,
           learn.d))                       stop("d.hyper' must be a vector of length 2", call.=FALSE)
    if(learn.d           &&
       any(d.hyper       <= 0))            stop("'Discount Beta prior hyperparameters must both be strictly positive", call.=FALSE)
    if(any(!is.logical(ind.slice),
           length(ind.slice)      != 1))   stop("'ind.slice' must be a single logical indicator", call.=FALSE)
    if(length(rho)        > 1     ||
      (rho >= 1  || rho   < 0))            stop("'rho' must be a single number in the interval [0, 1)", call.=FALSE)
    if(rho  < 0.5)                         warning("Are you sure 'rho' should be less than 0.5? This could adversely affect mixing\n", call.=FALSE, immediate.=TRUE)
    if(!missing(trunc.G) &&
       (length(trunc.G)   > 1     ||
        !is.numeric(trunc.G)      ||
        trunc.G  <= 0))                    stop("'trunc.G' must be a single strictly positive number", call.=FALSE)
    if(any(!is.numeric(kappa),
           length(kappa)          != 1))   stop("'kappa' must be a single number", call.=FALSE)
    if(kappa      <  0   || kappa  > 1)    stop("'kappa' must lie in the interval [0, 1]", call.=FALSE)
    if(isFALSE(learn.d)  &&
      !missing(discount) &&
      discount    > 0.5)                   warning("Fixing 'discount'>0.5 is discouraged\n", call.=FALSE, immediate.=TRUE)
    discount     <- ifelse(missing(discount), ifelse(learn.d, ifelse(kappa != 0 && stats::runif(1L) <= kappa, 0L, pmin(stats::rbeta(1L, d.hyper[1L], d.hyper[2L]), 1 - .Machine$double.eps)), 0L), discount)
    if(any(!is.numeric(discount),
           length(discount)       != 1))   stop("'discount' must be a single number", call.=FALSE)
    if(discount  >= 1)                     stop("'discount' must be less than 1",     call.=FALSE)
    if(discount   < 0)    {
      if(isTRUE(learn.d)) {
        if(isFALSE(miss.args$learn.d)) {   stop("'learn.d' must be FALSE when 'discount' is negative",      call.=FALSE)
        } else                             warning("'learn.d' forced to FALSE as 'discount' is negative\n", call.=FALSE, immediate.=TRUE)
      }
      if(!miss.args$learn.alpha   &&
         isFALSE(learn.alpha))             stop("'learn.alpha' must be TRUE when 'discount' is negative",   call.=FALSE)
     learn.d     <- FALSE
     learn.alpha <- TRUE
    }
    kappa        <- ifelse(all(!learn.d, discount == 0), 1L, kappa)
    if(all(kappa         == 0, !learn.d,
           discount      == 0))            stop("'kappa' is zero and yet 'discount' is fixed at zero:\neither learn the discount parameter or specify a non-zero value", call.=FALSE)
    if(all(kappa         == 1,  learn.d,
           discount      != 0))            stop(paste0("'kappa' is exactly 1 and yet", ifelse(learn.d, " 'discount' is being learned ", if(discount != 0) " the discount is fixed at a non-zero value"), ":\nthe discount should remain fixed at zero"), call.=FALSE)
    if(any(!is.logical(IM.lab.sw),
       length(IM.lab.sw) != 1))            stop("'IM.lab.sw' must be a single logical indicator", call.=FALSE)
    if(!missing(zeta)    &&
       any(!is.numeric(zeta),
       length(zeta)      != 1,
       zeta      <= 0))                    stop("'zeta' must be single strictly positive number", call.=FALSE)

    gibbs.def    <- all(kappa      < 1, learn.d,  learn.alpha)
    def.py       <- all(discount  != 0, !learn.d, learn.alpha)
    gibbs.may    <- gibbs.def     &&    kappa > 0
    zeta.names   <- c("heat", "lambda", "target", "start.zeta", "stop.zeta")
    zeta.null    <- stats::setNames(vapply(tune.zeta[zeta.names], is.null, logical(1L)), zeta.names)
    if(all(zeta.null)    ||
       !(tzdo    <- learn.alpha))  {
      tz         <- list(heat=0L, lambda=NULL, target=NULL, do=FALSE, start.zeta=100L, stop.zeta=Inf)
    } else  {
      tz         <- stats::setNames(tune.zeta[zeta.names], zeta.names)
      if(!inherits(tz, "list")    ||
         (!all(is.element(names(tz),
               zeta.names))       ||
          !all(lengths(tz)        == 1        |
               zeta.null)))                stop("'tune.zeta' must be a list containing named elements 'heat', 'lambda', 'target', 'start.zeta' and 'stop.zeta', all of length 1", call.=FALSE)
      if(is.null(tz$heat))   tz$heat         <- 1L
      if(is.null(tz$lambda)) tz$lambda       <- 1L
      if(is.null(tz$target)) tz$target       <- 0.441
      attr(tz, "startx")          <- is.null(tz$start.zeta)
      attr(tz, "stopx")           <- is.null(tz$stop.zeta)
      if(attr(tz, "startx")) tz$start.zeta   <- 100L
      if(attr(tz, "stopx"))  tz$stop.zeta    <- Inf
      if(!all(vapply(tz,
              is.numeric, logical(1L))))   stop("Not all elements of 'tune.zeta' are numeric", call.=FALSE)
      tz$do      <- tzdo
      if(tz$heat          < 0)             stop("Invalid 'heat': must be >= 0", call.=FALSE)
      if(tz$target        < 0     ||
         tz$target        > 1)             stop("Invalid 'target': must lie in the interval [0, 1]", call.=FALSE)
      if(tz$lambda       <= 0.5   ||
         tz$lambda        > 1)             stop("Invalid 'lambda': must lie in the interval (0.5, 1]", call.=FALSE)
      if(learn.alpha)     {
        if(tz$heat  > 0)  {
          if(isTRUE(gibbs.may))            warning("Are you sure you want to tune zeta?: Gibbs updates are possible as 'kappa' is between zero and 1\n", call.=FALSE, immediate.=TRUE)
        } else if(tz$do)                   warning("'heat' of 0 corresponds to no tuning: are you sure?\n", call.=FALSE, immediate.=TRUE)
      }
      if(tz$do   <- tzdo && 0     != tz$heat) {
       if(length(c(tz$start.zeta,
                   tz$stop.zeta)) != 2       ||
          !is.numeric(c(tz$start.zeta,
                        tz$stop.zeta)))    stop("'start.zeta' and 'stop.zeta' must both be numeric and of length 1", call.=FALSE)
      }
    }
    attr(tz,  "IM.Need") <- any(gibbs.def, def.py)
    BNP          <- list(learn.alpha = learn.alpha, a.hyper = alpha.hyper, discount = discount, learn.d = learn.d, d.hyper = d.hyper,
                         rho = rho, ind.slice = ind.slice, trunc.G = trunc.G, kappa = kappa, IM.lab.sw = IM.lab.sw, zeta = zeta, tune.zeta = tz)
    attr(BNP, "Missing") <- miss.args
      BNP
  }

#' Control settings for the MGP prior and AGS for infinite factor models
#'
#' Supplies a list of arguments for use in \code{\link{mcmc_IMIFA}} pertaining to the use of the multiplicative gamma process (MGP) shrinkage prior and adaptive Gibbs sampler (AGS) for use with the infinite factor models \code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"}.
#' @param alpha.d1 Shape hyperparameter of the column shrinkage on the first column of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to \code{2.1}.
#' @param alpha.d2 Shape hyperparameter of the column shrinkage on the subsequent columns of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to \code{3.1}.
#' @param phi.hyper A vector of length 2 giving the shape and rate hyperparameters for the gamma prior on the local shrinkage parameters. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to \code{c(3, 2)}. It is suggested that the rate be <= shape minus 1 to induce local shrinkage, though the cumulative shrinkage property is unaffected by these hyperparameters. Excessively small values may lead to critical numerical issues and should thus be avoided; indeed it is \emph{suggested} that the shape be >=1.
#' @param sigma.hyper A vector of length 2 giving the shape and rate hyperparameters for the gamma prior on the cluster shrinkage parameters. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to \code{c(3, 2)}. Again, it is \emph{suggested} that the shape be >= 1. Only relevant for the \code{"IMIFA"}, \code{"OMIFA"}, and \code{"MIFA"} methods when \code{isTRUE(cluster.shrink)}.
#' @param prop Proportion of loadings elements within the neighbourhood \code{eps} of zero necessary to consider a loadings column redundant. Defaults to \code{floor(0.7 * P)/P}, where \code{P} is the number of variables in the data set. However, if the data set is univariate or bivariate, the default is \code{0.5} (see Note).
#' @param eps Neighbourhood epsilon of zero within which a loadings entry is considered negligible according to \code{prop}. Defaults to \code{0.1}. Must be positive.
#' @param adapt A logical value indicating whether adaptation of the number of cluster-specific factors is to take place when the MGP prior is employed. Defaults to \code{TRUE}. Specifying \code{FALSE} and supplying \code{range.Q} within \code{\link{mcmc_IMIFA}} provides a means to either approximate the infinite factor model with a fixed high truncation level, or to use the MGP prior in a finite factor context, however this is NOT recommended for the \code{"OMIFA"} and \code{"IMIFA"} methods.
#' @param forceQg A logical indicating whether the upper limit on the number of cluster-specific factors \code{Q} is also cluster-specific. Defaults to \code{FALSE}: when \code{TRUE}, the number of factors in each cluster is kept below the number of observations in each cluster, in addition to the bound defined by \code{range.Q}. Only relevant for the \code{"IMIFA"}, \code{"OMIFA"}, and \code{"MIFA"} methods, and only invoked when \code{adapt} is \code{TRUE}. May be useful for low-dimensional data sets for which identifiable solutions are desired.
#' @param cluster.shrink A logical value indicating whether to place the prior specified by \code{sigma.hyper} on the cluster shrinkage parameters. Defaults to \code{TRUE}. Specifying \code{FALSE} is equivalent to fixing all cluster shrinkage parameters to 1. Only relevant for the \code{"IMIFA"}, \code{"OMIFA"}, and \code{"MIFA"} methods. If invoked, the posterior mean cluster shrinkage factors will be reported.
#' @param truncated A logical value indicating whether the version of the MGP prior based on left-truncated gamma distributions is invoked (see Zhang et al. reference below and additional relevant documentation in \code{\link{ltrgamma}} and \code{\link{MGP_check}}). Defaults to \code{FALSE}. Note that, when \code{TRUE}, only the shrinkage parameters for the first loadings column are affected and the conditions needed to pass \code{\link{MGP_check}} are much less strict. Moreover, more desirable shrinkage properties are easily obtained, at the expense of slightly longer run times.
#' @param b0,b1 Intercept & slope parameters for the exponentially decaying adaptation probability:
#'
#' \code{p(iter) = 1/exp(b0 + b1 * (iter - start.AGS))}.
#'
#' Defaults to \code{0.1} & \code{0.00005}, respectively. Must be non-negative and strictly positive, respectively, to ensure diminishing adaptation.
#' @param beta.d1 Rate hyperparameter of the column shrinkage on the first column of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to 1.
#' @param beta.d2 Rate hyperparameter of the column shrinkage on the subsequent columns of the loadings according to the MGP shrinkage prior. Passed to \code{\link{MGP_check}} to ensure validity. Defaults to 1.
#' @param start.AGS The iteration at which adaptation under the AGS is to begin. Defaults to \code{burnin} for the \code{"IFA"} and \code{"MIFA"} methods, defaults to \code{2} for the \code{"OMIFA"} and \code{"IMIFA"} methods, and defaults to \code{2} for all methods if the data set is univariate or bivariate. Cannot exceed \code{burnin}; thus defaults to the same value as \code{burnin} if necessary.
#' @param stop.AGS The iteration at which adaptation under the AGS is to stop completely. Defaults to \code{Inf}, such that the AGS is never explicitly forced to stop (thereby overriding the diminishing adaptation probability after \code{stop.AGS}). Must be greater than \code{start.AGS}. The diminishing adaptation probability prior to \code{stop.AGS} is still governed by the arguments \code{b0} and \code{b1}.
#' @param delta0g Logical indicating whether the \code{alpha.d1} and \code{alpha.d2} hyperparameters can be cluster-specific. Defaults to \code{FALSE}. Only relevant for the \code{"MIFA"} method and only allowed when \code{z.list} is supplied within \code{\link{mcmc_IMIFA}}.
#' @param ... Catches unused arguments.
#'
#' @return A named list in which the names are the names of the arguments related to the MGP and AGS and the values are the values supplied to the arguments.
#' @export
#' @keywords control
#'
#' @note Certain supplied arguments will be subject to further checks by \code{\link{MGP_check}} to ensure the cumulative shrinkage property of the MGP prior holds according to the given parameterisation.
#'
#' The adaptive Gibbs sampler (AGS) monitors the \code{prop} of loadings elements within the neighbourhood \code{eps} of 0 and discards columns or simulates new columns on this basis. However, if at any stage the number of group-specific latent factors reaches zero, the decision to add columns is instead based on a simple binary trial with probability \code{1-prop}, as there are no loadings entries to monitor.
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link{MGP_check}}, \code{\link{ltrgamma}}, \code{\link{mixfaControl}}, \code{\link{bnpControl}}, \code{\link{storeControl}}
#' @references Murphy, K., Viroli, C., and Gormley, I. C. (2020) Infinite mixtures of infinite factor analysers, \emph{Bayesian Analysis}, 15(3): 937-963. <\href{https://projecteuclid.org/euclid.ba/1570586978}{doi:10.1214/19-BA1179}>.
#'
#' Durante, D. (2017). A note on the multiplicative gamma process, \emph{Statistics & Probability Letters}, 122: 198-204.
#'
#' Bhattacharya, A. and Dunson, D. B. (2011) Sparse Bayesian infinite factor models, \emph{Biometrika}, 98(2): 291-306.
#'
#' Zhang, X., Dunson, D. B., and Carin, L. (2011) Tree-structured infinite sparse factor model. In Getoor, L. and Scheffer, T. (Eds.), \emph{Proceedings of the 28th International Conference on Machine Learning}, ICML'11, Madison, WI, USA, pp. 785-792. Omnipress.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' mgpControl(alpha.d1 = 2.1,
#'            alpha.d2 = 3.1,
#'            phi.hyper = c(3, 2),
#'            sigma.hyper = c(3, 2),
#'            prop = 0.7,
#'            eps = 0.1,
#'            adapt = TRUE,
#'            forceQg = FALSE,
#'            cluster.shrink = TRUE,
#'            truncated = FALSE,
#'            b0 = 0.1,
#'            b1 = 5e-05,
#'            beta.d1 = 1,
#'            beta.d2 = 1,
#'            start.AGS = 2L,
#'            stop.AGS = Inf,
#'            delta0g = FALSE,
#'            ...)
#' @examples
#' mgpctrl <- mgpControl(phi.hyper=c(2.5, 1), eps=1e-02, truncated=TRUE)
#'
#' # data(olive)
#' # sim   <- mcmc_IMIFA(olive, "IMIFA", n.iters=5000, MGP=mgpctrl)
#'
#' # Alternatively specify these arguments directly
#' # sim   <- mcmc_IMIFA(olive, "IMIFA", n.iters=5000,
#' #                     phi.hyper=c(2.5, 1), eps=1e-02, truncated=TRUE)
    mgpControl   <- function(alpha.d1 = 2.1, alpha.d2 = 3.1, phi.hyper = c(3, 2), sigma.hyper = c(3, 2), prop = 0.7, eps = 1e-01, adapt = TRUE, forceQg = FALSE,
                             cluster.shrink = TRUE, truncated = FALSE, b0 = 0.1, b1 = 5e-05, beta.d1 = 1, beta.d2 = 1, start.AGS = 2L, stop.AGS = Inf, delta0g = FALSE, ...) {
      miss.args  <- list(propx = missing(prop), startAGSx = missing(start.AGS), stopAGSx = missing(stop.AGS))
      if(any(!is.numeric(alpha.d1),
             !is.numeric(alpha.d2),
             c(alpha.d1, alpha.d2)   < 1)) stop("All column shrinkage shape hyperparameter values must be numeric and at least 1", call.=FALSE)
      if(prop     > 1          ||
         prop    <= 0)                     stop("'prop' must be lie in the interval (0, 1]",    call.=FALSE)
      if(eps     <= 0)                     stop("'eps' must be greater than 0", call.=FALSE)
      if(eps     >= 1)                     warning("'eps' should typically be smaller than 1, particularly for scaled data\n", call.=FALSE, immediate.=TRUE)
      if(any(length(adapt)     != 1,
             !is.logical(adapt)))          stop("'adapt' must be a single logical indicator",   call.=FALSE)
      if(any(length(forceQg)   != 1,
             !is.logical(forceQg)))        stop("'forceQg' must be a single logical indicator", call.=FALSE)
      forceQg   <- forceQg     && adapt
      if(any(length(cluster.shrink) != 1,
             !is.logical(cluster.shrink))) stop("'cluster.shrink' must be a single logical indicator", call.=FALSE)
      if(any(length(truncated) != 1,
             !is.logical(truncated)))      stop("'truncated' must be a single logical indicator",      call.=FALSE)
      if(any(length(b0)        != 1,
             !is.numeric(b0), b0     < 0)) stop("'b0' must be a non-negative scalar to ensure valid adaptation probability", call.=FALSE)
      if(any(length(b1)        != 1,
             !is.numeric(b1), b1    <= 0)) stop("'b1' must be a single strictly positive scalar to ensure adaptation probability decreases", call.=FALSE)
      if(length(phi.hyper)     != 2 ||
         !is.numeric(phi.hyper))           stop("'phi.hyper' must be a numeric vector of length 2", call.=FALSE)
      if(any(phi.hyper   < 1E-01))         stop("Excessively small values for the local shrinkage hyperparameters will lead to critical numerical issues & should thus be avoided", call.=FALSE)
      if(any(phi.hyper         <= 0))      stop("The shape and rate in 'phi.hyper' must both be strictly positive", call.=FALSE)
      if(length(sigma.hyper)   != 2 ||
         !is.numeric(sigma.hyper))         stop("'sigma.hyper' must be a numeric vector of length 2", call.=FALSE)
      if(any(sigma.hyper < 1E-01))         stop("Excessively small values for the cluster shrinkage hyperparameters will lead to critical numerical issues & should thus be avoided", call.=FALSE)
      if(any(sigma.hyper       <= 0))      stop("The shape and rate in 'sigma.hyper' must both be strictly positive", call.=FALSE)
      if(any(!is.numeric(beta.d1),
             !is.numeric(beta.d2),
             length(beta.d1)   != 1,
             length(beta.d2)   != 1,
             c(beta.d1,  beta.d2)   <= 0)) stop("'beta.d1' and 'beta.d2' must both be numeric, of length 1, and strictly positive", call.=FALSE)
      if(any(!is.numeric(prop),
             !is.numeric(start.AGS),
             !is.numeric(stop.AGS),
             !is.numeric(eps),
             length(prop)      != 1,
             length(start.AGS) != 1,
             length(stop.AGS)  != 1,
             length(eps)       != 1))      stop("'prop', 'start.AGS', 'stop.AGS', and 'eps' must all be numeric and of length 1", call.=FALSE)
      if(any(length(delta0g)   != 1,
             !is.logical(delta0g)))        stop("'delta0g' must be a single logical indicator", call.=FALSE)
      MGPAGS     <- list(alpha.d1 = alpha.d1, alpha.d2 = alpha.d2, delta0g = delta0g, phi.hyper = phi.hyper, sigma.hyper = sigma.hyper,
                         prop = prop, epsilon = eps, adapt = adapt, forceQg = forceQg, cluster.shrink = cluster.shrink, truncated = truncated,
                         b0 = b0, b1 = b1, beta.d1 = beta.d1, beta.d2 = beta.d2, start.AGS = start.AGS, stop.AGS = stop.AGS)
      attr(MGPAGS, "Missing")  <- miss.args
        MGPAGS
    }

#' Set storage values for use with the IMIFA family of models
#'
#' Supplies a list of values for logical switches indicating whether parameters of interest (means, scores, loadings, uniquenesses, and mixing proportions) should be stored when running models from the IMIFA family via \code{\link{mcmc_IMIFA}}. It may be useful not to store certain parameters if memory is an issue.
#' @param mu.switch Logical indicating whether the means are to be stored (defaults to \code{TRUE}).
#' @param score.switch Logical indicating whether the factor scores are to be stored.
#'
#' As the array containing each sampled scores matrix tends to be amongst the largest objects to be stored, this defaults to \code{FALSE} inside \code{\link{mcmc_IMIFA}} when \code{length(range.G) * length(range.Q) > 10}, otherwise the default is \code{TRUE}. For the \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"} methods, setting this switch to \code{FALSE} also offers a slight speed-up.
#'
#' Unlike other parameters, the scores need not be stored for posterior predictive checking (see Note below).
#' @param load.switch Logical indicating whether the factor loadings are to be stored (defaults to \code{TRUE}).
#' @param psi.switch Logical indicating whether the uniquenesses are to be stored (defaults to \code{TRUE}).
#' @param pi.switch Logical indicating whether the mixing proportions are to be stored (defaults to \code{TRUE}).
#' @param ... Catches unused arguments.
#'
#' @details \code{\link{storeControl}} is provided for assigning values for IMIFA models within \code{\link{mcmc_IMIFA}}. It may be useful not to store certain parameters if memory is an issue (e.g. for large data sets or for a large number of MCMC iterations after burnin and thinning).
#'
#' @note Posterior inference and plotting won't be possible for parameters not stored.
#'
#' Non-storage of parameters will almost surely prohibit the computation of posterior predictive checking error metrics within \code{\link{get_IMIFA_results}} also. In particular, if such error metrics are desired, \code{mu.switch} and \code{psi.switch} must be \code{TRUE} for all but the \code{"FA"} and \code{"IFA"} models, \code{load.switch} must be \code{TRUE} for all but the entirely zero-factor models, and \code{pi.switch} must be \code{TRUE} for models with clustering structure and unequal mixing proportions for all but the PPRE metric. \code{score.switch=TRUE} is not required for any posterior predictive checking.
#'
#' Finally, if loadings are not stored but scores are, caution is advised when examining posterior scores as Procrustes rotation will not occur within \code{\link{get_IMIFA_results}}.
#'
#' @return A named vector in which the names are the names of the storage switches and the values are logicals indicating whether that parameter is to be stored. The list also contains as an attribute a logical for each switch indicating whether it was actually supplied (\code{TRUE}) or the default was accepted (\code{FALSE}).
#' @export
#' @keywords control
#'
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link{get_IMIFA_results}}, \code{\link{mixfaControl}}, \code{\link{mgpControl}}, \code{\link{bnpControl}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @usage
#' storeControl(mu.switch = TRUE,
#'              score.switch = TRUE,
#'              load.switch = TRUE,
#'              psi.switch = TRUE,
#'              pi.switch = TRUE,
#'              ...)
#' @examples
#' stctrl <- storeControl(score.switch=FALSE)
#'
#' # data(olive)
#' # sim  <- mcmc_IMIFA(olive, "IMIFA", n.iters=5000, storage=stctrl)
#'
#' # Alternatively specify these arguments directly
#' # sim  <- mcmc_IMIFA(olive, "IMIFA", n.iters=5000, score.switch=FALSE)
   storeControl  <- function(mu.switch = TRUE, score.switch = TRUE, load.switch = TRUE, psi.switch = TRUE, pi.switch = TRUE, ...) {
      switches   <- c(mu.sw=mu.switch, s.sw=score.switch, l.sw=load.switch, psi.sw=psi.switch,pi.sw=pi.switch)
      attr(switches, "Missing") <- c(mu.sw=missing(mu.switch), s.sw=missing(score.switch), l.sw=missing(load.switch), psi.sw=missing(psi.switch), pi.sw=missing(pi.switch))
      if(any(!is.logical(switches)))       stop("All logical parameter storage switches must be single logical indicators", call.=FALSE)
        switches
   }

#' Decompose factor scores by cluster
#'
#' Takes posterior summaries of the overall factor scores matrix and returns lists of sub-matrices corresponding to the \code{G}-cluster MAP partition.
#' @param res An object of class \code{"Results_IMIFA"} generated by \code{\link{get_IMIFA_results}}.
#' @param dropQ A logical indicating whether columns of the factor scores matrix should be dropped such that the number of columns in each sub-matrix corresponds to the cluster-specific number of factors (if the number of factors is indeed cluster-specific). When \code{FALSE} (the default), the number of columns instead remains common to all sub-matrices - given by the largest of the cluster-specific numbers of latent factors.
#'
#' Note that this argument is irrelevant (i.e. always \code{FALSE}) for the finite factor methods (\code{"FA"}, \code{"MFA"}, \code{"OMFA"}, and \code{"IMFA"}).
#'
#' @details Under the models in the IMIFA family, there exists only one factor scores matrix. For the finite factor methods, this has dimensions \code{N * Q}.
#'
#' For the infinite factor methods (\code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"}), the factor scores matrix has dimensions \code{N * Qmax}, where \code{Qmax} is the largest of the cluster-specific numbers of latent factors \eqn{q_1,\ldots,q_g}{Q1,...,Qg}. Entries of this matrix thus may have been padded out with zero entries, as appropriate, prior to the Procrustes rotation-based correction applied within \code{\link{get_IMIFA_results}} (thus now these entries will be near-zero).
#'
#' In partitioning rows of the factor scores matrix into the same clusters the corresponding observations themselves belong to according to the MAP clustering, the number of columns \emph{may} vary according to the cluster-specific numbers of latent factors (depending on the value of \code{dropQ} and the method employed).
#' @return For models which achieve clustering, a list of lists (say \code{x}) decomposing the posterior mean scores (\code{x$post.eta}), the associated variance estimates (\code{x$var.eta}) and credible intervals (\code{x$ci.eta}), and the last valid sample of the scores (\code{x$last.eta}) into lists of length \code{G}, corresponding to the MAP clustering, with varying or common numbers of cluster-specific factors (depending on the value of \code{dropQ} and the method employed).
#'
#' For models with only one component, or the \code{"FA"} and \code{"IFA"} methods, scores cannot be decomposed, and posterior summaries of the scores will be returned unchanged.
#' @keywords utility
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @export
#' @usage
#' scores_MAP(res,
#'            dropQ = FALSE)
#' @seealso \code{\link{get_IMIFA_results}}
#'
#' @examples
#' \donttest{data(coffee)
#' sim <- mcmc_IMIFA(coffee, n.iters=1000)
#' res <- get_IMIFA_results(sim)
#'
#' # Examine the single posterior mean scores matrix
#' res$Scores$post.eta
#'
#' # Decompose into G matrices, common numbers of columns
#' eta <- scores_MAP(res)
#' eta$post.eta
#'
#' # Allow the number of columns be cluster-specific
#' scores_MAP(res, dropQ=TRUE)$post.eta}
   scores_MAP    <- function(res, dropQ = FALSE) {
      UseMethod("scores_MAP")
   }

#' @method scores_MAP Results_IMIFA
#' @export
   scores_MAP.Results_IMIFA     <- function(res, dropQ = FALSE) {
     if(!attr(res, "Switch")["s.sw"])       stop("Scores not stored!", call.=FALSE)
     resGQ       <- res$GQ.results
     scores      <- res$Scores
     if((G       <- resGQ$G)    == 1) {     message("Clustering has not taken place; cannot decompose scores\n")
       return(list(post.eta = scores$post.eta, var = scores$var.eta, ci.eta = scores$ci.eta, last.eta = scores$last.eta))
     }
     if(length(dropQ) != 1 ||
        !is.logical(dropQ))                 stop("'dropQ' must be a single logical indicator", call.=FALSE)
     varyQ       <- is.element(attr(res, "Method"), c("IFA", "MIFA", "OMIFA", "IMIFA"))
     if(all(!varyQ, dropQ))                 message("'dropQ' ignored as a finite factor method was employed\n")
     dropQ       <- dropQ  && varyQ
     Gseq        <- seq_len(G)
     Q           <- if(isTRUE(dropQ)) resGQ$Q else rep(max(resGQ$Q), G)
     MAP         <- res$Clust$MAP
       return(lapply(list(post.eta = lapply(Gseq, function(g) scores$post.eta[MAP == g, seq_len(Q[g]), drop=FALSE]),
                          var.eta  = lapply(Gseq, function(g) scores$var.eta[MAP  == g, seq_len(Q[g]), drop=FALSE]),
                          ci.eta   = lapply(Gseq, function(g) scores$ci.eta[,MAP  == g, seq_len(Q[g]), drop=FALSE]),
                          last.eta = lapply(Gseq, function(g) scores$last.eta[MAP == g, seq_len(Q[g]), drop=FALSE])),
                     stats::setNames, paste0("Cluster", Gseq)))
   }

#' Show the NEWS file
#'
#' Show the \code{NEWS} file of the \code{IMIFA} package.
#' @return The \code{IMIFA} \code{NEWS} file, provided the session is interactive.
#' @export
#' @keywords utility
#'
#' @usage IMIFA_news()
#' @examples
#' IMIFA_news()
  IMIFA_news     <- function() {
    news         <- file.path(system.file(package  = "IMIFA"), "NEWS.md")
    if(interactive()) file.show(news) else message("The session is not interactive\n")
  }

#' Pareto Scaling
#'
#' Pareto scaling of a numeric matrix, with or without centering. Observations are scaled by the square-root of their column-wise standard deviations.
#' @param x A numeric matrix-like object.
#' @param centering A logical vector indicating whether centering is to be applied (default=\code{TRUE}).
#'
#' @return The Pareto scaled version of the matrix \code{x}.
#' @importFrom matrixStats "colMeans2"
#' @importFrom Rfast "colVars"
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references van den Berg, R.A., Hoefsloot, H.C.J, Westerhuis, J.A. and Smilde, A.K. and van der Werf, M.J. (2006) Centering, scaling, and transformations: improving the biological information content of metabolomics data. \emph{BMC Genomics}, 7(142).
#' @keywords utility
#' @usage
#' pareto_scale(x,
#'              centering = TRUE)
#' @examples
#' dat  <- pareto_scale(olive[,-c(1,2)])
  pareto_scale <- function(x, centering=TRUE) {
    if(length(centering) != 1 ||
       !is.logical(centering))           stop("'centering' must be a single logical indicator", call.=FALSE)
    x          <- as.matrix(x)
      .scale2(x, centering, sqrt(colVars(x, std=TRUE)))
  }

#' Left Truncated Gamma Distributions
#'
#' Functions to draw pseudo-random numbers from, or calculate the expectation of, left-truncated gamma distributions (see Details below).
#' @param n Number of observations to generate.
#' @param shape Shape parameter for the desired gamma distribution. Must be strictly positive
#' @param rate Rate parameter for the desired gamma distribution. Must be strictly positive.
#' @param trunc The point of left truncation (corresponding to \eqn{\tau}{tau} below). Defaults to \code{1}. Must be non-negative. When \code{inverse} is \code{TRUE}, this becomes the point of \emph{right} truncation.
#'
#' @details The left-truncated gamma distribution has PDF:
#' \deqn{f(x|\alpha, \beta) = \frac{\beta^\alpha}{(\Gamma(\alpha)-\Gamma(\alpha, \tau\beta))}x^{\alpha-1}e^{-x\beta}}
#' for \eqn{0\le\tau\le x}, and \eqn{\min(\tau,\beta) > 0}{min(tau, beta) > 0}, where \eqn{\alpha}{alpha} and \eqn{\beta}{beta} are the \code{shape} and \code{rate} parameters, respectively, \eqn{\tau}{tau} is the cutoff point at which \code{trunc}ation occurs, and \eqn{\Gamma(\alpha, \tau\beta)} is the upper incomplete gamma function.
#'
#' @note \code{rltrgamma} is invoked internally for the \code{"IFA"}, \code{"MIFA"}, \code{"OMIFA"}, and \code{"IMIFA"} models to draw column shrinkage parameters for all but the first loadings column under the MGP prior when \code{truncated=TRUE} (which is \strong{not} the default) is supplied to \code{mgpControl}, at the expense of slightly longer run times. \code{exp_ltrgamma} is used within \code{MGP_check} to check the validity of the MGP hyperparameters when \code{truncated=TRUE} (which is again, \strong{not} the default). Both functions always assume \code{trunc=1} for these internal usages.
#'
#' Note also that no arguments are recycled, i.e. all arguments must be of length \code{1}.
#' @return For \code{rltrgamma}, a vector of length \code{n} giving draws from the left-truncated gamma distribution with the specified \code{shape} and \code{rate} parameters, and truncation point \code{trunc}.
#'
#' For \code{exp_ltrgamma}, the expected value of a left-truncated (inverse) gamma distribution.
#' @export
#' @seealso \code{\link{mgpControl}}, \code{\link{MGP_check}}
#' @name ltrgamma
#' @rdname ltrgamma
#' @references Dagpunar, J. S. (1978) Sampling of variates from a truncated gamma distribution, \emph{Statistical Computation and Simulation}, 8(1): 59-64.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' rltrgamma(n,
#'           shape,
#'           rate = 1,
#'           trunc = 1)
#'
#' @examples
#' # Generate left-truncated Ga(3.1, 2.1, 1) variates
#' rltrgamma(n=10, shape=3.1, rate=2.1)
#'
#' # Calculate the expectation of a Ga(3.1, 2.1, 1) distribution
#' exp_ltrgamma(shape=3.1, rate=2.1)
#'
#' # Calculate the expectation of an inverse gamma distribution right-truncated at 2
#' exp_ltrgamma(shape=3.1, rate=2.1, trunc=2, inverse=TRUE)
  rltrgamma    <- function(n, shape, rate = 1, trunc = 1) {
    EV         <- as.list(environment())
    if(!all(vapply(EV, is.numeric,
                   logical(1L))))        stop("All arguments must be numeric",              call.=FALSE)
    if(!all(lengths(EV) == 1))           stop("All arguments must be of length 1",          call.=FALSE)
    if(any(n   != floor(n)))             stop("'n' must be coercible to an integer",        call.=FALSE)
    if(shape   <= 1)                     stop("'shape' parameter must be greater than 1",   call.=FALSE)
    if(rate    <= 0)                     stop("'rate' parameter must be strictly positive", call.=FALSE)
    if(trunc    < 0)                     stop("'trunc' cannot be negative",                 call.=FALSE)
    U          <- stats::runif(n)
    G          <- stats::pgamma(trunc, shape, rate)
      stats::qgamma(G + U * (1 - G),   shape, rate)
  }

#' @keywords utility
#' @param inverse A logical indicating whether to calculate the expectation for a right-truncated \emph{inverse} gamma distribution instead of a left-truncated gamma distribution. Defaults to \code{FALSE}.
#' @rdname ltrgamma
#' @usage
#' exp_ltrgamma(shape,
#'              rate = 1,
#'              trunc = 1,
#'              inverse = FALSE)
#' @export
  exp_ltrgamma <- function(shape, rate = 1, trunc = 1, inverse = FALSE) {
    if(any(!is.numeric(shape),
           shape        <= 0))           stop("'shape' must be strictly positive and numeric",         call.=FALSE)
    if(any(!is.numeric(rate),
           rate         <= 0))           stop("'rate' must be strictly positive and numeric",          call.=FALSE)
    if(any(!is.numeric(trunc),
           trunc         < 0))           stop("'trunc' must be strictly non-negative and numeric",     call.=FALSE)
    if(isTRUE(inverse))  {
      if(any(shape      <= 1))           stop("'shape' must be > 1 when inverse=TRUE",                 call.=FALSE)
      stats::pgamma(rate / trunc, shape - 1, lower.tail=FALSE)/stats::pgamma(rate / trunc, shape, lower.tail=FALSE) * rate/(shape - 1)
    } else      {
      stats::pgamma(rate * trunc, shape + 1, lower.tail=FALSE)/stats::pgamma(rate * trunc, shape, lower.tail=FALSE) * shape/rate
    }
  }

  # Other Hidden Functions
    .a_drop      <- function(x, drop = TRUE, named.vector = TRUE, one.d.array = FALSE, ...) {
      if(is.null(dim(x)))                  stop("require an object with a dim attribute", call.=FALSE)
      x.dim      <- dim(x)
      if(is.logical(drop))  {
        if(length(drop) != length(x.dim))  stop("length of drop is not equal length of dim(x)", call.=FALSE)
        drop     <- which(drop)
      } else if(is.character(drop))  {
        if(any(is.na(i  <- match(drop,
                     names(x.dim)))))      stop("dimension names ", paste("'", drop[is.na(i)], "'", sep="", collapse=" "), " not found in x", call.=FALSE)
        drop     <- i
      } else if(is.null(drop))       {
        drop     <- numeric(0L)
      }
      if(!is.numeric(drop) ||
        any(is.na(drop))   ||
        any(drop  < 1 |
            drop  > length(x.dim)))        stop("drop must contain dimension numbers", call.=FALSE)
      if(!all(x.dim[drop]  == 1))          stop("dimensions to drop (", paste(drop, collapse = ", "), ") do not have length 1", call.=FALSE)
      x.dimnames        <- dimnames(x)
      dimnames(x)       <- NULL
      dim(x)     <- NULL
      keep       <- setdiff(seq_along(x.dim), drop)
      if(length(x.dim[keep]) > 1 || (length(x.dim[keep]) == 1 && one.d.array)) {
       dim(x)    <- x.dim[keep]
       if(!is.null(x.dimnames)) dimnames(x) <- x.dimnames[keep]
      } else if(length(x.dim[keep]) == 1 && named.vector) {
       names(x)  <- x.dimnames[keep][[1L]]
      }
        x
    }

    .as_numchar  <- function(x) {
       tryCatch(as.numeric(x), warning=function(w) as.numeric(factor(x, labels=seq_along(unique(x)))))
    }

    .chol        <- function(x, ...) tryCatch(chol(x, ...), error=function(e)    {
      d          <- nrow(x)
      eigs       <- eigen(x, symmetric = TRUE)
      eval       <- eigs$values
      evec       <- eigs$vectors
        return(chol(x + evec %*% tcrossprod(diag(pmax.int(0L, 2 * max(abs(eval)) * d * .Machine$double.eps - eval), d), evec), ...))
      }
    )

    #' @importFrom matrixStats "colSums2" "rowSums2"
    .class_agreement    <- function(tab, match.names = FALSE) {
      n          <- sum(tab)
      ni         <- rowSums2(tab)
      nj         <- colSums2(tab)
      if(match.names && !is.null(dimnames(tab))) {
        lev      <- intersect(colnames(tab), rownames(tab))
        p0       <- sum(diag(tab[lev, lev]))/n
        pc       <- sum(ni[lev]   * nj[lev])/n^2
      } else {
        m        <- seq_len(min(length(ni), length(nj)))
        p0       <- sum(diag(tab[m, m]))/n
        pc       <- sum((ni[m]/n) * (nj[m]/n))
      }
      n2         <- choose(n, 2)
      rand       <- 1 + (sum(tab^2)  - (sum(ni^2) + sum(nj^2))/2)/n2
      nis2       <- sum(choose(ni[ni > 1], 2))
      njs2       <- sum(choose(nj[nj > 1], 2))
      nijn2      <- (nis2 * njs2)/n2
      crand      <- (sum(choose(tab[tab > 1], 2)) - nijn2)/((nis2 + njs2)/2 - nijn2)
        list(diag = p0, kappa = (p0 - pc)/(1 - pc), rand = rand, crand = crand)
    }

    .clean_args  <- function(argstr, fn, exclude.repeats = FALSE, exclude.other = NULL, dots.ok = TRUE) {
      fnargs     <- names(formals(fn))
      if(length(argstr) > 0 && !("..." %in% fnargs && dots.ok)) {
        badargs  <- names(argstr)[!vapply(names(argstr), "%in%", logical(1L), c(fnargs, ""))]
        for(i in badargs) argstr[[i]]     <- NULL
      }
      if(exclude.repeats) {
        ntab     <- table(names(argstr))
        badargs  <- names(ntab)[ntab > 1 & names(ntab) != ""]
        for(i in badargs) argstr[[i]]     <- NULL
      }
      for(i in exclude.other) argstr[[i]] <- NULL
        argstr
    }

    #' @importFrom matrixStats "colMeans2"
    .col_vars    <- function(x, std = FALSE, avg = NULL) { # formerly replaced Rfast::colVars
      if(length(std) > 1 ||
         !is.logical(std))                 stop("'std' must be a single logical indicator")
      if(!is.matrix(x))                    stop("'x' must be a matrix")
      m          <- if(missing(avg)) colMeans2(x) else as.vector(avg)
      n          <- nrow(x)
      s          <- pmax((colMeans2(x^2) - m^2) * n/(n - 1), 0L)
        if(std) sqrt(s) else s
    }

    .detach_pkg  <- function(pkg, character.only = TRUE) {
      searches   <- paste("package", if(isTRUE(character.only)) pkg else deparse(substitute(pkg)), sep=":")
      while(searches %in% search()) {
        detach(searches, unload=TRUE, character.only=TRUE)
      }
    }

    .empty_mat   <- function(nr = 0L, nc = 0L) {
        base::matrix(0L, nrow=nr, ncol=nc)
    }

    .ent_exit    <- function(opts = options()) {
      ent        <- readline("Hit <Return> to see next plot or /hit <Esc> or type 'EXIT' to exit: ")
      options(show.error.messages=FALSE)
      on.exit(suppressWarnings(options(opts)), add=TRUE)
        if(ent  %in% c("exit", "EXIT"))    stop(call.=FALSE)
    }

    .gamma_inc   <- function(shape, trunc = 1, upper = TRUE) { # unused
        gamma(shape) * stats::pgamma(trunc, shape, rate=1, lower.tail=isFALSE(upper))
    }

    .geom_mean   <- function(x) { # unused
        return(if(any(x == 0, na.rm=TRUE)) 0L else exp(mean(log(x), na.rm=TRUE)))
    }

    .IntMult     <- Vectorize(function(a, b) {
      b          <- abs(b)
      x          <- a / b
        return(x  > 0 &
       (isTRUE(all.equal(x, round(x))) ||
        isTRUE(all.equal(0, a %% b))))
    })

    .logdensity     <- function(x, left = 0, ...) { # export and add ? ...
      d          <- tryCatch(stats::density(x, ...), error = function(e) stats::density(x))
      h          <- d$bw
      w          <- 1/stats::pnorm(left, mean = x, sd = h, lower.tail = FALSE)
        return(suppressWarnings(stats::density(x, bw = h, kernel = "gaussian", weights = w/length(x))))
    }

    .logitdensity   <- function(x, ...) { # export and add ? ...
      y          <- stats::qlogis(x[x > 0  & x < 1])
      g          <- tryCatch(stats::density(y, ...), error = function(e) stats::density(y))
      xgrid      <- stats::plogis(g$x)
      g$y        <- g$y/(xgrid * (1 - xgrid))
      g$x        <- xgrid
        return(g)
    }

    .match_classes  <- function(tab, method = "rowmax", iter = 1L, maxexact = 9L, verbose = TRUE) {
      methods    <- c("rowmax", "greedy", "exact")
      method     <- pmatch(method, methods)
      rmax       <- apply(tab, 1L, which.max)
      myseq      <- seq_len(ncol(tab))
      cn         <- colnames(tab)
      rn         <- rownames(tab)
      if(is.null(cn))  {
        cn       <- myseq
      }
      if(is.null(rn))  {
        rn       <- myseq
      }
      if(method  == 1) {
        retval   <- rmax
      }
      if(method  == 2 || method   == 3)  {
        if(ncol(tab)  != nrow(tab))        stop("Unique matching only for square tables.", call.=FALSE)
        dimnames(tab) <- list(myseq, myseq)
        cmax     <- apply(tab, 2L, which.max)
        retval   <- rep(NA, ncol(tab))
        names(retval) <- colnames(tab)
        baseok   <- cmax[rmax]     == myseq
        for(k in myseq[baseok])    {
          therow      <- (tab[k, ])[-rmax[k]]
          thecol      <- (tab[, rmax[k]])[-k]
          if(max(outer(therow, thecol, "+")) < tab[k, rmax[k]]) {
            retval[k] <- rmax[k]
          } else {
            baseok[k] <- FALSE
          }
        }
        if(isTRUE(verbose))                message("Direct agreement:", sum(baseok), "of", ncol(tab), "pairs\\n")
        if(!all(baseok))     {
          if(method   == 3)  {
            if(sum(!baseok)  > maxexact) { warning(paste("Would need permutation of", sum(!baseok), "numbers, resetting to greedy search\n"), call.=FALSE, immediate.=TRUE)
              method  <- 2
            } else     {
              iter    <- gamma(ncol(tab) - sum(baseok) + 1L)
              if(isTRUE(verbose))          message("Iterations for permutation matching:", iter, "\\n")
              perm    <- .permutations(ncol(tab) - sum(baseok))
            }
          }
          rest   <- if(any(baseok)) myseq[-retval[baseok]] else myseq
          for(l in 1L:iter)  {
            newretval <- retval
            if(method == 2)  {
              ok      <- baseok
              while(sum(!ok) > 1)  {
                rest  <- myseq[!ok]
                k     <- sample(rest, 1L)
                rmax  <- if(any(ok)) tab[k,-newretval[ok]] else tab[k,]
                ok[k] <- TRUE
                newretval[k]      <- as.numeric(names(rmax)[which.max(rmax)])
              }
              newretval[!ok]      <- myseq[-newretval[ok]]
            } else     {
              newretval[!baseok]  <- rest[perm[l,]]
            }
            if(l == 1) {
              retval  <- newretval
              agree   <- oldagree <- sum(diag(tab[,newretval]))/sum(tab)
            } else {
              agree   <- sum(diag(tab[,newretval]))/sum(tab)
              if(agree > oldagree) {
                retval            <- newretval
                oldagree          <- agree
              }
            }
          }
        }
      }
      if(isTRUE(verbose))                  message("Cases in matched pairs:", round(100 * sum(diag(tab[,retval]))/sum(tab), 2L), "%\\n")
      if(any(as.character(myseq) != cn)) {
        retval   <- cn[retval]
      }
      names(retval)   <- rn
        retval
    }

    .matnames    <- function(list, names, dim = 2L) {
        mapply(function(X, Y) { dimnames(X)[[dim]] <- Y; X }, list, names, SIMPLIFY=FALSE)
    }

    .ndeci       <- function(x, after.dot = TRUE) {
      scipen     <- options()$scipen
      digits     <- options()$digits
      options(scipen   = 999,
              digits   = 7)
      on.exit(options(scipen = scipen,
                      digits = digits))
      if(length(after.dot)  != 1   ||
         !is.logical(after.dot))           stop("'after.dot' must be a single logical indicator", call.=FALSE)
      if(!all(is.numeric(x)))              stop("'x' must be numeric", call.=FALSE)
      res        <- x
      na.ind     <- !is.na(x)
      x          <- abs(x[na.ind])
      if(all(icheck   <- (floor(x) == x), na.rm=TRUE)) {
        res[na.ind]   <- if(isTRUE(after.dot)) vector("integer", length(x)) else vapply(as.integer(x), format.info, integer(1L), digits=22)
      } else      {
        ipart    <- pmax(1L, floor(x))
        ichar    <- nchar(ipart)
        remain   <- x %% ipart
        res[na.ind]   <- (vapply(gsub("0+$", "", as.character(remain)), format.info, integer(1L), digits=22L) - !icheck) -
                      (if(isTRUE(after.dot)) pmin(1L, replace(ichar, remain == 0, 0L)) else ifelse(ichar > 1, -  ichar   + !icheck, 0L))
      }
        return(res)
    }

    .padding     <- function(x, nr) {
      length(x)  <- nr
        replace(x, is.na(x), 0L)
    }

    .permutations     <- function(n) {
      if(n == 1)       {                   return(matrix(1L))
      } else if(n < 2)                     stop("n must be a positive integer", call.=FALSE)
      z          <- matrix(1L)
      for(i in 2L:n)   {
       x         <- cbind(z, i)
       a         <- c(1L:i, 1L:(i - 1L))
       nr        <- nrow(x)
       z         <- matrix(0L, ncol=ncol(x), nrow=i * nr)
       z[1L:nr,] <- x
       for(j in (2L:i) - 1L) {
         z[j * nr + 1L:nr,] <- x[,a[(1L:i) + j]]
       }
      }
      dimnames(z)     <- NULL
        z
    }

    .plot_CI     <- function(x, y = NULL, uiw, liw = uiw, ui = NULL, li = NULL, err = "y", sfrac = 0.01,
                             gap = 0, slty = graphics::par("lty"), add = FALSE, scol = NULL, pt.bg = graphics::par("bg"), ...) {
      arglist    <- list(...)
      if(inherits(x, "list"))   {
        y        <- x$y
        x        <- x$x
      }
      if(is.null(y)) {
        if(is.null(x))                     stop("Both x and y are NULL", call.=FALSE)
        y        <- as.numeric(x)
        x        <- seq_along(x)
      }
      if(missing(uiw) &&
        (is.null(ui)  || is.null(li)))     stop("Must specify either relative limits or both lower and upper limits", call.=FALSE)
      if(!missing(uiw)) {
        z        <- if(err == "y") y else x
        ui       <- z  + uiw
        li       <- z  - liw
      }
      if(is.null(arglist$xlab)) arglist$xlab <- deparse(substitute(x))
      if(is.null(arglist$ylab)) arglist$ylab <- deparse(substitute(y))
      if(err == "y"   &&
         is.null(arglist$ylim)) arglist$ylim <- range(c(y, ui, li), na.rm = TRUE)
      if(err == "x"   &&
         is.null(arglist$xlim)) arglist$xlim <- range(c(x, ui, li), na.rm = TRUE)
      if(missing(scol)) {
        scol     <- if(!is.null(arglist$col)) arglist$col else graphics::par("col")
      }
      ppoints    <- TRUE
      if(!is.null(arglist$pch) && is.na(arglist$pch)) {
        arglist$pch       <- 1L
        ppoints           <- FALSE
      }
      if(isFALSE(add))       do.call(base::plot, c(list(x, y, type = "n"), .clean_args(arglist, base::plot)))
      if(isTRUE(gap)) gap <- 0.01
      ul         <- c(li, ui)
      pin        <- graphics::par("pin")
      usr        <- graphics::par("usr")
      x.to.in    <- pin[1L]/diff(usr[1L:2L])
      y.to.in    <- pin[2L]/diff(usr[3L:4L])
      if(err == "y") {
        gap      <- rep(gap, length(x)) * diff(graphics::par("usr")[3L:4L])
        smidge   <- graphics::par("fin")[1L] * sfrac
        nz       <- abs(li - pmax(y - gap, li)) * y.to.in > 0.001
        scols    <- rep(scol, length.out = length(x))[nz]
        arrow.args        <- c(list(lty = slty, angle = 90, length = smidge, code = 1, col = scols),
                               .clean_args(arglist, graphics::arrows, exclude.other = c("col", "lty", "axes")))
        do.call(graphics::arrows, c(list(x[nz], li[nz], x[nz], pmax(y - gap, li)[nz]), arrow.args))
        nz       <- abs(ui - pmin(y + gap, ui)) * y.to.in > 0.001
        scols    <- rep(scol, length.out = length(x))[nz]
        arrow.args$col    <- scols
        do.call(graphics::arrows, c(list(x[nz], ui[nz], x[nz], pmin(y + gap, ui)[nz]), arrow.args))
      } else if(err == "x") {
        gap      <- rep(gap, length(x)) * diff(graphics::par("usr")[seq_len(2L)])
        smidge   <- graphics::par("fin")[2L] * sfrac
        nz       <- abs(li - pmax(x - gap, li)) * x.to.in > 0.001
        scols    <- rep(scol, length.out = length(x))[nz]
        arrow.args        <- c(list(lty = slty, angle = 90, length = smidge, code = 1, col = scols),
                               .clean_args(arglist, graphics::arrows, exclude.other = c("col", "lty", "axes")))
        do.call(graphics::arrows, c(list(li[nz], y[nz], pmax(x - gap, li)[nz], y[nz]), arrow.args))
        nz       <- abs(ui - pmin(x + gap, ui)) * x.to.in > 0.001
        scols    <- rep(scol, length.out = length(x))[nz]
        arrow.args$col    <- scols
        do.call(graphics::arrows, c(list(ui[nz], y[nz], pmin(x + gap, ui)[nz], y[nz]), arrow.args))
      }
      if(ppoints) do.call(graphics::points, c(list(x, y, bg = pt.bg), .clean_args(arglist, graphics::points, exclude.other = c("xlab", "ylab", "xlim", "ylim", "axes"))))
        invisible(list(x = x, y = y))
    }

    .rgamma0     <- function(...) {
      tmp        <- stats::rgamma(...)
      tmp[tmp    == 0]     <- .Machine$double.eps
        tmp
    }

    #' @importFrom matrixStats "rowSums2"
    .row_vars    <- function(x, std = FALSE, suma = NULL) { # formerly replaced Rfast::rowVars
      if(length(std) > 1 ||
         !is.logical(std))                 stop("'std' must be a single logical indicator")
      if(!is.matrix(x))                    stop("'x' must be a matrix")
      m          <- if(missing(suma)) rowSums2(x) else suma
      n          <- ncol(x)
      s          <- (rowSums2(x^2) - m^2/n)/(n - 1L)
        if(std) sqrt(s) else s
    }

    #' @importFrom matrixStats "colMeans2" "rowSums2"
    .scale2      <- function(x, center = TRUE, scale = TRUE) { # replaces Rfast::standardise
      cmeans     <- if(isTRUE(center)) colMeans2(x)     else center
      center     <- if(is.logical(center))   center     else is.numeric(center)
      scaling    <- if(is.logical(scale))     scale     else is.numeric(scale)
      if(center  && scaling) {
        y        <- t(x) - cmeans
          if(isTRUE(scale)) t(y/sqrt(rowSums2(y^2)) * sqrt(nrow(x) - 1L)) else t(y/scale)
      } else if(center)      {
          t(t(x)  - cmeans)
      } else if(scaling)     {
          t(t(x)/if(isTRUE(scale)) colVars(x, std=TRUE) else scale)
      } else  x
    }

    #' @importFrom Rfast "colVars"
    .tune_beta0  <- function(dat, beta0 = 3, type = c("diag", "mse")) { # unused
      dat        <- as.matrix(dat)
      N          <- nrow(dat)
      P          <- ncol(dat)
      inv.cov    <- (beta0 + N/2L) * chol2inv(chol(diag(beta0, P) + 0.5 * crossprod(.scale2(dat)))) * tcrossprod(1/colVars(dat, std=TRUE))
      error      <- diag(P) - (stats::cov(dat) %*% inv.cov)
        switch(EXPR=match.arg(type), diag=sum(diag(error)^2)/P, mean(error^2))
    }

    #' @importFrom matrixStats "colSums2" "rowSums2"
    .vari_max    <- function(x, normalize = FALSE, eps = 1e-05, ...) {
      p          <- ncol(x)
      if(p        < 2)                     return(x)
      if(normalize)         {
        sc       <- sqrt(rowSums2(x^2))
        x        <- x/sc
      }
      n          <- nrow(x)
      TT         <- diag(p)
      d          <- 0
      ns         <- rep(1L, n)
      for(i in seq_len(1000L)) {
        z        <- x    %*% TT
        sB       <- La.svd(crossprod(x, z^3  - z %*% diag(colSums2(z^2)/n)))
        TT       <- sB$u %*% sB$vt
        dp       <- d
        d        <- sum(sB$d)
        if(d      < dp * (1 + eps))        break
      }
      z          <- if(normalize) (x %*% TT) * sc else x %*% TT
      dimnames(z)         <- dimnames(x)
      class(z)   <- "loadings"
        return(list(loadings = z, rotmat = TT))
    }

    .version_above        <- function(pkg, than) {
      pkg        <- as.character(utils::packageVersion(pkg))
        identical(pkg, than) || (utils::compareVersion(pkg, than) >= 0)
    }

    .which0      <- function(x) which(x == 0)
    #
