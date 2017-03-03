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
        mu.eta    + t(backsolve(u.eta, matrix(rnorm(Q * N), nrow=Q, ncol=N)))
    }

  # Loadings
    .sim_load    <- function(l.sigma, Q, c.data, eta, psi.inv, EtE, Q1)  {
      u.load     <- l.sigma + psi.inv * EtE
      u.load     <- if(Q1) sqrt(u.load) else .chol(u.load)
      mu.load    <- psi.inv * (if(Q1) 1/(u.load * u.load) else chol2inv(u.load)) %*% crossprod(eta, c.data)
        mu.load   + backsolve(u.load, rnorm(Q))
    }

    .sim_load_s  <- function(Q, c.data, eta, phi, tau, psi.inv, EtE, Q1) {
      u.load     <- diag(phi * tau, Q) + psi.inv * EtE
      u.load     <- if(Q1) sqrt(u.load) else .chol(u.load)
      mu.load    <- psi.inv  * (if(Q1) 1/(u.load * u.load) else chol2inv(u.load)) %*% crossprod(eta, c.data)
        mu.load   + backsolve(u.load, rnorm(Q))
    }

  # Uniquenesses
    .sim_psi_iu  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat) {
      S.mat      <- c.data - tcrossprod(eta, lmat)
        rgamma(P, shape=N/2 + psi.alpha, rate=colSums(S.mat * S.mat)/2 + psi.beta)
    }

    .sim_psi_ii  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat) {
      S.mat      <- c.data - tcrossprod(eta, lmat)
        rep(rgamma(1, shape=(N * P)/2 + psi.alpha, rate=sum(S.mat * S.mat)/2 + psi.beta), P)
    }

  # Local Shrinkage
    .sim_phi     <- function(Q, P, nu, tau, load.2, plus1) {
        matrix(rgamma(P * Q, shape=0.5 + nu + plus1, rate=(nu + sweep(load.2, 2, tau, FUN="*"))/2), nrow=P, ncol=Q)
    }

  # Global Shrinkage
    .sim_delta1  <- function(Q, P, alpha.d1, delta.1, beta.d1, tau, sum.term) {
        rgamma(1, shape=alpha.d1 + P * Q/2, rate=beta.d1 + 0.5/delta.1 * tau %*% sum.term)
    }

    .sim_deltak  <- function(Q, P, k, alpha.d2, beta.d2, delta.k, tau.kq, sum.term.kq) {
        rgamma(1, shape=alpha.d2 + P/2 * (Q - k + 1), rate=beta.d2 + 0.5/delta.k * tau.kq %*% sum.term.kq)
    }

  # Mixing Proportions
    .sim_pi      <- function(pi.alpha, nn, G) {
      tmp        <- rgamma(G, shape=pi.alpha + nn, rate=1)
        tmp/sum(tmp)
    }

    .sim_pi_inf  <- function(pi.alpha, nn, N = sum(nn), len, lseq, discount = 0L) {
      vs         <- if(discount == 0) rbeta(len, 1 + nn, pi.alpha + N - cumsum(nn)) else rbeta(len, 1 - discount + nn, pi.alpha + lseq * discount + N - cumsum(nn))
        return(list(Vs = vs, pi.prop = vapply(lseq, function(t) vs[t] * prod(1 - vs[seq_len(t - 1)]), numeric(1L))))
    }

  # Cluster Labels
#' Simulate Cluster Labels from Unnormalised Log-Probabilities using the Gumbel-Max Trick
#'
#' Samples cluster labels for N observations from G groups efficiently using log-probabilities and the so-called Gumbel-Max trick, without requiring that the log-probabilities need to be normalised; thus redunant computation can be avoided. Computation takes place on the log scale for stability/underflow reasons (to ensure negligible probabilities won't cause computational difficulties); in any case, many functions for calculating multivariate normal densities already output on the log scale. Please note that while the function is available for standalone use that no checks take place, in order to speed up repeated calls to the function inside \code{\link{mcmc_IMIFA}}.
#' @param probs An N x G matrix of unnormalised probabilities on the log scale.
#' @param N The number of observations that require labels to be sampled. Must be at least 1. Defaults to (and should not be changed from!) \code{nrow(log.probs)}. Only necessary if \code{slice} is \code{FALSE}.
#' @param G The number of 'active' clusters, s.t. sampled labels can take values in \code{1:G}. Must be at least 1. Defaults to (and should not be changed from!) \code{ncol(log.probs)}. Only necessary if \code{slice} is \code{FALSE}.
#' @param slice A logical indicating whether or not the indicator correction for slice sampling has been applied to \code{log.probs}. Defaults to \code{FALSE} but is \code{TRUE} for the "\code{IMIFA}" and "\code{IMFA}" methods under \code{\link{mcmc_IMIFA}}. Details of this correction are given in Murphy et. al. (2017).
#' @param log.like A logical indicating whether the normalising constant is to be computed. Defaults to \code{FALSE} but is \code{TRUE} for all methods under \code{\link{mcmc_IMIFA}} where it's necessary for computation of the log-likelihoods required for model choice.
#' @return Either a N-vector of sampled cluster labels, or if \code{isTRUE(log.like)}, a list with two elements:
#' \describe{
#' \item{z}{The numeric vector of \code{N} sampled cluster labels, with the largest label no greater than \code{G}.}
#' \item{log.like}{The log-likelihood(s), given by the normalising constant(s), computed with the aid of \code{\link[matrixStats]{rowLogSumExps}}.}
#' }
#' @seealso \code{\link{mcmc_IMIFA}}, \code{\link[matrixStats]{rowLogSumExps}}
#' @references Murphy, K., Gormley, I.C. and Viroli, C. (2017) Infinite Mixtures of Infinite Factor Analysers: Nonparametric Model-Based Clustering via Latent Gaussian Models, \code{https://arxiv.org/abs/1701.07010}
#'
#' Yellot, J. I. Jr. (1977) The relationship between Luce's choice axiom, Thurstone's theory of comparative judgment, and the double exponential distribution. \emph{Journal of Mathematical Psychology}, 15: 109-144.
#' @export
#'
#' @examples
#' # Set the dimensions & simulate a matrix of weights
#'   N         <- 1
#'   G         <- 3
#'   weights   <- matrix(c(1, 2, 3), nrow=N, ncol=G)
#'
#' # Call gumbel_max() repeatedly to obtain samples of the label zs
#'   iters     <- 10000
#'   zs        <- vapply(seq_len(iters), function(i)
#'                gumbel_max(probs=log(weights), N=N, G=G), numeric(1L))
#'
#' # Compare answer to the normalised weights
#'   tabulate(zs, nbins=G)/iters
#'   normalised <- as.numeric(weights/sum(weights))
#'
#' # Simulate a matrix of dirichlet weights & the associated vector of N labels
#'   N       <- 400
#'   G       <- 10
#'   tmp     <- matrix(rgamma(N * G, shape=1), nrow=N, ncol=G)
#'   weights <- sweep(tmp, 1, rowSums(tmp), FUN="/")
#'   zs      <- gumbel_max(probs=log(weights), N=N, G=G)
    gumbel_max   <- function(probs, N, G, log.like = FALSE, slice = FALSE) {
      if(isTRUE(slice))    {
        fp       <- is.finite(probs)
        zs       <- max.col(replace(probs, fp, probs[fp] - log(rexp(sum(fp)))))
      } else {
        if(missing(N)) N <- nrow(probs)
        if(missing(G)) G <- ncol(probs)
        zs       <- max.col(probs - log(matrix(rexp(N * G), nrow=N, ncol=G)))
      }
        return(if(isTRUE(log.like)) list(z = zs, log.like=rowLogSumExps(probs)) else zs)
    }

  # Alpha
    .sim_alpha_g <- function(alpha, shape, rate, G, N, discount = 0L) {
      shape2     <- shape + G - 1
      rate2      <- rate  - log(rbeta(1, alpha + discount + 1, N))
      weight     <- shape2/(shape2 + N * rate2)
        weight    * rgamma(1, shape=shape2 + 1, rate=rate2) + (1 - weight) * rgamma(1, shape=shape2, rate=rate2) - discount
    }

    .sim_alpha_m <- function(alpha, lower, upper, trunc.G, Vs, discount = 0L) {
      alpha.old  <- alpha   + discount
      alpha.new  <- runif(1, lower, upper) + discount
      a.prob     <- trunc.G * (log(alpha.new) - log(alpha.old))    + (alpha.new - alpha.old) * sum(log((1 - Vs[-trunc.G])))
      accept     <- a.prob >= 0 || - rexp(1)  < a.prob
        return(list(alpha   = ifelse(accept, alpha.new, alpha.old) - discount, rate = accept))
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
        sqrt(sigma.l) * matrix(rnorm(P * Q), nrow=P, ncol=Q)
    }

    .sim_load_ps <- function(Q, sigma.l, phi, tau) {
        sqrt(1/(phi * tau)) * rnorm(Q)
    }

  # Uniquenesses
    .sim_psi_ipu <- function(P, psi.alpha, psi.beta) {
        rgamma(n=P, shape=psi.alpha, rate=psi.beta)
    }

    .sim_psi_ipi <- function(P, psi.alpha, psi.beta) {
        rep(rgamma(1, shape=psi.alpha, rate=psi.beta), P)
    }

  # Local Shrinkage
    .sim_phi_p   <- function(Q, P, nu, plus1) {
        matrix(rgamma(n=P * Q, shape=nu + plus1, rate=nu), nrow=P, ncol=Q)
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
#' Takes a shape hyperparameter and covariance matrix, and finds data-driven rate hyperparameters in such a way that Heywood problems are avoided. Rates are allowed to be variable-specific. Used internally by \code{\link{mcmc_IMIFA}} when it's argument \code{psi_beta} is not supplied.
#' @param shape A positive shape hyperparameter.
#' @param covar A square, positive-semidefinite covariance matrix
#' @param type A switch indicating whether a single rate (\code{isotropic}) or variable-specific rates (\code{unconstrained}) are to be derived. uniquenesses are allowed to be variable specific
#'
#' @return Either a single rate hyperparameter or \code{ncol(covar)} variable specific hyperparameters.
#' @export
#' @importFrom matrixcalc "is.positive.semi.definite"
#'
#' @seealso \code{\link{mcmc_IMIFA}}
#' @references Fruwirth-Schnatter, S.  and Lopes, H.F. (2010). Parsimonious Bayesian factor analysis when the number of factors is unknown. \emph{Technical Report}. The University of Chicago Booth School of Business.
#'
#' @examples
#' data(olive)
#' olive2 <- olive[,-(1:2)]
#' rates  <- psi_hyper(shape=2.5, covar=cov(olive2), type="isotropic")
#' rates
#'
#' olive_scaled <- scale(olive2, center=TRUE, scale=TRUE)
#' rate   <- psi_hyper(shape=3, covar=cov(olive_scaled), type="unconstrained")
#' rate
    psi_hyper   <- function(shape, covar, type=c("unconstrained", "isotropic")) {
      if(!all(is.positive.semi.definite(covar),
              is.symmetric(covar),
              is.double(covar)))           stop("Invalid covariance matrix supplied")
      if(any(!is.numeric(shape),
             length(shape) != 1))          stop("'shape' must be a single digit")
      P         <- ncol(covar)
      inv.cov   <- try(base::solve(covar), silent=TRUE)
      if(inherits(inv.cov, "try-error"))  {
        inv.cov <- 1/covar
      }
        unname((shape - 1)/switch(match.arg(type), unconstrained=diag(inv.cov), isotropic=rep(sum(diag(inv.cov))/P, P)))
    }

  # Alpha/Discount Shifted Gamma Hyperparameters
    .shift_GA   <- function(shape, rate, shift = 0L, param = c("rate", "scale")) {
      var       <- shape/rate^2
      exp       <- var  * rate + shift
      rate      <- exp/var
      shape     <- rate * exp
        return(list(shape = shape, rate = switch(match.arg(param), rate=rate, 1/rate)))
    }

  # Check Shrinkage Hyperparemeters
#' Check the validity of MGP hyperparameters
#'
#' Checks the hyperparameters for the multiplicative gamma process (MGP) shrinkage prior in order to ensure that the property of cumulative shrinkage holds. This is called inside \code{\link{mcmc_IMIFA}} for the "\code{IFA}", "\code{MIFA}", "\code{OMIFA}" and "\code{IMIFA}" methods. The arguments \code{ad1, ad2, nu, bd1} and \code{bd2} are vectorised.
#' @param ad1 Shape hyperparameter for delta_1
#' @param ad2 Shape hyperparameter for delta_2
#' @param Q Number of latent factors
#' @param nu Hyperparameter for the local shrinkage parameters
#' @param bd1 Rate hyperparameter for delta_1. Defaults to 1.
#' @param bd2 Rate hyperparameter for delta_2. Defaults to 1.
#' @param plus1 Logical indicator for whether the Gamma prior on the local shrinkage parameters is of the form Ga(\code{nu + 1, nu}), the default, or Ga(\code{nu, nu}).
#' @param inverse Logical indicator for whether the increasing shrinkage property is assessed against the induced Inverse Gamma prior, the default. Always TRUE when used inside \code{\link{mcmc_IMIFA}}: the FALSE option exists only for demonstration purposes.
#'
#' @return A list of length 2 containing the following objects:
#' \describe{
#'   \item{expectation}{The vector of actual expected shrinkage factors}
#'   \item{valid}{A logical indicating whether the cumulative shrinkage property holds}
#' }
#' @export
#' @seealso \code{\link{mcmc_IMIFA}}
#' @references
#' Bhattacharya, A. and Dunson, D. B. (2011). Sparse Bayesian infinite factor models. \emph{Biometrika}, 98(2): 291–306.
#'
#' Durante, D. (2017). A note on the multiplicative gamma process. \emph{Statistics & Probability Letters}, 122: 198–204.
#'
#' @examples
#' # Check if the MGP global shrinkage parameters are increasing as the column index increases.
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, nu=2, inverse=FALSE)[[1]]$valid
#'
#' # Check if the induced IG prior on the MGP global shrinkage parameters is stochastically increasing.
#' MGP_check(ad1=1.5, ad2=1.8, Q=10, nu=2, inverse=TRUE)[[1]]$valid
    MGP_check   <- Vectorize(function(ad1, ad2, Q, nu, bd1 = 1L, bd2 = 1L, plus1 = TRUE, inverse = TRUE) {
      if(any(!is.logical(plus1),
             length(plus1)    != 1))       stop("'plus1' must be TRUE or FALSE")
      if(any(!is.logical(inverse),
             length(inverse)  != 1))       stop("'inverse' must be TRUE or FALSE")
      if(any(length(Q) > 1, Q  < 2))       stop("Q must be single value, greater than or equal to 2")
      if(any(nu <= !plus1,
             !is.numeric(nu)))             stop(paste0("'nu' must be a single ", ifelse(plus1,
                                                "strictly positive number for the Ga(nu + 1, nu) parameterisation",
                                                "number strictly greater than 1 for the Ga(nu, nu) parameterisation")))
      if(any(c(ad1, ad2)  < 1))            stop("All shrinkage shape hyperparameter values must be at least 1")
      if(any(c(bd1, bd2) <= 0))            stop("All shrinkage rate hyperparameter values must be strictly positive")
      if(any(ad1 < bd1, ad2 < bd2))        warning("Shrinkage shape hyperparameters should be greater than associated shrinkage rate hyperparameters", call.=FALSE)
      if(bd2/(ad2 - 1)   >= bd1/(ad1 - 1)) warning("Shrinkage in column k should be greater than shrinkage in column 1", call.=FALSE)
      rate      <- nu
      shape     <- ifelse(plus1, rate + 1, rate)
      if(inverse) {
        ad1     <- ifelse(ad1 == 1, ad1 + .Machine$double.eps, ad1)
        ad2     <- ifelse(ad2 == 1, ad2 + .Machine$double.eps, ad2)
        exp.seq <- rate/(shape - 1) * bd1/(ad1 - 1) * (bd2/(ad2 - 1))^(seq_len(Q) - 1)
        check   <- is.unsorted(exp.seq)
      } else {
        exp.seq <- shape/rate * ad1/bd1 * (ad2/bd2)^(seq_len(Q) - 1)
        check   <- !is.unsorted(exp.seq)
      }
        return(list(expectation = exp.seq, valid = check))
    }, vectorize.args = c("ad1", "ad2", "nu", "bd1", "bd2"), SIMPLIFY=FALSE)

  # Label Switching
    .lab_switch <- function(z.new, z.old, Gs, ng = tabulate(z.new)) {
      tab       <- table(z.new, z.old, dnn=NULL)
      tab.tmp   <- tab[rowsums(tab) != 0,colsums(tab) != 0, drop=FALSE]
      nc        <- ncol(tab.tmp)
      nr        <- nrow(tab.tmp)
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
        z.perm  <- suppressWarnings(matchClasses(tab.tmp, method="exact", verbose=FALSE))
        z.perm  <- setNames(as.numeric(z.perm), names(z.perm))
      }
      if(length(Gs) > length(z.perm)) {
        z.perm  <- c(z.perm, setNames(setdiff(Gs, z.perm), setdiff(Gs, names(z.perm))))
      }
      z.names   <- as.numeric(names(z.perm))
      z.perm    <- z.perm[Order(z.names)]
      z.sw      <- factor(z.new, labels=z.perm[seq_along(ng[ng > 0])])
        return(list(z = as.numeric(levels(z.sw))[z.sw], z.perm = z.perm))
    }

    # Move 1
    .lab_move1  <- function(nn.ind, pi.prop, nn) {
      sw        <- sample(nn.ind, 2L)
      pis       <- pi.prop[sw]
      nns       <- nn[sw]
      a.prob    <- (nns[1] - nns[2]) * (log(pis[1])  - log(pis[2]))
        return(list(rate1  = a.prob >= 0 || - exp(1) < a.prob, sw = sw))
    }

    # Move 2
    .lab_move2  <- function(G, Vs, nn) {
      sw        <- sample(G, 1L, prob=c(rep(1, G - 2), 0.5, 0.5))
      sw        <- if(is.element(sw, c(G, G - 1))) c(G - 1, G) else c(sw, sw + 1)
      nns       <- nn[sw]
      Vsw       <- Vs[sw]
      a.prob    <- nns[1] * log(1 - Vsw[2]) - nns[2]  * log(1 - Vsw[1])
        return(list(rate2 = a.prob >= 0  || - rexp(1) < a.prob, sw = sw))
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
#' Calculates the expected number of clusters under a Dirichlet process or Pitman-Yor process prior for a sample of size \code{N} at given values of the concentration parameter \code{alpha} and optionally also the \code{discount} parameter. Useful for soliciting sensible priors for \code{alpha} or suitable fixed values for \code{alpha} or \code{discount} under the "\code{IMFA}" and "\code{IMIFA}" methods for \code{\link{mcmc_IMIFA}}, All arguments are vectorised. Requires use of the \code{Rmpfr} and \code{gmp} libraries for non-zero \code{discount} values.
#' @param N The sample size.
#' @param alpha The concentration parameter. Must be specified and must be strictly greater than \code{-discount}.
#' @param discount The discount parameter for the Pitman-Yor process. Must lie in the interval [0, 1). Defaults to 0 (i.e. the Dirichlet process).
#'
#' @return The expected number of clusters under the specified prior conditions.
#' @export
#' @seealso \code{\link{G_variance}}, \code{\link{G_prior}}, \code{\link[Rmpfr]{Rmpfr}}
#'
#' @examples
#' G_expected(N=150, alpha=c(3, 10, 25))
    G_expected  <- Vectorize(function(N, alpha, discount = 0L) {
      if(!all(is.numeric(N), is.numeric(discount),
         is.numeric(alpha)))               stop("All inputs must be numeric")
      if(discount  < 0  || discount >= 1)  stop("Invalid discount value")
      if(alpha   < - discount)             stop("Invalid alpha value")
      if(discount == 0) {
          alpha * (digamma(alpha + N) - digamma(alpha))
      } else {
        if(suppressMessages(requireNamespace("Rmpfr", quietly=TRUE))) {
          on.exit(.detach_pkg("Rmpfr"))
          on.exit(.detach_pkg("gmp"), add=TRUE)
        } else                             stop("'Rmpfr' package not installed")
          gmp::asNumeric(alpha/discount * Rmpfr::pochMpfr(alpha + discount, N)/Rmpfr::pochMpfr(alpha, N) - alpha/discount)
      }
    })

#' 2nd Moment of Dirichlet / Pitman-Yor processes
#'
#' Calculates the variance in the number of clusters under a Dirichlet process or Pitman-Yor process prior for a sample of size \code{N} at given values of the concentration parameter \code{alpha} and optionally also the \code{discount} parameter. Useful for soliciting sensible priors for \code{alpha} or suitable fixed values for \code{alpha} or \code{discount} under the "\code{IMFA}" and "\code{IMIFA}" methods for \code{\link{mcmc_IMIFA}}, All arguments are vectorised. Requires use of the \code{Rmpfr} and \code{gmp} libraries for non-zero \code{discount} values.
#' @inheritParams G_expected
#' @return The variance of the number of clusters under the specified prior conditions.
#' @export
#' @seealso \code{\link{G_expected}}, \code{\link{G_prior}}, \code{\link[Rmpfr]{Rmpfr}}
#'
#' @examples
#' G_variance(N=150, alpha=c(3, 10, 25))
    G_variance  <- Vectorize(function(N, alpha, discount = 0L) {
      if(!all(is.numeric(N), is.numeric(discount),
         is.numeric(alpha)))               stop("All inputs must be numeric")
      if(discount  < 0  || discount >= 1)  stop("Invalid discount value")
      if(alpha   < - discount)             stop("Invalid alpha value")
      if(discount == 0) {
          alpha * (digamma(alpha + N) - digamma(alpha)) + (alpha^2) * (trigamma(alpha + N) - trigamma(alpha))
      } else {
        if(suppressMessages(requireNamespace("Rmpfr", quietly=TRUE))) {
          on.exit(.detach_pkg(Rmpfr))
          on.exit(.detach_pkg(gmp), add=TRUE)
        } else                             stop("'Rmpfr' package not installed")
        sum.ad  <- alpha + discount
        poch.a  <- Rmpfr::pochMpfr(alpha, N)
        poch.ad <- Rmpfr::pochMpfr(sum.ad, N)
        subterm <- alpha/discount * poch.ad/poch.a
          gmp::asNumeric((alpha * sum.ad)/discount^2 * Rmpfr::pochMpfr(sum.ad + discount, N)/poch.a - subterm - subterm^2)
      }
    })

  # Detach packages
    .detach_pkg <- function(pkg, character.only = FALSE) {
      if(!character.only) {
        pkg     <- deparse(substitute(pkg))
      }
      searches  <- paste("package", pkg, sep=":")
      while(searches %in% search()) {
        detach(searches, unload=TRUE, character.only=TRUE)
      }
    }

  # Print functions
#' @method print IMIFA
#' @export
    print.IMIFA <- function(x, ...) {
      meth      <- attr(x, "Method")
      name      <- attr(x, "Name")
      fac       <- attr(x, "Factors")
      grp       <- attr(x, "Groups")
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
      fac       <- attr(object, "Factors")
      grp       <- attr(object, "Groups")
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

#' @method print Results_IMIFA
#' @export
    print.Results_IMIFA  <- function(x, ...) {
      method    <- attr(x, "Method")
      G         <- x$GQ.results$G
      Q         <- x$GQ.results$Q
      if(is.element(method, c("FA", "IFA")))  {
        msg     <- paste0("The chosen ", method, " model has ", Q, " factor", ifelse(Q == 1, "\n", "s\n"))
      } else if(is.element(method, c("MFA", "OMFA", "IMFA"))) {
        msg     <- paste0("The chosen ", method, " model has ", G, " group",  ifelse(G == 1, " with ", "s, each with "), unique(Q), " factor", ifelse(unique(Q) == 1, "\n", "s\n"))
      } else {
        Q.msg   <- NULL
        for(i in seq_along(Q[-length(Q)])) {
          Q.msg <- c(Q.msg, (paste0(Q[i], ifelse(i + 1 < length(Q), ", ", " "))))
        }
        Q.msg   <- if(length(Q) > 1) paste(c(Q.msg, paste0("and ", Q[length(Q)])), sep="", collapse="") else Q
        msg     <- paste0("The chosen ", method, " model has ", G, " group",  ifelse(G == 1, " with ", "s, with "), Q.msg, " factor", ifelse(G == 1 && Q == 1, "\n", paste0("s", ifelse(G == 1, "\n", " respectively\n"))), sep="")
      }
        cat(msg)
    }

#' @method summary Results_IMIFA
#' @export
    summary.Results_IMIFA <- function(object, ...) {
      criterion <- unlist(strsplit(toupper(attr(object$GQ.results, "Criterion")), "[.]"))
      criterion <- ifelse(length(criterion) > 1, ifelse(criterion[1] != "LOG", paste0(criterion[1], ".", tolower(criterion[2])), "LogIntegratedLikelihood"), criterion)
      crit.mat  <- object$GQ.results[[paste0(criterion, "s")]]
      msg       <- NULL
      if(any(dim(crit.mat) > 1)) {
        msg     <- paste0(", and ", ifelse(substr(criterion, 1, 1) == "A", "an ", "a "),  criterion, " of ", round(max(crit.mat), 2), "\n")
      }
        cat(paste0(capture.output(print.Results_IMIFA(object)), msg))
    }

    .power2     <- function(x) x * x
    .which0     <- function(x) which(x == 0)
    .chol       <- function(x) tryCatch(chol(x), error=function(e) chol(make.positive.definite(x)))
    .ledermann  <- function(N, P) as.integer(min(N - 1, floor((2 * P + 1 - sqrt(8 * P + 1))/2)))
    .dim        <- Vectorize(function(Q, P) as.integer(P * Q - 0.5 * Q * (Q - 1) + 2 * P), vectorize.args = "Q", SIMPLIFY=TRUE)
    #
