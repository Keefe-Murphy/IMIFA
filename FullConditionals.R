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
        mu.eta    + t(backsolve(u.eta, matrix(rnorm(Q * N), nr=Q, nc=N)))
    }
      
  # Loadings
    .sim_load    <- function(l.sigma, Q, c.data, eta, psi.inv, EtE, Q1)  {
      u.load     <- l.sigma + psi.inv * EtE
      u.load     <- if(Q1) sqrt(u.load) else .chol(u.load)
      mu.load    <- psi.inv * (if(Q1) 1/(u.load * u.load) else chol2inv(u.load)) %*% crossprod(eta, c.data)
        mu.load   + backsolve(u.load, rnorm(Q))
    }
    
    .sim_load.s  <- function(Q, c.data, eta, phi, tau, psi.inv, EtE, Q1) {
      u.load     <- diag(phi * tau, Q) + psi.inv * EtE
      u.load     <- if(Q1) sqrt(u.load) else .chol(u.load)
      mu.load    <- psi.inv  * (if(Q1) 1/(u.load * u.load) else chol2inv(u.load)) %*% crossprod(eta, c.data)
        mu.load   + backsolve(u.load, rnorm(Q))
    }
    
  # Uniquenesses
    .sim_psi.iu  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat) { 
      S.mat      <- c.data    - tcrossprod(eta, lmat)
        rgamma(P,   shape=N/2 + psi.alpha, rate=colSums(S.mat * S.mat)/2 + psi.beta) 
    }
    
    .sim_psi.ii  <- function(N, P, psi.alpha, psi.beta, c.data, eta, lmat) { 
      S.mat      <- c.data - tcrossprod(eta, lmat)
        rep(rgamma(1, shape=(N * P)/2 + psi.alpha, rate=sum(S.mat * S.mat)/2 + psi.beta), P)
    }

  # Local Shrinkage
    .sim_phi     <- function(Q, P, nu, tau, load.2, plus1) {
        matrix(rgamma(P * Q, shape=1/2 + nu + plus1, rate=(nu + sweep(load.2, 2, tau, FUN="*"))/2), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    .sim_delta1  <- function(Q, P, alpha.d1, delta.1, beta.d1, tau, sum.term) {
        rgamma(1,   shape=alpha.d1 + P * Q/2, rate=beta.d1 + 0.5/delta.1 * tau %*% sum.term)
    }
    
    .sim_deltak  <- function(Q, P, k, alpha.d2, beta.d2, delta.k, tau.kq, sum.term.kq) {
        rgamma(1,   shape=alpha.d2 + P/2 * (Q - k + 1), rate=beta.d2 + 0.5/delta.k * tau.kq %*% sum.term.kq)
    }

  # Mixing Proportions
    .sim_pi      <- function(pi.alpha, nn) {
        rdirichlet(1, pi.alpha + nn)
    }
    
    .sim_pi.inf  <- function(pi.alpha, nn, N = sum(nn), len, lseq, discount = 0L) {
      vs         <- if(discount == 0) rbeta(len, 1 + nn, pi.alpha + N - cumsum(nn)) else rbeta(len, 1 - discount + nn, pi.alpha + lseq * discount + N - cumsum(nn))
        return(list(Vs = vs, pi.prop = vapply(lseq, function(t) vs[t] * prod(1 - vs[seq_len(t - 1)]), numeric(1))))
    }
  
  # Cluster Labels
    .sim_z       <- function(data, mu, sigma, Gseq, N, pi.prop, Q0) {
      log.num    <- vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
      log.denom  <- rowLogSumExps(log.num)
      lnp        <- sweep(log.num, 1, log.denom, FUN="-")
      for(g in Gseq[-1]) {
        lnp[,g]  <- rowLogSumExps(lnp[,g:(g - 1)])
      }
        return(list(z = rowsums(-rexp(N) > lnp) + 1, log.likes = log.denom))
    }
    
    .sim_z.inf   <- function(data, mu, sigma, Gseq, N, pi.prop, log.slice.ind, Q0) {
      log.num    <- vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N)) + log.slice.ind
      log.denom  <- rowLogSumExps(log.num)
      lnp        <- sweep(log.num, 1, log.denom, FUN="-")
      for(g in Gseq[-1]) {
        lnp[,g]  <- rowLogSumExps(lnp[,g:(g - 1)])
      }
        return(list(z = rowsums(-rexp(N) > lnp) + 1, log.likes = log.denom))
    }

  # Alpha
    .sim_alpha.g <- function(alpha, shape, rate, G, N, discount = 0L) {
      shape2     <- shape + G - 1
      rate2      <- rate  - log(rbeta(1, alpha + discount + 1, N))
      weight     <- shape2/(shape2 + N * rate2)
        weight    * rgamma(1, shape=shape2 + 1, rate=rate2) + (1 - weight) * rgamma(1, shape=shape2, rate=rate2) - discount
    }
    
    .sim_alpha.m <- function(alpha, lower, upper, trunc.G, Vs, discount = 0L) {
      alpha.old  <- alpha   + discount
      alpha.new  <- runif(1, lower, upper) + discount
      a.prob     <- trunc.G * (log(alpha.new) - log(alpha.old))    + (alpha.new - alpha.old) * sum(log((1 - Vs[-trunc.G])))
      accept     <- a.prob >= 0 || - rexp(1)  < a.prob
        return(list(alpha   = ifelse(accept, alpha.new, alpha.old) - discount, rate = accept))
    }

# Priors
  # Means
    .sim_mu.p    <- function(P, mu.zero, sig.mu.sqrt) {
      sig.mu.sqrt * rnorm(P) + mu.zero
    }
  
  # Scores
    .sim_eta.p   <- function(Q, N) {
        matrix(rnorm(N * Q), nr=N, nc=Q)
    }
  
  # Loadings
    .sim_load.p  <- function(Q, P, sigma.l) {
        sqrt(sigma.l) * matrix(rnorm(P * Q), nr=P, nc=Q)
    }
    
    .sim_load.ps <- function(Q, sigma.l, phi, tau) {
        sqrt(1/(phi * tau)) * rnorm(Q)
    }
  
  # Uniquenesses
    .sim_psi.ipu <- function(P, psi.alpha, psi.beta) {
        rgamma(n=P, shape=psi.alpha, rate=psi.beta) 
    }
    
    .sim_psi.ipi <- function(P, psi.alpha, psi.beta) {
        rep(rgamma(1, shape=psi.alpha, rate=psi.beta), P)
    }

  # Local Shrinkage
    .sim_phi.p   <- function(Q, P, nu, plus1) {
        matrix(rgamma(n=P * Q, shape=nu + plus1, rate=nu), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    .sim_delta.p <- function(Q = 2L, alpha, beta) {
        rgamma(n=Q - 1, shape=alpha, rate=beta)
    }
        
  # Cluster Labels
    .sim_z.p     <- function(N, prob.z) {
        factor(which(rmultinom(N, size=1, prob=prob.z) != 0, arr.ind=TRUE)[,1], levels=seq_along(prob.z))
    }

# Other Functions
  # Uniqueness Hyperparameters
    psi_hyper   <- function(alpha, covar, type=c("unconstrained", "isotropic"), P) {
      inv.cov   <- try(solve(covar), silent=TRUE)
      if(inherits(inv.cov, "try-error"))  {
        inv.cov <- 1/covar
      }
        unname((alpha - 1)/switch(match.arg(type), unconstrained=diag(inv.cov), isotropic=rep(sum(diag(inv.cov))/P, P)))
    }

  # Alpha Shifted Gamma Hyperparameters
    shift_gamma <- function(shape, rate, shift = 0L, param = c("rate", "scale"))   {
      var       <- shape/rate^2
      exp       <- var  * rate + shift
      rate      <- exp/var
      shape     <- rate * exp
        c(shape,   switch(match.arg(param), rate=rate, 1/rate))
    }
    
  # Check Shrinkage Hyperparemeters
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
    .lab.switch <- function(z.new, z.old, Gs, ng = tabulate(z.new)) {
      tab       <- table(z.new, z.old, dnn=NULL)
      tab.tmp   <- tab[rowsums(tab) != 0,colsums(tab) != 0, drop=FALSE]
      nc        <- ncol(tab.tmp)
      nr        <- nrow(tab.tmp)
      if(nc > nr) {
        tmp.mat <- matrix(rep(0, nc), nr=nc - nr, nc=nc)
        rownames(tmp.mat) <- setdiff(as.numeric(colnames(tab.tmp)), as.numeric(rownames(tab.tmp)))[seq_len(nc - nr)]
        tab.tmp <- rbind(tab.tmp, tmp.mat)
      } else if(nr > nc) {
        tmp.mat <- matrix(rep(0, nr), nr=nr, nc=nr - nc)
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
    .lab.move1  <- function(nn.ind, pi.prop, nn) {
      sw        <- sample(nn.ind, 2L)
      pis       <- pi.prop[sw]
      nns       <- nn[sw]
      a.prob    <- (nns[1] - nns[2]) * (log(pis[1])  - log(pis[2])) 
        return(list(rate1  = a.prob >= 0 || -rexp(1) < a.prob, sw = sw))
    }
    
    # Move 2
    .lab.move2  <- function(nn.ind, Vs, nn) {
      sw        <- sample(nn.ind, 1L)
      nn.x      <- which(nn.ind == sw)
      sw        <- c(sw, max(nn.ind[nn.x + 1], nn.ind[nn.x - 1], na.rm=TRUE))
      nns       <- nn[sw]
      Vsw       <- Vs[sw]
      a.prob    <- nns[1] * log(1 - Vsw[1]) - nns[2]  * log(1 - Vsw[2])
        return(list(rate2 = a.prob >= 0 ||  - rexp(1) < a.prob, sw = sw))
    }

  # Length Checker
    .len.check  <- function(obj0g, switch0g, method, P, range.G, P.dim = TRUE) {
      V         <- ifelse(P.dim, P, 1L)
      rGseq     <- seq_along(range.G)
      obj.name  <- deparse(substitute(obj0g))
      sw.name   <- deparse(substitute(switch0g))
      if(!is.list(obj0g))        obj0g  <- list(obj0g)
      if(length(obj0g) != length(range.G))    {
        if(!P.dim) {
          obj0g <- replicate(length(range.G), obj0g)
        } else                             stop(paste0(obj.name, " must be a list of length ", length(range.G)))
      }
      len       <- lengths(obj0g)
      if(is.element(method, c("FA", "IFA")))  {
        if(any(!is.element(len, c(1, V)))) stop(paste0(obj.name, " must be list of length 1 containing a scalar", ifelse(P.dim, paste0(" or a vector of length P=", V), ""), " for a 1-group model"))
      } else {
        if(is.element(len, c(1, range.G, V))) {
          if(all(len == range.G)) obj0g <- if(switch0g) lapply(rGseq, function(g) matrix(obj0g[[g]], nr=1))    else stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"))
          if(all(len == V))       obj0g <- if(V == 1)   lapply(rGseq, function(g) rep(obj0g[[g]], range.G[g])) else lapply(rGseq, function(g) matrix(obj0g[[g]], nr=V, nc=range.G[g]))
        } else if(!all(vapply(rGseq, function(g) is.matrix(obj0g[[g]]) && any(identical(dim(obj0g[[g]]), c(1, range.G[g])), identical(dim(obj0g[[g]]), c(V, range.G[g]))), logical(1)))) {
                                           stop(paste0(ifelse(length(range.G) > 1, "Each element of ", ""), obj.name, " must be either of length 1, ", ifelse(P.dim, paste0("P=", V, ", or it's corresponding range.G, or a matrix with P rows and it's corresponding range.G columns"), paste0("or G=", range.G)))) 
        } else if(all(vapply(obj0g, is.matrix, logical(1)), !switch0g) && any(vapply(rGseq, function(g) any(dim(obj0g[[g]]) == range.G[g]), logical(1)))) {
                                           stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"))
        }
      }
      if(all(length(unique(unlist(obj0g))) > 1,
             !switch0g, !P.dim))           stop(paste0(obj.name, " must be a scalar if ", sw.name, " is TRUE"))
        obj0g
    }

  # Moments of Dirichlet / Pitman-Yor Processes
    G_expected  <- Vectorize(function(N, alpha, discount = 0L) {
      if(!all(is.numeric(N), is.numeric(discount), 
         is.numeric(alpha)))               stop("All inputs must be numeric")
      if(discount  < 0  || discount >= 1)  stop("Invalid discount value")
      if(alpha   < - discount)             stop("Invalid alpha value")
      if(discount == 0) {
          alpha * (digamma(alpha + N) - digamma(alpha))
      } else {
        if(suppressMessages(require(Rmpfr))) {
          on.exit(.detach_pkg(Rmpfr))
          on.exit(.detach_pkg(gmp), add=TRUE)  
        } else                             stop("'Rmpfr' package not installed")
          asNumeric(alpha/discount * pochMpfr(alpha + discount, N)/pochMpfr(alpha, N) - alpha/discount)
      }
    })

    G_variance  <- Vectorize(function(N, alpha, discount = 0L) {
      if(!all(is.numeric(N), is.numeric(discount), 
         is.numeric(alpha)))               stop("All inputs must be numeric")
      if(discount  < 0  || discount >= 1)  stop("Invalid discount value")
      if(alpha   < - discount)             stop("Invalid alpha value")
      if(discount == 0) {
          alpha * (digamma(alpha + N) - digamma(alpha)) + (alpha^2) * (trigamma(alpha + N) - trigamma(alpha))
      } else {
        if(suppressMessages(require(Rmpfr))) {
          on.exit(.detach_pkg(Rmpfr))
          on.exit(.detach_pkg(gmp), add=TRUE)  
        } else                             stop("'Rmpfr' package not installed")
        sum.ad  <- alpha + discount
        poch.a  <- pochMpfr(alpha, N)
        poch.ad <- pochMpfr(sum.ad, N)
        subterm <- alpha/discount * poch.ad/poch.a
          asNumeric((alpha * sum.ad)/discount^2 * pochMpfr(sum.ad + discount, N)/poch.a - subterm - subterm^2)
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
    print.IMIFA <- summary.IMIFA <- function(imifa) {
      meth      <- attr(imifa, "Method")
      name      <- attr(imifa, "Name")
      fac       <- attr(imifa, "Factors")
      grp       <- attr(imifa, "Groups")
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
        cat(paste0(meth, " simulations for '", name, "' dataset", msg, " to be passed to tune.IMIFA(...)\n"))
    }
    
    print.Tuned_IMIFA   <- function(res) {
      method    <- attr(res, "Method")
      G         <- res$GQ.results$G
      Q         <- res$GQ.results$Q
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
  
    summary.Tuned_IMIFA <- function(res) {
      criterion <- unlist(strsplit(toupper(attr(res$GQ.results, "Criterion")), "[.]"))
      criterion <- ifelse(length(criterion) > 1, ifelse(criterion[1] != "LOG", paste0(criterion[1], ".", tolower(criterion[2])), "LogIntegratedLikelihood"), criterion)
      crit.mat  <- res$GQ.results[[paste0(criterion, "s")]]
      msg       <- NULL
      if(any(dim(crit.mat) > 1)) {
        msg     <- paste0(", and ", ifelse(substr(criterion, 1, 1) == "A", "an ", "a "),  criterion, " of ", round(max(crit.mat), 2), "\n")  
      }
        cat(paste0(capture.output(print(res)), msg))
    }
    
    .power2     <- function(x) x * x
    .which0     <- function(x) which(x == 0)
    .chol       <- function(x) tryCatch(chol(x), error=function(e) chol(make.positive.definite(x)))