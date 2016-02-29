#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.sims     <- function(sims = NULL, burnin = 0, thinning = 1, 
                          Q = NULL, Q.meth = NULL, recomp = F, ...) {
  defpar      <- par(no.readonly = T)
  defop       <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(missing(sims))             stop("Simulations must be supplied")
  if(!exists(deparse(substitute(sims)),
             envir=.GlobalEnv)) stop(paste0("Object ", match.call()$sims, " not found"))
  if(class(sims) != "IMIFA")    stop(paste0("Simulations object of class 'IMIFA' must be supplied"))
  store       <- seq(from=burnin + 1, to=attr(sims, "Store"), by=thinning)
  n.store     <- length(store)
  if(n.store  <= 1)             stop(paste0("burnin must be less than the stored number of iterations"))
  method      <- attr(sims, "Method")
  n.fac       <- attr(sims, "Factors")
  n.obs       <- attr(sims, "Obs")
  n.var       <- attr(sims, "Vars")
  sw          <- attr(sims, "Switch")
  cov.emp     <- sims[[1]]$cov.mat
  if(!is.logical(recomp))       stop("recomp must be TRUE or FALSE")
  if(thinning  > 1       ||
     burnin    > 0)   recomp <- T
  
  if(!missing(Q)) {
    Q.x       <- Q
    Q.ind     <- which(n.fac == Q.x)
    if(method == "FA"    &&
       !is.element(Q, n.fac))   stop("This Q value was not used during simulation")
    if(method == "IFA"   &&
      (Q * (n.fac - Q))   < 0)  stop(paste0("Q cannot be greater than the number of factors in ", match.call()$sims))
  } 
  Q.T         <- exists("Q.x", envir=environment())
  if(method == "IFA") {
    if(missing(Q.meth)) {
      Q.meth  <- "Mode"
    } else {
      if(!is.element(Q.meth, 
       c("Mode", "Median")))    stop("Q.meth must be MODE or MEDIAN")
    }
    
  # Retrieve distribution of Q, tabulate & plot
    Q.ind     <- 1
    Q.store   <- sims[[Q.ind]]$Q.store[store]
    Q.tab     <- table(Q.store, dnn=NULL)
    Q.prob    <- prop.table(Q.tab)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode    <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median  <- ceiling(median(Q.store) * 2)/2
    if(Q.T) {
       Q      <- Q.x
    } else if(Q.meth     == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- Q.median
    if(all(!sw["l.sw"],
           !sw["p.sw"]))       warning("Loadings &/or Uniquenesses not stored: can't calculate proportion of variation explained", call.=F)
    Q.CI      <- round(quantile(Q.store, c(0.025, 0.975)))
    Q.res     <- list(Q = Q, Mode = Q.mode, Median = Q.median, 
                      CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
  } else if(any(sw["l.sw"], sw["p.sw"])) {     
  
  # Calculate Proportion of Variation Explained & BIC
    Q.range   <- 1:length(n.fac)
    cum.var   <- rep(NA, length(n.fac))
    if(sw["p.sw"]   && (sw["l.sw"] || (length(n.fac) == 1 && n.fac == 0))) {
      bic     <- rep(NA, length(n.fac))
      if(!sw["mu.sw"]) {
        post.mu     <- rep(0, n.var) 
                                warning(paste0("Zero mean vector used in BIC calculation as means weren't stored"), call.=F)
      }
    }
    temp.b    <- max(1, burnin)
    data      <- attr(sims, "Name")
    if(!exists(data,
       envir=.GlobalEnv)) {     warning(paste0("Object ", data, " not found: can't compute BIC"), call.=F)
    } else {
      data    <- as.data.frame(get(data))
      data    <- data[sapply(data, is.numeric)]
      cent    <- attr(sims, "Center")
      scaling <- attr(sims, "Scaling")
      data    <- scale(data, center=cent, scale=scaling)
    }
    for(q in Q.range) {
      Q.fac   <- n.fac[q]
      if(sw["l.sw"] && Q.fac > 0) {
        lmat        <- sims[[q]]$load[,1:Q.fac,store, drop=F]
        l.temp      <- as.matrix(sims[[q]]$load[,1:Q.fac,temp.b])
        for(p in 1:n.store) {
          rot       <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
          lmat[,,p] <- lmat[,,p] %*% rot
        }
        post.load   <- rowMeans(lmat, dims=2)
        cum.var[q]  <- sum(colSums(post.load * post.load))/n.var
      } else {
        post.load   <- matrix(, nr=n.var, nc=0)
      }
      if(sw["mu.sw"]) {
        mu    <- sims[[q]]$mu[,store]
        post.mu     <- rowMeans(mu, 1)
      } 
      if(sw["p.sw"])  {
        psi     <- sims[[q]]$psi[,store]
        post.psi    <- rowMeans(psi, 1)
        if(!sw["l.sw"]   ||
           (sw["l.sw"]   && Q.fac == 0)) {
          cum.var[q]     <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
        }  
        if(sw["l.sw"]    || Q.fac == 0) {
          Sigma <- tcrossprod(post.load) + diag(post.psi)
          K     <- n.var * Q.fac - 0.5 * Q.fac * (Q.fac - 1) + n.var
          if(Q.fac   > 0) {
            U.Sig   <- chol(Sigma)
          } else {
            U.Sig   <- sqrt(Sigma)
          }
          rooti <- backsolve(U.Sig, diag(n.var))
          quad  <- crossprod(rooti, t(data) - post.mu)
          quads <- colSums(quad * quad)
          log.lik   <- sum(log(exp(- n.var/2 * log(2 * pi) + sum(log(diag(rooti))) - 0.5 * quads)))
          bic[q]    <- 2 * log.lik - K * log(n.obs)
        }
      }
    }
    if(max(cum.var[!is.na(cum.var)]) > 1    ||
       (exists("bic", envir=environment())  &&
       any(!is.finite(bic))))   warning("Chain may not have converged", call.=F)
    if(Q.T) {
      cum.var <- cum.var[Q.ind]
      bic     <- bic[Q.ind]
      Q       <- Q.x
      n.fac   <- Q
    } else if(sw["p.sw"] && (sw["l.sw"]     || (length(n.fac) == 1 && n.fac == 0))) {
      Q.ind   <- which.max(bic)
      Q       <- n.fac[Q.ind]
    } else {
      Q.ind   <- which.max(cum.var)
      Q       <- n.fac[Q.ind]
    }
    if(exists("bic", envir=environment())) {
      Q.res     <- list(Q = Q, BIC = bic, cum.var = cum.var)
    } else {
      Q.res     <- list(Q = Q, cum.var = cum.var)
    }
  }  else if(length(n.fac) == 1) {
     Q.ind    <- 1
     Q        <- n.fac
     Q.res    <- list(Q = Q)
                                warning("Loadings &/or Uniquenesses not stored: can't calculate BIC or proportion of variation explained", call.=F)
  } 
  if(Q == 0) {
    sw[c("f.sw", "l.sw")]  <- F
  }
  if(method   == "FA") {
    if(sw["f.sw"]) {
      f       <- sims[[Q.ind]]$f[,1:Q,store, drop=F]
    }
    if(sw["l.sw"]) {
      lmat    <- sims[[Q.ind]]$load[,1:Q,store, drop=F]
      temp.b  <- max(1, burnin)
      l.temp  <- as.matrix(sims[[Q.ind]]$load[,1:Q,temp.b])
    }
  } else {
    store     <- store[which(Q.store >= Q)]
    n.store   <- length(store)
   #store     <- tail(store, 0.9 * n.store)
    temp.b    <- store[1]
    if(sw["f.sw"]) {
      f       <- as.array(sims[[Q.ind]]$f)[,1:Q,store, drop=F]
    }
    if(sw["l.sw"]) {
      lmat    <- as.array(sims[[Q.ind]]$load)[,1:Q,store, drop=F]
      l.temp  <- as.matrix(as.array(sims[[Q.ind]]$load)[,1:Q,temp.b])
    }
  }
  if(sw["mu.sw"])  {
    mu        <- sims[[Q.ind]]$mu[,store]                            
  }
  if(sw["p.sw"])   {
    psi       <- sims[[Q.ind]]$psi[,store]
  }  
  
# Loadings matrix / identifiability / # etc.  
  if(sw["l.sw"])   {
    for(p in 1:n.store) {
      rot       <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
      lmat[,,p] <- lmat[,,p] %*% rot
      if(sw["f.sw"]) {
        f[,,p]  <- f[,,p]    %*% rot
      }  
    }
  }

  if(sw["mu.sw"]) post.mu  <- rowMeans(mu, 1)
  if(sw["f.sw"])  post.f   <- rowMeans(f, dims=2)
  if(sw["p.sw"])  post.psi <- rowMeans(psi, 1)
  if(sw["l.sw"])  {
    post.load <- rowMeans(lmat, dims=2)
    SS.load   <- colSums(post.load * post.load)
    comm      <- sum(SS.load)
    prop.var  <- SS.load/n.var
    cum.var   <- cumsum(prop.var)          
    prop.exp  <- comm/n.var
    prop.uni  <- 1 - prop.exp
  } else if(sw["p.sw"]) {
    prop.exp  <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
    cum.var   <- max(prop.exp)
  }
  cov.est     <- sims[[Q.ind]]$post.Sigma
  if(recomp) {
    if(!all(c(sw["l.sw"], 
              sw["p.sw"]))) {   warning("Loadings and/or uniquenesses not stored: can't recompute Sigma", call.=F)
    } else {
      cov.est <- replace(cov.est, is.numeric(cov.est), 0)
      for(r in 1:n.store) {
        Sigma   <- tcrossprod(lmat[,,r]) + diag(psi[,r])
        cov.est <- cov.est + Sigma/n.store
      } 
    }
  }
  error       <- cov.emp - cov.est
  MSE         <- mean(error * error)
  MAD         <- mean(abs(error))
  error       <- list(MSE = MSE, MAD = MAD) 
  if(any(ifelse(all(isTRUE(attr(sim, "Scaling")), attr(sim, "Center")), 
            sum(round(diag(cov.est)) != 
            round(diag(cov.emp)))    != 0, F), 
     ifelse(sw["p.sw"], 
            sum(abs(post.psi - (1 - post.psi)) < 0) != 0, F),
     ifelse(sw["l.sw"], 
            prop.exp > 1, F)))  warning("Chain may not have converged", call.=F)

  if(sw["l.sw"]) {
    class(post.load)       <- "loadings"
  }
  results     <- list(if(sw["mu.sw"]) list(means = mu, post.mu = post.mu),
                      if(sw["f.sw"])  list(scores = f, post.f = post.f),
                      if(sw["p.sw"])  list(uniquenesses = psi, post.psi = post.psi),
                      if(sw["l.sw"])  list(loadings = lmat, post.load = post.load,  
                                           communality = comm, SS.load = SS.load,
                                           prop.var = prop.var, prop.uni = prop.uni),
                      if(any(sw[c("p.sw", "l.sw")])) list(prop.exp = prop.exp, cum.var = cum.var),
                      if(exists("cov.emp")) list(error = error, cov.mat = cov.emp), list(post.Sigma = cov.est))
  results     <- unlist(results, recursive=F)
  if(Q.T) {
    attr(Q.res,
         "Factors")        <- Q.x
  } else  {
    attr(Q.res, 
         "Factors")        <- n.fac
  } 
  attr(Q.res, "Supplied")  <- Q.T
  results     <- c(results, Q.results = list(Q.res))
  class(results)           <- "IMIFA"
  attr(results, "Method")  <- attr(sims, "Method")
  attr(results, "Obs")     <- n.obs
  attr(results, "Store")   <- n.store
  attr(results, "Switch")  <- sw
  attr(results, "Vars")    <- n.var
  return(results)
}