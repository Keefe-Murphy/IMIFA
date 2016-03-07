#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.sims     <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL,
                          Q = NULL, Q.meth = c("Mode", "Median"), recomp = F, ...) {
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
  n.grp       <- attr(sims, "Groups")
  n.obs       <- attr(sims, "Obs")
  n.var       <- attr(sims, "Vars")
  sw          <- attr(sims, "Switch")
  G.ind       <- 1
  Q.ind       <- 1
  cov.emp     <- sims[[G.ind]][[Q.ind]]$cov.mat
  if(!is.logical(recomp))       stop("recomp must be TRUE or FALSE")
  if(any(thinning > 1, 
         burnin > 0)) recomp <- T
  
  if(!missing(G)) {
    G.x       <- G
    G.xind    <- which(n.grp == G.x)
    if(all(is.element(method, c("MFA", "MIFA")),
       !is.element(G, n.grp)))  stop("This G value was not used during simulation")
    if(all(method == "IFA", 
      (G * (n.grp - G)) < 0))   stop(paste0("G cannot be greater than the number of groups in ", match.call()$sims))
  } 
  G.T         <- exists("G.x", envir=environment())
  if(!missing(Q)) {
    Q.x       <- Q
    Q.xind    <- which(n.fac == Q.x)
    if(all(is.element(method, c("FA", "MFA")),
       !is.element(Q, n.fac)))  stop("This Q value was not used during simulation")
    if(all(method == "IFA", 
      (Q * (n.fac - Q)) < 0))   stop(paste0("Q cannot be greater than the number of factors in ", match.call()$sims))
  } 
  Q.T         <- exists("Q.x", envir=environment())
  
  if(method   == "IFA") {
    if(missing(Q.meth)) {
      Q.meth  <- "Mode"
    } else {
      Q.meth  <- match.arg(Q.meth)
    }
    
  # Retrieve distribution of Q, tabulate & plot
    Q.store   <- sims[[G.ind]][[Q.ind]]$Q.store[store]
    Q.tab     <- table(Q.store, dnn=NULL)
    Q.prob    <- prop.table(Q.tab)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode    <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median  <- ceiling(median(Q.store) * 2)/2
    if(Q.T) {
       Q      <- Q.x
    } else if(Q.meth == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- Q.median
    Q.CI      <- round(quantile(Q.store, c(0.025, 0.975)))
    Q.res     <- list(Q = Q, Mode = Q.mode, Median = Q.median, 
                      CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
  } else {     
  
  # Calculate Proportion of Variation Explained & BIC
    G.range   <- 1:length(n.grp)
    Q.range   <- 1:length(n.fac)
    cumvar    <- matrix(NA, nr=length(n.grp), nc=length(n.fac))
    if(any(sw["l.sw"], all(length(n.fac) == 1, n.fac == 0))) {
      bic     <- matrix(NA, nr=length(n.grp), nc=length(n.fac))
    }
    bic.x     <- exists("bic", envir=environment())
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
    for(g in G.range) {
      for(q in Q.range) {
        Q.fac <- n.fac[q]
        if(all(sw["l.sw"], Q.fac > 0)) {
          lmat        <- sims[[g]][[q]]$load[,1:Q.fac,store, drop=F]
          l.temp      <- as.matrix(sims[[g]][[q]]$load[,1:Q.fac,temp.b])
          for(p in 1:n.store) {
            rot       <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
            lmat[,,p] <- lmat[,,p] %*% rot
          }
          post.load   <- rowMeans(lmat, dims=2)
          cumvar[q]   <- sum(colSums(post.load * post.load))/n.var
        } else {
          post.load   <- matrix(, nr=n.var, nc=0)
        }
        post.mu       <- sims[[g]][[q]]$post.mu
        post.psi      <- sims[[g]][[q]]$post.psi
        if(sw["mu.sw"]) {
          mu  <- sims[[g]][[q]]$mu[,store]
          post.mu     <- rowMeans(mu, 1)
        } 
        if(sw["p.sw"])  {
          psi <- sims[[g]][[q]]$psi[,store]
          post.psi    <- rowMeans(psi, 1)
        }
        if(any(!sw["l.sw"], all(sw["l.sw"], Q.fac == 0))) {
          cumvar[g,q] <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
        }  
        if(any(sw["l.sw"], Q.fac == 0)) {
          Sigma       <- tcrossprod(post.load) + diag(post.psi)
          K           <- n.var * Q.fac - 0.5 * Q.fac * (Q.fac - 1) + n.var
          if(Q.fac > 0) {
            U.Sig     <- chol(Sigma)
          } else {
            U.Sig     <- sqrt(Sigma)
          }
          rooti       <- backsolve(U.Sig, diag(n.var))
          quad        <- crossprod(rooti, t(data) - post.mu)
          quads       <- colSums(quad * quad)
          log.lik     <- sum(log(exp(- n.var/2 * log(2 * pi) + sum(log(diag(rooti))) - 0.5 * quads)))
          bic[g,q]    <- 2 * log.lik - K * log(n.obs)
        }
      }
    }
    if(any(max(cumvar[!is.na(cumvar)]) > 1, bic.x &&
       any(!is.finite(bic))))   warning("Chain may not have converged", call.=F)
    if(bic.x && !any(is.nan(bic), !is.finite(bic))) {
      bic.max <- which(bic == max(bic), arr.ind = T)
      G.ind   <- bic.max[1]
      Q.ind   <- bic.max[2]
      G       <- n.grp[G.ind]
      Q       <- n.fac[Q.ind]
    } else {
      var.max <- which(cumvar == max(cumvar), arr.ind = TRUE)
      G.ind   <- var.max[1]
      Q.ind   <- var.max[2]
      G       <- n.grp[G.ind]
      Q       <- n.fac[Q.ind]
    }
    if(all(Q.T, G.T)) {
      cumvar  <- cumvar[G.xind,Q.xind]
      if(bic.x) {
        bic   <- bic[G.xind,Q.xind]
      }
      G       <- G.x
      Q       <- Q.x
      n.grp   <- G
      n.fac   <- Q  
      G.ind   <- G.xind
      Q.ind   <- Q.xind
    } else if(Q.T) {
      cumvar  <- cumvar[G.ind,Q.xind]
      if(bic.x) {
        bic   <- bic[G.ind,Q.xind]
      }
      Q       <- Q.x
      n.fac   <- Q
      Q.ind   <- Q.xind
    } else if(G.T) {
      cumvar  <- cumvar[G.xind,Q.ind]
      if(bic.x) {
        bic   <- bic[G.xind,Q.ind]
      }
      G       <- G.x
      n.grp   <- G
      G.ind   <- G.xind
    }
    if(bic.x) {
      Q.res   <- list(Q = Q, G = G, BIC = bic, cum.var = cumvar)
    } else {
      Q.res   <- list(Q = Q, G = G, cum.var = cumvar)
    }
  } 
  
  if(Q == 0) {
    sw[c("f.sw", "l.sw")]    <- F
  }
  if(method   == "FA") {
    if(sw["f.sw"]) {
      f       <- sims[[G.ind]][[Q.ind]]$f[,1:Q,store, drop=F]
    }
    if(sw["l.sw"]) {
      lmat    <- sims[[G.ind]][[Q.ind]]$load[,1:Q,store, drop=F]
      temp.b  <- max(1, burnin)
      l.temp  <- as.matrix(sims[[G.ind]][[Q.ind]]$load[,1:Q,temp.b])
    }
  } else {
    store     <- store[which(Q.store >= Q)]
    n.store   <- length(store)
   #store     <- tail(store, 0.9 * n.store)
    temp.b    <- store[1]
    if(sw["f.sw"]) {
      f       <- as.array(sims[[G.ind]][[Q.ind]]$f)[,1:Q,store, drop=F]
    }
    if(sw["l.sw"]) {
      lmat    <- as.array(sims[[G.ind]][[Q.ind]]$load)[,1:Q,store, drop=F]
      l.temp  <- as.matrix(as.array(sims[[G.ind]][[Q.ind]]$load)[,1:Q,temp.b])
    }
  }
  post.mu     <- sims[[G.ind]][[Q.ind]]$post.mu
  post.psi    <- sims[[G.ind]][[Q.ind]]$post.psi
  if(sw["mu.sw"])  {
    mu        <- sims[[G.ind]][[Q.ind]]$mu[,store]                            
    post.mu   <- rowMeans(mu, 1)
  }
  if(sw["p.sw"])   {
    psi       <- sims[[G.ind]][[Q.ind]]$psi[,store]
    post.psi  <- rowMeans(psi, 1)
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

  if(sw["f.sw"])  post.f     <-  rowMeans(f, dims=2)
  if(sw["l.sw"])  {
    post.load <- rowMeans(lmat, dims=2)
    SS.load   <- colSums(post.load * post.load)
    comm      <- sum(SS.load)
    prop.var  <- SS.load/n.var
    cum.var   <- cumsum(prop.var)          
    prop.exp  <- comm/n.var
    prop.uni  <- 1 - prop.exp
  } else {
    prop.exp  <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
    cum.var   <- max(prop.exp)
  }
  cov.est     <- sims[[G.ind]][[Q.ind]]$post.Sigma
  if(all(recomp, sw["l.sw"])) {
    cov.est   <- replace(cov.est, is.numeric(cov.est), 0)
    for(r in 1:n.store) {
      Sigma   <- tcrossprod(lmat[,,r]) + diag(psi[,r])
      cov.est <- cov.est + Sigma/n.store
    }
  } else if(recomp)             warning("Loadings not stored: can't recompute Sigma", call.=F)
  error       <- cov.emp - cov.est
  MSE         <- mean(error * error)
  RMSE        <- sqrt(MSE)
  MAD         <- mean(abs(error))
  error       <- list(MSE = MSE, RMSE = RMSE, MAD = MAD) 
  if(any(ifelse(all(isTRUE(attr(sim, "Scaling")), attr(sim, "Center")), 
            sum(round(diag(cov.est))  != 
            round(diag(cov.emp)))     != 0, F), 
     sum(abs(post.psi - (1 - post.psi))  <  0) != 0,
     ifelse(sw["l.sw"], 
            prop.exp > 1, F)))  warning("Chain may not have converged", call.=F)

  if(sw["l.sw"]) {
    class(post.load)         <- "loadings"
  }
  results     <- list(if(sw["mu.sw"])    list(means  = mu), 
                      list(post.mu    =  post.mu),
                      if(sw["f.sw"])     list(scores = f, 
                                              post.f = post.f),
                      if(sw["p.sw"])     list(uniquenesses = psi), 
                      list(post.psi   =  post.psi),
                      if(sw["l.sw"])     list(loadings     = lmat, 
                                              post.load    = post.load,  
                                              communality  = comm, 
                                              SS.load      = SS.load,
                                              prop.var     = prop.var, 
                                              prop.uni     = prop.uni),
                      list(prop.exp   =  prop.exp, 
                           cum.var    =  cum.var),
                      if(exists("cov.emp"))   list(error   = error, 
                                                   cov.mat = cov.emp), 
                      list(post.Sigma = cov.est))
  results     <- unlist(results, recursive=F)
  if(Q.T) {
    attr(Q.res,
         "Factors")          <- Q.x
  } else  {
    attr(Q.res, 
         "Factors")          <- n.fac
  } 
  if(G.T) {
    attr(Q.res,
         "Groups")           <- G.x
  } else  {
    attr(Q.res, 
         "Groups")           <- n.grp
  }
  attr(Q.res, "Supplied")    <- c(Q=Q.T, G=G.T)
  results     <- c(results, Q.results = list(Q.res))
  class(results)             <- "IMIFA"
  attr(results, "Method")    <- attr(sims, "Method")
  attr(results, "Obs")       <- n.obs
  attr(results, "Store")     <- n.store
  attr(results, "Switch")    <- sw
  attr(results, "Vars")      <- n.var
  return(results)
}