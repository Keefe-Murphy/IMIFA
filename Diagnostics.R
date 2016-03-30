#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.imifa       <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL, Q = NULL,
                             criterion = c("bic", "aic"), Q.meth = c("Mode", "Median"), recomp = F, ...) {
  defpar         <- par(no.readonly = T)
  defop          <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(missing(sims))               stop("Simulations must be supplied")
  if(!exists(deparse(substitute(sims)),
             envir=.GlobalEnv))   stop(paste0("Object ", match.call()$sims, " not found"))
  if(class(sims) != "IMIFA")      stop(paste0("Simulations object of class 'IMIFA' must be supplied"))
  store          <- seq(from=burnin + 1, to=attr(sims, "Store"), by=thinning)
  n.store        <- length(store)
  if(n.store     <= 1)            stop(paste0("burnin must be less than the stored number of iterations"))
  method         <- attr(sims, "Method")
  n.fac          <- attr(sims, "Factors")
  n.grp          <- attr(sims, "Groups")
  n.obs          <- attr(sims, "Obs")
  n.var          <- attr(sims, "Vars")
  sw             <- attr(sims, "Switch")
  criterion      <- match.arg(criterion)
  if(!is.logical(recomp))         stop("recomp must be TRUE or FALSE")
  if(any(burnin   > 0, 
     thinning     > 1)) recomp <- T
  
  G.ind          <- 1
  Q.ind          <- 1
  G.T            <- !missing(G)
  Q.T            <- !missing(Q)
  if(G.T) {
    if(!is.element(method, c("FA", "IFA"))) {
      G.ind      <- which(n.grp == G)
    } else if(G > 1)              warning(paste0("G must be equal to 1 for the ", method, " method"), call.=F)
    if(all(is.element(method, c("MFA", "MIFA")),
       !is.element(G, n.grp)))    stop("This G value was not used during simulation")
  } 
  if(Q.T) {
    if(!is.element(method, c("IFA", "MIFA", "IMIFA"))) {
      Q.ind      <- which(n.fac == Q)
    }
    if(all(is.element(method, c("FA", "MFA")),
       !is.element(Q, n.fac)))    stop("This Q value was not used during simulation")
    if(all(is.element(method, c("IFA", "classify")), 
      (Q * (n.fac - Q)) < 0))     stop(paste0("Q cannot be greater than the number of factors in ", match.call()$sims))
  } 
  G              <- ifelse(all(G.T, !is.element(method, c("FA", "IFA"))), G, 1)
  
  if(method == "IFA") {
    if(missing(Q.meth)) {
      Q.meth     <- "Mode"
    } else   {
      Q.meth     <- match.arg(Q.meth)
    }
    
  # Retrieve distribution of Q, tabulate & plot
    Q.store      <- sims[[G.ind]][[Q.ind]]$Q.store[store]
    Q.tab        <- table(Q.store, dnn=NULL)
    Q.prob       <- prop.table(Q.tab)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode       <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.med        <- ceiling(median(Q.store) * 2)/2
    if(!Q.T) {
      if(Q.meth == 'Mode') { 
        Q        <- min(Q.mode)
      } else {
        Q        <- Q.med
      }
    }
    Q.CI         <- round(quantile(Q.store, c(0.025, 0.975)))
    GQ.res       <- list(G = G, Q = Q, Mode = Q.mode, Median = Q.med, 
                         CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
  }
    
  if(is.element(method, c("FA", "MFA"))) {
    G.range      <- length(n.grp)
    Q.range      <- length(n.fac)
    
  # Retrieve AIC & BIC to tune G & Q   
    aic  <- bic  <- matrix(NA, nr=G.range, nc=Q.range, dimnames=list(paste0("G", n.grp), paste0("Q", n.fac)))
    for(g in seq_len(G.range)) { 
      for(q in seq_len(Q.range)) {
        aic[g,q] <- sim[[g]][[q]]$aic
        bic[g,q] <- sim[[g]][[q]]$bic  
      }  
    }
    crit         <- get(criterion)
    crit.max     <- which(crit == max(crit), arr.ind = T)
    if(all(Q.T, G.T)) {
      aic        <- aic[G.ind,Q.ind]
      bic        <- bic[G.ind,Q.ind]
    } else if(Q.T) {
      aic        <- aic[,Q.ind, drop=F]
      bic        <- bic[,Q.ind, drop=F]
      crit       <- crit[,Q.ind]
      G.ind      <- which(crit == max(crit))
      G          <- n.grp[G.ind]
    } else if(G.T) {
      aic        <- aic[G.ind,, drop=F]
      bic        <- bic[G.ind,, drop=F]
      crit       <- crit[G.ind,]
      Q.ind      <- which(crit == max(crit))
      Q          <- n.fac[Q.ind]
    } else {
      G.ind      <- crit.max[1]
      Q.ind      <- crit.max[2]
      G          <- n.grp[G.ind]
      Q          <- n.fac[Q.ind]
    }
    Q            <- setNames(rep(Q, G), paste0("Qg", seq_len(G)))
    GQ.res       <- list(G = G, Q = Q, AIC = aic, BIC = bic)
  }
  
  if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
    z            <- sims[[G.ind]][[Q.ind]]$z[,store]
    post.z       <- sims[[G.ind]][[Q.ind]]$post.z
    if(sw["pi.sw"])    {
      pi.prop    <- sims[[G.ind]][[Q.ind]]$pi.prop[,store]
      post.pi    <- rowMeans(pi.prop, dims=1)
    } else {
      post.pi    <- sims[[G.ind]][[Q.ind]]$post.pi
    }
    if(recomp) {
      post.z     <- replace(post.z, post.z, apply(z, 1, function(x) factor(which.max(tabulate(x)), levels=seq_len(G))))
      if(!sw["pi.sw"]) {
        post.pi  <- replace(post.pi, post.pi, prop.table(tabulate(post.z, nbins=G)))
      }
    }
    cluster      <- list(post.z = post.z, post.pi = post.pi, z = z)
    cluster      <- c(cluster, if(sw["pi.sw"]) list(pi.prop  = pi.prop))
  }
  
  if(all(Q == 0)) {
    if(sw["f.sw"])                warning("Scores not stored as model has zero factors", call.=F)
    sw["f.sw"]   <- F
  }
  if(sw["f.sw"])  {
    Qms          <- seq_len(max(Q))
    if(is.element(method, c("FA", "MFA"))) {
      f          <- sims[[G.ind]][[Q.ind]]$f[,Qms,store, drop=F]
    }
    if(is.element(method, c("IFA", "MIFA", "IMIFA"))) {
      f          <- as.array(sims[[G.ind]][[Q.ind]]$f)[,Qms,store, drop=F]
    }
  }
  
  result         <- list(list())
  temp.b         <- max(1, burnin)
  MSE  <- RMSE   <- NRMSE   <- CVRMSE   <- MAD   <- rep(NA, G)
  for(g in seq_len(G)) {
    Qg           <- Q[g]
    Qgs          <- seq_len(Qg)
    if(Qg == 0) {
      if(sw["l.sw"])              warning(paste0("Loadings not stored as", ifelse(G > 1, paste0(" group ", g), " model"), " has zero factors"), call.=F)
      sw["l.sw"] <- F
    } else {
      sw["l.sw"] <- attr(sims, "Switch")["l.sw"]
    }
    
    if(sw["l.sw"]) {
      if(all(method == "MFA", G > 1)) {
        lmat     <- adrop(sims[[G.ind]][[Q.ind]]$load[,Qgs,g,store, drop=F], drop=3)
        l.temp   <- adrop(sims[[G.ind]][[Q.ind]]$load[,Qgs,g,temp.b, drop=F], drop=3:4)
      }
      if(any(method == "FA",  all(method == "MFA",  G == 1))) {
        lmat     <- sims[[G.ind]][[Q.ind]]$load[,Qgs,store, drop=F]
        l.temp   <- adrop(sims[[G.ind]][[Q.ind]]$load[,Qgs,temp.b, drop=F], drop=3)
      }
    }
    
    if(any(method == "IFA", all(method == "MIFA", G == 1))) {
      store      <- store[which(Q.store >= Qg)]
      n.store    <- length(store)
      temp.b     <- store[1]
      if(sw["l.sw"]) {
        lmat     <- as.array(sims[[G.ind]][[Q.ind]]$load)[,Qgs,store, drop=F]
        l.temp   <- adrop(as.array(sims[[G.ind]][[Q.ind]]$load)[,Qgs,temp.b, drop=F], drop=3)
      }
    }
    
    if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
      cov.emp    <- sims[[G.ind]][[Q.ind]]$cov.mat[,,g]
      cov.est    <- sims[[G.ind]][[Q.ind]]$post.Sigma[,,g]
      post.mu    <- sims[[G.ind]][[Q.ind]]$post.mu[,g]
      post.psi   <- sims[[G.ind]][[Q.ind]]$post.psi[,g]
      if(sw["mu.sw"])  {
        mu       <- sims[[G.ind]][[Q.ind]]$mu[,g,store]                            
      }
      if(sw["psi.sw"]) {
        psi      <- sims[[G.ind]][[Q.ind]]$psi[,g,store]
      }
    } else {
      cov.emp    <- sims[[G.ind]][[Q.ind]]$cov.mat
      cov.est    <- sims[[G.ind]][[Q.ind]]$post.Sigma
      post.mu    <- sims[[G.ind]][[Q.ind]]$post.mu
      post.psi   <- sims[[G.ind]][[Q.ind]]$post.psi
      if(sw["mu.sw"])  {
        mu       <- sims[[G.ind]][[Q.ind]]$mu[,store]                            
      }
      if(sw["psi.sw"]) {
        psi      <- sims[[G.ind]][[Q.ind]]$psi[,store]
      }
    }
    
  # Loadings matrix / identifiability / error metrics / etc.  
    if(sw["l.sw"])     {
      for(p in seq_len(n.store)) {
        rot          <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
        lmat[,,p]    <- lmat[,,p] %*% rot
        if(sw["f.sw"]) {
          if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
            f[post.z == g,,p]  <- f[post.z == g,,p]   %*% rot
          } else {
            f[,,p]   <- f[,,p]    %*% rot
          }
        }  
      }
    }
    
    if(all(recomp, sw[c("l.sw", "psi.sw")])) {
      cov.est    <- replace(cov.est, is.numeric(cov.est), 0)
      for(r in seq_len(n.store)) {
        Sigma    <- tcrossprod(lmat[,,r]) + diag(psi[,r])
        cov.est  <- cov.est + Sigma/n.store
      }
    } else if(recomp)  {
      if(!sw["l.sw"])  {
                                  warning("Loadings not stored: can't re-estimate Sigma", call.=F)
      } else {
                                  warning("Uniquenesses not stored: can't re-estimate Sigma", call.=F)
      }
    }
    if(sw["mu.sw"])  post.mu   <- rowMeans(mu, dims=1)
    if(sw["psi.sw"]) post.psi  <- rowMeans(psi, dims=1)
    if(sw["l.sw"]) { post.load <- rowMeans(lmat, dims=2)
      var.exp    <- sum(colSums(post.load * post.load))/n.var
    } else   {
      var.exp    <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
    }
    error        <- cov.emp - cov.est
    MSE[g]       <- mean(error * error)
    RMSE[g]      <- sqrt(MSE)
    NRMSE[g]     <- RMSE/(max(cov.emp) - min(cov.emp))
    CVRMSE[g]    <- RMSE/mean(cov.emp)
    MAD[g]       <- mean(abs(error))
    if(any(all(isTRUE(attr(sims, "Scaling")), attr(sims, "Center")) && 
               sum(round(diag(cov.est))  != 
               round(diag(cov.emp)))     != 0,
       sum(abs(post.psi - (1 - post.psi)) < 0) != 0,
       var.exp    > 1))           warning(paste0(ifelse(G == 1, "C", paste0("Group ", g, "'s c")), "hain may not have converged"), call.=F)
  
    if(sw["l.sw"]) {
      class(post.load)         <- "loadings"
    }  
    results      <- list(if(sw["mu.sw"])   list(means = mu), 
                         list(post.mu    = post.mu),
                         if(sw["psi.sw"])  list(uniquenesses = psi), 
                         list(post.psi   = post.psi),
                         if(sw["l.sw"])    list(loadings     = lmat, 
                                                post.load    = post.load),
                         list(var.exp    = var.exp,
                              cov.mat    = cov.emp, 
                              post.Sigma = cov.est))
    result[[g]]  <- unlist(results, recursive=F)
    attr(result[[g]], "Store") <- n.store
  }
  names(result)  <- paste0("Group", seq_len(G))
  
  errors         <- list(MSE = mean(MSE), RMSE = mean(RMSE), NRMSE = mean(NRMSE),
                         CVRMSE = mean(CVRMSE), MAD = mean(MAD))

  attr(GQ.res, "Factors")      <- n.fac
  attr(GQ.res, "Groups")       <- n.grp
  attr(GQ.res, "Supplied")     <- c(Q=Q.T, G=G.T)
  
  if(sw["f.sw"])   {
    scores       <- list(f = f, post.f = rowMeans(f, dims=2))
  }    
  
  result         <- c(result, if(exists("cluster", envir=environment())) list(Clust = cluster), 
                      list(Error = errors, GQ.results = GQ.res, Scores = scores))
  class(result)                <- "IMIFA"
  attr(result, "Method")       <- attr(sims, "Method")
  attr(result, "Obs")          <- n.obs
  attr(result, "Switch")       <- sw
  attr(result, "Vars")         <- n.var
  return(result)
}