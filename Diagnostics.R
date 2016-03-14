#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.imifa       <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL,
                             Q = NULL, Q.meth = c("Mode", "Median"), recomp = F, ...) {
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
  if(!is.logical(recomp))         stop("recomp must be TRUE or FALSE")
  if(any(burnin   > 0, 
     thinning     > 1)) recomp <- T
  
  if(!missing(G)) {
    G.x          <- G
    G.xind       <- which(n.grp == G.x)
    if(all(is.element(method, c("MFA", "MIFA")),
       !is.element(G, n.grp)))    stop("This G value was not used during simulation")
  } 
  G.T            <- exists("G.x", envir=environment())
  if(!missing(Q)) {
    Q.x          <- Q
    Q.xind       <- which(n.fac == Q.x)
    if(all(is.element(method, c("FA", "MFA")),
       !is.element(Q, n.fac)))    stop("This Q value was not used during simulation")
    if(all(is.element(method, c("IFA", "classify")), 
      (Q * (n.fac - Q)) < 0))     stop(paste0("Q cannot be greater than the number of factors in ", match.call()$sims))
  } 
  Q.T            <- exists("Q.x", envir=environment())
  G.ind          <-     G      <- 1
  Q.ind          <- 1
  
  if(method == "IFA") {
    if(missing(Q.meth)) {
      Q.meth     <- "Mode"
    } else {
      Q.meth     <- match.arg(Q.meth)
    }
    
  # Retrieve distribution of Q, tabulate & plot
    Q.store      <- sims[[G.ind]][[Q.ind]]$Q.store[store]
    Q.tab        <- table(Q.store, dnn=NULL)
    Q.prob       <- prop.table(Q.tab)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode       <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.med        <- ceiling(median(Q.store) * 2)/2
    if(Q.T) {
      Q          <- Q.x
    } else if(Q.meth == 'Mode') { 
      Q          <- min(Q.mode)
    } else {
      Q          <- Q.med
    }
    Q.CI         <- round(quantile(Q.store, c(0.025, 0.975)))
    GQ.res       <- list(G = G, Q = Q, Mode = Q.mode, Median = Q.med, 
                         CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
  } else {
    
  # Retrieve BIC to tune G & Q   
    G.range      <- length(n.grp)
    Q.range      <- length(n.fac)
    bic          <- matrix(NA, nr=G.range, nc=Q.range, dimnames=list(paste0("G", n.grp), paste0("Q", n.fac)))
    for(g in 1:G.range) { 
      for(q in 1:Q.range) {
        bic[g,q] <- sim[[g]][[q]]$bic  
      }  
    }
    bic.max      <- which(bic == max(bic), arr.ind = T)
    G.ind        <- bic.max[1]
    Q.ind        <- bic.max[2]
    G            <- n.grp[G.ind]
    Q            <- n.fac[Q.ind]
    if(all(Q.T, G.T)) {
      bic        <- bic[G.xind,Q.xind]
      G          <- G.x
      Q          <- Q.x
      n.grp      <- G
      n.fac      <- Q  
      G.ind      <- G.xind
      Q.ind      <- Q.xind
    } else if(Q.T) {
      bic        <- bic[G.ind,Q.xind]
      Q          <- Q.x
      n.fac      <- Q
      Q.ind      <- Q.xind
    } else if(G.T) {
      bic        <- bic[G.xind,Q.ind]
      G          <- G.x
      n.grp      <- G
      G.ind      <- G.xind
    }
    Q            <- setNames(rep(Q, G), paste0("Qg", 1:G))
    GQ.res       <- list(G = G, Q = Q, BIC = bic)
  } 
  
  result         <- list(list())
  G.ind          <- which(n.grp == G)
  temp.b         <- max(1, burnin)
  MSE  <-  RMSE  <-  NRMSE  <-  CVRMSE  <-  MAD  <- rep(NA, G)
  for(g in 1:G) {
    Qg           <- Q[g]
    if(Qg == 0) {
      sw[c("f.sw", "l.sw")]    <- F
    } else {
      sw[c("f.sw", "l.sw")]    <- attr(sims, "Switch")[c("f.sw", "l.sw")]
    }
    
    if(method == "MFA") {
      if(sw["f.sw"]) {
        f        <- adrop(sims[[G.ind]][[Q.ind]]$f[,1:Qg,g,store, drop=F], drop=3)
      }
      if(sw["l.sw"]) {
        lmat     <- adrop(sims[[G.ind]][[Q.ind]]$load[,1:Qg,g,store, drop=F], drop=3)
        l.temp   <- adrop(sims[[G.ind]][[Q.ind]]$load[,1:Qg,g,temp.b, drop=F], drop=3:4)
      }
    }
    
    if(method == "FA")  {
      if(sw["f.sw"]) {
        f        <- sims[[G.ind]][[Q.ind]]$f[,1:Qg,store, drop=F]
      }
      if(sw["l.sw"]) {
        lmat     <- sims[[G.ind]][[Q.ind]]$load[,1:Qg,store, drop=F]
        l.temp   <- adrop(sims[[G.ind]][[Q.ind]]$load[,1:Qg,temp.b, drop=F], drop=3)
      }
    } 
    
    if(method == "IFA") {
      store      <- store[which(Q.store >= Qg)]
      n.store    <- length(store)
      temp.b     <- store[1]
      if(sw["f.sw"]) {
        f        <- as.array(sims[[G.ind]][[Q.ind]]$f)[,1:Qg,store, drop=F]
      }
      if(sw["l.sw"]) {
        lmat     <- as.array(sims[[G.ind]][[Q.ind]]$load)[,1:Qg,store, drop=F]
        l.temp   <- adrop(as.array(sims[[G.ind]][[Q.ind]]$load)[,1:Qg,temp.b, drop=F], drop=3)
      }
    }
    
    if(is.element(method, c("MFA", "MIFA", "IMIFA"))) {
      cov.emp    <- sims[[G.ind]][[Q.ind]]$cov.mat[,,g]
      cov.est    <- sims[[G.ind]][[Q.ind]]$post.Sigma[,,g]
      post.mu    <- sims[[G.ind]][[Q.ind]]$post.mu[,g]
      post.psi   <- sims[[G.ind]][[Q.ind]]$post.psi[,g]
      if(sw["mu.sw"])  {
        mu       <- sims[[G.ind]][[Q.ind]]$mu[,g,store]                            
      }
      if(sw["p.sw"])   {
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
      if(sw["p.sw"])   {
        psi      <- sims[[G.ind]][[Q.ind]]$psi[,store]
      }
    }
    
  # Loadings matrix / identifiability / error metrics / etc.  
    if(sw["l.sw"])   {
      for(p in 1:n.store) {
        rot           <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
        lmat[,,p]     <- lmat[,,p] %*% rot
        if(sw["f.sw"]) {
          f[,,p]      <- f[,,p]    %*% rot
        }  
      }
    }
    
    if(sw["mu.sw"])  post.mu   <- rowMeans(mu, dims=1)
    if(sw["f.sw"])   post.f    <- rowMeans(f, dims=2)
    if(sw["p.sw"])   post.psi  <- rowMeans(psi, dims=1)
    if(sw["l.sw"]) { post.load <- rowMeans(lmat, dims=2)
      var.exp    <- sum(colSums(post.load * post.load))/n.var
    } else {
      var.exp    <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
    }
    if(all(recomp, sw[c("l.sw", "p.sw")])) {
      cov.est    <- replace(cov.est, is.numeric(cov.est), 0)
      for(r in 1:n.store) {
        Sigma    <- tcrossprod(lmat[,,r]) + diag(psi[,r])
        cov.est  <- cov.est + Sigma/n.store
      }
    } else if(recomp) {
      if(!sw["l.sw"]) {
                                  warning("Loadings not stored: can't re-estimate Sigma", call.=F)
      } else {
                                  warning("Uniquenesses not stored: can't re-estimate Sigma", call.=F)
      }
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
    results      <- list(if(sw["mu.sw"])   list(means  = mu), 
                         list(post.mu    = post.mu),
                         if(sw["f.sw"])    list(scores = f, 
                                                post.f = post.f),
                         if(sw["p.sw"])    list(uniquenesses = psi), 
                         list(post.psi   = post.psi),
                         if(sw["l.sw"])    list(loadings     = lmat, 
                                                post.load    = post.load),
                         list(var.exp    = var.exp,
                              cov.mat    = cov.emp, 
                              post.Sigma = cov.est))
    result[[g]]  <- unlist(results, recursive=F)
  }
  names(result)  <- paste0("Group", 1:G)
  if(Q.T) {
    attr(GQ.res,
         "Factors")   <- Q.x
  } else  {
    attr(GQ.res, 
         "Factors")   <- n.fac
  } 
  if(G.T) {
    attr(GQ.res,
         "Groups")    <- G.x
  } else  {
    attr(GQ.res, 
         "Groups")    <- n.grp
  }
  attr(GQ.res, "Supplied")     <- c(Q=Q.T, G=G.T)
  errors         <- list(MSE = mean(MSE), RMSE = mean(RMSE), NRMSE = mean(NRMSE),
                         CVRMSE = mean(CVRMSE), MAD = mean(MAD))
  result         <- c(result, Error = list(errors), GQ.results = list(GQ.res))
  class(result)                <- "IMIFA"
  attr(result, "Method")       <- attr(sims, "Method")
  attr(result, "Obs")          <- n.obs
  attr(result, "Store")        <- n.store
  attr(result, "Switch")       <- sw
  attr(result, "Vars")         <- n.var
  return(result)
}