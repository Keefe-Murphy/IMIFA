#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.sims     <- function(sims=NULL, burnin=0, thinning=1, 
                          Q=NULL, Q.meth=NULL, recomp=F, ...) {
  if(missing(sims))             stop("Simulations must be supplied")
  if(!exists(as.character(match.call()$sims),
             envir=.GlobalEnv)) stop(paste0("Object ", match.call()$sims, " not found"))
  if(class(sims) != "IMIFA")    stop(paste0("Simulations object of class 'IMIFA' must be supplied"))
  store       <- seq(from=burnin + 1, to=attr(sims, "Store"), by=thinning)
  n.store     <- length(store)
  if(n.store  <= 1)             stop(paste0("burnin must be less than the stored number of iterations"))
  method      <- attr(sims, "Method")
  n.fac       <- attr(sims, "Factors")
  sw          <- attr(sims, "Switch")
  if(!is.logical(recomp))       stop("recomp must be TRUE or FALSE")
  if(thinning  > 1    ||
     burnin    > 0)   recomp <- T
  
  if(!missing(Q)) {
    if(Q > max(n.fac))          stop("Q cannot be greater than the number of factors in sim")
    if(method == 'FA' && length(n.fac) > 1) { 
      Q.ind   <- which(n.fac == Q) 
    }
  } else {
    Q.ind     <- 1
  }
  if(method   == 'FA' && length(n.fac) == 1) {
    Q         <- n.fac
    Q.res     <- list(Q = Q)
  } else if(method == 'IFA') {
    if(missing(Q.meth)) {
      Q.meth  <- "Mode"
    } else {
      if(!is.element(Q.meth, 
       c("Mode", "Median")))    stop("Q.meth must be MODE or MEDIAN")
    }
    
  # Retrieve distribution of Q, tabulate & plot
    Q.store   <- sims[[Q.ind]]$Q.store[store]
    Q.tab     <- table(Q.store, dnn=NULL)
    Q.prob    <- prop.table(Q.tab)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode    <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median  <- ceiling(median(Q.store) * 2)/2
    if(Q.meth == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- Q.median
    Q.CI      <- round(quantile(Q.store, c(0.025, 0.975)))
    Q.res     <- list(Q = Q, Mode = Q.mode, Median = Q.median, 
                      CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
  } else {    
      
  # Initialise
    Q.range   <- n.fac - min(n.fac) + 1
    P         <- length(sims[[1]]$psi[,1])
    prop.exp  <- rep(NA, length(n.fac))
  
  # Calculate Proportion of Variation Explained
    if(sw["l.sw"]) {
      temp.b  <- max(1, burnin)
      for(Q in Q.range) {
        lmat  <- sims[[Q]]$load[,1:n.fac[Q],store, drop=F]
        l.temp      <- as.matrix(sims[[Q]]$load[,1:n.fac[Q],temp.b])
        for(b in 1:n.store) {
          rot       <- procrustes(X=as.matrix(lmat[,,b]), Xstar=l.temp)$R
          lmat[,,b] <- lmat[,,b] %*% rot
        }
        post.load   <- rowMeans(lmat, dims=2)
        prop.exp[Q] <- sum(colSums(post.load * post.load))/nrow(lmat)
      }  
      if(max(prop.exp) > 1)       warning("chain may not have converged")
      Q.ind   <- which.max(prop.exp)
      Q       <- n.fac[Q.ind]
      Q.res   <- list(Q = Q, prop.exp = prop.exp)
    }
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

  if(sw["mu.sw"]) post.mu     <- rowMeans(mu, 1)
  if(sw["f.sw"])  post.f      <- rowMeans(f, dims=2)
  if(sw["p.sw"])  post.psi    <- rowMeans(psi, 1)
  if(sw["l.sw"])  {
    post.load <- rowMeans(lmat, dims=2)
    SS.load   <- colSums(post.load * post.load)
    comm      <- sum(SS.load)
    prop.var  <- SS.load/nrow(post.load)
    cum.var   <- cumsum(prop.var)          
    prop.exp  <- comm/nrow(post.load)
    prop.uni  <- 1 - prop.exp
  }
  data        <- as.data.frame(get(attr(sims, "Name")))
  data        <- data[sapply(data, is.numeric)]
  data        <- scale(data, center=attr(sims, "Center"), scale=attr(sims, "Scaling"))
  cov.empir   <- cov(data)
  cov.estim   <- sims[[Q.ind]]$post.Sigma
  if(recomp) {
    if(!all(c(sw["l.sw"], 
              sw["p.sw"]))) {    cat(paste0("Loadings and/or uniquenesses not stored: can't recompute Sigma", "\n"))
    } else {
    cov.estim <- replace(cov.estim, is.numeric(cov.estim), 0)
      for(r in 1:n.store) {
        Sigma       <- tcrossprod(lmat[,,r]) + diag(psi[,r])
        cov.estim   <- cov.estim + Sigma/n.store
      } 
    }
  }
  error       <- cov.empir - cov.estim
  MSE         <- mean(error * error)
  MAD         <- mean(abs(error))
  error       <- list(MSE=MSE, MAD=MAD)
  if(sum(round(diag(cov.estim))  !=
         round(diag(cov.empir))) != 0
  || sum(abs(post.psi - (1 - post.psi)) < 0) != 0
  || prop.exp  > 1)             warning("chain may not have converged")

  class(post.load)        <- "loadings"
  results     <- list(if(sw["mu.sw"]) list(means=mu, post.mu=post.mu),
                      if(sw["f.sw"])  list(scores=f, post.f=post.f),
                      if(sw["p.sw"])  list(uniquenesses=psi, post.psi=post.psi),
                      if(sw["l.sw"])  list(loadings=lmat, post.load=post.load,  
                                           communality=comm, SS.load=SS.load,
                                           prop.exp=prop.exp, prop.uni=prop.uni,
                                           prop.var=prop.var, cum.var=cum.var),
                      list(cov.mat=cov.empir, post.Sigma = cov.estim, error=error))
  results     <- unlist(results, recursive=F)
  attr(Q.res, "Factors")  <- n.fac 
  results     <- c(results, Q.results = list(Q.res))
  class(results)          <- "IMIFA"
  attr(results, "Method") <- attr(sims, "Method")
  attr(results, "Store")  <- n.store
  attr(results, "Switch") <- attr(sims, "Switch")
  return(results)
}