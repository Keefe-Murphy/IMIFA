#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.sims     <- function(sims=NULL, burnin=1, thinning=1, Q=NULL, Q.meth=NULL, ...) {
  if(missing(sims))             stop("Simulations must be supplied")
  if(!exists(as.character(match.call()$sims),
             envir=.GlobalEnv)) stop(paste0("Object ", match.call()$sims, " not found"))
  if(class(sims) != "IMIFA")    stop(paste0("Simulations object of class 'IMIFA' must be supplied"))
  if(!missing(burnin))          burnin <- burnin + 1
  store       <- seq(from=burnin, to=attr(sims, "Store"), by=thinning)
  method      <- attr(sims, "Method")
  n.fac       <- attr(sims, "Factors")
  
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
    Q.plot    <- barplot(Q.tab, main="Posterior Distribution of Q", 
                 ylab="Frequency", xlab="Q", xaxt="n")
    axis(1, at=Q.plot, labels=names(Q.tab), tick=F)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode    <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median  <- ceiling(median(Q.store) * 2)/2
    if(Q.meth == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- Q.median
    Q.CI      <- round(quantile(Q.store, c(0.025, 0.975)))
    res.bar   <- list(Mode = Q.mode, Median = Q.median, 
                      CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
    print(unlist(list(Q=Q, res.bar), recursive=F))
    cat(paste0("Warning: the user should choose Q based on the attached bar plot!", "\n"))
  } else {    
      
  # Initialise
    Q.range   <- n.fac - min(n.fac) + 1
    P         <- length(sims[[1]]$psi[,1])
    prop.exp  <- rep(NA, length(n.fac))
  
  # Calculate Proportion of Variation Explained
    for(Q in Q.range) {
      lmat    <- sims[[Q]]$load[,1:n.fac[Q],store, drop=F]
      l.temp  <- as.matrix(sims[[Q]]$load[,1:n.fac[Q],burnin])
      for(b in 1:length(store)) {
        rot       <- procrustes(X=as.matrix(lmat[,,b]), Xstar=l.temp)$R
        lmat[,,b] <- lmat[,,b] %*% rot
      }
      post.load   <- rowMeans(lmat, dims=2)
      prop.exp[Q] <- sum(colSums(post.load * post.load))/nrow(lmat)
    }  
    if(max(prop.exp) > 1)       cat(paste0("Warning: chain may not have converged", "\n"))
      
  # Produce Scree Plot & choose optimum Q
    plot(prop.exp, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
         ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
    axis(1, at=1:length(prop.exp), labels=n.fac)
    axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8)
    Q.ind     <- which.max(prop.exp)
    Q         <- n.fac[Q.ind]
    points(x=Q.ind, y=prop.exp[Q.ind], col="red", bg="red", pch=21)
    cat(paste0("Q = ", Q, "\n"))
    cat(paste0("Warning: the user should choose Q based on the attached scree plot!", "\n"))
  }
  
  if(method   == "FA") {
    store     <- store[-1]
    f         <- sims[[Q.ind]]$f[,1:Q,store, drop=F]
    lmat      <- sims[[Q.ind]]$load[,1:Q,store, drop=F]
    l.temp    <- as.matrix(sims[[Q.ind]]$load[,1:Q,burnin])
  } else {
    store     <- which(Q.store[store] >= Q)
    burnin    <- store[1]
    store     <- store[-1]
    f         <- as.array(sims[[Q.ind]]$f)[,1:Q,store, drop=F]
    lmat      <- as.array(sims[[Q.ind]]$load)[,1:Q,store, drop=F]
    l.temp    <- as.matrix(as.array(sims[[Q.ind]]$load)[,1:Q,burnin])
  }
  mu          <- sims[[Q.ind]]$mu[,store]                            
  psi         <- sims[[Q.ind]]$psi[,store]
  
# Loadings matrix / identifiability / # etc.  
  for(b in 1:length(store)) {
    rot       <- procrustes(X=as.matrix(lmat[,,b]), Xstar=l.temp)$R
    lmat[,,b] <- lmat[,,b] %*% rot
    f[,,b]    <- f[,,b]    %*% rot
  }

  post.mu     <- rowMeans(mu, 1)
  post.f      <- rowMeans(f, dims=2)
  post.load   <- rowMeans(lmat, dims=2)
  post.psi    <- rowMeans(psi, 1)       
        
  SS.load     <- colSums(post.load * post.load)
  communality <- sum(SS.load)
  prop.var    <- SS.load/nrow(post.load)
  cum.var     <- cumsum(prop.var)          
  prop.exp    <- communality/nrow(post.load)
  prop.uni    <- 1 - prop.exp
  if(sum(round(diag(tcrossprod(lmat[,,length(store)]) 
                  + psi[,length(store)])) != 1) != 0
  || prop.exp  > 1)         cat(paste0("Warning: chain may not have converged", "\n"))

  results     <- list(means = mu, scores = f, loadings = lmat, uniquenesses = psi,
                      post.mu = post.mu, post.f = post.f, post.load = post.load, post.psi = post.psi,
                      store = store, SS.load = SS.load, communality = communality, 
                      prop.var = prop.var, cum.var = cum.var, 
                      prop.exp = prop.exp, prop.uni = prop.uni, Q = Q)
  if(method   == "IFA") {
    results   <- unlist(list(results, res.bar), recursive=F)
  }
  class(results)          <- "IMIFA"
  attr(results, "Method") <- attr(sims, "Method")
  return(results)
}