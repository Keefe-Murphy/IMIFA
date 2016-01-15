#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.params   <- function(burnin=1, thinning=1, Q=NULL, ...) {
  
  if(!missing("burnin"))    burnin <- burnin + 1
  store       <- seq(from=burnin + 1, to=sim[[1]]$n.store, by=thinning)
  #method     <- "IFA"
  
  if(!missing(Q)) {
    if(method == 'FA' && length(range.Q) > 1) { 
      
      Q.ind   <- which(range.Q == Q) 
      
    } else {
      
      Q.ind   <- 1
      
    }
  } else if(method == 'FA' && length(range.Q) == 1) {
    
    Q         <- range.Q
      
  } else if(method == 'IFA') {
    
    post.Q    <- "Mode"
    Q.ind     <- 1 
    
  # Retrieve distribution of Q, tabulate & plot
    Q.store   <- sim[[Q.ind]]$Q.store[store]
    Q.tab     <- table(Q.store, dnn=NULL)
    Q.prob    <- prop.table(Q.tab)
    Q.plot    <- barplot(Q.tab, main="Posterior Distribution of Q", 
                 ylab="Frequency", xlab="Q", xaxt="n")
    axis(1, at=Q.plot, labels=names(Q.tab), tick=F)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode    <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median  <- ceiling(median(Q.store))
    if(post.Q == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- Q.median
      
    Q.CI      <- quantile(Q.store, c(0.025, 0.975))
    res.bar   <- list(Mode = Q.mode, Median = Q.median, 
                      CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
    print(list(Q=Q, Warning="But the user should choose Q based on the attached scree plot!"))
    } else {    
      
  # Initialise
    if(!exists("range.Q")) {
      range.Q <- attr(sim, "Factors")
    }
    Q.range   <- range.Q - min(range.Q) + 1
    P         <- length(sim[[1]]$psi[,1])
    prop.exp  <- rep(NA, length(range.Q))
  
  # Calculate Proportion of Variation Explained
    for(Q in Q.range) {
      load    <- sim[[Q]]$load[,1:range.Q[Q],store, drop=F]
      l.temp  <- as.matrix(sim[[Q]]$load[,1:range.Q[Q],burnin])
      for(b in 1:length(store)) {
        rot       <- procrustes(X=as.matrix(load[,,b]), Xstar=l.temp)$R
        load[,,b] <- load[,,b] %*% rot
      }
      post.load   <- apply(load, c(1,2), mean)
      prop.exp[Q] <- sum(colSums(post.load * post.load))/P
    }
      
  # Produce Scree Plot & choose optimum Q
    plot(prop.exp, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
         ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
    axis(1, at=1:length(prop.exp), labels=range.Q)
    axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8)
    assign("Q.ind", which.max(prop.exp))
    assign("Q", range.Q[Q.ind])
    points(x=Q.ind, y=prop.exp[Q.ind], col="red", bg="red", pch=21)
    print(list(Q=Q, Warning="But the user should choose Q based on the attached scree plot!"))
  }
  
  mu     <- sim[[Q.ind]]$mu[,store]                            
  f      <- sim[[Q.ind]]$f[,1:Q,store, drop=F]
  load   <- sim[[Q.ind]]$load[,1:Q,store, drop=F]
  psi    <- sim[[Q.ind]]$psi[,store]
    
# Loadings matrix / identifiability / # etc.
  l.temp <- as.matrix(sim[[Q.ind]]$load[,1:Q,burnin])
  for(b in 1:length(store)) {
    rot       <- procrustes(X=as.matrix(load[,,b]), Xstar=l.temp)$R
    load[,,b] <- load[,,b] %*% rot
    f[,,b]    <- t(f[,,b]  %*% rot)
  }

  assign("mu",        mu,   envir=.GlobalEnv)
  assign("f",         f,    envir=.GlobalEnv)
  assign("load",      load, envir=.GlobalEnv)
  assign("psi",       psi,  envir=.GlobalEnv)
  
  assign("post.mu",   apply(mu, 1, mean),        envir=.GlobalEnv)
  assign("post.f",    apply(f, c(1,2), mean),    envir=.GlobalEnv)
  assign("post.load", apply(load, c(1,2), mean), envir=.GlobalEnv)
  assign("post.psi",  apply(psi, 1, mean),       envir=.GlobalEnv)
        
  SS.load     <- colSums(post.load * post.load)
  communality <- sum(SS.load)
  prop.var    <- SS.load/P
  cum.var     <- cumsum(prop.var)
  prop.exp    <- communality/P
  prop.uni    <- 1 - prop.exp

  results     <- list(SS.load = SS.load, communality = communality, 
                      prop.var = prop.var, cum.var = cum.var, 
                      prop.exp = prop.exp, prop.uni = prop.uni, Q = Q)
  if(method   == "IFA") {
    results   <- lappend(results, res.bar)
  }
  return(results)
};tune.params <- cmpfun(tune.params)

plot.cum.var  <- function(n.fac=res$Q, cum.var=res$cum.var) {
  prop.exp    <- cum.var[n.fac]
  if(n.fac     > 1) {
    plot(cum.var, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
         ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
    axis(1, at=1:length(cum.var), labels=1:n.fac)
    axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8) 
    points(x=n.fac, y=prop.exp, col="red", bg="red", pch=21)
  }
 cat(paste0("Proportion of Variation Explained = ",
            round(prop.exp[length(prop.exp)], 3), "\n"))
}