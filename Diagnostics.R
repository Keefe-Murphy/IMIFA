#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.params   <- function(burnin=1, thinning=1, Q=NULL, ...) {
  
  if(!missing("burnin"))    burnin <- burnin + 1
  store       <- seq(from=burnin + 1, to=sim[[1]]$n.store, by=thinning)
  #method     <- "IFA"
  
  if(!missing(Q)) {
    if(method == 'FA') { 
      
      Q.ind   <- which(range.Q == Q) 
      
    } else {
      
      Q.ind   <- 1
      
    }
  } else if(method == 'FA' && length(range.Q) == 1) {
    
    Q         <- range.Q
    rm("range.Q", envir=.GlobalEnv)
    
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
    Q.median  <- median(Q.store)
    if(post.Q == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- min(Q.median)
      
    Q.CI      <- quantile(Q.store, c(0.025, 0.975))
    print(list(Q=Q, Mode = Q.mode, Median = Q.median, 
               Credible_Interval = Q.CI, Probs = Q.prob, Counts = Q.tab,
               Warning="But the user should choose Q based on the attached bar plot!"))
    return(list(Mode = Q.mode, Median = Q.median, CI = Q.CI, Probs= Q.prob, Counts = Q.tab))
  } else {    
      
  # Initialise
    if(!exists("range.Q")) {
      range.Q <- attr(sim, "Factors")
    }
    Q.star    <- range.Q - min(range.Q) + 1
    P         <- length(sim[[1]]$psi[,1])
    prop.exp  <- rep(NA, length(range.Q))
  
  # Calculate Proportion of Variation Explained
    for(Q in Q.star) {
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

  assign("mu",   mu,   envir=.GlobalEnv)
  assign("f",    f,    envir=.GlobalEnv)
  assign("load", load, envir=.GlobalEnv)
  assign("psi",  psi,  envir=.GlobalEnv)
};tune.params <- cmpfun(tune.params)