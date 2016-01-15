#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

invisible(ifelse(!exists("burnin"), 
                  assign("burnin", 1), 
                  assign("burnin", burnin + 1)))
if(!exists("thinning")) assign("thinning", 1)
if(!exists("store"))    assign("store", seq(from=burnin + 1, to=sim[[1]]$n.store, by=thinning))

if(method == 'FA' && length(range.Q) == 1) {
  assign("Q", range.Q)
  rm("range.Q")
} else if(method == 'IFA') {
  assign("post.Q", 'Mode')
  assign("Q.ind", 1) 
  tune.parameters  <- function(...) { 
  # Retrieve distribution of Q, tabulate & plot
    assign("Q.store", sim[[Q.ind]]$Q.store[store])
    assign("Q.tab",   table(Q.store, dnn=NULL))
    assign("Q.prob",  prop.table(Q.tab))
    Q.plot     <- barplot(Q.tab, main="Posterior Distribution of Q", 
                          ylab="Frequency", xlab="Q", xaxt="n")
    axis(1, at=Q.plot, labels=names(Q.tab), tick=F)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode     <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median   <- median(Q.store)
    if(post.Q == 'Mode') { 
       assign("Q", min(Q.mode), envir=.GlobalEnv)
    } else { assign("Q", Q.median) }
      Q.CI     <- quantile(Q.store, c(0.025, 0.975))
      print(list(Q=Q, Mode = Q.mode, Median = Q.median, 
                 Credible_Interval = Q.CI, Probabilities = Q.prob, Counts = Q.tab,
                 Warning="But the user should choose Q based on the attached bar plot!"))
    }
  } else {    
    tune.parameters  <- function(...) {
  # Initialise
    if(!exists("range.Q")) {
      range.Q  <- attr(sim, "Factors")
    }
    Q.star     <- range.Q - min(range.Q) + 1
    P          <- length(sim[[1]]$psi[,1])
    prop.exp   <- rep(NA, length(range.Q))
  
  # Calculate Proportion of Variation Explained
    for(Q in Q.star) {
      load        <- sim[[Q]]$load[,1:range.Q[Q],store, drop=F]
      l.temp      <- as.matrix(sim[[Q]]$load[,1:range.Q[Q],burnin])
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
};tune.parameters <- cmpfun(tune.parameters)

if(method == 'FA') { 
  new.Q  <- function(Q, ...) {
    assign("Q", Q)
    assign("Q.ind", which(range.Q == Q)) 
  }
} else {
  new.Q  <- function(Q, ...) {
    assign("Q", Q)
    assign("Q.ind", 1)
  }
}; new.Q <- cmpfun(new.Q)

extract.results <- function(Q) {
  mu       <- sim[[Q.ind]]$mu[,store]                            
  f        <- sim[[Q.ind]]$f[,1:Q,store, drop=F]
  load     <- sim[[Q.ind]]$load[,1:Q,store, drop=F]
  psi      <- sim[[Q.ind]]$psi[,store]
  
  # Loadings matrix / identifiability / # etc.
    l.temp <- as.matrix(sim[[Q.ind]]$load[,1:Q,burnin])
    for(b in 1:length(store)) {
      rot       <- procrustes(X=as.matrix(load[,,b]), Xstar=l.temp)$R
      load[,,b] <- load[,,b] %*% rot
      f[,,b]    <- t(f[,,b]  %*% rot)
    }
  
  assign("mu",   mu)
  assign("f",    f)
  assign("load", load)
  assign("psi",  psi)
}; extract.results <- cmpfun(extract.results)