#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.sims     <- function(sims=NULL, burnin=1, thinning=1, Q=NULL, Q.meth=NULL, ...) {
  if(missing(sims))               stop("Simulations must be supplied")
  if(!exists(as.character(match.call()$sims),
             envir=.GlobalEnv))   stop(paste0("Object ", match.call()$sims, " not found"))
  if(!missing(burnin))    burnin <- burnin + 1
  store       <- seq(from=burnin + 1, to=sims[[1]]$n.store, by=thinning)
  method      <- attr(sim, "Method")
  
  if(!missing(Q)) {
    if(Q > attr(sims, "Factors")) stop("Q cannot be greater than the number of factors in sim")
    if(method == 'FA' && length(range.Q) > 1) { 
      
      Q.ind   <- which(range.Q == Q) 
      
    } else {
      
      Q.ind   <- 1
      
    }
  } else if(method == 'FA' && length(range.Q) == 1) {
    
    Q         <- range.Q
    Q.ind     <- 1
      
  } else if(method == 'IFA') {
    
    if(missing(Q.meth)) {
      Q.meth  <- "Mode"
    } else {
      if(!is.element(Q.meth, 
       c("Mode", "Median")))      stop("Q.meth must be MODE or MEDIAN")
    }
    Q.ind     <- 1 
    
  # Retrieve distribution of Q, tabulate & plot
    Q.store   <- sims[[Q.ind]]$Q.store[store]
    Q.tab     <- table(Q.store, dnn=NULL)
    Q.prob    <- prop.table(Q.tab)
    Q.plot    <- barplot(Q.tab, main="Posterior Distribution of Q", 
                 ylab="Frequency", xlab="Q", xaxt="n")
    axis(1, at=Q.plot, labels=names(Q.tab), tick=F)
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode    <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median  <- ceiling(median(Q.store))
    if(Q.meth == 'Mode') { 
       Q      <- min(Q.mode)
    } else Q  <- Q.median
      
    Q.CI      <- quantile(Q.store, c(0.025, 0.975))
    res.bar   <- list(Mode = Q.mode, Median = Q.median, 
                      CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
    print(unlist(list(Q=Q, res.bar, 
                      Warning="But the user should choose Q based on the attached bar plot!"), recursive=F))
    } else {    
      
  # Initialise
    if(!exists("range.Q")) {
      range.Q <- attr(sims, "Factors")
    }
    Q.range   <- range.Q - min(range.Q) + 1
    P         <- length(sims[[1]]$psi[,1])
    prop.exp  <- rep(NA, length(range.Q))
  
  # Calculate Proportion of Variation Explained
    for(Q in Q.range) {
      lmat    <- sims[[Q]]$load[,1:range.Q[Q],store, drop=F]
      l.temp  <- as.matrix(sims[[Q]]$load[,1:range.Q[Q],burnin])
      for(b in 1:length(store)) {
        rot       <- procrustes(X=as.matrix(lmat[,,b]), Xstar=l.temp)$R
        lmat[,,b] <- lmat[,,b] %*% rot
      }
      post.load   <- apply(lmat, c(1,2), mean)
      prop.exp[Q] <- sum(colSums(post.load * post.load))/P
    }
      
  # Produce Scree Plot & choose optimum Q
    plot(prop.exp, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
         ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
    axis(1, at=1:length(prop.exp), labels=range.Q)
    axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8)
    Q.ind <- which.max(prop.exp)
    Q     <- range.Q[Q.ind]
    points(x=Q.ind, y=prop.exp[Q.ind], col="red", bg="red", pch=21)
    print(list(Q=Q, Warning="But the user should choose Q based on the attached scree plot!"))
  }
  
  mu      <- sims[[Q.ind]]$mu[,store]                            
  f       <- sims[[Q.ind]]$f[,1:Q,store, drop=F]
  lmat    <- sims[[Q.ind]]$load[,1:Q,store, drop=F]
  psi     <- sims[[Q.ind]]$psi[,store]
    
# Loadings matrix / identifiability / # etc.
  l.temp  <- as.matrix(sims[[Q.ind]]$load[,1:Q,burnin])
  for(b in 1:length(store)) {
    rot       <- procrustes(X=as.matrix(lmat[,,b]), Xstar=l.temp)$R
    lmat[,,b] <- lmat[,,b] %*% rot
    f[,,b]    <- t(f[,,b]  %*% rot)
  }

  post.mu     <- apply(mu, 1, mean)
  post.f      <- apply(f, c(1,2), mean)
  post.load   <- apply(lmat, c(1,2), mean)
  post.psi    <- apply(psi, 1, mean)
        
  SS.load     <- colSums(post.load * post.load)
  communality <- sum(SS.load)
  prop.var    <- SS.load/P
  cum.var     <- cumsum(prop.var)
  prop.exp    <- communality/P
  prop.uni    <- 1 - prop.exp

  results     <- list(means = mu, scores = f, loadings = lmat, uniquenesses = psi,
                      post.mu = post.mu, post.f = post.f, post.load = post.load, post.psi = post.psi,
                      store = store, SS.load = SS.load, communality = communality, 
                      prop.var = prop.var, cum.var = cum.var, 
                      prop.exp = prop.exp, prop.uni = prop.uni, Q = Q)
  if(method   == "IFA") {
    results   <- unlist(list(results, res.bar), recursive=F)
  }
  return(results)
};tune.sims   <- cmpfun(tune.sims)

source(paste(dataDirectory, "/IMIFA-GIT/PlottingFunctions.R", sep=""))