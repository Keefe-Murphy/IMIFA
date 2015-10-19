# BAYESIAN FACTOR ANALYSIS

# preamble
  dataDirectory <- "C:/Users/Windows/Documents/Claire IMIFA/IMIFA-GIT"
  setwd(dataDirectory)
  set.seed(21092015)
  library(MASS)
  library(compiler)
  library(pgmm)
  library(MCMCpack)

# read in the data
  data(wine); labels <- wine[,1]; wine <- wine[,-1]
  #subjectmarks <- read.csv(paste(dataDirectory,"/Data/","SubjectMarks.csv",sep=""))
  #cereal       <- read.csv(paste(dataDirectory,"/Data/","Cereal.csv",sep=""))

# simulate data
  P <- 50
  Q <- 10
  N <- 1500
  mu.true   <- mvrnorm(mu=rep(0,P), Sigma=diag(P));      names(mu.true)      <- c(1:P)
  f.true    <- mvrnorm(n=N, mu=rep(0,Q), Sigma=diag(Q)); colnames(f.true)    <- paste("Factor", 1:Q); rownames(f.true)    <- c(1:N)
  load.true <- mvrnorm(n=P, mu=rep(0,Q), Sigma=diag(Q)); colnames(load.true) <- paste("Factor", 1:Q); rownames(load.true) <- c(1:P)
  psi.true  <- abs(rnorm(P)*P); names(psi.true) <- c(1:P)
  eps.true  <- mvrnorm(n=N, mu=rep(0,P), Sigma=diag(psi.true))
  data      <- matrix(0,nr=N,nc=P)
  for (i in 1:N) {  
    data[i, ]      <- mu.true + load.true%*%f.true[i,] + eps.true[i,]
  }; rownames(data) <- c(1:N); colnames(data) <- c(1:P)

  save(mu.true, f.true, load.true, psi.true, eps.true, file="/home/kmurphy/Simulated Data.Rdata")

# define full conditional functions
  # means
    sim.mu      <- function(mu.sigma, N, P, psi, data, f, load, iter, ...) {
      mu.omega  <- solve(solve(mu.sigma) + N * solve(diag(psi[,iter-1])))
      mvrnorm(mu=rep(0, P), Sigma=mu.omega) +
        mu.omega %*% solve(diag(psi[,iter-1])) %*% t(t(apply(data, 2, sum)) - t(apply(f[,,iter-1], 2, sum)) %*% t(load[,,iter-1])) }
    sim.mu      <- cmpfun(sim.mu)
  
  # scores
    sim.omega.f <- function(Q, load, psi, iter, ...) {
      solve(diag(Q) + t(load[,,iter-1]) %*% solve(diag(psi[,iter-1])) %*% load[,,iter-1]) }
    sim.omega.f <- cmpfun(sim.omega.f)
    
    sim.scores  <- function(Q, f.omega, load, psi, data, mu, i, iter, ...) {
      mvrnorm(mu=rep(0, Q), Sigma=f.omega) + 
        (f.omega %*% t(load[,,iter-1]) %*% solve(diag(psi[,iter-1]))) %*% (data[i,] - mu[,iter]) }
    sim.scores  <- cmpfun(sim.scores)
  
  # loadings
    sim.load    <- function(l.sigma, Q, f, psi, data, mu, j, iter, ...) {
      l.omega   <- solve(solve(l.sigma) + (1/psi[j,iter-1]) * t(f[,,iter]) %*% f[,,iter])
      t(mvrnorm(mu=rep(0, Q), Sigma=l.omega) + 
        (l.omega %*% t(f[,,iter]) * (1/psi[j,iter-1])) %*% (data[,j] - mu[j,iter])) }
    sim.load    <- cmpfun(sim.load)
    
  # uniquenesses
    sim.psi     <- function(N, psi.alpha, psi.beta, data, mu, load, f, j, iter, ...) {
      rinvgamma(1, shape=(N+psi.alpha)/2, 
        scale=(sum(data[,j] - mu[j,iter] - load[j,,iter]%*%t(f[,,iter]))^2 + psi.beta)/2) }
    sim.psi     <- cmpfun(sim.psi)
  
# gibbs sampler function
  gibbs  <- function(data, n.iters=100000, Q, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5, scaling=T, ...) {
    
  # centre the data (optional)
    if (scaling) {data <- scale(data, center=T, scale=F)} else  {data <- as.matrix(data)}
    
  # define & initialise variables
    N         <- nrow(data)
    P         <- ncol(data)
    if (Q>=P) stop ("Number of factors must be less than the number of variables")
    mu        <- matrix(NA, nr=P, nc=n.iters);    rownames(mu)   <- colnames(data)
    f         <- array(NA, dim=c(N, Q, n.iters)); colnames(f)    <- paste("Factor",1:Q)
    load      <- array(NA, dim=c(P, Q, n.iters)); rownames(load) <- colnames(data); colnames(load) <- paste("Factor",1:Q)
    psi       <- matrix(NA, nr=P, nc=n.iters);    rownames(psi)  <- colnames(data)
    mu.sigma  <- sigma.mu * diag(P)
    l.sigma   <- sigma.l  * diag(Q)
    mu.omega  <- matrix(NA, P, P)
    f.omega   <- matrix(NA, Q, Q)
    l.omega   <- matrix(NA, Q, Q)
    mu[,1]    <- mvrnorm(mu=rep(0, P), Sigma=mu.sigma)
    f[,,1]    <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
    load[,,1] <- mvrnorm(n=P, mu=rep(0, Q), Sigma=l.sigma)
    psi[,1]   <- rinvgamma(P, shape=psi.alpha/2, scale=psi.beta/2)
    
  # iterate
    for(iter in 2:n.iters) { 
      mu[,iter]       <- sim.mu(mu.sigma, N, P, psi, data, f, load, iter, ...)
      f.omega         <- sim.omega.f(Q, load, psi, iter)
      for (i in 1:N)    {
        f[i,,iter]    <- sim.scores(Q, f.omega, load, psi, data, mu, i, iter)
      }
      for (j in 1:P) {
        load[j,,iter] <- sim.load(l.sigma, Q, f, psi, data, mu, j, iter, ...)
        psi[j,iter]   <- sim.psi(N, psi.alpha, psi.beta, data, mu, load, f, j, iter)
      }
    }
    return(list(mu   = mu,
                f    = f, 
                load = load, 
                psi  = psi))
  }; gibbs.comp <- cmpfun(gibbs)

# run the gibbs sampler
  n.iters <- 100000
  sim     <- gibbs.comp(data=wine, n.iters, Q=3, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)
  #system.time(gibbs(data=data, n.iters, Q=2, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)
  #system.time(gibbs.comp(data=data, n.iters, Q=2, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)

# save / load results
  save(sim,file="/home/kmurphy/Wine Simulations.Rdata") # in server, tick box, export
  load(file="Simulations/Wine Simulations.Rdata",envir=.GlobalEnv)

# convergence diagnostics
  burnin <- 10000
  thin   <- 3
  mu     <- sim$mu[,seq(from=burnin+1, to=n.iters, by=thin)]
  f      <- sim$f[,,seq(from=burnin+1, to=n.iters, by=thin)]
  load   <- sim$load[,,seq(from=burnin+1, to=n.iters, by=thin)]
  psi    <- sim$psi[,seq(from=burnin+1, to=n.iters, by=thin)]

# NB: You can check your answer by plotting
#     the f of a 2-factor model to the
#     wine dataset. Expect to see a horseshoe.

# loadings matrix / identifiability / # etc.
  l.temp <- sim$load[,,burnin]
  for(b in 1:dim(load)[3]){
    rot       <- procrustes(X=load[,,b], Xstar=l.temp)$R
    load[,,b] <- load[,,b]%*%rot
    f[,,b]    <- t(t(rot)%*%t(f[,,b]))
  }

# plots & posterior summaries etc.
  # means
    plot(mu[1,], type="l")
    matplot(t(mu[,]), type="l")
    post.mu <- apply(mu, 1, mean)
    plot(post.mu, type="n")
    text(x=1:length(post.mu), y=post.mu, names(post.mu))
    acf(mu[1,])
  
  # scores
    plot(f[1,1,], type="l")
    matplot(t(f[1,,]), type="l")
    post.f <- apply(f, c(1,2), mean)
    plot(post.f, type="n")
    text(post.f[,1], post.f[,2], 1:nrow(post.f), col=labels)
    plot(f[,,dim(f)[3]], type="n")
    text(f[,1,dim(f)[3]], f[,2,dim(f)[3]], 1:nrow(post.f), col=labels)
    acf(f[1,1,])
  
  # uniquenesses
    plot(psi[1,], type="l")
    matplot(t(psi[,]), type="l")
    post.psi <- apply(psi, 1, mean)
    plot(post.psi, type="n")
    text(1:length(post.psi), post.psi, names(post.psi))
    acf(psi[1,])
  
  # loadings
    plot(load[1,1,], type="l")
    matplot(t(load[1,,]), type="l")
    post.load <- apply(load, c(1,2), mean)
    plot(post.load, type="n")
    text(post.load[,1], post.load[,2], rownames(post.load))
    acf(load[1,1,])
  
  # summaries 
    P <- nrow(post.load)
    sum(post.psi)/P # % of variance which is unique
    communality <- P - sum(post.psi)
    communality/P   # % of variance 
    #sum(post.load[,]^2)