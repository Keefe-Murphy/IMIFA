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
data(wine); wine.true=wine[,1]; wine=wine[,-1]
#subjectmarks <- read.csv(paste(dataDirectory,"/Data/","SubjectMarks.csv",sep=""))
#cereal <- read.csv(paste(dataDirectory,"/Data/","Cereal.csv",sep=""))

# simulate data
P=8
Q=3
N=150
mu.true           <- mvrnorm(mu=rep(0, P), Sigma=diag(P)); names(mu.true) <- c(1:P)
scores.true       <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q)); colnames(scores.true)   <- paste("Factor",1:Q); rownames(scores.true) <- c(1:N)
loadings.true     <- mvrnorm(n=P, mu=rep(0, Q), Sigma=diag(Q)); colnames(loadings.true) <- paste("Factor",1:Q); rownames(loadings.true) <- c(1:P)
uniquenesses.true <- abs(rnorm(P)*P); names(uniquenesses.true) <- c(1:P)
errors.true       <- mvrnorm(n=N, mu=rep(0, P), Sigma=diag(uniquenesses.true))
data              <- matrix(0,nr=N,nc=P)
for (i in 1:N) {  
  data[i, ]      <- mu.true + loadings.true%*%scores.true[i,] + errors.true[i,]
};rownames(data)<-c(1:N);colnames(data)<-c(1:P)

save(mu.true, scores.true, loadings.true, uniquenesses.true, errors.true, file="/home/kmurphy/Simulated Data.Rdata")

# gibbs sampler function
gibbs  <- function(data, n.iters=10000, Q, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5, scaling=T) {
  
  # centre the data (optional)
  if (scaling) {data <- scale(data, center=T, scale=F)} else  {data <- as.matrix(data)}
  
  # define & initialise variables
  N                <- nrow(data)
  P                <- ncol(data)
  if (Q>=P) stop ("Number of factors must be less than the number of variables")
  mu               <- matrix(NA, nr=P, nc=n.iters);           rownames(mu) <- colnames(data)
  scores           <- array(NA, dim=c(N, Q, n.iters));    colnames(scores) <- paste("Factor",1:Q)
  loadings         <- array(NA, dim=c(P, Q, n.iters));  rownames(loadings) <- colnames(data); colnames(loadings) <- paste("Factor",1:Q)
  uniquenesses     <- matrix(NA, nr=P, nc=n.iters); rownames(uniquenesses) <- colnames(data)
  mu.sigma         <- sigma.mu * diag(P)
  l.sigma          <- sigma.l * diag(Q)
  mu.omega         <- array(NA, dim=c(P, P, n.iters))
  f.omega          <- array(NA, dim=c(Q, Q, n.iters))
  l.omega          <- array(NA, dim=c(Q, Q, P, n.iters))
  mu[,1]           <- mvrnorm(mu=rep(0, P), Sigma=mu.sigma)
  scores[,,1]      <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
  loadings[,,1]    <- mvrnorm(n=P, mu=rep(0, Q), Sigma=l.sigma)
  uniquenesses[,1] <- rinvgamma(P, shape=psi.alpha/2, scale=psi.beta/2)
  
  # iterate
  for(iter in 2:n.iters) { 
    mu.omega[,,iter]        <- solve(solve(mu.sigma) + N*solve(diag(uniquenesses[,iter-1])))
    f.omega[,,iter]         <- solve(diag(Q) + t(loadings[,,iter-1])%*%solve(diag(uniquenesses[,iter-1]))%*%loadings[,,iter-1])
    mu[,iter]               <- (mvrnorm(mu=rep(0, P), Sigma=mu.omega[,,iter])
                                + mu.omega[,,iter]%*%solve(diag(uniquenesses[,iter-1]))%*%t(t(apply(data, 2, sum)) - t(apply(scores[,,iter-1], 2, sum))%*%t(loadings[,,iter-1])))
    for (i in 1:N)    {
      scores[i,,iter]       <- (mvrnorm(mu=rep(0, Q), Sigma=f.omega[,,iter])
                                + (f.omega[,,iter]%*%t(loadings[,,iter-1])%*%solve(diag(uniquenesses[,iter-1])))%*%(data[i,]-mu[,iter]))
    }
    for (j in 1:P) {
      l.omega[,,j,iter]     <- solve(solve(l.sigma) + (1/uniquenesses[j,iter-1])*t(scores[,,iter])%*%scores[,,iter])
      loadings[j,,iter]     <- t(mvrnorm(mu=rep(0, Q), Sigma=l.omega[,,j,iter])
                                 + (l.omega[,,j,iter]%*%t(scores[,,iter])*(1/uniquenesses[j,iter-1]))%*%(data[,j]-mu[j,iter]))
      uniquenesses[j,iter]  <- rinvgamma(1, shape=(N+psi.alpha)/2, scale=(sum(data[,j] - mu[j,iter] - loadings[j,,iter]%*%t(scores[,,iter]))^2+psi.beta)/2)
    }
  }
  return(list(mu             = mu,
              scores         = scores, 
              loadings       = loadings, 
              uniquenesses   = uniquenesses))
}; gibbs.comp <- cmpfun(gibbs)

# run the gibbs sampler
n.iters <- 100000
sim <- gibbs.comp(data=wine, n.iters, Q=3, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)
#system.time(gibbs(data=data, n.iters, Q=2, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)
#system.time(gibbs.comp(data=data, n.iters, Q=2, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)

# save / load results
save(sim,file="/home/kmurphy/Wine Simulations.Rdata") # in server, tick box, export
load(file="Simulations/Wine Simulations.Rdata",envir=.GlobalEnv)