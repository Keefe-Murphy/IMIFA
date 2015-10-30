######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
  dataDirectory <- "C:/Users/Windows/Documents/Claire IMIFA"
  setwd(dataDirectory)
  set.seed(21092015)
  library(MASS)
  library(compiler)
  library(pgmm)
  library(MCMCpack)

# Read in the data
  data(wine); labels <- wine[,1]; wine <- wine[,-1]
  #subjectmarks <- read.csv(paste(dataDirectory,"/Data/","SubjectMarks.csv",sep=""))
  #cereal       <- read.csv(paste(dataDirectory,"/Data/","Cereal.csv",sep=""))

  # Simulate data
    #source(paste(dataDirectory,"/IMIFA-GIT/Simulate_Data.R",sep=""))
    #save(data, mu.true, f.true, load.true, psi.true, eps.true, file=paste(dataDirectory,"/Data/Simulated_Data.Rdata",sep=""))
    load(file=paste(dataDirectory,"/Data/Simulated_Data.Rdata",sep=""),envir=.GlobalEnv)

# Define full conditional functions
  source(paste(dataDirectory,"/IMIFA-GIT/FullConditionals_BFA_Single.R",sep=""))

# Gibbs Sampler function
  source(paste(dataDirectory,"/IMIFA-GIT/Gibbs_BFA_Single.R",sep=""))

# Run the Gibbs Sampler
  data    <- data
  n.iters <- 50000
  Q       <- 2
  #sim    <- gibbs.single(data, n.iters, Q, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5, burnin=(n.iters/5) - 1, thin=2)
  sim     <- gibbs.single(data, n.iters, Q)

# Save / Load results
  save(sim,file=paste(dataDirectory,"/Simulations/Wine Simulations.Rdata",sep="")) # in server, tick box, export
  load(file=paste(dataDirectory,"/Simulations/Wine Simulations.Rdata",sep=""),envir=.GlobalEnv)

# NB: You can check your answer by plotting
#     the scores of a 2-factor model to the
#     wine dataset. Expect to see a horseshoe.

# Convergence diagnostics (optional additional burnin & thinning)
  burnin  <- 1
  thin    <- 1
  n.store <- seq(from=burnin + 1, to=tail(dim(sim$mu), n=1), by=thin)
  mu      <- sim$mu[,n.store]
  f       <- sim$f[,,n.store]
  load    <- sim$load[,,n.store]
  psi     <- sim$psi[,n.store]

# Loadings matrix / identifiability / # etc.
  l.temp  <- sim$load[,,burnin]
  for(b in 1:tail(dim(load), n=1)) {
    rot       <- procrustes(X=load[,,b], Xstar=l.temp)$R
    load[,,b] <- load[,,b] %*% rot
    f[,,b]    <- t(t(rot) %*% t(f[,,b]))
  }

# Plots & posterior summaries etc.
  # Means
  plot(mu[1,], type="l")
  matplot(t(mu[,]), type="l")
  post.mu <- apply(mu, 1, mean)
  plot(post.mu, type="n")
  text(x=1:length(post.mu), y=post.mu, names(post.mu))
  acf(mu[1,])
  
  # Scores
  plot(f[1,1,], type="l")
  matplot(t(f[1,,]), type="l")
  post.f <- apply(f, c(1,2), mean)
  plot(post.f, type="n")
  text(post.f[,1], post.f[,2], 1:nrow(post.f), col=if(exists("labels")) labels else 1)
  plot(f[,,tail(dim(f), n=1)], type="n")
  text(f[,1,tail(dim(f), n=1)], f[,2,tail(dim(f), n=1)], 1:nrow(post.f), col=if(exists("labels")) labels else 1)
  acf(f[1,1,])
  
  # Uniquenesses
  plot(psi[1,], type="l")
  matplot(t(psi[,]), type="l")
  post.psi <- apply(psi, 1, mean)
  plot(post.psi, type="n")
  text(1:length(post.psi), post.psi, names(post.psi))
  acf(psi[1,])
  
  # Loadings
  plot(load[1,1,], type="l")
  matplot(t(load[1,,]), type="l")
  post.load <- apply(load, c(1,2), mean)
  plot(post.load, type="n")
  text(post.load[,1], post.load[,2], rownames(post.load))
  acf(load[1,1,])
  
  # Summaries 
  P <- nrow(post.load)
  sum(post.psi)/P # % of variance which is unique
  communality <- P - sum(post.psi)
  communality/P   # % of variance which is explained