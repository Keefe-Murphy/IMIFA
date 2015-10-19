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
  source(paste(dataDirectory,"/IMIFA-GIT/IMIFA_FullConditionals.R",sep=""))

# Gibbs Sampler function
  source(paste(dataDirectory,"/IMIFA-GIT/Gibbs_BFA_Single.R",sep=""))

# Run the gibbs sampler
  n.iters <- 100000
  sim     <- gibbs.single(data=wine, n.iters, Q=3, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5)

# Save / Load results
  save(sim,file=paste(dataDirectory,"/Simulations/Wine Simulations.Rdata",sep="")) # in server, tick box, export
  load(file=paste(dataDirectory,"Simulations/Wine Simulations.Rdata",sep""),envir=.GlobalEnv)

# NB: You can check your answer by plotting
#     the scores of a 2-factor model to the
#     wine dataset. Expect to see a horseshoe.

# Convergence diagnostics
  burnin <- 10000
  thin   <- 3
  mu     <- sim$mu[,seq(from=burnin+1, to=n.iters, by=thin)]
  f      <- sim$f[,,seq(from=burnin+1, to=n.iters, by=thin)]
  load   <- sim$load[,,seq(from=burnin+1, to=n.iters, by=thin)]
  psi    <- sim$psi[,seq(from=burnin+1, to=n.iters, by=thin)]

# Loadings matrix / identifiability / # etc.
  l.temp <- sim$load[,,burnin]
  for(b in 1:dim(load)[3]) {
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
  text(post.f[,1], post.f[,2], 1:nrow(post.f), col=labels)
  plot(f[,,dim(f)[3]], type="n")
  text(f[,1,dim(f)[3]], f[,2,dim(f)[3]], 1:nrow(post.f), col=labels)
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
  communality/P   # % of variance 
  #sum(post.load[,]^2)