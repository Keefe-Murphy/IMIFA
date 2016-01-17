######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
  if(getwd() == "/home/kmurphy") {
    dataDirectory <- getwd()
  } else {
    dataDirectory <- "C:/Users/Windows/Documents/Claire IMIFA"
    setwd(dataDirectory)
  }
  source(paste(dataDirectory, "/IMIFA-GIT/Preamble.R", sep=""))
    
# Read in the data (& call it data)
  data(wine); Label <- as.factor(wine[,1]); wine <- wine[,-1]; data <- wine; rm("wine")
  #subjectmarks     <- read.csv(paste(dataDirectory, "/Data/", "SubjectMarks.csv", sep="")); data <- subjectmarks; rm("subjectmarks")
  #cereal           <- read.csv(paste(dataDirectory, "/Data/", "Cereal.csv", sep=""));       data <- cereal; rm("cereal")

# Simulate data
  #source(paste(dataDirectory, "/IMIFA-GIT/Simulate_Data.R", sep=""))
  #save(data, mu.true, load.true, psi.true, file=paste(dataDirectory,"/Data/Simulated_Data.Rdata", sep=""))
  load(file=paste(dataDirectory, "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv); if(exists("Label")) rm("Label")

# Initialise the Gibbs Sampler & override default hyperparameters if desired
  imifa     <- imifa.gibbs(data, method="IFA")

# Run the Gibbs Sampler
{ Rprof()
  start.time   <- proc.time()
  if(method == 'FA') {
    if(length(range.Q) == 1) {
      sim[[1]] <- do.call(paste0("gibbs.", method), 
                          args=list(dat=data, n.iters=ifelse(exists("n.iters"), n.iters, 50000), Q=range.Q))
    } else {
      for(q in range.Q) { 
        Q.ind  <- q - min(range.Q) + 1
        sim[[Q.ind]] <- do.call(paste0("gibbs.", method),
                                args=list(dat=data, n.iters=ifelse(exists("n.iters"), n.iters, 50000), Q=q))
        cat(paste0(round(Q.ind/length(range.Q) * 100, 2), "% Complete", "\n"))
      }
    }
  } else if (method == 'IFA') {
      sim[[1]] <- do.call(paste0("gibbs.", method),
                          args=list(dat=data, n.iters=ifelse(exists("n.iters"), n.iters, 50000), Q=Q.star))
  }
  total.time   <- proc.time() - start.time
  avg.time     <- total.time/ifelse(exists('range.Q'), length(range.Q), length(Q.star))
  attr(sim, "Date")    <- Sys.time()
  attr(sim, "Factors") <- if(method == 'FA') range.Q else Q.star
  #attr(sim, "Method") <- 
  #attr(sim, "Name")   <- 
  attr(sim, "Time")    <- list(Total = total.time, Average = avg.time); print(attr(sim, "Time"))  
  Rprof(NULL)
}
  summaryRprof()
  invisible(file.remove("Rprof.out"))

# Save / Load results
  sim.name     <- "Wine"
  save(sim, file=paste(dataDirectory, "/Simulations/", sim.name, "_Simulations_", method, ".Rdata", sep="")) # in server, tick box, export
  load(file=paste(dataDirectory, "/Simulations/", sim.name, "_Simulations_", method, ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional: additional 'burnin' & 'thinning' & user-defined Q)
  source(paste(dataDirectory, "/IMIFA-GIT/Diagnostics.R", sep=""))
  res          <- tune.params()
  
# Posterior Summaries & Plots, etc.  
  plot.cum.var()
  
# Means
  scatterplot(x=store, y=mu[1,])
  matplot(t(mu[,]), type="l")
  plot.posterior("m")
  acf(mu[1,])
  
# Scores
  scatterplot(x=store, y=f[1,1,])
  matplot(t(f[1,,]), type="l")
  plot.posterior("s")
  acf(f[1,1,])
      
# Loadings
  scatterplot(x=store, y=lmat[1,1,])
  matplot(t(lmat[1,,]), type="l")
  plot.posterior("l")
  acf(lmat[1,1,])
  load.heat()

# Uniquenesses
  scatterplot(x=store, y=psi[1,])
  matplot(t(psi[,]), type="l")
  plot.posterior("u")
  acf(psi[1,])
####