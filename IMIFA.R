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
    
# Read in the data
  data(wine); Lab <- as.factor(wine[,1]); wine <- wine[,-1]
  #subjectmarks   <- read.csv(paste(dataDirectory, "/Data/", "SubjectMarks.csv", sep=""))
  #cereal         <- read.csv(paste(dataDirectory, "/Data/", "Cereal.csv", sep=""))

# Simulate data
  #source(paste(dataDirectory, "/IMIFA-GIT/Simulate_Data.R", sep=""))
  #save(data, mu.true, load.true, psi.true, file=paste(dataDirectory,"/Data/Simulated_Data.Rdata", sep=""))
  load(file=paste(dataDirectory, "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)

# Run the Gibbs Sampler
  sim          <- imifa.gibbs(wine, n.iters=50, method="IFA")

# Run the Gibbs Sampler
{ Rprof()
  start.time   <- proc.time()
  total.time   <- proc.time() - start.time
  avg.time     <- total.time/ifelse(exists('range.Q'), length(range.Q), length(Q.star))
  attr(sim, "Time")    <- list(Total = total.time, Average = avg.time); print(attr(sim, "Time"))  
  Rprof(NULL)
}
  summaryRprof()
  invisible(file.remove("Rprof.out"))

# Save / Load results
  save(sim, file=paste(dataDirectory, "/Simulations/", 
                       attr(sim, "Name"), "_Simulations_", attr(sim, "Method"), "_", 
                       attr(sim, "Date"), ".Rdata", sep=""))
  load(file=paste(dataDirectory, "/Simulations/", 
                  attr(sim, "Name"), "_Simulations_", attr(sim, "Method"), "_", 
                  attr(sim, "Date"), ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional: additional 'burnin' & 'thinning' & user-defined Q)
  source(paste(dataDirectory, "/IMIFA-GIT/Diagnostics.R", sep=""))
  res          <- tune.sims(sim)
  
# Posterior Summaries & Plots, etc.  
  plot.cum.var(res)
  
# Means
  plot.trace(res, "m", F)
  plot.trace(res, "m")
  plot.posterior(res, "m")
  plot.acf(res, "m")
  
# Scores
  plot.trace(res, "s", F)
  plot.trace(res, "s")
  plot.posterior(res, "s", Lab)
  plot.acf(res, "s")
      
# Loadings
  plot.trace(res, "l", F)
  plot.trace(res, "l")
  plot.posterior(res, "l")
  plot.acf(res, "l")
  plot.load.heat(res)

# Uniquenesses
  plot.trace(res, "u", F)
  plot.trace(res, "u")
  plot.posterior(res, "u")
  plot.acf(res, "u")
####