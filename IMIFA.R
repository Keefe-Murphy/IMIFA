######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
  if(getwd() != "/home/kmurphy") {
    setwd("C:/Users/Windows/Documents/Claire IMIFA")
  }
  source(paste(getwd(), "/IMIFA-GIT/PackageSetup.R", sep=""))
    
# Read in the data
  data(wine)
  wine[,1]   <- as.factor(wine[,1])
  Lab        <- wine[,1]
  #subjects  <- read.csv(paste(getwd(), "/Data/", "SubjectMarks.csv", sep=""))
  #cereal    <- read.csv(paste(getwd(), "/Data/", "Cereal.csv", sep=""))

# Simulate data
  SimData    <- sim.imifa()
  #save(SimData, file=paste(getwd(),"/Data/Simulated_Data.Rdata", sep=""))
  load(file=paste(getwd(), "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)

# Run the Gibbs Sampler
  sim <- imifa.gibbs(wine, method="IFA")

# Save / Load results
  save(sim, file=paste(getwd(), "/Simulations/", 
                       attr(sim, "Name"), "__Simulations_", attr(sim, "Method"), "_", 
                       attr(sim, "Date"), ".Rdata", sep=""))
  load(file=paste(getwd(), "/Simulations/", 
                  attr(sim, "Name"), "__Simulations_", attr(sim, "Method"), "_", 
                  attr(sim, "Date"), ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional: additional 'burnin' & 'thinning' & user-defined Q)
  res <- tune.sims(sim)
  
# Posterior Summaries & Plots, etc.  
  plot(res, "c")
  
# Means
  plot(res, "t", "m")
  plot(res, "t", "m", mat=F)
  plot(res, "p", "m")
  plot(res, "a", "m")
  
# Scores
  plot(res, "t", "s")
  plot(res, "t", "s", mat=F)
  plot(res, "p", "s", Lab)
  plot(res, "a", "s")
      
# Loadings
  plot(res, "t", "l")
  plot(res, "t", "l", mat=F)
  plot(res, "p", "l", heat=F)
  plot(res, "p", "l")
  plot(res, "a", "l")
  
# Uniquenesses
  plot(res, "t", "u")
  plot(res, "t", "u", mat=F)
  plot(res, "p", "u")
  plot(res, "a", "u")
####