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
  load(file=paste(dataDirectory, "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)

# Initialise the Gibbs Sampler & override default hyperparameters if desired
  sim          <- imifa.gibbs(data, method="IFA")

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
  attr(sim, "Date")    <- format(Sys.Date(), "%d-%b-%Y")
  attr(sim, "Factors") <- if(method == 'FA') range.Q else Q.star
  attr(sim, "Method")  <- method
  #attr(sim, "Name")   <- 
  attr(sim, "Time")    <- list(Total = total.time, Average = avg.time); print(attr(sim, "Time"))  
  Rprof(NULL)
}
  summaryRprof()
  invisible(file.remove("Rprof.out"))

# Save / Load results
  sim.name     <- "Wine"
  save(sim, file=paste(dataDirectory, "/Simulations/", 
                       sim.name, "_Simulations_", attr(sim, "Method"), "_", 
                       attr(sim, "Date"), ".Rdata", sep=""))
  load(file=paste(dataDirectory, "/Simulations/", 
                  sim.name, "_Simulations_", attr(sim, "Method"), "_", 
                  attr(sim, "Date"), ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional: additional 'burnin' & 'thinning' & user-defined Q)
  source(paste(dataDirectory, "/IMIFA-GIT/Diagnostics.R", sep=""))
  res          <- tune.sims(sim)
  
# Posterior Summaries & Plots, etc.  
  plot.cum.var(res)
  
# Means
  scatterplot(x=res$store, y=res$means[1,])
  matplot(t(res$means[,]), type="l")
  plot.posterior(res, "m")
  plot.acf(res, "m")
  
# Scores
  scatterplot(x=res$store, y=res$scores[1,1,])
  matplot(t(res$scores[1,,]), type="l")
  plot.posterior(res, "s", Label)
  plot.acf(res, "s")
      
# Loadings
  scatterplot(x=res$store, y=res$loadings[1,1,])
  matplot(t(res$loadings[1,,]), type="l")
  plot.posterior(res, "l")
  plot.acf(res, "l")
  load.heat(res)

# Uniquenesses
  scatterplot(x=res$store, y=res$uniquenesses[1,])
  matplot(t(res$uniquenesses[,]), type="l")
  plot.posterior(res, "u")
  plot.acf(res, "u")
####