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
  seed <- 21092015
  set.seed(seed)
  def.par         <- par()
  methods         <- c("FA", "IFA", "MIFA", "IMIFA")
  method          <- "IFA"
  if(!is.element(method, methods)) stop("'method' must be one of 'FA', 'IFA', 'MIFA', or 'IMIFA'.")
  pkgs            <- c("pgmm", "car")
  invisible(lapply(pkgs, library, ch=T))
  # WARNING: Remove everything
    # update.packages()
    # search()
    # searchpaths()
    # rm(list = ls(all = TRUE))
  # WARNING: Remove loaded libraries
    # pkgs <- names(sessionInfo()$otherPkgs)
    # pkgs <- paste('package:', pkgs, sep = "")
    # invisible(lapply(pkgs, detach, ch = T, unload = T, force= T))
    
# Read in the data (& call it data)
  data(wine); Label <- as.factor(wine[,1]); wine <- wine[,-1]; data <- wine; rm("wine")
  #subjectmarks     <- read.csv(paste(dataDirectory, "/Data/", "SubjectMarks.csv", sep="")); data <- subjectmarks; rm("subjectmarks")
  #cereal           <- read.csv(paste(dataDirectory, "/Data/", "Cereal.csv", sep="")); data <- cereal; rm("cereal")

  # Simulate data
    #source(paste(dataDirectory, "/IMIFA-GIT/Simulate_Data.R", sep=""))
    #save(data, mu.true, load.true, psi.true, file=paste(dataDirectory,"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste(dataDirectory, "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv); if(exists("Label")) rm("Label")

# Vanilla 'factanal' for comparison purposes
  if(!exists("Label")) stop("Should the data be labelled?")
  N        <- nrow(data)
  P        <- sum(sapply(data, is.numeric))
  (res     <- factanal(data[,sapply(data, is.numeric)], 
                       factors=round(sqrt(sum(sapply(data, is.numeric)))), control=list(nstart=50)))

# Initialise the Gibbs Sampler & set hyperparameters
  range.Q  <- 1:3   # can be SCALAR or VECTOR; scalar preferred!
  n.iters  <- 50000
  sigma.mu <- 0.5; sigma.l <- 0.5; psi.alpha <- 2; psi.beta <- 0.6
  if(method != "FA") { phi.nu <- 3; delta.a1 <- 2.1; delta.a2 <- 12.1; rm('sigma.l', 'range.Q') }

# Define full conditional & Gibbs Sampler functions for desired method
  source(paste(dataDirectory, "/IMIFA-GIT/Gibbs_", method, ".R", sep=""))

# Run the Gibbs Sampler
{ Rprof()
  start.time   <- proc.time()
  if(method == 'FA') {
    if(length(range.Q) == 1) {
      Q.ind    <- 1
      sim[[1]] <- gibbs.single(data=data, n.iters=ifelse(exists("n.iters"), n.iters, 50000), Q=range.Q)
    } else {
      for(q in range.Q) { 
        Q.ind  <- q - min(range.Q) + 1
        sim[[Q.ind]] <- gibbs.single(data=data, n.iters=ifelse(exists("n.iters"), n.iters, 50000), Q=q)
        cat(paste0(round(Q.ind/length(range.Q) * 100, 2), "% Complete", "\n"))
      }
    }
  } else if (method == 'IFA') {
      sim[[1]] <- gibbs.shrink(data=data, n.iters=ifelse(exists("n.iters"), n.iters, 50000), Q=Q.star)
  }
  total.time   <- proc.time() - start.time
  average.time <- total.time/ifelse(exists('range.Q'), length(range.Q), length(Q.star))
  sim$time     <- list(Total = total.time, Average = average.time); print(sim$time)  
  attr(sim, "Factors") <- if(method == 'FA') range.Q else Q.star
  attr(sim, "Date")    <- Sys.time()
  Rprof(NULL)
}
  summaryRprof()
  invisible(file.remove("Rprof.out"))

# Save / Load results
  sim.name <- "Wine"
  save(sim,file=paste(dataDirectory, "/Simulations/", sim.name, "_Simulations_", method, ".Rdata", sep="")) # in server, tick box, export
  load(file=paste(dataDirectory, "/Simulations/", sim.name, "_Simulations_", method, ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional additional burnin & thinning)
  burnin   <- 1
  thinning <- 1
  source(paste(dataDirectory, "/IMIFA-GIT/Diagnostics.R", sep=""))
  tune.parameters()
  
  # For user defined Q based on scree plot or bar plot
    # new.q(2)

# Posterior Summaries & Plots, etc.
  extract.results(Q)
  post.mu     <- apply(mu, 1, mean)
  post.f      <- apply(f, c(1,2), mean)
  post.load   <- apply(load, c(1,2), mean)
  post.psi    <- apply(psi, 1, mean)
  
  P           <- nrow(post.load)
  SS.load     <- colSums(post.load * post.load)
  communality <- sum(SS.load)
  prop.var    <- SS.load/P
  cum.var     <- cumsum(prop.var)
  prop.exp    <- communality/P
  prop.uni    <- 1 - prop.exp

  if(Q > 1) {
  plot(cum.var, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
       ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
  axis(1, at=1:length(cum.var), labels=1:Q)
  axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8) 
  points(x=Q, y=prop.exp, col="red", bg="red", pch=21)
  }; print(prop.exp[length(prop.exp)])

  # Means
    scatterplot(x=store, y=mu[1,])
    matplot(t(mu[,]), type="l")
    plot(post.mu, type="n", main="Posterior Means")
   #scatterplot(x=1:P, y=post.mu, main="Posterior Means", pch=NA, col=c("blue", "brown"))
    text(x=1:P, y=post.mu, names(post.mu), cex=0.5)
    acf(mu[1,])
  
  # Scores
    scatterplot(x=store, y=f[1,1,])
    matplot(t(f[1,,]), type="l")
    plot(post.f, type="n", main="Posterior Scores")
    text(post.f[,1], post.f[,2], 1:nrow(post.f), col=if(exists("Label")) as.numeric(Label) else 1, cex=0.5)
    plot(f[,,length(store)], type="n")
    text(f[,1,length(store)], f[,2,length(store)], 1:nrow(post.f), col=if(exists("Label")) as.numeric(Label) else 1, cex=0.5)
    acf(f[1,1,])
      
  # Loadings
    scatterplot(x=store, y=load[1,1,])
    matplot(t(load[1,,]), type="l")
    plot(post.load, type="n")
    text(post.load[,1], post.load[,2], rownames(post.load), cex=0.5)
    acf(load[1,1,])
    
    # Heatmaps
      par(mfrow=c(1, 1), mar=c(5.1, 7.1, 4.1, 2.1), xpd=F)
      image(z=t(post.load), xlab="", ylab="", 
            main="Posterior Loadings", xaxt="n", yaxt="n")
      axis(1, cex.axis=0.8, line=-0.5, tick=F, 
           at=if(Q != 1) seq(0, 1, 1/(Q-1)) else 0, labels=1:Q)
      axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
           at=seq(0, 1, 1/(nrow(post.load)-1)), labels=rownames(post.load))
      box(lwd=2)
      mtext("Factors", side=1, line=2)
      if(Q != 1) abline(v=seq(1/(2*(Q-1)), 1-1/(2*(Q-1)), 1/(Q-1)), lty=2, lwd=1)
      invisible(par(def.par))

  # Uniquenesses
    scatterplot(x=store, y=psi[1,])
    matplot(t(psi[,]), type="l")
    plot(post.psi, type="n")
   #scatterplot(x=1:P, y=post.psi, main="Posterior Uniquenesses", pch=NA, col=c("blue", "brown"))
    text(1:length(post.psi), post.psi, names(post.psi), cex=0.5)
    acf(psi[1,])
####