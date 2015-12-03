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
  set.seed(21092015)
  def.par         <- par()
  cases           <- c("Single", "Shrinkage", "Grouped", "IMIFA")
  case            <- 'Shrinkage'
  if(!is.element(case, cases)) stop("'case' must be one of 'Single', 'Shrinkage', 'Grouped', or 'IMIFA'.")
  pkgs            <- c("pgmm")
  invisible(lapply(pkgs, library, ch=T))
  # WARNING: Remove everything
    # rm(list = ls(all = TRUE))
  # WARNING: Remove loaded libraries
    # pkgs <- names(sessionInfo()$otherPkgs)
    # pkgs <- paste('package:', pkgs, sep = "")
    # invisible(lapply(pkgs, detach, ch = T, unload = T, force= T))
    
# Read in the data (& call it data)
  data(wine); Label <- as.factor(wine[,1]); wine <- wine[,-1]; data <- wine; rm("wine")
  #subjectmarks   <- read.csv(paste(dataDirectory, "/Data/", "SubjectMarks.csv", sep="")); data <- subjectmarks; rm("subjectmarks")
  #cereal         <- read.csv(paste(dataDirectory, "/Data/", "Cereal.csv", sep="")); data <- cereal; rm("cereal")

  # Simulate data
    #source(paste(dataDirectory, "/IMIFA-GIT/Simulate_Data.R", sep=""))
    #save(data, mu.true, load.true, psi.true, file=paste(dataDirectory,"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste(dataDirectory, "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv); if(exists("Label")) rm("Label")

# Vanilla 'factanal' for comparison purposes
  if(!exists("Label")) stop("Should the data be labelled?")
  N        <- nrow(data)
  P        <- sum(sapply(data, is.numeric))
  res      <- factanal(data[,sapply(data, is.numeric)], 
                       factors=round(sqrt(sum(sapply(data, is.numeric)))), control=list(nstart=50))
  res

# Initialise the Gibbs Sampler & set hyperparameters
  n.iters  <- 50000
  range.Q  <- 1:3   # can be SCALAR or VECTOR; scalar preferred!
  sigma.mu <- 0.5; sigma.l <- 0.5; psi.alpha <- 2; psi.beta <- 0.6
  if(case != 'Single') { phi.nu <- 3; delta.a1 <- 2.1; delta.a2 <- 3.1; rm('sigma.l', 'range.Q') }

  # Define full conditional & Gibbs Sampler functions for desired case
    source(paste(dataDirectory, "/IMIFA-GIT/Gibbs_BFA_", case, ".R", sep=""))
    if(case == 'Single') {
      sim    <- vector("list", length(range.Q))
    } else if(case == 'Shrinkage') {
      Q.ind  <- 1
      Q.star <- min(round(5 * log(P)), P)
      sim    <- vector("list", length(Q.star))
    } else {
      stop("Not yet implented for other cases.")
    }

# Run the Gibbs Sampler
{ start.time   <- proc.time()
  if(case == 'Single') {
    if(length(range.Q) == 1) {
      Q.ind    <- 1
      sim[[1]] <- gibbs.single(data, n.iters, Q=range.Q)
    } else {
      for(q in range.Q) { 
        Q.ind  <- q - min(range.Q) + 1
        sim[[Q.ind]] <- gibbs.single(data, n.iters, Q=q)
        cat(paste0(round(Q.ind/length(range.Q) * 100, 2), "% Complete", "\n"))
      }
    }
  } else if (case == 'Shrinkage') {
      sim[[1]] <- gibbs.shrink(data, n.iters, Q=Q.star)
  }
  total.time   <- proc.time() - start.time
  average.time <- total.time/if(exists('range.Q')) length(range.Q) else length(Q.star)
  sim$time     <- list(Total = total.time, Average = average.time); sim$time
}

# Save / Load results
  save(sim,file=paste(dataDirectory, "/Simulations/Wine_Simulations_", case, ".Rdata", sep="")) # in server, tick box, export
  load(file=paste(dataDirectory, "/Simulations/Wine_Simulations_", case, ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional additional burnin & thinning)
  burnin  <- 1
  thin    <- 1
  store   <- seq(from=burnin + 1, to=sim[[1]]$n.store, by=thin)
  
  if(case == 'Single' && length(range.Q) == 1) {
    Q     <- range.Q
    rm("range.Q")
  } else {
    if(case == 'Shrinkage') { post.Q <- 'Mode' }
    source(paste(dataDirectory, "/IMIFA-GIT/Tune_Parameters.R", sep=""))
  }
  
  # For user defined Q based on scree plot or bar plot
  # Rather than Q.ind <- which.max(prop.exp); Q <- range.Q[Q.ind] for 'Single' case
  # Rather than Q     <- names(Q.store[Q.store == max(Q.store)]) for 'Shrinkage' case
    # Q   <- 2
    # if(case == 'Single') { Q.ind <- which(range.Q == Q) } else Q.ind <- 1

  mu      <- sim[[Q.ind]]$mu[,store]                            
  f       <- array(sim[[Q.ind]]$f[,1:Q,store], dim=c(N, Q, length(store)))
  load    <- array(sim[[Q.ind]]$load[,1:Q,store], dim=c(P, Q, length(store)))
  psi     <- sim[[Q.ind]]$psi[,store]
  colnames(f) <- paste("Factor", 1:Q); rownames(load) <- colnames(data); colnames(load) <- paste("Factor", 1:Q)

# Loadings matrix / identifiability / # etc.
  l.temp  <- matrix(sim[[Q.ind]]$load[,1:Q,burnin], nr=P, nc=Q)
  for(b in 1:length(store)) {
    rot       <- procrustes(X=as.matrix(load[,,b]), Xstar=as.matrix(l.temp))$R
    load[,,b] <- load[,,b] %*% rot
    f[,,b]    <- t(f[,,b]  %*% rot)
  }

# Posterior Summaries & Plots, etc.
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

  plot(cum.var, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
       ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
  axis(1, at=1:length(cum.var), labels=1:Q)
  axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8) 
  points(x=Q, y=prop.exp, col="red", bg="red", pch=21)
    
  # Means
    plot(mu[1,], type="l")
    matplot(t(mu[,]), type="l")
    plot(post.mu, type="n")
    text(x=1:length(post.mu), y=post.mu, names(post.mu))
    acf(mu[1,])
  
  # Scores
    plot(f[1,1,], type="l")
    matplot(t(f[1,,]), type="l")
    plot(post.f, type="n")
    text(post.f[,1], post.f[,2], 1:nrow(post.f), col=if(exists("Label")) as.numeric(Label) else 1)
    plot(f[,,length(store)], type="n")
    text(f[,1,length(store)], f[,2,length(store)], 1:nrow(post.f), col=if(exists("Label")) as.numeric(Label) else 1)
    acf(f[1,1,])
      
  # Loadings
    plot(load[1,1,], type="l")
    matplot(t(load[1,,]), type="l")
    plot(post.load, type="n")
    text(post.load[,1], post.load[,2], rownames(post.load))
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
    plot(psi[1,], type="l")
    matplot(t(psi[,]), type="l")
    plot(post.psi, type="n")
    text(1:length(post.psi), post.psi, names(post.psi))
    acf(psi[1,])
####