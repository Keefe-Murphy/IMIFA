################################
### IMIFA Plotting Functions ###
################################

# Plot the cumulative proportion of variation explained 
  plot.cum.var   <- function(results=NULL, n.fac=results$Q, cum.var=results$cum.var, ...) {
    if(missing(results))          stop("Results must be supplied")
    if(!exists(as.character(match.call()$results),
               envir=.GlobalEnv)) stop(paste0("Object ", match.call()$results, " not found"))
    prop.exp     <- cum.var[n.fac]
    if(n.fac      > 1) {
      plot(cum.var, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
           ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
      axis(1, at=1:length(cum.var), labels=1:n.fac)
      axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8) 
      points(x=n.fac, y=prop.exp, col="red", bg="red", pch=21)
    }
    cat(paste0("Proportion of Variation Explained = ",
               round(prop.exp[length(prop.exp)], 3), "\n"))
  }

# Plot heatmap of Posterior Loadings
  load.heat      <- function(results=NULL, n.fac=results$Q, ...) {
    if(missing(results))          stop("Results must be supplied")
    if(!exists(as.character(match.call()$results),
               envir=.GlobalEnv)) stop(paste0("Object ", match.call()$results, " not found"))
    if(n.fac > results$Q)         stop("Cannot plot this many loadings columns")
    par(mfrow=c(1, 1), mar=c(5.1, 7.1, 4.1, 2.1))
    image(z=t(results$post.load[,1:n.fac]), xlab="", ylab="", 
          main="Posterior Loadings", xaxt="n", yaxt="n")
    axis(1, cex.axis=0.8, line=-0.5, tick=F, 
         at=if(n.fac != 1) seq(0, 1, 1/(n.fac - 1)) else 0, labels=1:n.fac)
    axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
         at=seq(0, 1, 1/(nrow(results$post.load) - 1)), 
         labels=rownames(results$post.load))
    box(lwd=2)
    mtext("Factors", side=1, line=2)
    if(n.fac != 1) abline(v=seq(1/(2*(n.fac - 1)), 
                            1-1/(2*(n.fac - 1)), 
                            1/(n.fac - 1)), lty=2, lwd=1)
    par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))
  }

# Plot posteriors
  plot.posterior <- function(results=NULL, x=NULL, Labels=NULL, 
                             type=c("n", "p"), facs=NULL, ...) {
    if(missing(results))          stop("Results must be supplied")
    if(!exists(as.character(match.call()$results),
               envir=.GlobalEnv)) stop(paste0("Object ", match.call()$results, " not found"))
    x    <- match.arg(x, c("means", "scores", "loadings", "uniquenesses"))
    type <- match.arg(type)
    if(x == "scores" || x == "loadings") {
      if(missing(facs)) facs <- c(1, 2)
      if(length(facs) > 2)        stop("Only two columns can be plotted")
      if(max(facs) > results$Q)   stop(paste0("Only the first ", res$Q, " columns can be plotted"))
    }
    if(x == "means") {
      plot.x  <- results$post.mu
      plot(plot.x, type=type, main="Posterior Means", ylab="Means")
      if(type == "n") text(x=1:length(plot.x), y=plot.x, names(plot.x), cex=0.5)
    }
    if(x == "scores") {
      plot.x  <- results$post.f
      if(!missing(Labels)) {
        if((length(Labels) != N) ||
             !is.factor(Labels))   stop("Labels must be a factor of length N")
      } else Labels <- 1
      print(Labels)
      plot(plot.x[,facs[1]], plot.x[,facs[2]], type=type, main="Posterior Scores", 
           xlab=paste("Factor ", facs[1]), ylab=paste("Factor ", facs[2]), col=as.numeric(Labels))
      if(type == "n") text(plot.x[,facs[1]], plot.x[,facs[2]], 1:nrow(plot.x), 
                           col=as.numeric(Labels), cex=0.5)
    }
    if(x == "loadings") {
      plot.x  <- results$post.load
      plot(plot.x[,facs[1]], plot.x[,facs[2]], type=type, main="Posterior Loadings", 
           xlab=paste("Factor ", facs[1]), ylab=paste("Factor ", facs[2]))
      if(type == "n") text(plot.x[,facs[1]], plot.x[,facs[2]], rownames(plot.x), cex=0.5)
    }
    if(x == "uniquenesses") {
      plot.x  <- results$post.psi
      plot(plot.x, type=type, main="Posterior Uniquenesses", ylab="Uniquenesses")
      if(type == "n") text(1:length(plot.x), plot.x, names(plot.x), cex=0.5)
    }
  } 

# Autocorrelation Functions
  plot.acf    <- function(results=NULL, x=NULL, ind=NULL, ...) {
    if(missing(results))          stop("Results must be supplied")
    if(!exists(as.character(match.call()$results),
               envir=.GlobalEnv)) stop(paste0("Object ", match.call()$results, " not found"))
    x    <- match.arg(x, c("means", "scores", "loadings", "uniquenesses"))
    if(x == "scores" || x == "loadings") {
      if(missing(ind)) ind <- c(1, 1)
      if(length(ind) > 2)         stop("Length of indexes for plotting cannot be greater than 2")
      if(ind[1] > N)              stop(paste0("Length of first index cannot be greater than ", N))
      if(ind[2] > res$Q)          stop(paste0("Length of second index cannot be greater than ", res$Q))
    } else {
      if(missing(ind)) ind <- 1
      if(length(ind) > 1)         stop("Length of indexes for plotting cannot be greater than 1")
      if(ind > P)                 stop(paste0("Length of second index cannot be greater than ", P))
    }
    if(x == "means") {
      acf(results$means[ind,], main="Means")
    }
    if(x == "scores") { 
      acf(results$scores[ind[1],ind[2],], main="Scores")
    }
    if(x == "loadings") { 
      acf(results$loadings[ind[1],ind[2],], main="Loadings")
    }
    if(x == "uniquenesses") { 
      acf(results$uniquenesses[ind,], main="Uniquenesses")
    }
  }