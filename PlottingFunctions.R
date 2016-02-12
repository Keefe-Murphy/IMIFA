################################
### IMIFA Plotting Functions ###
################################

plot.IMIFA  <- function(results=NULL, plot.meth=NULL, var=NULL, Label=NULL, fac=NULL,
                        ind=NULL, heat=T, n.fac=NULL, type=c("n", "p", "l"), mat=T, ... ) {
 
  if(missing(results))                stop("Results must be supplied")
  if(!exists(as.character(match.call()$results),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "IMIFA")       stop(paste0("Results object of class 'IMIFA' must be supplied"))
  if(missing(n.fac))     n.fac <- results$Q
  if(n.fac   > results$Q)             stop("Cannot plot this many factors")
  n.var     <- nrow(results$post.load)
  n.obs     <- nrow(results$post.f)
  if(!is.logical(heat))               stop("Arg. must be TRUE or FALSE")
  if(missing(plot.meth))              stop("What type of plot would you like to produce?")
  plot.meth <- match.arg(plot.meth, c("acf", "cum.var", "posterior", "trace"))
  type      <- match.arg(type)
  if(plot.meth != "cum.var") {
    if(missing(var))                  stop("What variable would you like to plot?")               
    var     <- match.arg(var, c("means", "scores", "loadings", "uniquenesses"))
  }
  if(!is.logical(mat))                stop("Arg. must be TRUE or FALSE")
  if(n.fac  == 1 ||
     !missing(ind))        mat <- F
  
  if(plot.meth == "acf") {
    if(var == "scores" || 
       var == "loadings") {
      if(missing(ind))     ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)             stop("Length of indexes for plotting cannot be greater than 2")
      if(var == "scores"   && 
         ind[1] >  n.obs)             stop(paste0("Length of first index cannot be greater than ", n.obs))
      if(var == "loadings" && 
         ind[1] >  n.var)             stop(paste0("Length of first index cannot be greater than ", n.var))
      if(ind[2] >  n.fac)             stop(paste0("Length of second index cannot be greater than ", n.fac))
    } else {
      if(missing(ind))     ind <- 1
      if(length(ind) > 1)             stop("Length of indexes for plotting cannot be greater than 1")
      if(ind    > n.var)              stop(paste0("Length of index cannot be greater than ", n.var))
    }
    if(var == "means") {
      acf(results$means[ind,], main="Means")
    }
    if(var == "scores") { 
      acf(results$scores[ind[1],ind[2],], main="Scores")
    }
    if(var == "loadings") { 
      acf(results$loadings[ind[1],ind[2],], main="Loadings")
    }
    if(var == "uniquenesses") { 
      acf(results$uniquenesses[ind,], main="Uniquenesses")
    }
  }
  
  if(plot.meth == "cum.var") {
    plot.x     <- results$cum.var
    prop.exp   <- plot.x[n.fac]
    if(n.fac    > 1) {
       plot(plot.x, type="l", main=paste0("Cumulative Variance:\n", n.fac, " Factors"), xlab="# Factors", 
            ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
       axis(1, at=1:length(plot.x), labels=1:n.fac)
       axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8) 
       points(x=n.fac, y=prop.exp, col="red", bg="red", pch=21)
     }
    cat(paste0("Proportion of Variation Explained = ",
                round(prop.exp[length(prop.exp)]*100, 2), "%", "\n"))
    if(max(prop.exp) > 1)             cat(paste0("Warning: chain may not have converged", "\n"))
  }
  
  if(plot.meth == "posterior") {
    if(!missing(ind))     heat <- F
    if(var == "scores" || 
       var == "loadings") {
      if(missing(ind))     ind <- c(1, 2)
      if(!missing(fac)) ind[2] <- max(fac, 2)
      if(n.fac == 1) {
        ind    <- 1
      } else {
        if(length(ind) > 2 ||
           ind[1] == ind[2])          stop("Only two columns can be plotted")
        if(ind[2]  > n.fac)           stop(paste0("Only the first ", n.fac, " columns can be plotted"))
      }
    }
    if(var  == "means") {
      plot.x   <- results$post.mu
      plot(plot.x, type=type, main="Posterior Means", ylab="Means", xlab="Variable", ylim=if(is.element(attr(results, "Method"), c("FA", "IFA"))) c(-1,1))
      if(type  == "n") text(x=1:length(plot.x), y=plot.x, names(plot.x), cex=0.5)
    }
    if(var == "scores") {
      plot.x   <- results$post.f
      if(ind[1] > n.obs)              stop(paste0("Only the first ", n.obs, " scores can be plotted"))
      if(!missing(Label)) {
        if(!exists(as.character(match.call()$Label),
                   envir=.GlobalEnv)) stop(paste0("Object ", match.call()$Label, " not found"))
        Label  <- as.factor(Label)
        if(length(Label) != n.obs)    stop(paste0("Labels must be a factor of length N=",  n.obs))
      } else {
        Label  <- 1
        cat(paste0("Should the data be labelled?", "\n"))
      }
      if(n.fac != 1) {
        plot(plot.x[,ind[1]], plot.x[,ind[2]], type=type, main="Posterior Scores", 
             xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind[2]), col=as.numeric(Label))
        if(type == "n") text(plot.x[,ind[1]], plot.x[,ind[2]], 1:n.obs, 
                             col=as.numeric(Label), cex=0.5)
      } else {
        plot(plot.x[,ind], type=type, main="Posterior Scores", 
             xlab="Observation", ylab=paste0("Factor ", ind), col=as.numeric(Label))
        if(type == "n") text(plot.x[,ind], col=as.numeric(Label), cex=0.5)
      }
    }
    if(var == "loadings") {
      plot.x   <- results$post.load
      if(heat) {
        par(mfrow=c(1, 1), mar=c(5.1, 7.1, 4.1, 2.1))
        image(z=t(plot.x[n.var:1,1:n.fac]), xlab="", ylab="", 
              main="Posterior Loadings", xaxt="n", yaxt="n")
        axis(1, cex.axis=0.8, line=-0.5, tick=F, 
             at=if(n.fac != 1) seq(0, 1, 1/(n.fac - 1)) else 0, labels=1:n.fac)
        if(n.var < 100) {
          axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
               at=seq(0, 1, 1/(n.var - 1)), 
               labels=rownames(plot.x)[n.var:1])
        }
        box(lwd=2)
        mtext("Factors", side=1, line=2)
        if(n.fac != 1) abline(v=seq(1/(2*(n.fac - 1)), 
                                    1-1/(2*(n.fac - 1)), 
                                    1/(n.fac - 1)), lty=2, lwd=1)
        par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))
      } else {
        if(ind[1] > n.var)            stop(paste0("Only the first ", n.var, " variables can be plotted"))
        if(n.fac != 1) {
          plot(plot.x[,ind[1]], plot.x[,ind[2]], type=type, main="Posterior Loadings", 
               xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind[2]))
          if(type == "n") text(plot.x[,ind[1]], plot.x[,ind[2]], rownames(plot.x), cex=0.5)
        } else { 
          plot(plot.x[,ind], type=type, main="Posterior Loadings",
               xlab="Variable", ylab=paste0("Factor ", ind))
          if(type == "n") text(plot.x[,ind], rownames(plot.x), cex=0.5)
        }
      }
    }
    if(var == "uniquenesses") {
      plot.x     <- results$post.psi
      plot(plot.x, type=type, main="Posterior Uniquenesses", ylab="Uniquenesses", xlab="Variable")
      if(type    == "n") text(1:length(plot.x), plot.x, names(plot.x), cex=0.5)
    }
  }  
  
  if(plot.meth == "trace") {
    if(var == "scores" || 
       var  == "loadings") {
      if(missing(ind))     ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)             stop("Length of indexes for plotting cannot be greater than 2")
      if(ind[1] >  n.var)             stop(paste0("Length of first index cannot be greater than ", n.var))
      if(ind[2] >  n.fac)             stop(paste0("Length of second index cannot be greater than ", n.fac))
    } else {
      if(missing(ind))     ind <- 1
      if(length(ind) > 1)             stop("Length of indexes for plotting cannot be greater than 1")
      if(ind    >  n.var)             stop(paste0("Length of second index cannot be greater than ", n.var))
    }
    if(!mat) iter <- 1:length(results$store)
    
    if(var == "means") {
      plot.x   <- results$means
      if(mat) {
        matplot(t(plot.x[,]), type="l", main="Trace Plot:\nMeans", ylab="Means", xlab="Iteration", ylim=if(is.element(attr(results, "Method"), c("FA", "IFA"))) c(-1,1))
      } else plot(x=iter, y=plot.x[ind,], type="l", ylab="Mean", xlab="Iteration", 
                         main=paste0("Trace Plot:\nMean of ", rownames(plot.x)[ind]), ylim=if(is.element(attr(results, "Method"), c("FA", "IFA"))) c(-1,1))
    }
    if(var == "scores") {
      plot.x   <- results$scores
      if(mat) {
        matplot(t(plot.x[ind[1],,]), type="l", main="Trace Plot:\nScores", ylab="Scores", xlab="Iteration")
      } else plot(x=iter, y=plot.x[ind[1],ind[2],], type="l", ylab="Scores", xlab="Iteration",
                         main=paste0("Trace Plot:\nScores - Observation ", ind[1], ", Factor", ind[2]))
    }
    if(var == "loadings") {
      plot.x   <- results$loadings
      if(mat) {
        matplot(t(plot.x[ind[1],,]), type="l", main="Trace Plot:\nLoadings", ylab="Loadings", xlab="Iteration")
      } else plot(x=iter, y=plot.x[ind[1],ind[2],], type="l", ylab="Loadings", xlab="Iteration",
                         main=paste0("Trace Plot:\nLoadings - Variable ", ind[1], ", Factor", ind[2]))
    }
    if(var == "uniquenesses") {
      plot.x   <- results$uniquenesses
      if(mat) {
        matplot(t(plot.x[,]), type="l", main="Trace Plot\nUniquenesses", ylab="Uniquenesses", xlab="Iteration")
      } else plot(x=iter, y=plot.x[ind,], ylab="Uniquenesses", type="l", xlab="Iteration",
                         main=paste0("Trace Plot:\nUniquenesses of ", rownames(plot.x)[ind]))
    }
  }
}