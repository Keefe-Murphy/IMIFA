################################
### IMIFA Plotting Functions ###
################################

plot.IMIFA  <- function(results=NULL, plot.meth=c("all", "correlation", "density", "posterior", "Q", "trace"), 
                        var=c("means", "scores", "loadings", "uniquenesses"), Label=NULL, 
                        fac=NULL, ind=NULL, heat=T, n.fac=NULL, type=c("n", "p", "l"), mat=T, ... ) {
 
  if(missing(results))                stop("Results must be supplied")
  if(!exists(as.character(match.call()$results),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "IMIFA")       stop(paste0("Results object of class 'IMIFA' must be supplied"))
  if(missing(n.fac))     n.fac <- results$Q.results$Q
  if(n.fac   > results$Q.results$Q)   stop("Cannot plot this many factors")
  n.var     <- attr(results, "Vars")
  n.obs     <- attr(results, "Obs")
  if(!is.logical(heat))               stop("heat must be TRUE or FALSE")
  if(missing(plot.meth))              stop("What type of plot would you like to produce?")
  plot.meth <- match.arg(plot.meth)
  type      <- match.arg(type)
  pardef    <- par(no.readonly = T)
  pardef$mfrow <- c(1, 1)
  m.sw      <- c(Q.sw = F, cor.sw = F, den.sw = F, pos.sw = F, tra.sw = F)
  if(plot.meth == "all")  {
    m.sw[-1]   <- !m.sw[-1]
    layout(matrix(c(1, 2, 3, 4), nr=2, nc=2, byrow = TRUE))
    cex.t   <- 0.75
  } else {
    sw.n    <- paste0(substring(plot.meth, 1, 3), ".sw")
    m.sw[sw.n] <- T
    cex.t   <- 1
  }
  if(!m.sw["Q.sw"]      &&
     missing(var))                  { stop("What variable would you like to plot?"); par(pardef) }              
  v.sw      <- attr(results, "Switch")
  names(v.sw)  <- formals(sys.function(sys.parent()))$var
  var       <- match.arg(var)
  method    <- attr(results, "Method")
  if(!v.sw[var]         && 
     !m.sw["Q.sw"])                 { stop(paste0(var, " were not stored")); par(pardef) }
  if(!is.logical(mat))              { stop("mat must be TRUE or FALSE"); par(pardef) }
  if(!missing(ind))      x.ind <- ind
  ind.x     <- !exists("x.ind", envir=environment())
  if(n.fac  == 1 ||
     !missing(ind))        mat <- F

  if(m.sw["tra.sw"]) {
    if(var == "scores" || 
       var == "loadings") {
      if(ind.x)            ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)           { stop("Length of indexes for plotting cannot be greater than 2"); par(pardef) }
      if(var == "scores") {
        if(ind[1] >  n.obs)         { stop(paste0("First index cannot be greater than ",  n.obs)); par(pardef) }
      } else if(ind[1] > n.var)     { stop(paste0("First index cannot be greater than ",  n.var)); par(pardef) }
      if(ind[2]   >  n.fac)         { stop(paste0("Second index cannot be greater than ", n.fac)); par(pardef) }
    } else {
      if(ind.x)            ind <- 1
      if(length(ind) >   1)         { stop("Length of indexes for plotting cannot be greater than 1"); par(pardef) }
      if(ind      >  n.var)         { stop(paste0("Index cannot be greater than ", n.var)); par(pardef) }
    }
    if(!mat) iter <- 1:attr(results, "Store")
    
    if(var == "means") {
      plot.x   <- results$means
      if(mat) {
        matplot(t(plot.x[,]), type="l", ylab="Means", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
        title(main=list("Trace Plot:\nMeans", cex=cex.t))
      } else {
        plot(x=iter, y=plot.x[ind,], type="l", ylab="Mean", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
        title(main=list(paste0("Trace Plot:\nMean of ", rownames(plot.x)[ind], " Variable"), cex=cex.t))
      }
    }
    if(var == "scores") {
      plot.x   <- results$scores
      if(mat) {
        matplot(t(plot.x[ind[1],,]), type="l", ylab="Scores", xlab="Iteration")
        title(main=list("Trace Plot:\nScores", cex=cex.t))
      } else {
        plot(x=iter, y=plot.x[ind[1],ind[2],], type="l", ylab="Scores", xlab="Iteration")
        title(main=list(paste0("Trace Plot:\nScores - Observation ", ind[1], ", Factor ", ind[2]), cex=cex.t))
      }
    }
    if(var == "loadings") {
      plot.x   <- results$loadings
      if(mat) {
        matplot(t(plot.x[ind[1],,]), type="l", ylab="Loadings", xlab="Iteration")
        title(main=list("Trace Plot:\nLoadings", cex=cex.t))
      } else {
        plot(x=iter, y=plot.x[ind[1],ind[2],], type="l", ylab="Loadings", xlab="Iteration")
        title(main=list(paste0("Trace Plot:\nLoadings - ", rownames(plot.x)[ind[1]], " Variable, Factor ", ind[2]), cex=cex.t))
      }
    }
    if(var == "uniquenesses") {
      plot.x   <- results$uniquenesses
      if(mat) {
        matplot(t(plot.x[,]), type="l", ylab="Uniquenesses", xlab="Iteration")
        title(main=list("Trace Plot\nUniquenesses", cex=cex.t))
      } else {
        plot(x=iter, y=plot.x[ind,], ylab="Uniquenesses", type="l", xlab="Iteration")
        title(main=list(paste0("Trace Plot:\nUniqueness of ", rownames(plot.x)[ind], " Variable"), cex=cex.t))
      }
    }
    if(!ind.x)             ind <- x.ind
  }
  
  if(m.sw["den.sw"]) {
    if(var == "scores" || 
       var  == "loadings") {
      if(ind.x)            ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)           { stop("Length of indexes for plotting cannot be greater than 2"); par(pardef) }
      if(var == "scores") {
        if(ind[1] >  n.obs)         { stop(paste0("First index cannot be greater than ",  n.obs)); par(pardef) }
      } else if(ind[1] > n.var)     { stop(paste0("First index cannot be greater than ",  n.var)); par(pardef) }
      if(ind[2]   >  n.fac)         { stop(paste0("Second index cannot be greater than ", n.fac)); par(pardef) }
    } else {
      if(ind.x)            ind <- 1
      if(length(ind) >   1)         { stop("Length of indexes for plotting cannot be greater than 1"); par(pardef) }
      if(ind      >  n.var)         { stop(paste0("Index cannot be greater than ", n.var)); par(pardef) }
    }
    if(var == "means") {
      plot.x   <- results$means
      if(mat) {
        plot.x <- apply(plot.x, 1, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        title(main=list(paste0("Density Plot:\n Means"), cex=cex.t))
      } else {
        plot.d <- density(plot.x[ind,])
        plot(plot.d, main="")
        title(main=list(paste0("Density Plot:\n Mean of ", rownames(plot.x)[ind], " Variable"), cex=cex.t))
        polygon(plot.d, col="black")
      }
    }
    if(var == "scores") {
      plot.x   <- results$scores
      if(mat) {
        plot.x <- apply(plot.x, 2, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        title(main=list(paste0("Density Plot:\n Scores"), cex=cex.t))
      } else {
        plot.d <- density(plot.x)
        plot(plot.d, main="")
        title(main=list(paste0("Density Plot:\n Scores - Observation ", ind[1], ", Factor ", ind[2]), cex=cex.t))
        polygon(plot.d, col="black")
      }
    }
    if(var == "loadings") {
      plot.x   <- results$loadings
      if(mat) {
        plot.x <- apply(plot.x, 2, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        title(main=list(paste0("Density Plot:\n Loadings"), cex=cex.t))
      } else {
        plot.d <- density(plot.x)
        plot(plot.d, main="")
        title(main=list(paste0("Density Plot:\n Loadings - ", rownames(plot.x)[ind[1]], " Variable, Factor ", ind[2]), cex=cex.t))
        polygon(plot.d, col="black")
      }
    }
    if(var == "uniquenesses") {
      plot.x   <- results$uniquenesses
      if(mat) {
        plot.x <- apply(plot.x, 1, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        title(main=list(paste0("Density Plot:\n Uniquenesses"), cex=cex.t))
      } else {
        plot.d <- density(plot.x[ind,])
        plot(plot.d, main="")
        title(main=list(paste0("Density Plot:\n Uniqueness of ", rownames(plot.x)[ind], " Variable"), cex=cex.t))
        polygon(plot.d, col="black")
      }
    }
    if(!ind.x)             ind <- x.ind
  }
  
  if(m.sw["pos.sw"]) {
    if(var == "scores" || 
       var == "loadings") {
      if(ind.x)            ind <- c(1, 2)
      if(!missing(fac)) ind[2] <- max(fac, 2)
      if(n.fac == 1) {
        ind    <- 1
      } else {
        if(length(ind) >  2)        { stop("Only two columns can be plotted"); par(pardef) }
        if(ind[1]  >  n.fac)        { stop(paste0("Only the first ", n.fac, " columns can be plotted")); par(pardef) }
        if(ind[2]  >  n.fac)        { stop(paste0("Only the first ", n.fac, " columns can be plotted")); par(pardef) }
      }
    }
    if(var  == "means") {
      plot.x   <- results$post.mu
      plot(plot.x, type=type, ylab="Means", xlab="Variable", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
      title(main=list("Posterior Means", cex=cex.t))
      if(type  == "n") text(x=1:length(plot.x), y=plot.x, names(plot.x), cex=0.5)
    }
    if(var == "scores") {
      plot.x   <- results$post.f
      if(ind[1] > n.obs)            { stop(paste0("Only the first ", n.obs, " scores can be plotted")); par(pardef) }
      if(!missing(Label)) {
        if(!exists(as.character(match.call()$Label),
           envir=.GlobalEnv))       { stop(paste0("Object ", match.call()$Label, " not found")); par(pardef) }
        Label  <- as.factor(Label)
        if(length(Label) != n.obs)  { stop(paste0("Labels must be a factor of length N=",  n.obs)); par(pardef) }
      } else {
        Label  <- 1
      }
      if(n.fac != 1) {
        plot(plot.x[,ind[1]], plot.x[,ind[2]], type=type, col=as.numeric(Label),
             xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind[2]))
        title(main=list("Posterior Scores", cex=cex.t))
        if(type == "n") text(plot.x[,ind[1]], plot.x[,ind[2]], 1:n.obs, 
                             col=as.numeric(Label), cex=0.5)
      } else {
        plot(plot.x[,ind], type=type, col=as.numeric(Label),
             xlab="Observation", ylab=paste0("Factor ", ind))
        title(main=list("Posterior Scores", cex=cex.t))
        if(type == "n") text(plot.x[,ind], col=as.numeric(Label), cex=0.5)
      }
    }
    if(var == "loadings") {
      plot.x   <- results$post.load
      if(heat) {
        image(z=t(plot.x[n.var:1,1:n.fac]), xlab="", 
              ylab="", xaxt="n", yaxt="n")
        title(main=list("Posterior Loadings", cex=cex.t))
        axis(1, cex.axis=0.8, line=-0.5, tick=F, 
             at=if(n.fac != 1) seq(0, 1, 1/(n.fac - 1)) else 0, labels=1:n.fac)
        if(n.var < 100) {
          axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
               at=seq(0, 1, 1/(n.var - 1)), 
               labels=substring(rownames(plot.x)[n.var:1], 1, 10))
        }
        box(lwd=2)
        mtext("Factors", side=1, line=2)
        if(n.fac != 1) abline(v=seq(1/(2 * (n.fac - 1)), 
                                    1 - 1/(2 * (n.fac - 1)), 
                                    1/(n.fac - 1)), lty=2, lwd=1)
    } else {
      if(ind[1] > n.var)            { stop(paste0("Only the first ", n.var, " variables can be plotted")); par(pardef) }
      if(n.fac != 1) {
        plot(plot.x[,ind[1]], plot.x[,ind[2]], type=type, 
             xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind[2]))
        title(main=list("Posterior Loadings", cex=cex.t))
        if(type == "n") text(plot.x[,ind[1]], plot.x[,ind[2]], rownames(plot.x), cex=0.5)
      } else { 
        plot(plot.x[,ind], type=type,
             xlab="Variable", ylab=paste0("Factor ", ind))
        title(main=list("Posterior Loadings", cex=cex.t))
        if(type == "n") text(plot.x[,ind], rownames(plot.x), cex=0.5)
        }
      }
    }
    if(var == "uniquenesses") {
      plot.x     <- results$post.psi
      plot(plot.x, type=type, ylab="Uniquenesses", xlab="Variable")
      title(main=list("Posterior Uniquenesses", cex=cex.t))
      if(type    == "n") text(1:length(plot.x), plot.x, names(plot.x), cex=0.5)
    }
    if(!ind.x)             ind <- x.ind
  } 
  
  if(m.sw["Q.sw"]) {
    Q.res    <- results$Q.results
    print(Q.res[-length(Q.res)])
    range.Q  <- attr(Q.res, "Factors")
    if(method  == "IFA") {
      par(mfrow =  c(1, 2))
      plot.Q     <- Q.res$Counts
      col.Q      <- c("black", "red")[(names(plot.Q) == n.fac) + 1]
      Q.plot     <- barplot(plot.Q, ylab="Frequency", xlab="Q", xaxt="n", col=col.Q)
      title(main=list("Posterior Distribution of Q", cex=cex.t))
      axis(1, at=Q.plot, labels=names(plot.Q), tick=F)
    } else if(length(range.Q) > 1) {
      par(mfrow = c(1, 2))
      plot.Q     <- Q.res$prop.exp
      plot(plot.Q, type="l", xlab="# Factors", ylim=c(0,1),
           ylab="% Variation Explained", xaxt="n", yaxt="n")
      title(main=list("Scree Plot to Choose Q", cex=cex.t))
      axis(1, at=1:length(plot.Q), labels=range.Q)
      axis(2, at=seq(0, 1, 0.1), labels=seq(0, 100, 10), cex.axis=0.8)
    }
    plot.x       <- results$cum.var
    prop.exp     <- plot.x[n.fac]
    if(n.fac      > 1) {
      plot(plot.x, type="l", xlab="# Factors", ylim=c(0,1),
           ylab="% Variation Explained", xaxt="n", yaxt="n")
      title(main=list(paste0("Cumulative Variance:\n", n.fac, " Factors"), cex=cex.t))
      axis(1, at=1:length(plot.x), labels=1:n.fac)
      axis(2, at=seq(0, 1, 0.1), labels=seq(0, 100, 10), cex.axis=0.8) 
    }
    cat(paste0("Proportion of Variation Explained = ",
               round(prop.exp[length(prop.exp)]*100, 2), "%", "\n"))
    if(max(prop.exp) > 1)             warning("chain may not have converged")
    par(mfrow = c(1, 1))
  }

  if(m.sw["cor.sw"]) {
    if(var == "scores" || 
         var == "loadings") {
      if(ind.x)            ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)           { stop("Length of indexes for plotting cannot be greater than 2"); par(pardef) }
      if(var == "scores") {
        if(ind[1] >  n.obs)         { stop(paste0("First index cannot be greater than ",  n.obs)); par(pardef) }
      } else if(ind[1] > n.var)     { stop(paste0("First index cannot be greater than ",  n.var)); par(pardef) }
      if(ind[2]   >  n.fac)         { stop(paste0("Second index cannot be greater than ", n.fac)); par(pardef) }
    } else {
      if(ind.x)            ind <- 1
      if(length(ind) >   1)         { stop("Length of indexes for plotting cannot be greater than 1"); par(pardef) }
      if(ind      >  n.var)         { stop(paste0("Index cannot be greater than ", n.var)); par(pardef) }
    }
    if(var == "means") {
      plot.x   <- results$means 
      acf(plot.x[ind,], main="")
      title(main=list(paste0("ACF:\n Mean of ", rownames(plot.x)[ind], " Variable"), cex=cex.t))
    }
    if(var == "scores") { 
      plot.x   <- results$scores
      acf(plot.x[ind[1],ind[2],], main="")
      title(main=list(paste0("ACF:\n Scores - Observation ", ind[1], ", Factor ", ind[2]), cex=cex.t))
    }
    if(var == "loadings") { 
      plot.x   <- results$loadings
      acf(plot.x[ind[1],ind[2],], main="")
      title(main=list(paste0("ACF:\n Loadings - ", rownames(plot.x)[ind[1]], " Variable, Factor ", ind[2]), cex=cex.t))
    }
    if(var == "uniquenesses") { 
      plot.x   <- results$uniquenesses
      acf(plot.x[ind,], main="")
      title(main=list(paste0("ACF:\n Uniqueness of ", rownames(plot.x)[ind], " Variable"), cex=cex.t))
    }
    if(!ind.x)             ind <- x.ind
  }
  par(pardef)
}