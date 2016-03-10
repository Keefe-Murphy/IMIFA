################################
### IMIFA Plotting Functions ###
################################

plot.IMIFA  <- function(results = NULL, plot.meth = c("all", "correlation", "density", "posterior", "Q", "trace"), 
                        vars = c("means", "scores", "loadings", "uniquenesses"), Label = NULL, G = 1, Q = NULL,
                        fac = NULL, by.fac = T, ind = NULL, type = c("h", "n", "p", "l"), mat = T, ... ) {
 
  defpar    <- par(no.readonly = T)
  defop     <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  par(cex.axis=0.8, new=F)
  if(missing(results))                stop("Results must be supplied")
  if(!exists(deparse(substitute(results)),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "IMIFA")       stop(paste0("Results object of class 'IMIFA' must be supplied"))
  GQ.res    <- results$GQ.results
  n.grp     <- GQ.res$G
  n.fac     <- attr(GQ.res, "Factors")
  G.supp    <- attr(GQ.res, "Supplied")["G"]
  Q.supp    <- attr(GQ.res, "Supplied")["Q"]
  method    <- attr(results, "Method")
  if(all(is.element(method, c("MFA", "MIFA")),
         !is.element(G, 1:n.grp)))    stop(paste0("G cannot be greater than the number of groups in ", match.call()$results))
  G.ind     <- which(n.grp == G)
  result    <- results[[G.ind]]
  if(missing(Q))             Q <- GQ.res$Q[G]
  if(Q       > GQ.res$Q[G])           stop("can't plot this many factors")
  n.var     <- attr(results, "Vars")
  n.obs     <- attr(results, "Obs")
  if(missing(plot.meth))              stop("What type of plot would you like to produce?")
  plot.meth <- match.arg(plot.meth)
  type.x    <- missing(type)
  type      <- match.arg(type)
  m.sw      <- c(Q.sw = F, C.sw = F, D.sw = F, P.sw = F, T.sw = F)
  v.sw      <- attr(results, "Switch")
  names(v.sw)  <- formals(sys.function(sys.parent()))$vars
  missing   <- missing(vars)
  vars      <- match.arg(vars)
  var.names <- rownames(result$post.load)
  obs.names <- rownames(result$post.f)
  all.ind   <- plot.meth == "all"
  grp.ind   <- !is.element(method, c("FA", "IFA"))
  if(all.ind)   {
    if(v.sw[vars]) {
      m.sw[-1]  <- !m.sw[-1]
      if(vars   == "loadings") {
        layout(matrix(c(1, 2, 3, 4, 3, 5), nr=3, nc=2, byrow = TRUE))
      } else {
        layout(matrix(c(1, 2, 3, 4), nr=2, nc=2, byrow = TRUE))
      }
      par(oma=c(0, 0, 2, 0), mai=c(0.7, 0.7, 0.5, 0.2), mgp=c(2, 1, 0), cex=0.8)
    }
  } else {
    sw.n    <- paste0(toupper(substring(plot.meth, 1, 1)), ".sw")
    m.sw[sw.n]  <- T
  }
  if(all(!m.sw["Q.sw"],
     missing(vars)))                  stop("What variable would you like to plot?")
  if(all(any(vars == "scores",
     vars   == "loadings"),
     Q      == 0))                    stop(paste0("Can't plot ", vars, " as they contain no columns/factors"))
  if(all(any(m.sw["P.sw"], all.ind),
     any(vars   == "means",
         vars   == "uniquenesses"),
     !v.sw[vars])) {  
    if(all.ind)                       warning(paste0("Can only plot posterior mean, as ", vars, " weren't stored"), call.=F)
   v.sw[vars]   <- !v.sw[vars]
   all.ind      <- F
   m.sw["P.sw"] <- T
  } 
  if(all(!v.sw[vars],
     !m.sw["Q.sw"]))                  stop(paste0(vars, " weren't stored"))
  if(!is.logical(mat))                stop("mat must be TRUE or FALSE")
  if(!missing(ind))      x.ind <- ind
  ind.x     <- !exists("x.ind", envir=environment())
  if(any(all(any(vars == "means", 
     vars   == "uniquenesses"),
     !missing(ind)),
     Q      == 1))         mat <- F
  
  if(m.sw["T.sw"]) {
    if(any(vars == "scores",
           vars == "loadings")) {
      if(ind.x)            ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)             stop("Length of plotting indices can't be greater than 2")
      if(vars == "scores") {
        if(ind[1] > n.obs)            stop(paste0("First index can't be greater than the number of observations - ",  n.obs))
      } else if(ind[1] > n.var)       stop(paste0("First index can't be greater than the number of variables - ",  n.var))
      if(ind[2]   > Q)                stop(paste0("Second index can't be greater than the number of factors - ", Q))
    } else {
      if(ind.x)            ind <- 1
      if(length(ind) >  1)            stop("Length of plotting indices can't be greater than 1")
      if(ind      > n.var)            stop(paste0("Index can't be greater than the number of variables - ", n.var))
    }
    if(!mat)              iter <- 1:attr(results, "Store")
    
    if(vars == "means") {
      plot.x   <- result$means
      if(mat) {
        matplot(t(plot.x[,]), type="l", ylab="Means", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
        title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
      } else {
        plot(x=iter, y=plot.x[ind,], type="l", ylab="Mean", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
        title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind], " Variable")))
      }
    }
    if(vars == "scores") {
      plot.X   <- result$scores
      if(by.fac) {
        plot.x <- plot.X[ind[1],,]
      } else {
        plot.x <- plot.X[,ind[2],]
      }
      if(mat) {
        matplot(t(plot.x), type="l", ylab="Scores", xlab="Iteration")    
        if(by.fac) {
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Observation ", obs.names[ind[1]])))))
        } else {
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Factor ", ind[2])))))
        }
      } else {
        plot(x=iter, y=plot.X[ind[1],ind[2],], type="l", ylab="Scores", xlab="Iteration")
        title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))))
      }
    }
    if(vars == "loadings") {
      plot.X   <- result$loadings
      if(by.fac) {
        plot.x <- plot.X[ind[1],,]
      } else {
        plot.x <- plot.X[,ind[2],]
      }
      if(mat) {
        matplot(t(plot.x), type="l", ylab="Loadings", xlab="Iteration")
        if(by.fac) {
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), var.names[ind[1]], " Variable")))))
        } else {
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Factor ", ind[2])))))
        }
      } else {
        plot(x=iter, y=plot.X[ind[1],ind[2],], type="l", ylab="Loadings", xlab="Iteration")
        title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), var.names[ind[1]], " Variable, Factor ", ind[2])))))
      }
    }
    if(vars == "uniquenesses") {
      plot.x   <- result$uniquenesses
      if(mat) {
        matplot(t(plot.x[,]), type="l", ylab="Uniquenesses", xlab="Iteration")
        title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
      } else {
        plot(x=iter, y=plot.x[ind,], ylab="Uniquenesses", type="l", xlab="Iteration")
        title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind], " Variable")))
      }
    }
    if(!ind.x)             ind <- x.ind
    if(all.ind) title(paste0(toupper(substr(vars, 1, 1)),
                      substr(vars, 2, nchar(vars)), 
                      ifelse(grp.ind, paste0(" - Group ", G), "")), outer=T)
  }
  
  if(m.sw["D.sw"]) {
    if(any(vars == "scores", 
           vars == "loadings")) {
      if(ind.x)            ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)             stop("Length of plotting indices can't be greater than 2")
      if(vars == "scores") {
        if(ind[1] > n.obs)            stop(paste0("First index can't be greater than the number of observations - ",  n.obs))
      } else if(ind[1] > n.var)       stop(paste0("First index can't be greater than the number of variables - ",  n.var))
      if(ind[2]   > Q)                stop(paste0("Second index can't be greater than the number of factors - ", Q))
    } else {
      if(ind.x)            ind <- 1
      if(length(ind) >  1)            stop("Length of plotting indices can't be greater than 1")
      if(ind      > n.var)            stop(paste0("Index can't be greater than the number of variables - ", n.var))
    }
    if(vars == "means") {
      plot.X   <- result$means
      if(mat) {
        plot.x <- apply(plot.X, 1, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
      } else {
        plot.d <- density(plot.X[ind,])
        plot(plot.d, main="")
        title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind], " Variable")))
        polygon(plot.d, col="black")
      }
    }
    if(vars == "scores") {
      plot.X   <- result$scores
      if(by.fac) {
        plot.x <- plot.X[ind[1],,]
      } else {
        plot.x <- plot.X[,ind[2],]
      }
      if(mat) {
        plot.x <- apply(plot.x, 1, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        if(by.fac) {
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Observation ", obs.names[ind[1]])))))
        } else {
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Factor ", ind[2])))))
        }
      } else {
        plot.d <- density(plot.X[ind[1],ind[2],])
        plot(plot.d, main="")
        title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))))
        polygon(plot.d, col="black")
      }
    }
    if(vars == "loadings") {
      plot.X   <- result$loadings
      if(by.fac) {
        plot.x <- plot.X[ind[1],,]
      } else {
        plot.x <- plot.X[,ind[2],]
      }
      if(mat) {
        plot.x <- apply(plot.x, 1, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        if(by.fac) {
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), var.names[ind[1]], " Variable")))))
        } else {
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), "Factor ", ind[2])))))
        }
      } else {
        plot.d <- density(plot.X[ind[1],ind[2],])
        plot(plot.d, main="")
        title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""), var.names[ind[1]], " Variable, Factor ", ind[2])))))
        polygon(plot.d, col="black")
      }
    }
    if(vars == "uniquenesses") {
      plot.X   <- result$uniquenesses
      if(mat) {
        plot.x <- apply(plot.X, 1, density)
        plot.x <- sapply(plot.x, "[[", "y")
        matplot(plot.x, type="l", ylab="Density")
        title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
      } else {
        plot.d <- density(plot.X[ind,])
        plot(plot.d, main="")
        title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind], " Variable")))
        polygon(plot.d, col="black")
      }
    }
    if(!ind.x)             ind <- x.ind
  }
  
  if(m.sw["P.sw"]) {
    if(any(vars == "scores", 
           vars == "loadings")) {
      if(ind.x) {
       if(vars == "scores") {
           ind <- c(1, 2)
       } else {
           ind <- c(1, 1)
       } 
      }
      if(!missing(fac)) {
       if(vars == "loadings") {
        ind[2] <- fac
       }
        ind[2] <- max(fac, 2)
      }
      if(Q  == 1) {
        ind    <- 1
      } else {
        if(length(ind) > 2)           stop("Only two columns can be plotted")
        if(ind[1] > Q)                stop(paste0("Only the first ", Q, " columns can be plotted"))
        if(ind[2] > Q)                stop(paste0("Only the first ", Q, " columns can be plotted"))
      }
    }
    if(vars == "means") {
      plot.x   <- result$post.mu
      plot(plot.x, type=type, ylab="Means", xlab="Variable", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
      title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
      if(type  == "n") text(x=1:length(plot.x), y=plot.x, var.names, cex=0.5)
    }
    if(vars == "scores") {
      plot.x   <- result$post.f
      if(ind[1] > n.obs)              stop(paste0("Only the first ", n.obs, " scores can be plotted"))
      Labs     <- 1
      if(!missing(Label)) {
        if(!exists(as.character(match.call()$Label),
            envir=.GlobalEnv)) {      warning(paste0("Object ", match.call()$Label, " not found"), call.=F)
        } else {
          Labs <- as.factor(Label)
          if(length(Labs) != n.obs)   stop(paste0("Labels must be a factor of length N=",  n.obs))
        }
      }
      type.f   <- ifelse(any(type.x, type == "l"), "p", type)
      if(Q     != 1) {
        plot(plot.x[,ind[1]], plot.x[,ind[2]], type=type.f, col=as.numeric(Labs), 
             xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind[2]))
        title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nScores", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
        if(type.f == "n") text(plot.x[,ind[1]], plot.x[,ind[2]], obs.names, 
                             col=as.numeric(Labs), cex=0.5)
      } else {
        plot(plot.x[,ind], type=type.f, col=as.numeric(Labs), 
             xlab="Observation", ylab=paste0("Factor ", ind))
        title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nScores", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
        if(type.f == "n") text(plot.x[,ind], col=as.numeric(Labs), cex=0.5)
      }
    }
    if(vars == "loadings") {
      if(!all.ind) {
        par(mfrow=c(1, 2), oma=c(0, 0, 1, 0), mai=c(1.25, 1, 0.75, 0.5))
      }
      plot.x   <- result$post.load
      image(z=t(plot.x[n.var:1,1:Q]), xlab="", 
            ylab="", xaxt="n", yaxt="n")
      title(main=list(paste0("Posterior Mean Heatmap")))
      axis(1, line=-0.5, tick=F, 
           at=if(Q != 1) seq(0, 1, 1/(Q - 1)) else 0, labels=1:Q)
      if(n.var < 100) {
        axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
             at=seq(0, 1, 1/(n.var - 1)), 
             labels=substring(var.names[n.var:1], 1, 10))
      }
      box(lwd=2)
      mtext("Factors", side=1, line=2)
      if(Q != 1) abline(v=seq(1/(2 * (Q - 1)), 
                                  1 - 1/(2 * (Q - 1)), 
                                  1/(Q - 1)), lty=2, lwd=1)
      if(by.fac) {
        if(Q != 1) {
          plot(plot.x[,ind[2]], type=type, xaxt="n", xlab="", ylab="Loading")
          axis(1, line=-0.5, tick=F, at=1:n.var, labels=1:n.var)
          mtext("Variable #", side=1, line=2)
          title(main=list(paste0("Factor ", ind[2])))
          if(type == "n") text(x=plot.x, var.names, cex=0.5)
        } else {
          plot(plot.x[,ind], type=type, xaxt="n", xlab="", ylab="Loading")
          axis(1, line=-0.5, tick=F, at=1:n.var, labels=1:n.var)
          mtext("Variable #", side=1, line=2)
          title(main=list(paste0("Factor ", ind)))
          if(type == "n") text(x=plot.x[,ind], var.names, cex=0.5)
        }
      } else {
        plot(plot.x[ind[1],], type=type, xaxt="n", xlab="", ylab="Loading")
        axis(1, line=-0.5, tick=F, at=1:Q, labels=1:Q)
        mtext("Factors", side=1, line=2)
        title(main=list(paste0(var.names[ind[1]], " Variable")))
        if(type == "n") text(x=plot.x[ind[1],], paste0("Factor ", 1:Q), cex=0.5)
      }
      if(!all.ind) {
        title(paste0("Loadings", ifelse(grp.ind, paste0(" - Group ", G), "")), outer=T)
      }
    }
    if(vars == "uniquenesses") {
      plot.x   <- result$post.psi
      plot(plot.x, type=type, ylab="Uniquenesses", xlab="Variable")
      title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", G), ""))))))
      if(type  == "n") text(1:length(plot.x), plot.x, var.names, cex=0.5)
    }
    if(!ind.x)             ind <- x.ind
  } 
  
  if(m.sw["Q.sw"]) {
    no.fac   <- all(Q == 0)
    if(all(method == "FA", any(v.sw["loadings"], no.fac))) {
      bic    <- round(GQ.res$BIC, 2)
    }
    if(method  == "IFA") {
      if(all(!Q.supp, Q > 1, v.sw["loadings"], result$prop.exp <= 1)) {
        par(mfrow = c(1, 2))
      }
      if(!Q.supp) {
        plot.Q <- GQ.res$Counts
        col.Q  <- c("black", "red")[(names(plot.Q) == Q) + 1]
        Q.plot <- barplot(plot.Q, ylab="Frequency", xlab="Q", xaxt="n", col=col.Q)
        title(main=list("Posterior Distribution of Q"))
        axis(1, at=Q.plot, labels=names(plot.Q), tick=F) 
      }
    } else {
      if(all(Q > 1, length(n.fac) > 1, v.sw["loadings"])) {
        par(mfrow = c(1, 2))
      }     
      plot.Q     <- GQ.res$cum.var
      if(length(n.fac) > 1) {
        plot.Q   <- plot.Q[!is.na(plot.Q)]  
        if(all(length(plot.Q) > 1, !any(plot.Q > 1))) {
          plot(plot.Q, type="l", xlab="# Factors", ylim=c(0,1),
               ylab="% Variation Explained", xaxt="n", yaxt="n")
          title(main=list("Scree Plot to Choose Q"))
          axis(1, at=1:length(plot.Q), labels=n.fac)
          axis(2, at=seq(0, 1, 0.1), labels=seq(0, 100, 10), las=1)
        }
      } 
    }
    plot.x       <- result$cum.var
    if(length(plot.x) == 1) {
      prop.exp   <- result$prop.exp
    } else {
      prop.exp   <- plot.x[max(1, Q)]
      if(all(Q > 1, !any(plot.x > 1))) {
        plot(plot.x, type="l", xlab="# Factors", ylim=c(0,1),
             ylab="% Variation Explained", xaxt="n", yaxt="n")
        title(main=list(paste0("Cumulative Variance:\n", Q, " Factors")))
        axis(1, at=1:length(plot.x), labels=1:Q)
        axis(2, at=seq(0, 1, 0.1), labels=seq(0, 100, 10), las=1) 
      }
    } 
    if(all(!exists("Q.plot",  envir=environment()),
           length(plot.x) == 1,
           length(Q)      == 1)) {
                                      warning("Nothing to plot", call.=F)
    }
    if(method == "IFA") {
        print(GQ.res[2:length(GQ.res)])
    } else {
        cat(paste0("Q = ", Q, "\n"))
    }
    if(all(method == "FA", any(v.sw["loadings"], no.fac))) {
        cat(paste0("BIC = ", bic[which.max(bic)], "\n"))
    }
        cat(paste0("Proportion of Variation Explained = ",
            round(prop.exp[length(prop.exp)] * 100, 2), "%\n"))
    if(max(prop.exp) > 1)             warning("Chain may not have converged", call.=F)
  }

  if(m.sw["C.sw"]) {
    if(any(vars == "scores",
           vars == "loadings")) {
      if(ind.x)            ind <- c(1, 1)
      if(!missing(fac)) ind[2] <- fac
      if(length(ind) > 2)             stop("Length of plotting indices can't be greater than 2")
      if(vars == "scores") {
        if(ind[1] > n.obs)            stop(paste0("First index can't be greater than the number of observations -",  n.obs))
      } else if(ind[1] > n.var)       stop(paste0("First index can't be greater than the number of variables - ",  n.var))
      if(ind[2] > Q)                  stop(paste0("Second index can't be greater than the number of factors - ", Q))
    } else {
      if(ind.x)            ind <- 1
      if(length(ind) > 1)             stop("Length of plotting indices can't be greater than 1")
      if(ind    >  n.var)             stop(paste0("Index can't be greater than the number of variables - ", n.var))
    }
    if(vars == "means") {
      plot.x   <- result$means 
      acf(plot.x[ind,], main="")
      title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind], " Variable")))
    }
    if(vars == "scores") { 
      plot.x   <- result$scores
      acf(plot.x[ind[1],ind[2],], main="")
      title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\nScores - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
    }
    if(vars == "loadings") { 
      plot.x   <- result$loadings
      acf(plot.x[ind[1],ind[2],], main="")
      title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\n Loadings - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
    }
    if(vars == "uniquenesses") { 
      plot.x   <- result$uniquenesses
      acf(plot.x[ind,], main="")
      title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", G, " - "), ""))), var.names[ind], " Variable")))
    }
    if(!ind.x)             ind <- x.ind
  }
}