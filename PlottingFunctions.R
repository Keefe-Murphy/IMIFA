################################
### IMIFA Plotting Functions ###
################################

plot.IMIFA     <- function(results = NULL, plot.meth = c("all", "correlation", "density", "posterior", "GQ", "trace"), 
                           vars = c("means", "scores", "loadings", "uniquenesses"), Label = NULL, fac = NULL, g = NULL,
                           by.fac = T, ind = NULL, type = c("h", "n", "p", "l"), mat = T, ... ) {
 
  defpar  <- par(no.readonly = T)
  defop   <- options()
  options(warn=1)
  par(cex.axis=0.8, new=F)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(missing(results))                stop("Results must be supplied")
  if(!exists(deparse(substitute(results)),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "IMIFA")       stop(paste0("Results object of class 'IMIFA' must be supplied"))
  GQ.res  <- results$GQ.results
  G       <- GQ.res$G
  Qs      <- GQ.res$Q
  n.grp   <- attr(GQ.res, "Groups")
  n.fac   <- attr(GQ.res, "Factors")
  G.supp  <- attr(GQ.res, "Supplied")["G"]
  Q.supp  <- attr(GQ.res, "Supplied")["Q"]
  method  <- attr(results, "Method")
  vars    <- match.arg(vars)
  n.var   <- attr(results, "Vars")
  n.obs   <- attr(results, "Obs")
  if(missing(plot.meth))              stop("What type of plot would you like to produce?")
  if(is.element(plot.meth, 
     c("G", "Q", 
       "QG"))) {       plot.meth <- "GQ"
  }
  plot.meth    <- match.arg(plot.meth)
  type.x       <- missing(type)
  type         <- match.arg(type)
  m.sw         <- c(G.sw = F, C.sw = F, D.sw = F, P.sw = F, T.sw = F)
  v.sw         <- attr(results, "Switch")
  names(v.sw)  <- formals(sys.function(sys.parent()))$vars
  var.names    <- rownames(results[[1]]$post.load)
  obs.names    <- rownames(results$Scores$post.f)
  all.ind <- plot.meth == "all"
  grp.ind <- all(G != 1, !is.element(method, c("FA", "IFA")))
  if(grp.ind) {
    clust <- results$Clust
  }
  if(all.ind)   {
    if(v.sw[vars]) {
      m.sw[-1]    <- !m.sw[-1]
      if(vars  == "loadings") {
        layout(matrix(c(1, 2, 3, 4, 3, 5), nr=3, nc=2, byrow = TRUE))
      } else {
        layout(matrix(c(1, 2, 3, 4), nr=2, nc=2, byrow = TRUE))
      }
      par(cex=0.8, mai=c(0.7, 0.7, 0.5, 0.2), mgp=c(2, 1, 0), oma=c(0, 0, 2, 0))
    }
  } else {
    sw.n  <- paste0(toupper(substring(plot.meth, 1, 1)), ".sw")
    m.sw[sw.n]    <- T
  }
  if(all(!m.sw["G.sw"],
     missing(vars)))                  stop("What variable would you like to plot?")
  if(all(any(m.sw["P.sw"], all.ind),
     is.element(vars, c("means", "uniquenesses")),
     !v.sw[vars])) {  
    if(all.ind)                       warning(paste0("Can only plot posterior mean, as ", vars, " weren't stored"), call.=F)
   v.sw[vars]     <- !v.sw[vars]
   all.ind        <- F
   m.sw["P.sw"]   <- T
  } 
  if(all(!v.sw[vars],
     !m.sw["G.sw"]))                  stop(paste0(vars, " weren't stored"))
  if(!is.logical(mat))                stop("mat must be TRUE or FALSE")
  indx    <- missing(ind)
  if(!indx)                 xind <- ind
  if(any(vars  == "scores",
         m.sw["G.sw"])) {
    Gs    <- 1
  } else if(!missing(g)) {
    if(g > 1)                         warning(paste0("g must be equal to 1 for the ", method, " method"), call.=F)
    if(all(is.element(method, c("MFA", "MIFA")),
      !is.element(g, seq_len(G))))    stop("This g value was not used during simulation")
    Gs    <- g
    rm(g)
  } else if(!interactive()) {         stop("g must be supplied for non-interactive sessions")
  } else {
    Gs    <- seq_len(G)
  }
  
  for(g in Gs) {
    Q     <- Qs[g]
    msg   <- "Hit <Return> to see next plot: "
    msgx  <- all(interactive(), g != max(Gs))
    result     <- results[[g]]
    if(any(all(Q  == 0,  vars == "loadings"),
       all(all(Qs == 0), vars == "scores"))) {            
                                       warning(paste0("Can't plot ", vars, " as they contain no columns/factors"), call.=F)
      if(length(unique(tail(Qs, - g))) == 1) {
        break
      } else {
        ifelse(msgx, readline(msg), "")
        next
      }
    }
    if(any(all(is.element(vars, c("means", "uniquenesses")),
               !indx),
           all(is.element(vars, c("scores", "loadings")),
               Q  == 1))) { matx <- F
    } else {
      matx     <- mat
    }                           
    
    if(m.sw["T.sw"]) {
      if(is.element(vars, c("scores", "loadings"))) {
        if(indx)             ind <- c(1, 1)
        if(!missing(fac)) ind[2] <- fac
        if(length(ind) > 2)           stop("Length of plotting indices can't be greater than 2")
        if(vars   == "scores") {
          if(ind[1] > n.obs)          stop(paste0("First index can't be greater than the number of observations - ",  n.obs))
        } else if(ind[1] > n.var)     stop(paste0("First index can't be greater than the number of variables - ",  n.var))
        if(ind[2]   > Q)              stop(paste0("Second index can't be greater than the number of factors - ", Q))
      } else {
        if(indx)             ind <- 1
        if(length(ind) >  1)          stop("Length of plotting indices can't be greater than 1")
        if(ind      > n.var)          stop(paste0("Index can't be greater than the number of variables - ", n.var))
      }
      if(!matx)             iter <- seq_len(attr(result, "Store"))
      
      if(vars  == "means") {
        plot.x <- result$means
        if(matx) {
          matplot(t(plot.x[,]), type="l", ylab="Means", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot(x=iter, y=plot.x[ind,], type="l", ylab="Mean", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(vars  == "scores") {
        plot.X <- results$Scores$f
        if(by.fac) {
          plot.x  <- plot.X[ind[1],,]
        } else {
          plot.x  <- plot.X[,ind[2],]
        }
        if(matx) {
          matplot(t(plot.x), type="l", ylab="Scores", xlab="Iteration")    
          if(by.fac) {
            title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else {
          plot(x=iter, y=plot.X[ind[1],ind[2],], type="l", ylab="Scores", xlab="Iteration")
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
        }
      }
      if(vars  == "loadings") {
        plot.X <- result$loadings
        if(by.fac) {
          plot.x  <- plot.X[ind[1],,]
        } else {
          plot.x  <- plot.X[,ind[2],]
        }
        if(matx) {
          matplot(t(plot.x), type="l", ylab="Loadings", xlab="Iteration")
          if(by.fac) {
            title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else {
          plot(x=iter, y=plot.X[ind[1],ind[2],], type="l", ylab="Loadings", xlab="Iteration")
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$uniquenesses
        if(matx) {
          matplot(t(plot.x[,]), type="l", ylab="Uniquenesses", xlab="Iteration")
          title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot(x=iter, y=plot.x[ind,], ylab="Uniquenesses", type="l", xlab="Iteration")
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(!indx)              ind <- xind
      if(all.ind) title(paste0(toupper(substr(vars, 1, 1)),
                        substr(vars, 2, nchar(vars)), 
                        ifelse(all(grp.ind, vars != "scores"), 
                               paste0(" - Group ", g), "")), outer=T)
    }
    
    if(m.sw["D.sw"]) {
      if(is.element(vars, c("scores", "loadings"))) {
        if(indx)             ind <- c(1, 1)
        if(!missing(fac)) ind[2] <- fac
        if(length(ind) > 2)           stop("Length of plotting indices can't be greater than 2")
        if(vars   == "scores") {
          if(ind[1] > n.obs)          stop(paste0("First index can't be greater than the number of observations - ",  n.obs))
        } else if(ind[1] > n.var)     stop(paste0("First index can't be greater than the number of variables - ",  n.var))
        if(ind[2]   > Q)              stop(paste0("Second index can't be greater than the number of factors - ", Q))
      } else {
        if(indx)             ind <- 1
        if(length(ind) >  1)          stop("Length of plotting indices can't be greater than 1")
        if(ind      > n.var)          stop(paste0("Index can't be greater than the number of variables - ", n.var))
      }
      if(vars  == "means") {
        plot.X <- result$means
        if(matx) {
          plot.x  <- apply(plot.X, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot.d  <- density(plot.X[ind,])
          plot(plot.d, main="")
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col="black")
        }
      }
      if(vars  == "scores") {
        plot.X <- results$Scores$f
        if(by.fac) {
          plot.x  <- plot.X[ind[1],,]
        } else {
          plot.x  <- plot.X[,ind[2],]
        }
        if(matx) {
          plot.x  <- apply(plot.x, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(by.fac) {
            title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else {
          plot.d  <- density(plot.X[ind[1],ind[2],])
          plot(plot.d, main="")
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
          polygon(plot.d, col="black")
        }
      }
      if(vars  == "loadings") {
        plot.X <- result$loadings
        if(by.fac) {
          plot.x  <- plot.X[ind[1],,]
        } else {
          plot.x  <- plot.X[,ind[2],]
        }
        if(matx) {
          plot.x  <- apply(plot.x, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(by.fac) {
            title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else {
          plot.d  <- density(plot.X[ind[1],ind[2],])
          plot(plot.d, main="")
          title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
          polygon(plot.d, col="black")
        }
      }
      if(vars  == "uniquenesses") {
        plot.X <- result$uniquenesses
        if(matx) {
          plot.x  <- apply(plot.X, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot.d  <- density(plot.X[ind,])
          plot(plot.d, main="")
          title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col="black")
        }
      }
      if(!indx)              ind <- xind
    }
    
    if(m.sw["P.sw"]) {
      if(is.element(vars, c("scores", "loadings"))) {
        if(indx) {
         if(vars  == "scores") {
             ind  <- c(1, 2)
         } else {
             ind  <- c(1, 1)
         } 
        }
        if(!missing(fac)) {
         if(vars  == "loadings") {
          ind[2]  <- fac
         }
          ind[2]  <- max(fac, 2)
        }
        if(Q   == 1) {
          ind  <- 1
        } else {
          if(length(ind) > 2)         stop("Only two columns can be plotted")
          if(ind[1] > Q)              stop(paste0("Only the first ", Q, " columns can be plotted"))
          if(ind[2] > Q)              stop(paste0("Only the first ", Q, " columns can be plotted"))
        }
      }
      if(vars  == "means") {
        plot.x <- result$post.mu
        plot(plot.x, type=type, ylab="Means", xlab="Variable", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
        title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(x=seq_along(plot.x), y=plot.x, var.names, cex=0.5)
      }
      if(vars  == "scores") {
        plot.x <- results$Scores$post.f
        if(ind[1] > n.obs)            stop(paste0("Only the first ", n.obs, " scores can be plotted"))
        if(grp.ind) {
          Labs <- clust$post.z
        } else {
          Labs <- 1
        }
        if(!missing(Label)) {
          if(!exists(as.character(match.call()$Label),
              envir=.GlobalEnv)) {    warning(paste0("Object ", match.call()$Label, " not found"), call.=F)
          } else {
            Labs  <- as.factor(Label)
            if(length(Labs) != n.obs) stop(paste0("Labels must be a factor of length N=",  n.obs))
          }
        }
        ind2   <- ifelse(Q > 1, max(2, ind[2]), ind[2])
        type.f <- ifelse(any(type.x, type == "l"), "p", type)
        if(Q   != 1) {
          plot(plot.x[,ind[1]], plot.x[,ind2], type=type.f, col=as.numeric(Labs), 
               xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
          title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"))))
          if(type.f == "n") text(plot.x[,ind[1]], plot.x[,ind2], obs.names, 
                               col=as.numeric(Labs), cex=0.5)
        } else {
          plot(plot.x[,ind], type=type.f, col=as.numeric(Labs), 
               xlab="Observation", ylab=paste0("Factor ", ind))
          title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"))))
          if(type.f == "n") text(plot.x[,ind], col=as.numeric(Labs), cex=0.5)
        }
      }
      if(vars  == "loadings") {
        if(!all.ind) {
          par(mai=c(1.25, 1, 0.75, 0.5), mfrow=c(1, 2), oma=c(0, 0, 1, 0))
        }
        plot.x <- result$post.load
        image(z=t(plot.x[seq(n.var, 1),seq_len(Q)]), xlab="", 
              ylab="", xaxt="n", yaxt="n")
        title(main=list(paste0("Posterior Mean Heatmap")))
        axis(1, line=-0.5, tick=F, 
             at=if(Q != 1) seq(0, 1, 1/(Q - 1)) else 0, labels=seq_len(Q))
        if(n.var < 100) {
          axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
               at=seq(0, 1, 1/(n.var - 1)), 
               labels=substring(var.names[n.var:1], 1, 10))
        }
        box(lwd=2)
        mtext("Factors", side=1, line=2)
        if(Q   != 1) abline(v=seq(1/(2 * (Q - 1)), 
                                  1 - 1/(2 * (Q - 1)), 
                                  1/(Q - 1)), lty=2, lwd=1)
        if(by.fac) {
          if(Q != 1) {
            plot(plot.x[,ind[2]], type=type, xaxt="n", xlab="", ylab="Loading")
            axis(1, line=-0.5, tick=F, at=seq_len(n.var), labels=seq_len(n.var))
            mtext("Variable #", side=1, line=2)
            title(main=list(paste0("Factor ", ind[2])))
            if(type == "n") text(x=plot.x, var.names, cex=0.5)
          } else {
            plot(plot.x[,ind], type=type, xaxt="n", xlab="", ylab="Loading")
            axis(1, line=-0.5, tick=F, at=seq_len(n.var), labels=seq_len(n.var))
            mtext("Variable #", side=1, line=2)
            title(main=list(paste0("Factor ", ind)))
            if(type == "n") text(x=plot.x[,ind], var.names, cex=0.5)
          }
        } else {
          plot(plot.x[ind[1],], type=type, xaxt="n", xlab="", ylab="Loading")
          axis(1, line=-0.5, tick=F, at=seq_len(Q), labels=seq_len(Q))
          mtext("Factors", side=1, line=2)
          title(main=list(paste0(var.names[ind[1]], " Variable")))
          if(type == "n") text(x=plot.x[ind[1],], paste0("Factor ", seq_len(Q)), cex=0.5)
        }
        if(!all.ind) {
          title(paste0("Loadings", ifelse(grp.ind, paste0(" - Group ", g), "")), outer=T)
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$post.psi
        plot(plot.x, type=type, ylab="Uniquenesses", xlab="Variable")
        title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(seq_along(plot.x), plot.x, var.names, cex=0.5)
      }
      if(!indx)              ind <- xind
    } 
    
    if(m.sw["G.sw"]) {
      if(is.element(method, c("FA", "MFA"))) {
        aic    <- round(GQ.res$AIC, 2)
        bic    <- round(GQ.res$BIC, 2)
      }
      if(all(method == "IFA", !Q.supp)) {
        plot.Q <- GQ.res$Counts
        col.Q  <- c("black", "red")[(names(plot.Q) == Q) + 1]
        Q.plot <- barplot(plot.Q, ylab="Frequency", xlab="Q", xaxt="n", col=col.Q)
        title(main=list("Posterior Distribution of Q"))
        axis(1, at=Q.plot, labels=names(plot.Q), tick=F) 
      }  
      if(!exists("Q.plot",  envir=environment())) 
                                      warning("Nothing to plot", call.=F)
      if(method  == "MIFA") {
          print(GQ.res)
      } else if(method == "IFA") {
          print(tail(GQ.res, -1))
      } else if(method == "MFA") {
          print(head(GQ.res, -2))
      } else {
          cat(paste0("Q = ", Q, "\n"))
      }
      if(is.element(method, c("FA", "MFA"))) {
          cat(paste0("AIC = ", aic[which.max(aic)], "\n"))
          cat(paste0("BIC = ", bic[which.max(bic)], "\n"))
      }
    }
  
    if(m.sw["C.sw"]) {
      if(is.element(vars, c("scores", "loadings"))) {
        if(indx)             ind <- c(1, 1)
        if(!missing(fac)) ind[2] <- fac
        if(length(ind) > 2)           stop("Length of plotting indices can't be greater than 2")
        if(vars   == "scores") {
          if(ind[1] > n.obs)          stop(paste0("First index can't be greater than the number of observations -",  n.obs))
        } else if(ind[1] > n.var)     stop(paste0("First index can't be greater than the number of variables - ",  n.var))
        if(ind[2] > Q)                stop(paste0("Second index can't be greater than the number of factors - ", Q))
      } else {
        if(indx)             ind <- 1
        if(length(ind) > 1)           stop("Length of plotting indices can't be greater than 1")
        if(ind    >  n.var)           stop(paste0("Index can't be greater than the number of variables - ", n.var))
      }
      if(vars  == "means") {
        plot.x <- result$means 
        acf(plot.x[ind,], main="")
        title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
      }
      if(vars  == "scores") { 
        plot.x <- results$Scores$f
        acf(plot.x[ind[1],ind[2],], main="")
        title(main=list(paste0("ACF", ifelse(all.ind, ":\n", ":\n Scores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
      }
      if(vars  == "loadings") { 
        plot.x <- result$loadings
        acf(plot.x[ind[1],ind[2],], main="")
        title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
      }
      if(vars  == "uniquenesses") { 
        plot.x <- result$uniquenesses
        acf(plot.x[ind,], main="")
        title(main=list(paste0("ACF", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
      }
      if(!indx)              ind <- xind
    }
    ifelse(msgx, readline(msg), "")
  }
}