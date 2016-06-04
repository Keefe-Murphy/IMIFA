################################
### IMIFA Plotting Functions ###
################################

plot.IMIFA     <- function(results = NULL, plot.meth = c("all", "correlation", "density", "posterior", "GQ", "trace"), 
                           vars = c("means", "scores", "loadings", "uniquenesses"), Labels = NULL, load.meth = c("all", "heatmap", "raw"),
                           g = NULL, fac = NULL, by.fac = T, ind = NULL, type = c("h", "n", "p", "l"), intervals = T, mat = T, partial = F, titles = T) {

  defpar  <- suppressWarnings(par(no.readonly = T))
  defpar$new   <- F
  defop   <- options()
  options(warn=1)
  suppressWarnings(par(cex.axis=0.8, new=F))
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(missing(results))                stop("Results must be supplied")
  if(!exists(deparse(substitute(results)),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "IMIFA")       stop(paste0("Results object of class 'IMIFA' must be supplied"))
  GQ.res  <- results$GQ.results
  G       <- GQ.res$G
  Qs      <- GQ.res$Q
  Q.max   <- max(Qs)
  n.grp   <- attr(GQ.res, "Groups")
  n.fac   <- attr(GQ.res, "Factors")
  G.supp  <- attr(GQ.res, "Supplied")["G"]
  Q.supp  <- attr(GQ.res, "Supplied")["Q"]
  method  <- attr(results, "Method")
  store   <- attr(results, "Store")
  vars    <- match.arg(vars)
  n.var   <- attr(results, "Vars")
  n.obs   <- attr(results, "Obs")
  if(missing(plot.meth))              stop("What type of plot would you like to produce?")
  if(is.element(plot.meth, 
     c("G", "Q", 
       "QG")))  {      plot.meth <- "GQ"
  }
  plot.meth    <- match.arg(plot.meth)
  load.meth    <- match.arg(load.meth)
  type.x       <- missing(type)
  type         <- match.arg(type)
  m.sw         <- c(G.sw = F, C.sw = F, D.sw = F, P.sw = F, T.sw = F)
  v.sw         <- attr(results, "Switch")
  names(v.sw)  <- formals(sys.function(sys.parent()))$vars
  ci.sw        <- v.sw
  var.names    <- rownames(results[[1]]$post.load)
  obs.names    <- rownames(results$Scores$post.f)
  all.ind      <- plot.meth == "all"
  grp.ind      <- all(G != 1, !is.element(method, c("FA", "IFA")))
  load.all     <- load.meth == "all"
  if(grp.ind)   {
    clust <- results$Clust
  }
  if(all.ind)   {
    if(v.sw[vars]) {
      m.sw[-1]    <- !m.sw[-1]
      if(all(vars  == "loadings", load.all)) {
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
     !v.sw[vars],
     is.element(method, c("FA", "IFA")))) {  
    if(all.ind)                       warning(paste0("Can only plot posterior mean, as ", vars, " weren't stored"), call.=F)
   v.sw[vars]     <- !v.sw[vars]
   all.ind        <- F
   m.sw["P.sw"]   <- T
  } 
  if(all(!v.sw[vars],
     !m.sw["G.sw"]))                  stop(paste0(vars, " weren't stored"))
  if(!is.logical(intervals))          stop("'intervals' must be TRUE or FALSE")
  if(!is.logical(mat))                stop("'mat' must be TRUE or FALSE")
  if(!is.logical(partial))            stop("'partial' must be TRUE or FALSE")
  if(!is.logical(titles))             stop("'titles' must be TRUE or FALSE")
  indx    <- missing(ind)
  facx    <- missing(fac)
  if(!indx)                 xind <- ind
  if(!facx) {
    if(length(fac) == 1)     fac <- rep(fac, G)
    if(length(fac) != G)              stop(paste0("'fac' must be supplied for each of the ", G, " groups"))
  }
  g.score <- all(grp.ind, !all.ind, vars == "scores")
  if(any(all(vars  == "scores", any(all.ind, !m.sw["P.sw"])), 
         m.sw["G.sw"]))  {
    Gs    <- 1
  } else if(!missing(g)) {
    if(is.element(method, c("MFA", "MIFA"))) {
      if(!is.element(g, seq_len(G)))  stop("This g value was not used during simulation")
      Gs  <- g
    } else if(g > 1)     {            message(paste0("Forced g=1 for the ", method, " method"))
      Gs  <- 1
    }
  } else if(!interactive()) {         stop("g must be supplied for non-interactive sessions")
  } else {
    Gs    <- seq_len(G)
  }
  
  for(g in Gs) {
    Q     <- Qs[g]
    msg   <- "Hit <Return> to see next plot: "
    msgx  <- all(interactive(), g != max(Gs))
    result     <- results[[g]]
    if(any(all(Q  == 0, vars == "loadings"),
           all(Qs == 0, vars == "scores")))  {            
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
    } else   {
      matx     <- mat
    }  
    if(!matx) {
      if(vars  == "scores")   {
        iter   <- store
      } else {
        iter   <- seq_len(attr(result, "Store"))
      }
    }               
    if(is.element(vars, c("scores", "loadings"))) {
      if(indx)               ind <- c(1, 1)
      if(!facx)           ind[2] <- fac[g]
      if(length(ind) > 2)             stop("Length of plotting indices can't be greater than 2")
      if(vars  == "scores")  {
        if(ind[1] >  n.obs)           stop(paste0("First index can't be greater than the number of observations - ",  n.obs))
        if(ind[2] >  Q.max)  {        warning(paste0("Second index can't be greater than ", Q.max, ", the total number of factors", if(grp.ind) paste0(" across groups"), ".\n Try specifying a vector of fac values with maximum entries ", paste0(Qs, collapse = ", "), "."), call.=F)
        ifelse(msgx, readline(msg), "")
        next
        }
      } else {
        if(ind[1] > n.var)            stop(paste0("First index can't be greater than the number of variables - ",  n.var))
        if(ind[2] > Q) {              warning(paste0("Second index can't be greater than ", Q, ", the number of factors", if(grp.ind) paste0(" in group ", g), ".\n Try specifying a vector of fac values with maximum entries ", paste0(Qs, collapse = ", "), "."), call.=F)
        ifelse(msgx, readline(msg), "")
        next
        }
      }
    } else   {
      if(indx)               ind <- 1
      if(length(ind) >  1)            stop("Length of plotting indices can't be greater than 1")
      if(ind      > n.var)            stop(paste0("Index can't be greater than the number of variables - ", n.var))
    }
    
    if(m.sw["T.sw"]) {
      if(vars  == "means")  {
        plot.x <- result$means
        if(matx) {
          matplot(t(plot.x[,]), type="l", ylab="Means", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot(x=iter, y=plot.x[ind,], type="l", ylab="Mean", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
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
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else {
          plot(x=iter, y=plot.X[ind[1],ind[2],], type="l", ylab="Scores", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
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
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else   {
          plot(x=iter, y=plot.X[ind[1],ind[2],], type="l", ylab="Loadings", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$psi
        if(matx) {
          matplot(t(plot.x[,]), type="l", ylab="Uniquenesses", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot(x=iter, y=plot.x[ind,], ylab="Uniquenesses", type="l", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(!indx) {         ind[1] <- xind[1]
        if(facx)          ind[2] <- xind[2]
      }
      if(all.ind)          xxind <- ind                            
    }
    
    if(m.sw["D.sw"]) {
      if(vars  == "means") {
        plot.X <- result$means
        if(matx) {
          plot.x  <- apply(plot.X, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot.d  <- density(plot.X[ind,])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col="black")
        }
      }
      if(vars  == "scores") {
        plot.X <- results$Scores$f
        if(by.fac) {
          plot.x  <- plot.X[ind[1],,]
        } else   {
          plot.x  <- plot.X[,ind[2],]
        }
        if(matx) {
          plot.x  <- apply(plot.x, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(by.fac) {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else   {
          plot.d  <- density(plot.X[ind[1],ind[2],])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
          polygon(plot.d, col="black")
        }
      }
      if(vars  == "loadings") {
        plot.X <- result$loadings
        if(by.fac) {
          plot.x  <- plot.X[ind[1],,]
        } else   {
          plot.x  <- plot.X[,ind[2],]
        }
        if(matx) {
          plot.x  <- apply(plot.x, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(by.fac) {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else   {
          plot.d  <- density(plot.X[ind[1],ind[2],])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
          polygon(plot.d, col="black")
        }
      }
      if(vars  == "uniquenesses") {
        plot.X <- result$psi
        if(matx) {
          plot.x  <- apply(plot.X, 1, density)
          plot.x  <- sapply(plot.x, "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot.d  <- density(plot.X[ind,])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col="black")
        }
      }
    }
    
    if(m.sw["P.sw"])  {
      if(is.element(vars, c("scores", "loadings"))) {
        if(indx)  {
         if(vars  == "scores")   {
            ind   <- c(1, min(Q.max, 2))
         } else   {
            ind   <- c(1, 1)
         } 
        }
        if(!facx) {
          ind[2]  <- fac[g]
        }
        if(vars == "scores") {
          if(any(ind[1]  > Q.max,
                 ind[2]  > Q.max))    stop(paste0("Only the first ", Q.max, " columns can be plotted"))  
        } else if(ind[2] > Q)         stop(paste0("Only the first ", Q, " columns can be plotted"))  
      }
      if(vars  == "means") {
        plot.x <- result$post.mu
        if(ci.sw[vars]) ci.x   <- result$CI.mu
        plot(plot.x, type=type, ylab="Means", xlab="Variable", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1) else if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
        if(all(intervals, ci.sw[vars])) plotCI(plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol="grey", add=T, gap=T, pch=ifelse(type == "n", NA, 16))
        if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(x=seq_along(plot.x), y=plot.x, var.names, cex=0.5)
      }
      if(vars  == "scores") {
        if(grp.ind)   {
          Labs <- clust$post.z
        } else   {
          Labs <- 1
        }
        if(!missing(Labels)) {
          if(!exists(as.character(match.call()$Labels),
              envir=.GlobalEnv)) {    warning(paste0("Object ", match.call()$Labels, " not found"), call.=F)
          } else {
            Labs  <- as.factor(Labels)
            if(length(Labs) != n.obs) stop(paste0("Labels must be a factor of length N=",  n.obs))
          }
        }
        if(g.score)  {
          if(g == 1) temp.labs  <- Labs
          z.ind  <- as.numeric(temp.labs) == g
          plot.x <- results$Scores$post.f[z.ind,]
          ind2   <- ifelse(any(!facx, Q <= 1), ind[2], if(Q > 1) max(2, ind[2]))
          if(ci.sw[vars]) ci.x  <- results$Scores$CI.f[,z.ind,]
          Labs   <- g
        } else      {
          plot.x <- results$Scores$post.f
          ind2   <- ifelse(any(!facx, Q.max <= 1), ind[2], if(Q.max > 1) max(2, ind[2]))
          if(ci.sw[vars]) ci.x  <- results$Scores$CI.f
        }
        type.f <- ifelse(any(type.x, type == "l"), "p", type)
        if(ind2 != 1)  {
          if(all(intervals, ci.sw[vars])) {
            plotCI(plot.x[,ind[1]], plot.x[,ind2], li=ci.x[1,,ind2], ui=ci.x[2,,ind2], gap=T, pch=NA, scol="grey", slty=3, xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
            plotCI(plot.x[,ind[1]], plot.x[,ind2], li=ci.x[1,,ind[1]], ui=ci.x[2,,ind[1]], add=T, gap=T, pch=NA, scol="grey", slty=3, err="x")
            if(type.f != "n") points(plot.x[,ind[1]], plot.x[,ind2], type=type.f, col=as.numeric(Labs))
          } else {
            plot(plot.x[,ind[1]], plot.x[,ind2], type=type.f, col=as.numeric(Labs),
                 xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"), ifelse(g.score, paste0(" - Group ", g), ""))))
          if(type.f == "n") text(plot.x[,ind[1]], plot.x[,ind2], obs.names, 
                               col=as.numeric(Labs), cex=0.5)
        } else   {
          if(all(intervals, ci.sw[vars])) {
            plotCI(if(!g.score) seq_len(n.obs) else seq_len(sum(z.ind)), plot.x[,ind[1]], li=ci.x[1,,ind[1]], ui=ci.x[2,,ind[1]], gap=T, pch=NA, scol="grey", slty=3, xlab="Observation", ylab=paste0("Factor ", ind[1]))
            points(plot.x[,ind[1]], type=type.f, col=as.numeric(Labs))
          } else {
            plot(plot.x[,ind[1]], type=type.f, col=as.numeric(Labs), xlab="Observation", ylab=paste0("Factor ", ind[1]))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"), ifelse(g.score, paste0(" - Group ", g), ""))))
          if(type.f == "n") text(plot.x[,ind[1]], col=as.numeric(Labs), cex=0.5)
        }
      }
      if(vars  == "loadings") {
        if(all(!all.ind, load.all)) {
          par(mai=c(1.25, 1, 0.75, 0.5), mfrow=c(1, 2), oma=c(0, 0, 1, 0))
        }
        plot.x <- result$post.load
        if(is.element(load.meth, c("all", "heatmap"))) {
          source(paste(getwd(), "/IMIFA-GIT/FullConditionals.R", sep=""), local=T)
          mcol <- mat2color(plot.x)
          plotcolors(mcol)
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all(!all.ind, !load.all), " Loadings ", " "), "Heatmap", ifelse(all(!all.ind, grp.ind, !load.all), paste0(" - Group ", g), ""))))
          axis(1, line=-0.5, tick=F, at=if(Q != 1) seq_len(Q) else 0, labels=seq_len(Q))
          if(n.var < 100) {
            axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1, at=seq_len(n.var), labels=substring(var.names[n.var:1], 1, 10))
          }
          box(lwd=2)
          mtext("Factors", side=1, line=2)
          if(Q   != 1) abline(v=seq(1, Q - 1, 1) + 0.5, lty=2, lwd=1)
        }
        if(is.element(load.meth, c("all", "raw"))) {
          if(ci.sw[vars]) ci.x   <- result$CI.load  
          if(by.fac) {
            if(ci.sw[vars]) ci.x <- ci.x[,,ind[2]]
            plot(plot.x[,ind[2]], type=type, xaxt="n", xlab="", ylab="Loading", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
            if(all(intervals, ci.sw[vars])) plotCI(plot.x[,ind[2]], li=ci.x[1,], ui=ci.x[2,], slty=3, scol="grey", add=T, gap=T, pch=ifelse(type == "n", NA, 16))
            axis(1, line=-0.5, tick=F, at=seq_len(n.var), labels=seq_len(n.var))
            mtext("Variable #", side=1, line=2)
            if(titles) title(main=list(paste0(ifelse(all(!all.ind, !load.all), paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), "")), ""), "Factor ", ind[2])))
            if(type == "n") text(x=plot.x, var.names, cex=0.5)
          } else     {
            if(ci.sw[vars]) ci.x <- ci.x[,ind[1],]
            plot(plot.x[ind[1],], type=type, xaxt="n", xlab="", ylab="Loading", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
            if(all(intervals, ci.sw[vars])) plotCI(plot.x[ind[1],], li=ci.x[1,], ui=ci.x[2,], slty=3, scol="grey", add=T, gap=T, pch=ifelse(type == "n", NA, 16))
            axis(1, line=-0.5, tick=F, at=seq_len(Q), labels=seq_len(Q))
            mtext("Factors", side=1, line=2)
            if(titles) title(main=list(paste0(ifelse(all(!all.ind, !load.all), paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), "")), ""), var.names[ind[1]], " Variable")))
            if(type == "n") text(x=plot.x[ind[1],], paste0("Factor ", seq_len(Q)), cex=0.5)
          }
        }
        if(all(!all.ind, load.all)) {
          if(titles) title(paste0("Loadings", ifelse(grp.ind, paste0(" - Group ", g), "")), outer=T)
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$post.psi
        if(ci.sw[vars]) ci.x   <- result$CI.psi
        plot(plot.x, type=type, ylab="Uniquenesses", xlab="Variable", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
        if(all(intervals, ci.sw[vars])) plotCI(plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol="grey", add=T, gap=T, pch=ifelse(type == "n", NA, 16))
        if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(seq_along(plot.x), plot.x, var.names, cex=0.5)
      }
      if(!indx) {         ind[1] <- xind[1]
        if(facx)          ind[2] <- xind[2]
      }
      if(all.ind)          ind   <- xxind
    } 
    
    if(m.sw["G.sw"]) {
      if(is.element(method, c("FA", "MFA", "MIFA"))) {
        aicm        <- round(GQ.res$AICM, 2)
        bicm        <- round(GQ.res$BICM, 2)
        if(method   != "MIFA") {
          aic.mcmc  <- round(GQ.res$AIC.mcmc, 2)
          bic.mcmc  <- round(GQ.res$BIC.mcmc, 2)
        }
      }
      if(method == "IFA")  {
        plot.Q <- GQ.res$Q.Counts
        Q.name <- names(plot.Q)
        rangeq <- as.numeric(Q.name)
        rangeq <- seq(from=min(rangeq), to=max(rangeq), by=1)
        miss   <- setdiff(rangeq, Q.name)
        miss   <- setNames(rep(0, length(miss)), as.character(miss))
        plot.Q <- c(plot.Q, miss)
        plot.Q <- plot.Q[order(as.numeric(names(plot.Q)))]
        col.Q  <- c("black", "red")[(rangeq == Q) + 1]
        Q.plot <- barplot(plot.Q, ylab="Frequency", xaxt="n", col=col.Q)
        if(titles) title(main=list("Posterior Distribution of Q"))
        axis(1, at=Q.plot, labels=Q.name, tick=F) 
        axis(1, at=median(Q.plot), labels="Q", tick=F, line=1.5) 
      }  
      if(method == "MIFA") {
        plot.Q <- GQ.res$Q.Counts
        Q.name <- lapply(plot.Q, names)
        rangeq <- as.numeric(unique(unlist(Q.name, use.names=F)))
        rangeq <- seq(from=min(rangeq), to=max(rangeq), by=1)
        miss   <- lapply(seq_len(G), function(g) setdiff(rangeq, as.numeric(Q.name[[g]])))
        miss   <- lapply(seq_len(G), function(g) setNames(rep(0, length(miss[[g]])), as.character(miss[[g]])))
        plot.Q <- lapply(seq_len(G), function(g) c(plot.Q[[g]], miss[[g]]))
        plot.Q <- do.call(rbind, lapply(seq_len(G), function(g) plot.Q[[g]][order(as.numeric(names(plot.Q[[g]])))]))
        Q.plot <- barplot(plot.Q, beside=T, ylab="Frequency", xaxt="n", col=seq_len(G + 1)[-1])
        if(titles) title(main=list(expression('Posterior Distribution of Q'["g"])))
        axis(1, at=apply(Q.plot, 2, median), labels=colnames(plot.Q), tick=F)
        axis(1, at=median(Q.plot), labels="Q", tick=F, line=1.5)
        if(titles) legend("topright", legend=paste0("Group ", seq_len(G)), bty="n", pch=15, col=seq_len(G + 1)[-1])
      }
      if(!exists("Q.plot",  envir=environment())) 
                                      message("Nothing to plot")
      class(GQ.res)    <- "listof"
      if(is.element(method, c("MFA", "MIFA"))) {
          print(GQ.res)
      } else if(method == "IFA") {
          print(tail(GQ.res, -1))
      } else {
          cat(paste0("Q = ", Q, "\n"))
      }
      if(is.element(method, c("FA", "MFA", "MIFA"))) {
        G.ind  <- ifelse(G.supp, 1, which(n.grp == G))
        Q.ind  <- ifelse(any(Q.supp, method == "MIFA"), 1, which(n.fac == Q))
        if(any(nrow(bicm) > 1, ncol(bicm) > 1)) {
          cat(paste0("AICM = ", aicm[G.ind,Q.ind], "\n"))
          cat(paste0("BICM = ", bicm[G.ind,Q.ind], "\n"))
          if(method != "MIFA") {
            cat(paste0("AIC.mcmc = ", aic.mcmc[G.ind,Q.ind], "\n"))
            cat(paste0("BIC.mcmc = ", bic.mcmc[G.ind,Q.ind], "\n"))
          }
        }
      }
    }
  
    if(m.sw["C.sw"]) {
      if(!all.ind) partial <- F
      if(!all.ind) {
        par(mai=c(1.25, 1, 0.75, 0.5), mfrow=c(1, 2), oma=c(0, 0, 2, 0))
      }
      if(vars  == "means")    {
        plot.x <- result$means 
        if(!partial) { 
          acf(plot.x[ind,], main="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind], " Variable"), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind], " Variable"), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Means - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind], " Variable")), outer=T)
        }
      }
      if(vars  == "scores")   { 
        plot.x <- results$Scores$f
        if(!partial) {
          acf(plot.x[ind[1],ind[2],], main="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", "Observation ", obs.names[ind[1]], ", Factor ", ind[2]), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind[1],ind[2],], main="", type="partial")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", "Observation ", obs.names[ind[1]], ", Factor ", ind[2]), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Scores - ", "Observation ", obs.names[ind[1]], ", Factor ", ind[2])), outer=T)
        }
      }
      if(vars  == "loadings") { 
        plot.x <- result$loadings
        if(!partial) {
          acf(plot.x[ind[1],ind[2],], main="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind[1]], " Variable, Factor ", ind[2]), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind[1],ind[2],], main="", type="partial")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind[1]], " Variable, Factor ", ind[2]), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind[1]], " Variable, Factor ", ind[2])), outer=T)
        }
      }
      if(vars  == "uniquenesses")  { 
        plot.x <- result$psi
        if(!partial) {
          acf(plot.x[ind,], main="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind], " Variable"), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind], " Variable"), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Uniquenesses - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind], " Variable")), outer=T)
        }
      }
    }
    if(all(all.ind, titles)) title(paste0(toupper(substr(vars, 1, 1)),
                             substr(vars, 2, nchar(vars)), 
                             ifelse(all(grp.ind, vars != "scores"), 
                                    paste0(" - Group ", g), "")), outer=T)
    ifelse(msgx, readline(msg), "")
  }
}