################################
### IMIFA Plotting Functions ###
################################

plot.Tuned_IMIFA    <- function(results = NULL, plot.meth = c("all", "correlation", "density", "errors", "GQ", "means", "parallel.coords", "trace", "zlabels"), 
                                vars = c("means", "scores", "loadings", "uniquenesses", "pis", "alpha"), zlabels = NULL, load.meth = c("heatmap", "raw"), palette = NULL, g = NULL, 
                                fac = NULL, by.fac = TRUE, ind = NULL, type = c("h", "n", "p", "l"), intervals = TRUE, mat = TRUE, partial = FALSE, titles = TRUE, transparency = NULL) {

  defpar  <- suppressWarnings(par(no.readonly=TRUE))
  defpar$new        <- FALSE
  if(missing(palette))   palette <- c("#999999", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if(!all(are.cols(cols=palette)))    stop("Supplied colour palette contains invalid colours")
  if(length(palette) < 3)             stop("Palette must contain 3 or more colours")
  if(missing(transparency)) {
    transparency    <- 0.75
  }
  if(transparency    < 0   || 
     transparency    > 1)             stop("'transparency' must be a single number in [0, 1]")
  tmp.pal <- palette
  palette <- adjustcolor(palette, alpha.f=transparency)
  palette(palette)
  grey    <- adjustcolor("#999999", alpha.f=0.3)
  defopt  <- options()
  options(warn=1)
  suppressWarnings(par(cex.axis=0.8, new=FALSE))
  on.exit(suppressWarnings(par(defpar)))
  on.exit(do.call("clip", as.list(defpar$usr)), add=TRUE)
  on.exit(palette("default"), add=TRUE)
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
  if(missing(results))                stop("Results must be supplied")
  if(!exists(deparse(substitute(results)),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "Tuned_IMIFA") stop(paste0("Results object of class 'Tuned_IMIFA' must be supplied"))
  GQ.res  <- results$GQ.results
  G       <- GQ.res$G
  Gseq    <- seq_len(G)
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
  m.sw         <- c(G.sw = FALSE, Z.sw = FALSE, E.sw = FALSE, P.sw = FALSE, C.sw = FALSE, D.sw = FALSE, M.sw = FALSE, T.sw = FALSE)
  v.sw         <- attr(results, "Switch")
  names(v.sw)  <- formals(sys.function(sys.parent()))$vars
  ci.sw        <- v.sw
  var.names    <- rownames(results[[1]]$post.load)
  obs.names    <- rownames(results$Scores$post.eta)
  all.ind      <- plot.meth == "all"
  grp.ind      <- !is.element(method, c("FA", "IFA"))
  if(grp.ind)   {
    clust      <- results$Clust
    labelmiss  <- !attr(clust, "Label.Sup")
  }
  grp.ind      <- all(G != 1, grp.ind)
  if(all.ind)   {
    if(v.sw[vars]) {
      m.sw[-(1:4)]  <- !m.sw[-(1:4)]
      layout(matrix(c(1, 2, 3, 4), nr=2, nc=2, byrow=TRUE))
      par(cex=0.8, mai=c(0.7, 0.7, 0.5, 0.2), mgp=c(2, 1, 0), oma=c(0, 0, 2, 0))
    }
  } else {
    sw.n  <- paste0(toupper(substring(plot.meth, 1, 1)), ".sw")
    m.sw[sw.n] <- TRUE
  }
  if(m.sw["P.sw"]) {
    if(!is.element(vars, c("means",
       "loadings", "uniquenesses")))  stop("Can only plot parallel coordinates for means, loadings or uniquenesses")
  }
  if(!grp.ind)  {
    if(m.sw["Z.sw"])                  stop("Can't use 'Z' for 'plot.meth' as no clustering has taken place")
    if(vars == "pis")                 stop("Can't plot mixing proportions as no clustering has taken place")
  }
  if(all(m.sw["E.sw"], 
         !attr(results, "Errors")))   stop("Can't plot error metrics as they were not calculated due to storage switches")
  if(all(!m.sw["G.sw"], !m.sw["Z.sw"], !m.sw["E.sw"],
     missing(vars)))                  stop("What variable would you like to plot?")
  if(all(any(m.sw["M.sw"], all.ind),
     is.element(vars, c("means", "uniquenesses")),
     !v.sw[vars],
     is.element(method, c("FA", "IFA")))) {  
    if(all.ind)                       warning(paste0("Can only plot posterior mean, as ", vars, ifelse(vars == "alpha", " wasn't", " weren't"), " stored"), call.=FALSE)
    v.sw[vars]     <- !v.sw[vars]
    all.ind        <- FALSE
    m.sw["M.sw"]   <- TRUE
  } 
  if(all(!v.sw[vars], !m.sw["G.sw"], 
     !m.sw["Z.sw"],   !m.sw["E.sw"])) stop(paste0("Nothing to plot: ", vars, ifelse(vars == "alpha", paste0(" was fixed at ", attr(results, "Alpha")), " weren't stored")))
  if(!is.logical(intervals))          stop("'intervals' must be TRUE or FALSE")
  if(!is.logical(mat))                stop("'mat' must be TRUE or FALSE")
  if(!is.logical(partial))            stop("'partial' must be TRUE or FALSE")
  if(!is.logical(titles))             stop("'titles' must be TRUE or FALSE")
  indx    <- missing(ind)
  facx    <- missing(fac)
  gx      <- missing(g)
  if(!indx)                 xind <- ind
  if(!facx) {
    flen  <- length(fac)
    if(flen == 1 && gx)      fac <- rep(fac, G)
    if(flen != G && gx)               stop(paste0("'fac' must be supplied for each of the ", G, " groups"))
  }
  g.score <- all(grp.ind, !all.ind, vars == "scores")
  if(any(all(is.element(method, c("IMIFA", "OMIFA")), m.sw["G.sw"]), m.sw["Z.sw"])) {
    Gs    <- if(gx) seq_len(2) else ifelse(g <= 2, g, 
                                      stop("Invalid 'g' value"))
  } else if(any(all(is.element(vars, c("scores", "pis", "alpha")), any(all.ind, vars != "scores", !m.sw["M.sw"])), 
            m.sw["G.sw"], all(m.sw["P.sw"], vars != "loadings"), m.sw["E.sw"])) {
    Gs    <- 1
  } else if(!gx) {
    if(!is.element(method, c("FA", "IFA"))) {
      if(!is.element(g, Gseq))        stop("This g value was not used during simulation")
      Gs  <- g
    } else if(g > 1)     {            message(paste0("Forced g=1 for the ", method, " method"))
      Gs  <- 1
    }
  } else if(!interactive())  {        stop("g must be supplied for non-interactive sessions")
  } else {
    Gs    <- Gseq
  }
  
  for(g in Gs) {
    Q     <- Qs[g]
    g.ind <- which(Gs == g)
    msg   <- "Hit <Return> to see next plot or type 'EXIT'/hit <Esc> to exit: "
    msgx  <- all(interactive(), g != max(Gs))
    ent.exit   <- function() {
      if(msgx) {
        ent  <- readline(msg)
        options(show.error.messages=FALSE)
        on.exit(suppressWarnings(options(defopt)), add=TRUE)
        if(is.element(ent, 
           c("exit", "EXIT")))        stop()
      }
    }
    result     <- results[[g]]
    if(any(all(Q  == 0, vars == "loadings"),
           all(Qs == 0, vars == "scores")))  {            
                                      warning(paste0("Can't plot ", vars, paste0(ifelse(all(vars == "loadings", G > 1), paste0(" for group ", g), "")), " as they contain no columns/factors"), call.=FALSE)
      if(g == max(Gs)) {
        break
      } else {
        ent.exit()
      }
    }
    if(any(vars   == "alpha",
           all(is.element(vars, c("means", "uniquenesses")), !indx),
           all(is.element(vars, c("scores", "loadings")),
               Q  == 1))) { matx <- FALSE
    } else   {
      matx     <- mat
    }  
    if(!matx) {
      iter     <- if(vars == "scores") store else seq_len(attr(result, "Store"))
    }               
    if(is.element(vars, c("scores", "loadings"))) {
      if(indx)               ind <- c(1, 1)
      if(!facx)           ind[2] <- fac[g]
      if(length(ind) > 2)             stop("Length of plotting indices can't be greater than 2")
      if(vars  == "scores")  {
        if(ind[1] >  n.obs)           stop(paste0("First index can't be greater than the number of observations: ",  n.obs))
        if(ind[2] >  Q.max)  {        warning(paste0("Second index can't be greater than ", Q.max, ", the total number of factors", if(grp.ind) paste0(" across groups"), ".\n Try specifying a vector of fac values with maximum entries ", paste0(Qs, collapse=", "), "."), call.=FALSE)
        ent.exit()
        next
        }
      } else {
        if(ind[1] > n.var)            stop(paste0("First index can't be greater than the number of variables: ",  n.var))
        if(ind[2] > Q) {              warning(paste0("Second index can't be greater than ", Q, ", the number of factors", if(grp.ind) paste0(" in group ", g), ".\n Try specifying a vector of fac values with maximum entries ", paste0(Qs, collapse=", "), "."), call.=FALSE)
        ent.exit()
        next
        }
      }
    } else   {
      if(any(vars == "alpha",
             indx))       ind    <- 1
      if(length(ind) >  1)            stop("Length of plotting indices can't be greater than 1")
      if(vars  == "pis")    {
        if(ind       >  G)            stop(paste0("Index can't be greater than the number of groups: ", G))
      } else {
        if(ind       > n.var)         stop(paste0("Index can't be greater than the number of variables: ", n.var))
      }
    }
    
    if(m.sw["T.sw"]) {
      if(vars  == "means")  {
        plot.x <- result$means
        if(matx) {
          matplot(t(plot.x), type="l", ylab="", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot(x=iter, y=plot.x[ind,], type="l", ylab="", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nMean - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(vars  == "scores") {
        x.plot <- results$Scores$eta
        if(by.fac) {
          plot.x  <- if(Q > 1) x.plot[ind[1],,] else t(x.plot[ind[1],,])
        } else {
          plot.x  <- x.plot[,ind[2],]
        }
        if(matx) {
          matplot(t(plot.x), type="l", ylab="", xlab="Iteration")    
          if(by.fac) {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else {
          plot(x=iter, y=x.plot[ind[1],ind[2],], type="l", ylab="", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
        }
      }
      if(vars  == "loadings") {
        x.plot <- result$loadings
        plot.x <- if(by.fac) x.plot[ind[1],,] else x.plot[,ind[2],]
        if(matx) {
          matplot(t(plot.x), type="l", ylab="", xlab="Iteration")
          if(by.fac) {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else   {
          plot(x=iter, y=x.plot[ind[1],ind[2],], type="l", ylab="", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$psi
        if(matx) {
          matplot(t(plot.x), type="l", ylab="", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot(x=iter, y=plot.x[ind,], ylab="", type="l", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nUniqueness - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(vars  == "pis") {
        plot.x <- clust$pi.prop
        if(matx) {
          matplot(t(plot.x), type="l", ylab="", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMixing Proportions")))))
        } else   {
          plot(x=iter, y=plot.x[ind,], ylab="", type="l", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMixing Proportion - Group ", ind)))))
        }
      }
      if(vars  == "alpha") {
        plot.x <- clust$DP.alpha
        plot(plot.x$alpha, ylab="", type="l", xlab="Iteration", main="")
        if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nAlpha")))))
        if(all(intervals, ci.sw[vars])) {
          ci.x <- plot.x$ci.alpha  
          abline(h=plot.x$post.alpha,  col=2,    lty=2)
          abline(h=plot.x$ci.alpha[1], col=grey, lty=2)
          abline(h=plot.x$ci.alpha[2], col=grey, lty=2)
        }
      }
      if(!indx) {         ind[1] <- xind[1]
        if(facx)          ind[2] <- xind[2]
      }
      if(all.ind)          xxind <- ind                            
    }
    
    if(m.sw["D.sw"]) {
      if(vars  == "means") {
        x.plot <- result$means
        if(matx) {
          plot.x  <- sapply(apply(x.plot, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot.d  <- density(x.plot[ind,])
          plot(plot.d, main="", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "scores") {
        x.plot <- results$Scores$eta
        if(by.fac) {
          plot.x  <- if(Q > 1) x.plot[ind[1],,] else t(x.plot[ind[1],,])
        } else   {
          plot.x  <- x.plot[,ind[2],]
        }
        if(matx) {
          plot.x  <- sapply(apply(plot.x, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="")
          if(by.fac) {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else   {
          plot.d  <- density(x.plot[ind[1],ind[2],])
          plot(plot.d, main="", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "loadings") {
        x.plot <- result$loadings
        plot.x    <- if(by.fac) x.plot[ind[1],,] else x.plot[,ind[2],]
        if(matx) {
          plot.x  <- sapply(apply(plot.x, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="")
          if(by.fac) {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else   {
          plot.d  <- density(x.plot[ind[1],ind[2],])
          plot(plot.d, main="", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "uniquenesses") {
        x.plot <- result$psi
        if(matx) {
          plot.x  <- sapply(apply(x.plot, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot.d  <- density(x.plot[ind,])
          plot(plot.d, main="", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "pis") {
        x.plot <- clust$pi.prop
        if(matx) {
          plot.x  <- sapply(apply(x.plot, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMixing Proportions")))))
        } else   {
          plot.d  <- density(x.plot[ind,])
          plot(plot.d, main="", ylab="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMixing Proportions - Group ", ind)))))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "alpha") {
        plot.x <- clust$DP.alpha
        plot.d <- density(plot.x$alpha)
        plot(plot.d, main="", ylab="")
        if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nAlpha")))))
        polygon(plot.d, col=grey)
        if(intervals) {
          avg  <- plot.x$post.alpha
          clip(avg, avg, 0, plot.d$y[which.min(abs(plot.d$x - avg))])
          abline(v=avg, col=2, lty=2)
        }
      }
    }
    
    if(m.sw["M.sw"])  {
      if(is.element(vars, c("scores", "loadings"))) {
        if(indx)  {
          ind     <- if(vars == "scores") c(1, min(Q.max, 2)) else c(1, 1)
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
        if(ci.sw[vars])   ci.x   <- result$ci.mu
        plot(plot.x, type=type, ylab="", xlab="Variable", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1) else if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
        if(all(intervals, ci.sw[vars])) plotCI(plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
        if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(x=seq_along(plot.x), y=plot.x, var.names, cex=0.5)
      }
      if(vars  == "scores") {
        labs   <- if(grp.ind) clust$map else 1
        if(!missing(zlabels)) {
          if(!exists(as.character(match.call()$zlabels),
              envir=.GlobalEnv)) {    warning(paste0("Object ", match.call()$zlabels, " not found"), call.=FALSE)
          } else {
            labs  <- as.numeric(as.factor(zlabels))
            if(length(labs) != n.obs) stop(paste0("'zlabels' must be a factor of length N=",  n.obs))
          }
        }
        if(g.score)  { 
          if(g.ind == 1)  tmplab <- labs
          z.ind  <- as.numeric(levels(tmplab))[tmplab] == g
          plot.x <- results$Scores$post.eta[z.ind,,drop=FALSE]
          ind2   <- ifelse(any(!facx, Q <= 1), ind[2], if(Q > 1) max(2, ind[2]))
          if(ci.sw[vars])  ci.x  <- results$Scores$ci.eta[,z.ind,, drop=FALSE]
          labs   <- g
        } else       {
          plot.x <- results$Scores$post.eta
          ind2   <- ifelse(any(!facx, Q.max <= 1), ind[2], if(Q.max > 1) max(2, ind[2]))
          if(ci.sw[vars])  ci.x  <- results$Scores$ci.eta
        }
        col.s  <- if(is.factor(labs)) as.numeric(levels(labs))[labs] + 1 else labs + 1
        type.s <- ifelse(any(type.x, type == "l"), "p", type)
        if(ind2 != 1)  {
          if(all(intervals, ci.sw[vars])) {
            plotCI(plot.x[,ind[1]], plot.x[,ind2], li=ci.x[1,,ind2], ui=ci.x[2,,ind2], gap=TRUE, pch=NA, scol=grey, slty=3, xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
            plotCI(plot.x[,ind[1]], plot.x[,ind2], li=ci.x[1,,ind[1]], ui=ci.x[2,,ind[1]], add=TRUE, gap=TRUE, pch=NA, scol=grey, slty=3, err="x")
            if(type.s != "n") points(plot.x[,ind[1]], plot.x[,ind2], type=type.s, col=col.s, pch=20)
          } else {
            plot(plot.x[,ind[1]], plot.x[,ind2], type=type.s, col=col.s, pch=20,
                 xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"), ifelse(g.score, paste0(" - Group ", g), ""))))
          if(type.s == "n") text(plot.x[,ind[1]], plot.x[,ind2], obs.names, col=col.s, cex=0.5)
        } else   {
          if(all(intervals, ci.sw[vars])) {
            plotCI(if(!g.score) seq_len(n.obs) else seq_len(sum(z.ind)), plot.x[,ind[1]], li=ci.x[1,,ind[1]], ui=ci.x[2,,ind[1]], gap=TRUE, pch=NA, scol=grey, slty=3, xlab="Observation", ylab=paste0("Factor ", ind[1]))
            points(plot.x[,ind[1]], type=type.s, col=col.s, pch=20)
          } else {
            plot(plot.x[,ind[1]], type=type.s, col=col.s, xlab="Observation", ylab=paste0("Factor ", ind[1]), pch=20)
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"), ifelse(g.score, paste0(" - Group ", g), ""))))
          if(type.s == "n") text(plot.x[,ind[1]], col=col.s, cex=0.5)
        }
      }
      if(vars  == "loadings") {
        plot.x <- result$post.load
        if(load.meth == "heatmap") {
          if(Q  > 1) {
            plotcolors(mat2cols(plot.x))
          } else {
            image(z=t(plot.x[seq(n.var, 1),seq_len(Q)]), xlab="", ylab="", xaxt="n", yaxt="n", col=dichromat(heat.colors(30)))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(!all.ind, " Loadings ", " "), "Heatmap", ifelse(all(!all.ind, grp.ind), paste0(" - Group ", g), ""))))
          axis(1, line=-0.5, tick=FALSE, at=if(Q != 1) seq_len(Q) else 0, labels=seq_len(Q))
          if(n.var < 100) {
            axis(2, cex.axis=0.5, line=-0.5, tick=FALSE, las=1, at=if(Q > 1) seq_len(n.var) else seq(from=0, to=1, by=1/(n.var - 1)), labels=substring(var.names[n.var:1], 1, 10))
          }
          box(lwd=2)
          mtext(ifelse(Q > 1, "Factors", "Factor"), side=1, line=2)
          if(Q != 1) abline(v=seq(1, Q - 1, 1) + 0.5, lty=2, lwd=1)
        } else {
          if(ci.sw[vars])  ci.x  <- result$ci.load  
          if(by.fac) {
            if(ci.sw[vars]) ci.x <- ci.x[,,ind[2]]
            plot(plot.x[,ind[2]], type=type, xaxt="n", xlab="", ylab="Loading", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
            if(all(intervals, ci.sw[vars])) plotCI(plot.x[,ind[2]], li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
            axis(1, line=-0.5, tick=FALSE, at=seq_len(n.var), labels=seq_len(n.var))
            mtext("Variable #", side=1, line=2, cex=0.8)
            if(titles) title(main=list(paste0(ifelse(!all.ind, paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), "")), ""), "Factor ", ind[2])))
            if(type == "n") text(x=plot.x, var.names, cex=0.5)
          } else     {
            if(ci.sw[vars]) ci.x <- ci.x[,ind[1],]
            plot(plot.x[ind[1],], type=type, xaxt="n", xlab="", ylab="Loading", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
            if(all(intervals, ci.sw[vars])) plotCI(plot.x[ind[1],], li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
            axis(1, line=-0.5, tick=FALSE, at=seq_len(Q), labels=seq_len(Q))
            mtext("Factors", side=1, line=2)
            if(titles) title(main=list(paste0(ifelse(!all.ind, paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), "")), ""), var.names[ind[1]], " Variable")))
            if(type == "n") text(x=plot.x[ind[1],], paste0("Factor ", seq_len(Q)), cex=0.5)
          }  
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$post.psi
        if(ci.sw[vars])   ci.x   <- result$ci.psi
        plot(plot.x, type=type, ylab="", xlab="Variable", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
        if(all(intervals, ci.sw[vars])) plotCI(plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
        if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(seq_along(plot.x), plot.x, var.names, cex=0.5)
      }
      if(vars  == "pis") {
        plot.x <- clust$post.pi
        if(ci.sw[vars])   ci.x   <- clust$ci.pi
        if(matx) {
          if(all(intervals, ci.sw[vars])) {
            plotCI(barplot(plot.x, ylab="", xlab="", col=grey, ylim=c(0, 1), cex.names=0.7),
                   plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol=2, add=TRUE, gap=TRUE, pch=20)
          } else {
            barplot(plot.x, ylab="", xlab="", ylim=c(0, 1), cex.names=0.7)
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMixing Proportions")))))  
        } else {
          if(all(intervals, ci.sw[vars])) {
            plotCI(barplot(plot.x[ind], ylab="", xlab="", ylim=c(0, 1), cex.names=0.7),
                   plot.x[ind], li=ci.x[1,ind], ui=ci.x[2,ind], slty=3, scol=2, add=TRUE, gap=TRUE, pch=20)
          } else {
            barplot(plot.x[ind], ylab="", xlab="Variable", ylim=c(0, 1), cex.names=0.7)
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMixing Proportions - Group ", ind)))))
        }
      }
      if(vars  == "alpha") {
        plot(c(0, 1), c(0, 1), ann=FALSE, bty='n', type='n', xaxt='n', yaxt='n')
        if(titles) title(main=list(paste0("Summary Statistics", ifelse(all.ind, "", ":\nAlpha"))))
        plot.x <- clust$DP.alpha[-1]
        a.step <- attr(results, "Alpha.step")
        conf   <- attr(results, "Conf.Level")
        digits <- options()$digits
        a.adj  <- rep(0.5, 2)
        a.cex  <- par()$fin[2]/ifelse(a.step != "metropolis", 4, 5)
        pen    <- ifelse(a.step != "metropolis", 0.125, 0)
        text(x=0.5, y=0.85 - pen, cex=a.cex, col="black", adj=a.adj, expression(bold("Posterior Mean:\n")))
        text(x=0.5, y=0.85 - pen, cex=a.cex, col="black", adj=a.adj, bquote(.(round(plot.x$post.alpha, digits))))
        text(x=0.5, y=0.57 - pen, cex=a.cex, col="black", adj=a.adj, expression(bold("\nVariance:\n")))
        text(x=0.5, y=0.57 - pen, cex=a.cex, col="black", adj=a.adj, bquote(.(round(plot.x$var.alpha, digits))))
        text(x=0.5, y=0.4  - pen, cex=a.cex, col="black", adj=a.adj, bquote(bold(.(100 * conf))~bold("% Confidence Interval:")))
        text(x=0.5, y=0.28 - pen, cex=a.cex, col="black", adj=a.adj, bquote(paste("[", .(round(plot.x$ci.alpha[1], digits)), ", ", .(round(plot.x$ci.alpha[2], digits)), "]")))
        if(a.step == "metropolis") {
          text(x=0.5, y=0.17, cex=a.cex, col="black", adj=a.adj, expression(bold("Acceptance Rate:")))
          text(x=0.5, y=0.1,  cex=a.cex, col="black", adj=a.adj, bquote(paste(.(round(100 * plot.x$acceptance.rate, 2)), "%")))
        }
      }
      if(!indx) {         ind[1] <- xind[1]
        if(facx)          ind[2] <- xind[2]
      }
      if(all.ind)          ind   <- xxind
    } 
    
    if(m.sw["G.sw"]) {
      plotG.ind  <- is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))
      plotQ.ind  <- any(any(g > 1, is.element(method, c("IFA", "MIFA"))), all(is.element(method, c("IMIFA", "OMIFA")), g != 1))
      aicm       <- round(GQ.res$AICM, 2)
      bicm       <- round(GQ.res$BICM, 2)
      log.iLLH   <- round(GQ.res$LogIntegratedLikelihoods, 2)
      if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
        aic.mcmc <- round(GQ.res$AIC.mcmc, 2)
        bic.mcmc <- round(GQ.res$BIC.mcmc, 2)
      }
      if(all(plotG.ind, g == 1))  {
        layout(1)
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        plot.G <- GQ.res$G.Counts
        G.name <- names(plot.G)
        rangeG <- as.numeric(G.name)
        rangeG <- seq(from=min(rangeG), to=max(rangeG), by=1)
        missG  <- setdiff(rangeG, G.name)
        missG  <- setNames(rep(0, length(missG)), as.character(missG))
        plot.G <- c(plot.G, missG)
        plot.G <- plot.G[order(as.numeric(names(plot.G)))]
        col.G  <- c(1, 2)[(rangeG == G) + 1]
        G.plot <- barplot(plot.G, ylab="Frequency", xaxt="n", col=col.G)
        if(titles) title(main=list("Posterior Distribution of G"))
        axis(1, at=G.plot, labels=names(plot.G), tick=FALSE) 
        axis(1, at=median(G.plot), labels="G", tick=FALSE, line=1.5) 
      }
      if(all(method == "IFA", plotQ.ind)) {
        layout(1)
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        plot.Q <- GQ.res$Q.Counts
        Q.name <- names(plot.Q)
        rangeQ <- as.numeric(Q.name)
        rangeQ <- seq(from=min(rangeQ), to=max(rangeQ), by=1)
        missQ  <- setdiff(rangeQ, Q.name)
        missQ  <- setNames(rep(0, length(missQ)), as.character(missQ))
        plot.Q <- c(plot.Q, missQ)
        plot.Q <- plot.Q[order(as.numeric(names(plot.Q)))]
        col.Q  <- c(1, 2)[(rangeQ == Q) + 1]
        Q.plot <- barplot(plot.Q, ylab="Frequency", xaxt="n", col=col.Q)
        if(titles) title(main=list("Posterior Distribution of Q"))
        axis(1, at=Q.plot, labels=Q.name, tick=FALSE) 
        axis(1, at=median(Q.plot), labels="Q", tick=FALSE, line=1.5) 
      }  
      if(all(method != "IFA", plotQ.ind)) {
        plot.Q <- GQ.res$Q.Counts
        plot.Q <- if(is.list(plot.Q)) plot.Q else list(plot.Q)
        Q.name <- lapply(plot.Q, names)
        rangeQ <- as.numeric(unique(unlist(Q.name, use.names=FALSE)))
        rangeQ <- seq(from=min(rangeQ), to=max(rangeQ), by=1)
        missQ  <- lapply(Gseq, function(g) setdiff(rangeQ, as.numeric(Q.name[[g]])))
        missQ  <- lapply(Gseq, function(g) setNames(rep(0, length(missQ[[g]])), as.character(missQ[[g]])))
        plot.Q <- lapply(Gseq, function(g) c(plot.Q[[g]], missQ[[g]]))
        plot.Q <- do.call(rbind, lapply(Gseq, function(g) plot.Q[[g]][order(as.numeric(names(plot.Q[[g]])))]))
        if(titles) {
          layout(rbind(1, 2), heights=c(9, 1))
          par(mar=c(3.1, 4.1, 4.1, 2.1))
        }
        Q.plot <- barplot(plot.Q, beside=TRUE, ylab="Frequency", xaxt="n", col=Gseq, space=c(0, 2))
        if(titles) title(main=list(expression('Posterior Distribution of Q'["g"])))
        axis(1, at=apply(Q.plot, 2, median), labels=colnames(plot.Q), tick=FALSE)
        axis(1, at=median(Q.plot), labels="Q", tick=FALSE, line=1)
        if(titles) {
          par(mar=c(0, 0, 0, 0))
          plot.new()
          tmp  <- if(G > 5) unlist(lapply(Gseq, function(g) c(Gseq[g], Gseq[g + ceiling(G/2)])))[Gseq] else Gseq
          ltxt <- paste0("Group ", tmp)
          lcol <- Gseq[tmp]
          legend("center", legend=ltxt, ncol=if(G > 5) ceiling(G/2) else G, bty="n", pch=15, col=lcol, cex=max(0.7, 1 - 0.03 * G))
        }
      }
      if(!any(plotQ.ind, plotG.ind))  message("Nothing to plot")      
      gq.nam <- substring(names(GQ.res), 1, 1)
      if(is.element(method, c("IMIFA", "OMIFA"))) {
        if(g == 1) {
          print(GQ.res[gq.nam == "G"])
        } else {
          print(GQ.res[gq.nam == "Q"])
        }
        if(g == max(Gs)) {
          print(GQ.res[gq.nam != "G" & gq.nam != "Q" & gq.nam != "S"])
        }
      } else if(is.element(method, c("MFA", "MIFA", "OMFA", "IMFA"))) {
          print(GQ.res[gq.nam != "S"])
      } else if(method == "IFA") {
          print(tail(GQ.res[gq.nam != "S"], -1))
      } else   {
          cat(paste0("Q = ", Q, "\n"))
      }
      if(any(dim(bicm) > 1)) {
        G.ind  <- ifelse(any(G.supp, !is.element(method, c("MFA", "MIFA"))), 1, which(n.grp == G))
        Q.ind  <- ifelse(any(Q.supp, !is.element(method, c("FA", "MFA"))),   1, which(n.fac == Q))
          cat(paste0("AICM = ", aicm[G.ind,Q.ind], "\n"))
          cat(paste0("BICM = ", bicm[G.ind,Q.ind], "\n"))
        if(!is.element(method, c("IFA", "MIFA"))) {
          cat(paste0("AIC.mcmc = ", aic.mcmc[G.ind,Q.ind], "\n"))
          cat(paste0("BIC.mcmc = ", bic.mcmc[G.ind,Q.ind], "\n"))
        }
        cat(paste0("Log Integrated Likelihood = ", log.iLLH[G.ind,Q.ind], "\n"))
      }
      if(all(plotG.ind, 
             attr(GQ.res, "G.big")))  warning("G had to be prevented from exceeding the maximum allowable number of groups.\n Consider re-running the model with a higher value for 'trunc.G'", call.=FALSE)
      if(all(plotQ.ind, 
             attr(GQ.res, "Q.big")))  warning("Q had to be prevented from exceeding its initial value.\n Consider re-running the model with a higher value for 'range.Q'", call.=FALSE)
    }
    
    if(m.sw["Z.sw"]) {
      if(type == "l")                 stop("'type' cannot be 'l' for clustering uncertainty plots")
      plot.x <- clust$uncertainty
      if(g == 1) {
        col.x  <- c(1, 4)[(plot.x >= 1/G) + 1]
        if(type != "h") col.x[plot.x == 0] <- NA
        if(titles) {
          layout(rbind(1, 2), heights=c(1, 6))
          par(mar=c(0, 4.1, 0.5, 2.1))
          plot.new()
          legend("center", legend=bquote(1/G == 1/.(G)), title="", lty=2, col=2, bty="n", y.intersp=par()$fin[2] * 7/5)
          legend("center", legend=c(" "," "), title=expression(bold("Clustering Uncertainty")), bty='n', y.intersp=par()$fin[2] * 2/5, cex=par()$cex.main)
          par(mar=c(5.1, 4.1, 0.5, 2.1))
        }
        plot(plot.x, type=type, ylim=c(0, 1 - 1/G), col=col.x, axes=FALSE, ylab="Uncertainty", xlab="Observation", pch=ifelse(type == "n", NA, 16))
        rect(0, 0, n.obs, 1 - 1/G) 
        if(G == 2) {
          abline(h=0.5, col=par()$bg)
          abline(v=0,   col=par()$bg)
        }
        lines(x=c(0, n.obs), y=c(1/G, 1/G), lty=2, col=2)  
        axis(1, las=1, pos=0, cex.axis=0.9)
        axis(2, at=c(seq(from=0, to=min(1 - 1/G - 1/1000, 0.8), by=0.1), 1 - 1/G), labels=c(seq(from=0, to=min(1 - 1/G - 1/1000, 0.8), by=0.1), "1 - 1/G"), las=2, pos=0, cex.axis=0.9)
        if(type == "n")  {
          znam  <- obs.names
          znam[plot.x == 0] <- ""
          text(x=seq_along(plot.x), y=plot.x, znam, col=col.x, cex=0.5)
        }
      } else {
        if(titles) {
          layout(rbind(1, 2), heights=c(1, 6))
          par(mar=c(0, 4.1, 0.5, 2.1))
          plot.new()
          legend("center", legend=bquote({NA >= 1/G} == 1/.(G)), title="", pch=15, col=4, bty="n", y.intersp=par()$fin[2] * 7/5)
          legend("center", legend=c(" "," "), title=expression(bold("Clustering Uncertainty")), bty='n', y.intersp=par()$fin[2] * 2/5, cex=par()$cex.main)
          par(mar=c(5.1, 4.1, 0.5, 2.1))
        }
        x.plot  <- hist(plot.x, plot=FALSE)
        breaks  <- x.plot$breaks
        cols    <- 3 + (x.plot$breaks >= 1/G)
        cols[cols == 3] <- grey
        plot(x.plot, main="", xlab="Uncertainties", xlim=c(0, 1 - 1/G), col=cols, xaxt="n")
        axis(1, at=c(breaks[round(breaks, 1) < min(0.8, 1 - 1/G)], 1 - 1/G), labels=c(breaks[round(breaks, 1) < min(0.8, 1 - 1/G)], "1 - 1/G"), las=2, pos=0, cex.axis=0.9)
      }
      if(g == min(Gs)) {
        if(any(!labelmiss,  !missing(zlabels))) {
          if(all(!labelmiss, missing(zlabels))) {
           prf  <- clust$perf
          } else   {
           labs <- as.factor(zlabels)
           if(length(labs) != n.obs)  stop(paste0("'zlabels' must be a factor of length N=",  n.obs))
           pzs  <- clust$map
           if(nlevels(pzs) == nlevels(labs)) {
            lsw <- lab.switch(z.new=pzs, z.old=labs, Gs=seq_len(G))
            pzs <- factor(lsw$z)
           }
           tab  <- table(pzs, labs, dnn=list("Predicted", "Observed"))
           prf  <- c(classAgreement(tab), classError(pzs, labs))
           if(nrow(tab) != ncol(tab))   {
            prf <- prf[-seq_len(2)]
            names(prf)[4]        <- "error.rate"
           } else {
            names(prf)[6]        <- "error.rate"
           }
           if(prf$error.rate     == 0) {
            prf$misclassified    <- NULL
           }
           prf  <- c(list(confusion.matrix = tab), prf)
           if(nlevels(pzs)  == nlevels(labs)) {
            names(prf)[1]  <- "matched.confusion.matrix"
           }
           class(prf)      <- "listof"
          }
          ucert <- attr(plot.x, "Obs")
          if(!is.null(ucert)) {
           prf  <- c(prf, list(uncertain = ucert))  
          }
          prf$error.rate   <- paste0(round(100 * prf$error.rate, 2), "%")
          print(prf)
        } else                        message("Nothing to print: try supplying known cluster labels")
      }
    }
    
    if(m.sw["P.sw"]) {
      plot.x <- if(vars == "means") results$Means$post.mu else if(vars == "uniquenesses") results$Uniquenesses$post.psi else results$Loadings$post.load[[g]]
      x.plot <- apply(plot.x, 1L, range, na.rm=TRUE)
      plot.x <- apply(plot.x, 2L, function(x) (x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))
      varnam <- paste0(toupper(substr(vars, 1, 1)), substr(vars, 2, nchar(vars)))
      if(grp.ind) {
        layout(rbind(1, 2), heights=c(9, 1))
        par(mar=c(3.1, 4.1, 4.1, 2.1))
      }
      matplot(seq_len(n.var), plot.x, type="p", col=if(vars == "loadings") seq_len(Q) + 1 else seq_len(G) + 1, pch=15, xlab="Variable", ylab=paste0("Standardised ", varnam), axes=FALSE, main=paste0("Parallel Coordinates: ", varnam, ifelse(all(grp.ind, vars == "loadings"), paste0("\n Group ", g), "")))
      axis(1, at=seq_len(n.var), labels=if(titles && n.var < 100) rownames(plot.x) else rep("", n.var), cex.axis=0.5, tick=FALSE)
      for(i in seq_len(n.var))    {
        lines(c(i, i), c(0, 1), col=grey)
        if(titles && n.var < 100) { 
          text(c(i, i), c(0, 1), labels=format(x.plot[,i], digits=3), xpd=NA, offset=0.3, pos=c(1, 3), cex=0.5)
        }
      }
      if(grp.ind) {
        par(mar=c(0, 0, 0, 0))
        plot.new()
        Xp   <- ifelse(vars == "loadings", Q, G)
        Xseq <- seq_len(Xp)
        tmp  <- if(Xp > 5) unlist(lapply(Xseq, function(x) c(Xseq[x], Xseq[x + ceiling(Xp/2)])))[Xseq] else Xseq
        ltxt <- paste0(ifelse(vars == "loadings", "Factor ", "Group "), tmp)
        lcol <- Xseq[tmp]
        legend("center", pch=15, col=lcol + 1, legend=ltxt, ncol=if(Xp > 5) ceiling(Xp/2) else Xp, bty="n", cex=max(0.7, 1 - 0.03 * Xp))
      }
    } 
    
    if(m.sw["E.sw"]) {
      palette(tmp.pal)
      x.plot <- results$Error
      plot.x      <- if(G > 1) cbind(do.call(rbind, x.plot[-length(x.plot)]), Averages = x.plot$Averages) else x.plot
      if(titles) {
        layout(rbind(1, 2), heights=c(9, 1))
        par(mar=c(3.1, 4.1, 4.1, 2.1))
      }
      col.e  <- if(G > 1) seq_len(nrow(plot.x)) else seq_along(plot.x)
      if(G > 1)  {
        dens <- matrix(-1, nr=5, nc=G + 1)
        dens[,G + 1]       <- 30
      } else {
        dens <- NULL
      }
      pl.x   <- barplot(plot.x, beside=TRUE, col=col.e, main="", ylab="Deviation", density=dens)
      na.x   <- if(G > 1) is.na(results$Error[[1]]) else FALSE
      if(G > 1) points(x=apply(as.matrix(pl.x[,which(na.x)]), 2, median), y=rep(0, sum(na.x)), pch=16, col=6)
      if(titles) title(main=list("Error Metrics"))
      if(titles) {
        par(mar=c(0, 0, 0, 0))
        plot.new()
        ltxt <- c("MSE", "RMSE", "NRMSE", "CVRMSE", "MAD")
        lnc  <- length(col.e)
        lcol <- col.e
        xna  <- sum(na.x) > 0
        lpch <- rep(15, 5)
        temp <- legend("center", legend=if(xna) c(ltxt, "Missing") else ltxt, ncol=ifelse(xna, lnc + 1, lnc), bty="n",
                       pch=if(xna) c(lpch, max(lpch) + 1) else lpch, col=if(xna) c(lcol, length(lcol) + 1) else lcol, cex=0.8)
        if(xna) text(x=temp$text$x[6] - 0.025, y=temp$text$y[6] + 0.015, "__")
      }  
      if(G > 1) {
        avg  <- setNames(list(x.plot$Averages), "Average Error Metrics") 
        class(avg)         <- "listof"
      } else {
        avg  <- x.plot
      }
      print(avg)  
    }
  
    if(m.sw["C.sw"]) {
      if(!all.ind)   {
       partial <- FALSE
       par(mai=c(1.25, 1, 0.75, 0.5), mfrow=c(1, 2), oma=c(0, 0, 2, 0))
      }
      if(vars  == "means")    {
        plot.x <- result$means 
        if(!partial) { 
          acf(plot.x[ind,], main="", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind], " Variable"), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind], " Variable"), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Means - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind], " Variable")), outer=TRUE)
        }
      }
      if(vars  == "scores")   { 
        plot.x <- results$Scores$eta
        if(!partial) {
          acf(plot.x[ind[1],ind[2],], main="", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", "Observation ", obs.names[ind[1]], ", Factor ", ind[2]), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind[1],ind[2],], main="", type="partial", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", "Observation ", obs.names[ind[1]], ", Factor ", ind[2]), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Scores - ", "Observation ", obs.names[ind[1]], ", Factor ", ind[2])), outer=TRUE)
        }
      }
      if(vars  == "loadings") { 
        plot.x <- result$loadings
        if(!partial) {
          acf(plot.x[ind[1],ind[2],], main="", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind[1]], " Variable, Factor ", ind[2]), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind[1],ind[2],], main="", type="partial", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind[1]], " Variable, Factor ", ind[2]), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind[1]], " Variable, Factor ", ind[2])), outer=TRUE)
        }
      }
      if(vars  == "uniquenesses")  { 
        plot.x <- result$psi
        if(!partial) {
          acf(plot.x[ind,], main="", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind], " Variable"), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind], " Variable"), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Uniquenesses - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind], " Variable")), outer=TRUE)
        }
      }
      if(vars  == "pis")  { 
        plot.x <- clust$pi.prop
        if(!partial) {
          acf(plot.x[ind,], main="", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("ACF", ifelse(all(all.ind, matx), paste0(" - Group ", ind), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("PACF", ifelse(all(all.ind, matx), paste0(" - Group ", ind), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Mixing Proportions - Group ", ind)), outer=TRUE)
        }
      }
      if(vars  == "alpha") {          
        plot.x <- clust$DP.alpha$alpha
        if(clust$DP.alpha$acceptance.rate == 0)  {
                                      warning(paste0("0% acceptance rate: can't plot ", ifelse(all.ind, ifelse(partial, "partial-", "auto-"), ""), "correlation function", ifelse(all.ind, "", "s")), call.=FALSE)
          next
        }
        if(!partial) {
          acf(plot.x, main="", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("ACF")))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x, main="", type="partial", ci.col=4, ylab="")
          if(titles) title(main=list(paste0("PACF")))
          if(all(!all.ind, titles)) title(main=list(paste0("Alpha")), outer=TRUE)
        }
      }
    }
    if(all(all.ind, titles)) title(ifelse(vars != "pis", paste0(toupper(substr(vars, 1, 1)), substr(vars, 2, nchar(vars)), 
                             ifelse(all(grp.ind, !is.element(vars, c("scores", "pis", "alpha"))), paste0(" - Group ", g), "")), 
                             paste0("Mixing Proportions", ifelse(matx, "", paste0(" - Group ", ind)))), outer=TRUE)
    ent.exit()
  }
}

# Loadings Heatmaps
  mat2cols     <- function(m, cols = dichromat(heat.colors(30)), 
                           byrank = FALSE, breaks = length(cols)) { 
    m1         <- if(isTRUE(byrank)) rank(m) else m
    facs       <- cut(m1, breaks, include.lowest=TRUE)
    answer     <- cols[as.numeric(facs)]
    if(is.matrix(m)) {
      answer   <- matrix(answer, nrow(m), ncol(m))
      rownames(answer)  <- rownames(m)
      colnames(answer)  <- colnames(m)
    }
    answer    
  }

# Colour Checker
  are.cols     <- function(cols) {
    vapply(cols, function(x) { tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) }, logical(1))
  }

# Prior No. Groups (DP & PY)
  G.prior      <- function(N, alpha, discount = 0, plot = TRUE, 
                           avg = FALSE, col = "black", ...) {
    if(length(setdiff("Rmpfr", rownames(installed.packages()))) > 0) {
                                      stop("'Rmpfr' package not installed")
    }
    if(length(setdiff("Rmpfr", (.packages()))) > 0) {
      suppressMessages(library(Rmpfr))
      on.exit(detach.pkg(Rmpfr))
      on.exit(detach.pkg(gmp), add=TRUE)  
    }
    on.exit(do.call("clip", as.list(par("usr"))), add=TRUE)
    if(any(c(length(N), length(plot),
             length(avg)) > 1))       stop("Arguments 'N', 'plot', 'add', and 'avg' must be strictly of length 1")
    max.len    <- max(length(alpha),  length(discount))
    if(!is.element(length(alpha),
       c(1, max.len)))                stop("'alpha' must be of length 1 or length(discount)")
    if(!is.element(length(discount),
       c(1, max.len)))                stop("'discount' must be of length 1 or length(alpha)")
    if(!all(is.numeric(discount), is.numeric(alpha), 
            is.numeric(N)))           stop("All inputs must be numeric")
    if(any(discount < 0 |
       discount >= 1))                stop("'discount' must lie in the interval [0,1)")
    if(any(alpha < -discount))        stop("'alpha' must be strictly greater than -discount")
    if(missing(col))    {
      col      <- c("#999999", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    }   
    if(!all(are.cols(col)))           stop("Supplied 'col' is not a valid colour")
    if(!is.logical(avg))              stop("'avg' must be TRUE or FALSE")
    if(!is.logical(plot))             stop("'plot' must be TRUE or FALSE")
    if(length(alpha)    != max.len) {
      alpha    <- rep(alpha,    max.len)
    }
    if(length(discount) != max.len) {
      discount <- rep(discount, max.len)
    }
    rx         <- matrix(0, nr=N + 1, nc=max.len)
    for(i in seq_len(max.len))      {
      tmp      <- if(discount[i] == 0) alpha[i]^seq_len(N)/pochMpfr(alpha[i], N) else {
                  exp(unlist(vapply(seq_len(N), function(Gs=seq_len(k), x=0) { for(g in Gs) {
                    x <- x + log(alpha[i] + g * discount[i]) }; x}, numeric(1))) - 
                           log(pochMpfr(alpha[i] + 1, N - 1))) / discount[i]^seq_len(N) }
      if(discount[i]  == 0) {
        rx[,i] <- c(0, asNumeric(abs(Stirling1.all(N) * tmp)))
      } else                          stop("Plotting with non-zero discount not yet implemented\nTry supplying the same arguments to G.expected() or G.variance()")
    }
    rx         <- scale(rx, center=FALSE, scale=colSums(rx))
    max.rx     <- apply(rx, 2, max)
    if(plot)   {
      matplot(x=c(0, seq_len(N)), y=rx, type="l", col=col, lty=1, xlab="Groups", 
              ylab="Density", main=paste0("Prior Distribution of G\nN=", N), ...)
    }
    if(avg)    {
      exp.g    <- G.expected(N, alpha, discount)
      cat("\t");  cat(paste("E(G) = ", round(exp.g, options()$digits), "\n"))  
      if(plot) {
        col    <- rep(col, max.len)
        for(i in seq_len(max.len))  {
          clip(exp.g[i], exp.g[i], 0, max.rx[i])
          abline(v=exp.g[i],  lty=2, col=col[i])  
        }
      }
    }
      invisible(rx)
  }
#