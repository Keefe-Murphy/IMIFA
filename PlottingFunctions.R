################################
### IMIFA Plotting Functions ###
################################

plot.IMIFA     <- function(results = NULL, plot.meth = c("all", "correlation", "density", "errors", "posterior", "GQ", "trace", "Z"), 
                           vars = c("means", "scores", "loadings", "uniquenesses", "pis", "alpha"), labels = NULL, load.meth = c("all", "heatmap", "raw"), palette = NULL, g = NULL,
                           fac = NULL, by.fac = TRUE, ind = NULL, type = c("h", "n", "p", "l"), intervals = TRUE, mat = TRUE, partial = FALSE, titles = TRUE, transparency = NULL) {

  source(paste(getwd(), "/IMIFA-GIT/FullConditionals.R", sep=""), local=TRUE)
  defpar  <- suppressWarnings(par(no.readonly=TRUE))
  defpar$new   <- FALSE
  if(missing(palette))   palette <- c("#999999", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if(!all(are.colours(cols=palette))) stop("Supplied colour palette contains invalid colours")
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
  on.exit(palette("default"), add=TRUE)
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
  if(missing(results))                stop("Results must be supplied")
  if(!exists(deparse(substitute(results)),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$results, " not found"))
  if(class(results) != "IMIFA")       stop(paste0("Results object of class 'IMIFA' must be supplied"))
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
  m.sw         <- c(G.sw = FALSE, Z.sw = FALSE, E.sw = FALSE, C.sw = FALSE, D.sw = FALSE, P.sw = FALSE, T.sw = FALSE)
  v.sw         <- attr(results, "Switch")
  names(v.sw)  <- formals(sys.function(sys.parent()))$vars
  ci.sw        <- v.sw
  var.names    <- rownames(results[[1]]$post.load)
  obs.names    <- rownames(results$Scores$post.f)
  all.ind      <- plot.meth == "all"
  grp.ind      <- all(G != 1, !is.element(method, c("FA", "IFA")))
  load.all     <- all(load.meth == "all", vars == "loadings")
  if(grp.ind)   {
    clust      <- results$Clust
    labelmiss  <- !attr(clust, "Label.Sup")
  }
  if(all.ind)   {
    if(v.sw[vars]) {
      m.sw[-(1:3)]  <- !m.sw[-(1:3)]
      if(all(vars   == "loadings", load.all)) {
        layout(matrix(c(1, 2, 3, 4, 3, 5), nr=3, nc=2, byrow=TRUE))
      } else {
        layout(matrix(c(1, 2, 3, 4), nr=2, nc=2, byrow=TRUE))
      }
      par(cex=0.8, mai=c(0.7, 0.7, 0.5, 0.2), mgp=c(2, 1, 0), oma=c(0, 0, 2, 0))
    }
  } else {
    sw.n  <- paste0(toupper(substring(plot.meth, 1, 1)), ".sw")
    m.sw[sw.n] <- TRUE
  }
  if(!grp.ind)  {
    if(m.sw["Z.sw"])                  stop("Can't use 'Z' for 'plot.meth' as no clustering has taken place")
    if(vars == "pis")                 stop("Can't plot mixing proportions as no clustering has taken place")
  }
  if(all(m.sw["E.sw"], 
         !attr(results, "Errors")))   stop("Can't plot error metrics as they were not calculated due to storage switches")
  if(all(!m.sw["G.sw"], !m.sw["Z.sw"], !m.sw["E.sw"],
     missing(vars)))                  stop("What variable would you like to plot?")
  if(all(any(m.sw["P.sw"], all.ind),
     is.element(vars, c("means", "uniquenesses")),
     !v.sw[vars],
     is.element(method, c("FA", "IFA")))) {  
    if(all.ind)                       warning(paste0("Can only plot posterior mean, as ", vars, " weren't stored"), call.=FALSE)
   v.sw[vars]     <- !v.sw[vars]
   all.ind        <- FALSE
   m.sw["P.sw"]   <- TRUE
  } 
  if(all(!v.sw[vars], !m.sw["G.sw"], 
     !m.sw["Z.sw"],   !m.sw["E.sw"])) stop(paste0("Nothing to plot: ", vars, " weren't stored"))
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
  if(all(is.element(method, c("IMIFA", "OMIFA")), m.sw["G.sw"])) {
    Gs    <- seq_len(2)
    if(!missing(g))                   warning(paste0("Removed 'g'=", g ," for the ", plot.meth, " plotting method"), call.=FALSE)
  } else if(any(all(is.element(vars, c("scores", "pis", "alpha")), any(all.ind, vars != "scores", !m.sw["P.sw"])), 
            m.sw["G.sw"], m.sw["Z.sw"], m.sw["E.sw"])) {
    Gs    <- 1
  } else if(!missing(g)) {
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
    msg   <- "Hit <Return> to see next plot or type 'EXIT' to exit: "
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
           all(is.element(vars, c("means", "uniquenesses")),
               !indx),
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
          matplot(t(plot.x), type="l", ylab="Means", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else {
          plot(x=iter, y=plot.x[ind,], type="l", ylab="Mean", xlab="Iteration", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1))
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(vars  == "scores") {
        x.plot <- results$Scores$f
        if(by.fac) {
          plot.x  <- if(Q > 1) x.plot[ind[1],,] else t(x.plot[ind[1],,])
        } else {
          plot.x  <- x.plot[,ind[2],]
        }
        if(matx) {
          matplot(t(plot.x), type="l", ylab="Scores", xlab="Iteration")    
          if(by.fac) {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else {
          plot(x=iter, y=x.plot[ind[1],ind[2],], type="l", ylab="Scores", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
        }
      }
      if(vars  == "loadings") {
        x.plot <- result$loadings
        plot.x <- if(by.fac) x.plot[ind[1],,] else x.plot[,ind[2],]
        if(matx) {
          matplot(t(plot.x), type="l", ylab="Loadings", xlab="Iteration")
          if(by.fac) {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else   {
          plot(x=iter, y=x.plot[ind[1],ind[2],], type="l", ylab="Loadings", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$psi
        if(matx) {
          matplot(t(plot.x), type="l", ylab="Uniquenesses", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot(x=iter, y=plot.x[ind,], ylab="Uniquenesses", type="l", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
        }
      }
      if(vars  == "pis") {
        plot.x <- clust$pi.prop
        if(matx) {
          matplot(t(plot.x), type="l", ylab="Mixing Proportions", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMixing Proportions")))))
        } else   {
          plot(x=iter, y=plot.x[ind,], ylab="Mixing Proportions", type="l", xlab="Iteration")
          if(titles) title(main=list(paste0("Trace", ifelse(all.ind, "", paste0(":\nMixing Proportions - Group ", ind)))))
        }
      }
      if(vars  == "alpha") {
        plot.x <- clust$MH.alpha
        plot(plot.x$alpha.pi, ylab="Alpha", type="l", xlab="Iteration", main="")
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
          matplot(plot.x, type="l", ylab="Density")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot.d  <- density(x.plot[ind,])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nMeans - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "scores") {
        x.plot <- results$Scores$f
        if(by.fac) {
          plot.x  <- if(Q > 1) x.plot[ind[1],,] else t(x.plot[ind[1],,])
        } else   {
          plot.x  <- x.plot[,ind[2],]
        }
        if(matx) {
          plot.x  <- sapply(apply(plot.x, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(by.fac) {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]])))
          } else {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Factor ", ind[2])))
          }
        } else   {
          plot.d  <- density(x.plot[ind[1],ind[2],])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", ":\nScores - "), "Observation ", obs.names[ind[1]], ", Factor ", ind[2])))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "loadings") {
        x.plot <- result$loadings
        plot.x    <- if(by.fac) x.plot[ind[1],,] else x.plot[,ind[2],]
        if(matx) {
          plot.x  <- sapply(apply(plot.x, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(by.fac) {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable")))
          } else {
            if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), "Factor ", ind[2])))
          }
        } else   {
          plot.d  <- density(x.plot[ind[1],ind[2],])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nLoadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind[1]], " Variable, Factor ", ind[2])))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "uniquenesses") {
        x.plot <- result$psi
        if(matx) {
          plot.x  <- sapply(apply(x.plot, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        } else   {
          plot.d  <- density(x.plot[ind,])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, ":\n", paste0(":\nUniquenesses - ", ifelse(grp.ind, paste0("Group ", g, " - "), ""))), var.names[ind], " Variable")))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "pis") {
        x.plot <- clust$pi.prop
        if(matx) {
          plot.x  <- sapply(apply(x.plot, 1, density), "[[", "y")
          matplot(plot.x, type="l", ylab="Density")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMixing Proportions")))))
        } else   {
          plot.d  <- density(x.plot[ind,])
          plot(plot.d, main="")
          if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nMixing Proportions - Group ", ind)))))
          polygon(plot.d, col=grey)
        }
      }
      if(vars  == "alpha") {
        plot.x <- clust$MH.alpha
        plot.d <- density(plot.x$alpha.pi)
        plot(plot.d, main="")
        if(titles) title(main=list(paste0("Density", ifelse(all.ind, "", paste0(":\nAlpha")))))
        polygon(plot.d, col=grey)
        if(intervals) abline(v=plot.x$post.alpha, col=2, lty=2)
      }
    }
    
    if(m.sw["P.sw"])  {
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
        plot(plot.x, type=type, ylab="Means", xlab="Variable", ylim=if(is.element(method, c("FA", "IFA"))) c(-1, 1) else if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
        if(all(intervals, ci.sw[vars])) plotCI(plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
        if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMeans", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(x=seq_along(plot.x), y=plot.x, var.names, cex=0.5)
      }
      if(vars  == "scores") {
        labs   <- if(grp.ind) clust$post.z else 1
        if(!missing(labels)) {
          if(!exists(as.character(match.call()$labels),
              envir=.GlobalEnv)) {    warning(paste0("Object ", match.call()$labels, " not found"), call.=FALSE)
          } else {
            labs  <- as.factor(labels)
            if(length(labs) != n.obs) stop(paste0("'labels' must be a factor of length N=",  n.obs))
          }
        }
        if(g.score)  { 
          if(g.ind == 1)  tmplab <- labs
          z.ind  <- as.numeric(levels(tmplab))[tmplab] == g
          plot.x <- results$Scores$post.f[z.ind,,drop=FALSE]
          ind2   <- ifelse(any(!facx, Q <= 1), ind[2], if(Q > 1) max(2, ind[2]))
          if(ci.sw[vars])  ci.x  <- results$Scores$ci.f[,z.ind,, drop=FALSE]
          labs   <- g
        } else       {
          plot.x <- results$Scores$post.f
          ind2   <- ifelse(any(!facx, Q.max <= 1), ind[2], if(Q.max > 1) max(2, ind[2]))
          if(ci.sw[vars])  ci.x  <- results$Scores$ci.f
        }
        col.f  <- if(is.factor(labs)) as.numeric(levels(labs))[labs] + 1 else labs + 1
        type.f <- ifelse(any(type.x, type == "l"), "p", type)
        if(ind2 != 1)  {
          if(all(intervals, ci.sw[vars])) {
            plotCI(plot.x[,ind[1]], plot.x[,ind2], li=ci.x[1,,ind2], ui=ci.x[2,,ind2], gap=TRUE, pch=NA, scol=grey, slty=3, xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
            plotCI(plot.x[,ind[1]], plot.x[,ind2], li=ci.x[1,,ind[1]], ui=ci.x[2,,ind[1]], add=TRUE, gap=TRUE, pch=NA, scol=grey, slty=3, err="x")
            if(type.f != "n") points(plot.x[,ind[1]], plot.x[,ind2], type=type.f, col=col.f, pch=20)
          } else {
            plot(plot.x[,ind[1]], plot.x[,ind2], type=type.f, col=col.f, pch=20,
                 xlab=paste0("Factor ", ind[1]), ylab=paste0("Factor ", ind2))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"), ifelse(g.score, paste0(" - Group ", g), ""))))
          if(type.f == "n") text(plot.x[,ind[1]], plot.x[,ind2], obs.names, 
                               col=col.f, cex=0.5)
        } else   {
          if(all(intervals, ci.sw[vars])) {
            plotCI(if(!g.score) seq_len(n.obs) else seq_len(sum(z.ind)), plot.x[,ind[1]], li=ci.x[1,,ind[1]], ui=ci.x[2,,ind[1]], gap=TRUE, pch=NA, scol=grey, slty=3, xlab="Observation", ylab=paste0("Factor ", ind[1]))
            points(plot.x[,ind[1]], type=type.f, col=col.f, pch=20)
          } else {
            plot(plot.x[,ind[1]], type=type.f, col=col.f, xlab="Observation", ylab=paste0("Factor ", ind[1]), pch=20)
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", ":\nScores"), ifelse(g.score, paste0(" - Group ", g), ""))))
          if(type.f == "n") text(plot.x[,ind[1]], col=col.f, cex=0.5)
        }
      }
      if(vars  == "loadings") {
        if(all(!all.ind, load.all)) {
          par(mai=c(1.25, 1, 0.75, 0.5), mfrow=c(1, 2), oma=c(0, 0, 1, 0))
        }
        plot.x <- result$post.load
        if(is.element(load.meth, c("all", "heatmap"))) {
          if(Q > 1) {
            plotcolors(mat2color(plot.x))
          } else {
            image(z=t(plot.x[seq(n.var, 1),seq_len(Q)]), xlab="", ylab="", xaxt="n", yaxt="n", col=dichromat(heat.colors(30)))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all(!all.ind, !load.all), " Loadings ", " "), "Heatmap", ifelse(all(!all.ind, grp.ind, !load.all), paste0(" - Group ", g), ""))))
          axis(1, line=-0.5, tick=FALSE, at=if(Q != 1) seq_len(Q) else 0, labels=seq_len(Q))
          if(n.var < 100) {
            axis(2, cex.axis=0.5, line=-0.5, tick=FALSE, las=1, at=if(Q > 1) seq_len(n.var) else seq(from=0, to=1, by=1/(n.var - 1)), labels=substring(var.names[n.var:1], 1, 10))
          }
          box(lwd=2)
          mtext(ifelse(Q > 1, "Factors", "Factor"), side=1, line=2)
          if(Q   != 1) abline(v=seq(1, Q - 1, 1) + 0.5, lty=2, lwd=1)
        }
        if(is.element(load.meth, c("all", "raw"))) {
          if(ci.sw[vars])  ci.x  <- result$ci.load  
          if(by.fac) {
            if(ci.sw[vars]) ci.x <- ci.x[,,ind[2]]
            plot(plot.x[,ind[2]], type=type, xaxt="n", xlab="", ylab="Loading", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
            if(all(intervals, ci.sw[vars])) plotCI(plot.x[,ind[2]], li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
            axis(1, line=-0.5, tick=FALSE, at=seq_len(n.var), labels=seq_len(n.var))
            mtext("Variable #", side=1, line=2, cex=0.8)
            if(titles) title(main=list(paste0(ifelse(all(!all.ind, !load.all), paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), "")), ""), "Factor ", ind[2])))
            if(type == "n") text(x=plot.x, var.names, cex=0.5)
          } else     {
            if(ci.sw[vars]) ci.x <- ci.x[,ind[1],]
            plot(plot.x[ind[1],], type=type, xaxt="n", xlab="", ylab="Loading", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
            if(all(intervals, ci.sw[vars])) plotCI(plot.x[ind[1],], li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
            axis(1, line=-0.5, tick=FALSE, at=seq_len(Q), labels=seq_len(Q))
            mtext("Factors", side=1, line=2)
            if(titles) title(main=list(paste0(ifelse(all(!all.ind, !load.all), paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, " - "), "")), ""), var.names[ind[1]], " Variable")))
            if(type == "n") text(x=plot.x[ind[1],], paste0("Factor ", seq_len(Q)), cex=0.5)
          }
        }
        if(all(!all.ind, load.all)) {
          if(titles) title(paste0("Loadings", ifelse(grp.ind, paste0(" - Group ", g), "")), outer=TRUE)
        }
      }
      if(vars  == "uniquenesses") {
        plot.x <- result$post.psi
        if(ci.sw[vars])   ci.x   <- result$ci.psi
        plot(plot.x, type=type, ylab="Uniquenesses", xlab="Variable", ylim=if(ci.sw[vars]) c(min(ci.x[1,]), max(ci.x[2,])))
        if(all(intervals, ci.sw[vars])) plotCI(plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol=grey, add=TRUE, gap=TRUE, pch=ifelse(type == "n", NA, 20))
        if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nUniquenesses", ifelse(grp.ind, paste0(" - Group ", g), ""))))))
        if(type  == "n") text(seq_along(plot.x), plot.x, var.names, cex=0.5)
      }
      if(vars  == "pis") {
        plot.x <- clust$post.pi
        if(ci.sw[vars])   ci.x   <- clust$ci.pi
        if(matx) {
          if(all(intervals, ci.sw[vars])) {
            plotCI(barplot(plot.x, ylab="Mixing Proportions", xlab="", col=grey, ylim=c(0, 1)),
                   plot.x, li=ci.x[1,], ui=ci.x[2,], slty=3, scol=2, add=TRUE, gap=TRUE, pch=20)
          } else {
            barplot(plot.x, ylab="Mixing Proportions", xlab="", ylim=c(0, 1))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMixing Proportions")))))  
        } else {
          if(all(intervals, ci.sw[vars])) {
            plotCI(barplot(plot.x[ind], ylab="Mixing Proportions", xlab="", ylim=c(0, 1)),
                   plot.x[ind], li=ci.x[1,ind], ui=ci.x[2,ind], slty=3, scol=2, add=TRUE, gap=TRUE, pch=20)
          } else {
            barplot(plot.x[ind], ylab="Mixing Proportions", xlab="Variable", ylim=c(0, 1))
          }
          if(titles) title(main=list(paste0("Posterior Mean", ifelse(all.ind, "", paste0(":\nMixing Proportions - Group ", ind)))))
        }
      }
      if(vars  == "alpha") {
        plot(c(0, 1), c(0, 1), ann=FALSE, bty='n', type='n', xaxt='n', yaxt='n')
        if(titles) title(main=list(paste0("Summary Statistics", ifelse(all.ind, "", ":\nAlpha"))))
        plot.x <- clust$MH.alpha[-1]
        conf   <- attr(results, "Conf.Level")
        digits <- options()$digits
        a.cex  <- par()$fin[2]/5
        a.adj  <- rep(0.5, 2)
        text(x=0.5, y=0.775, cex=a.cex, col="black", adj=a.adj, expression(bold("Posterior Mean:\n")))
        text(x=0.5, y=0.775, cex=a.cex, col="black", adj=a.adj, bquote(.(round(plot.x$post.alpha, digits))))
        text(x=0.5, y=0.55,  cex=a.cex, col="black", adj=a.adj, expression(bold("\nVariance:\n")))
        text(x=0.5, y=0.55,  cex=a.cex, col="black", adj=a.adj, bquote(.(round(plot.x$var.alpha, digits))))
        text(x=0.5, y=0.4,   cex=a.cex, col="black", adj=a.adj, bquote(bold(.(100 * conf))~bold("% Confidence Interval:")))
        text(x=0.5, y=0.3,   cex=a.cex, col="black", adj=a.adj, bquote(paste("[", .(round(plot.x$ci.alpha[1], digits)), ", ", .(round(plot.x$ci.alpha[2], digits)), "]")))
        text(x=0.5, y=0.175, cex=a.cex, col="black", adj=a.adj, expression(bold("Acceptance Rate:")))
        text(x=0.5, y=0.1,   cex=a.cex, col="black", adj=a.adj, bquote(paste(.(round(100 * plot.x$acceptance.rate, 2)), "%")))
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
      if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
        aic.mcmc <- round(GQ.res$AIC.mcmc, 2)
        bic.mcmc <- round(GQ.res$BIC.mcmc, 2)
      }
      if(all(plotG.ind, g == 1))  {
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
      if(is.element(method, c("IMIFA", "OMIFA"))) {
        if(g == 1) {
          print(GQ.res[substring(names(GQ.res), 1, 1) == "G"])
        } else {
          print(GQ.res[substring(names(GQ.res), 1, 1) != "G"])
        }
      } else if(is.element(method, c("MFA", "MIFA", "OMFA", "IMFA"))) {
          print(GQ.res)
      } else if(method == "IFA") {
          print(tail(GQ.res, -1))
      } else   {
          cat(paste0("Q = ", Q, "\n"))
      }
      if(any(nrow(bicm) > 1, ncol(bicm) > 1)) {
        G.ind  <- ifelse(any(G.supp, !is.element(method, c("MFA", "MIFA"))), 1, which(n.grp == G))
        Q.ind  <- ifelse(any(Q.supp, !is.element(method, c("FA", "MFA"))),   1, which(n.fac == Q))
          cat(paste0("AICM = ", aicm[G.ind,Q.ind], "\n"))
          cat(paste0("BICM = ", bicm[G.ind,Q.ind], "\n"))
        if(!is.element(method, c("IFA", "MIFA"))) {
          cat(paste0("AIC.mcmc = ", aic.mcmc[G.ind,Q.ind], "\n"))
          cat(paste0("BIC.mcmc = ", bic.mcmc[G.ind,Q.ind], "\n"))
        }
      }
    }
    
    if(m.sw["Z.sw"]) {
      if(type == "l")       stop("'type' cannot be 'l' for clustering uncertainty plots")
      plot.x <- clust$uncertainty
      col.x  <- c(1, 4)[(plot.x >= 1/G) + 1]
      if(type != "h") col.x[plot.x == 0] <- NA
      if(titles) {
        layout(rbind(1, 2), heights=c(1, 6))
        par(mar=c(0, 4.1, 0.5, 2.1))
        plot.new()
        legend("center", legend=paste0("1/G = 1/", G), title="", lty=2, col=2, bty="n", y.intersp=par()$fin[2] * 5/4)
        legend("center", legend=c(" "," "), title=expression(bold("Clustering Uncertainty")), bty='n', y.intersp=par()$fin[2] * 2/5, cex=par()$cex.main)
        par(mar=c(5.1, 4.1, 0, 2.1))
      }
      plot(plot.x, type=type, ylim=c(0, 1 - 1/G), col=col.x, axes=FALSE, ylab="Uncertainty", xlab="Observation", pch=ifelse(type == "n", NA, 16))
      rect(0, 0, n.obs, 1 - 1/G) 
      axis(1, las=1, pos=0)
      axis(2, las=2, pos=0)
      lines(x=c(0, n.obs), y=c(1/G, 1/G), lty=2, col=2)
      axis(2, at=1 - 1/G, label="1 - 1/G", las=2, line=-0.7, tick=T)
      if(type == "n")  {
        znam  <- obs.names
        znam[plot.x == 0] <- ""
        text(x=seq_along(plot.x), y=plot.x, znam, col=col.x, cex=0.5)
      }
      if(any(!labelmiss, !missing(labels))) {
        if(!labelmiss) {
          perf    <- clust$perf
        } else   {
          labs    <- as.factor(labels)
          if(length(labs) != n.obs)   stop(paste0("'labels' must be a factor of length N=",  n.obs))
          pzs     <- clust$post.z
          if(nlevels(pzs) == nlevels(labs)) {
            l.sw  <- lab.switch(z.new=pzs, z.old=labs, Gs=seq_len(G))
            pzs   <- factor(l.sw$z)
          }
          tab     <- table(pzs, labs, dnn=list("Predicted", "Observed"))
          perf    <- c(classAgreement(tab), classError(pzs, labs))
          if(nrow(tab) != ncol(tab))   {
            perf  <- perf[-seq_len(2)]
            names(perf)[4]       <- "error.rate"
          } else {
            names(perf)[6]       <- "error.rate"
          }
          if(perf$error.rate     == 0) {
            perf$misclassified   <- NULL
          }
          perf    <- c(list(confusion.matrix = tab), perf)
          if(nlevels(pzs)  == nlevels(labs)) {
            names(perf)[1] <- "matched.confusion.matrix"
          }
          class(perf)      <- "listof"
        }
        uncert    <- attr(plot.x, "Obs")
        if(!is.null(uncert)) {
          perf    <- c(perf, list(uncertain = uncert))  
        }
        perf$error.rate    <- paste0(round(100 * perf$error.rate, 2), "%")
        print(perf)
      } else                          message("Nothing to print: try supplying known cluster labels")
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
      na.x   <- if(G > 1) is.na(res$Error[[1]]) else F
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
        if(xna) text(x=temp$text$x[6] - 0.015, y=temp$text$y[6] + 0.015, "__")
      }  
      avg    <- if(G > 1) setNames(list(x.plot$Averages), "Average Error Metrics") else x.plot
      class(avg)           <- "listof"
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
          acf(plot.x[ind,], main="", ci.col=4)
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind], " Variable"), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial", ci.col=4)
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind], " Variable"), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Means - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind], " Variable")), outer=TRUE)
        }
      }
      if(vars  == "scores")   { 
        plot.x <- results$Scores$f
        if(!partial) {
          acf(plot.x[ind[1],ind[2],], main="", ci.col=4)
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", "Observation ", obs.names[ind[1]], ", Factor ", ind[2]), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind[1],ind[2],], main="", type="partial", ci.col=4)
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", "Observation ", obs.names[ind[1]], ", Factor ", ind[2]), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Scores - ", "Observation ", obs.names[ind[1]], ", Factor ", ind[2])), outer=TRUE)
        }
      }
      if(vars  == "loadings") { 
        plot.x <- result$loadings
        if(!partial) {
          acf(plot.x[ind[1],ind[2],], main="", ci.col=4)
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind[1]], " Variable, Factor ", ind[2]), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind[1],ind[2],], main="", type="partial", ci.col=4)
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind[1]], " Variable, Factor ", ind[2]), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Loadings - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind[1]], " Variable, Factor ", ind[2])), outer=TRUE)
        }
      }
      if(vars  == "uniquenesses")  { 
        plot.x <- result$psi
        if(!partial) {
          acf(plot.x[ind,], main="", ci.col=4)
          if(titles) title(main=list(paste0("ACF", ifelse(all.ind, paste0(":\n", var.names[ind], " Variable"), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial", ci.col=4)
          if(titles) title(main=list(paste0("PACF", ifelse(partial, paste0(":\n", var.names[ind], " Variable"), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Uniquenesses - ", ifelse(grp.ind, paste0("Group ", g, ":\n "), ""), var.names[ind], " Variable")), outer=TRUE)
        }
      }
      if(vars  == "pis")  { 
        plot.x <- clust$pi.prop
        if(!partial) {
          acf(plot.x[ind,], main="", ci.col=4)
          if(titles) title(main=list(paste0("ACF", ifelse(all(all.ind, matx), paste0(" - Group ", ind), ""))))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x[ind,], main="", type="partial", ci.col=4)
          if(titles) title(main=list(paste0("PACF", ifelse(all(all.ind, matx), paste0(" - Group ", ind), ""))))
          if(all(!all.ind, titles)) title(main=list(paste0("Mixing Proportions - Group ", ind)), outer=TRUE)
        }
      }
      if(vars  == "alpha") {          
        plot.x <- clust$MH.alpha$alpha.pi
        if(clust$MH.alpha$acceptance.rate == 0)  {
                                      warning(paste0("0% acceptance rate: can't plot ", ifelse(all.ind, ifelse(partial, "partial-", "auto-"), ""), "correlation function", ifelse(all.ind, "", "s")), call.=FALSE)
          next
        }
        if(!partial) {
          acf(plot.x, main="", ci.col=4)
          if(titles) title(main=list(paste0("ACF")))
        }
        if(any(!all.ind, partial)) {
          acf(plot.x, main="", type="partial", ci.col=4)
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