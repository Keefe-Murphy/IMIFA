################################
### IMIFA Plotting Functions ###
################################

# Plot the cumulative proportion of variation explained 
  plot.cum.var   <- function(n.fac=res$Q, cum.var=res$cum.var, ...) {
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
  load.heat      <- function(Q = res$Q, ...) {
    if(Q > res$Q)  stop("Cannot plot this many loadings columns")
    par(mfrow=c(1, 1), mar=c(5.1, 7.1, 4.1, 2.1))
    image(z=t(post.load[,1:Q]), xlab="", ylab="", 
          main="Posterior Loadings", xaxt="n", yaxt="n")
    axis(1, cex.axis=0.8, line=-0.5, tick=F, 
         at=if(Q != 1) seq(0, 1, 1/(Q-1)) else 0, labels=1:Q)
    axis(2, cex.axis=0.5, line=-0.5, tick=F, las=1,
         at=seq(0, 1, 1/(P-1)), labels=rownames(post.load))
    box(lwd=2)
    mtext("Factors", side=1, line=2)
    if(Q != 1) abline(v=seq(1/(2*(Q-1)), 1-1/(2*(Q-1)), 1/(Q-1)), lty=2, lwd=1)
    par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))
  }

# Plot posteriors
  plot.posterior <- function(x=NULL, type=c("n", "l", "p"), cols=NULL, ...) {
    x    <- match.arg(x, c("means", "scores", "loadings", "uniquenesses"))
    type <- match.arg(type)
    if(x == "scores" | x == "loadings") {
      if(missing(cols)) cols <- c(1, 2)
      if(length(cols) > 2)   stop("Only two columns can be plotted")
      if(max(cols) > res$Q)  stop(paste0("Only the first ", res$Q, " columns can be plotted"))
    }
    if(x == "means") {
      plot.x  <- post.mu
      plot(plot.x, type=type, main="Posterior Means", ylab="Means")
      if(type == "n") text(x=1:length(plot.x), y=plot.x, names(plot.x), cex=0.5)
    }
    if(x == "scores") {
      plot.x  <- post.f
      if(exists("Label")) plot.col <- as.numeric(Label) else plot.col <- 1
      plot(plot.x, type=type, main="Posterior Scores", col=plot.col)
      if(type == "n") text(plot.x[,cols[1]], plot.x[,cols[2]], 1:nrow(plot.x), 
                                col=plot.col, cex=0.5)
    }
    if(x == "loadings") {
      plot.x  <- post.load
      plot(plot.x, type=type, main="Posterior Loadings")
      if(type == "n") text(plot.x[,cols[1]], plot.x[,cols[2]], rownames(plot.x), cex=0.5)
    }
    if(x == "uniquenesses") {
      plot.x  <- post.psi
      plot(plot.x, type=type, main="Posterior Uniquenesses", ylab="Uniquenesses")
      if(type == "n") text(1:length(plot.x), plot.x, names(plot.x), cex=0.5)
    }
  } 