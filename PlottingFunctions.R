################################
### IMIFA Plotting Functions ###
################################

# Plot the cumulative proportion of variation explained 
  plot.cum.var  <- function(n.fac=res$Q, cum.var=res$cum.var) {
    prop.exp    <- cum.var[n.fac]
    if(n.fac     > 1) {
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
  load.heat     <- function(Q = res$Q) {
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