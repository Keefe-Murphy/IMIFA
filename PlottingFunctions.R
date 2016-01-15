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