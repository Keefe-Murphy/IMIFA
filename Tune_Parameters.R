#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

if(case == 'Shrinkage') {
  # Retrieve distribution of Q, tabulate & plot
    Q.store    <- sim[[Q.ind]]$Q.store
    Q.tab      <- table(Q.store)
    Q.plot     <- barplot(Q.tab, main="Posterior Distribution of Q", 
                          ylab="Frequency", xlab="Q", xaxt="n")
    axis(1, at=Q.plot, labels=names(Q.tab), tick=F)
  
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.mode     <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.median   <- median(Q.store)
    if(post.Q == 'Mode') { 
      Q        <- min(Q.mode)
    } else { Q <- Q.median }
    Q.CI       <- quantile(Q.store,c(0.025,0.975))
    print(list(Q=Q, Mode = Q.mode, Median = Q.median, 
               Credible_Interval = Q.CI, Warning="But the user should choose Q based on the attached bar plot!"))
  } else {
  
  # Initialise
    Q.star   <- range.Q - min(range.Q) + 1
    P        <- length(sim[[1]]$psi[,1])
    prop.var <- rep(NA, length(range.Q))

  # Calculate Proportion of Variation Explained
    for(Q in Q.star) {
      store       <- seq(from=burnin + 1, to=sim[[1]]$n.store, by=thin)
      psi         <- sim[[Q]]$psi[,store]
      post.psi    <- apply(psi, 1, mean)
      prop.var[Q] <- (P - sum(post.psi))/P
    }
  
  # Produce Scree Plot & choose optimum Q
    plot(prop.var, type="l", main="Scree Plot to Choose Q", xlab="# Factors", 
         ylab="% Variation Explained", xaxt="n", yaxt="n", ylim=c(0,1))
    axis(1, at=1:length(prop.var), labels=range.Q)
    axis(2, at=seq(0,1,0.1), labels=seq(0,100,10), cex.axis=0.8)
    Q.ind <- which.max(prop.var)
    Q     <- range.Q[Q.ind]
    points(x=Q.ind, y=prop.var[Q.ind], col="red", cex=1.5, bg="red", pch=21)
    print(list(Q=Q, Warning="But the user should choose Q based on the attached scree plot!"))
  }