###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.imifa      <- function(N = 1000, G = 2, P = 25, Q = rep(5, G)) {
  
  if(any(N  < 0, P  < 0, Q < 0, G <= 0)) stop("N, P, and Q must be strictly non-negative and G must be strictly positive")
  if(any(Q >= P, Q >= N - 1))            stop(paste0("Cannot generate this many factors relative to N=", N, " and P=", P))
  if(length(Q)  != G)                    
  if(!missing(Q)) {
   if(length(Q) == 1) {
    Q    <- rep(Q, G)
   } else if(length(Q != G)) {           stop(paste0("Q must supplied for each of the G=", G, " groups"))  
   }                                       
  }
  nn     <- rnorm(G, N/G, exp(G * G))
  if(abs(sum(nn)) < 0.01) nn <- nn + 1
  nn     <- round(nn / sum(nn) * N)
  dev    <- N - sum(nn)
  for(. in seq_len(abs(dev))) {
    nn[i]       <- nn[i <- sample(G, 1)] + sign(dev)
  }
  while(any(nn   < 0)) {
    neg  <- nn   < 0
    pos  <- nn   > 0
    nn[neg][i]  <- nn[neg][i <- sample(sum(neg), 1)] + 1
    nn[pos][i]  <- nn[pos][i <- sample(sum(pos), 1)] - 1
  }
  if(any(nn == 0)) {
    G    <- G - sum(nn == 0)
    Q    <- Q[1:length(G)]
    nn   <- nn[nn > 0]
  }
  
  gnames        <- paste0("Group ", 1:G)
  vnames        <- paste0("Var ", 1:P)
  SimData       <- matrix(0, nr=0, nc=P)
  true.mu       <- setNames(vector("list", G), gnames)
  true.l        <- setNames(vector("list", G), gnames)
  true.psi      <- setNames(vector("list", G), gnames)
  true.cov      <- setNames(vector("list", G), gnames)
  true.zlab     <- rep(0, N)
  
# Simulate true parameter values
  for(g in 1:G) {
    Q.g         <- Q[g]
    U.mu        <- sqrt(abs(rnorm(P, Q.g, g)) * diag(P))
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- U.mu %*% z.mu
    mu.mu       <- rnorm(P, 0, max(Q.g, 1))
    mu.true     <- as.vector(mu.mu + v.mu)
    
    if(Q.g > 0) {
      U.load    <- sqrt(abs(rnorm(Q.g, 0, g/2)) * diag(Q.g))
      z.load    <- matrix(rnorm(P * Q.g, 0, 1), nr=Q.g, ncol=P)
      v.load    <- t(U.load %*% z.load)
      mu.load   <- rnorm(Q.g, 0, Q.g/2)
      l.true    <- mu.load + v.load
    } else {
      l.true    <- matrix(, nr=P, nc=0)
    }
  
    shape.g     <- floor((g * max(Q.g, 2))/2)
    psi.true    <- 1/rgamma(n=P, shape=shape.g, rate=shape.g - 1)                            
    
  # Simulate data
    covmat      <- tcrossprod(l.true, l.true) + diag(psi.true)
    if(Q.g > 0) {
      U.om      <- chol(covmat)
    } else {
      U.om      <- sqrt(covmat)
    }
    z           <- matrix(rnorm(P * nn[g], 0, 1), nr=P, ncol=nn[g])
    var         <- U.om %*% z
    SimData     <- rbind(SimData, t(rep(mu.true, nn[g]) + var))
    names(mu.true)      <- vnames
    names(psi.true)     <- vnames
    dimnames(covmat)    <- list(vnames, vnames)
    dimnames(l.true)    <- list(vnames, if(Q.g > 0) paste0("Factor ", 1:Q.g))
    true.mu[[g]]        <- mu.true
    true.l[[g]]         <- l.true
    true.psi[[g]]       <- psi.true
    true.cov[[g]]       <- covmat
  }
  
# Post-process data
  true.zlab     <- factor(rep(nn, nn), labels=1:length(nn))
  permute       <- sample(1:N, N, replace=F)
  SimData       <- SimData[permute,]
  true.zlab     <- true.zlab[permute]
  dimnames(SimData)     <- list(paste0("Obs ", 1:N), vnames)
  SimData       <- as.data.frame(SimData)
  attr(SimData, 
       "Means")         <- true.mu
  attr(SimData, 
       "Loadings")      <- true.l
  attr(SimData, 
       "Uniquenesses")  <- true.psi
  attr(SimData, 
       "Labels")        <- true.zlab
  attr(SimData,
       "Covariance")    <- true.cov
  class(SimData)        <- c("data.frame", "IMIFA")
  return(SimData)
}