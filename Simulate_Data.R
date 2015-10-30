###################################
### Simulate Data (Single Case) ###
###################################

P <- 50
Q <- 10
N <- 1500
mu.true   <- mvrnorm(mu=rep(0, P), Sigma=diag(P));      names(mu.true)      <- c(1:P)
f.true    <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q)); colnames(f.true)    <- paste("Factor", 1:Q); rownames(f.true)    <- c(1:N)
load.true <- mvrnorm(n=P, mu=rep(0, Q), Sigma=diag(Q)); colnames(load.true) <- paste("Factor", 1:Q); rownames(load.true) <- c(1:P)
psi.true  <- rinvgamma(n=P, 1, 1);                      names(psi.true)     <- c(1:P)
eps.true  <- mvrnorm(n=N, mu=rep(0, P), Sigma=diag(psi.true))
data      <- matrix(0, nr=N, nc=P)

for (i in 1:N) {  
  data[i, ]    <- mu.true + load.true %*% f.true[i,] + eps.true[i,]
}

rownames(data) <- c(1:N)
colnames(data) <- c(1:P)
data           <- as.data.frame(data)