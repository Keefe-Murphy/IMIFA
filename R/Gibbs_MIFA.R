################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  .gibbs_MIFA      <- function(Q, data, iters, N, P, G, sw, mu, mu.zero, uni.type,
                               sigma.mu, burnin, thinning, verbose, nu, cluster, 
                               psi.alpha, psi.beta, adapt, adapt.at, prop, b0, 
                               b1, beta.d1, beta.d2, epsilon, nuplus1, ...) {
        
  # Define & initialise variables
    start.time     <- proc.time()
    total          <- max(iters)
    if(verbose)       pb    <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    Nseq           <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Gseq)
    iternames      <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store     <- provideDimnames(array(0, dim=c(P, G, n.store)), base=list(varnames, gnames, iternames))
    }
    if(sw["s.sw"])   {
      eta.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)   <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store   <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store)  <- list(varnames, if(Q > 0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store    <- provideDimnames(array(0, dim=c(P, G, n.store)), base=list(varnames, gnames, iternames))
    }
    if(sw["pi.sw"])  {
      pi.store     <- provideDimnames(matrix(0, nr=G, nc=n.store), base=list(gnames, iternames))
    }
    z.store        <- provideDimnames(matrix(0, nr=N, nc=n.store), base=list(obsnames, iternames))
    ll.store       <- stats::setNames(rep(0, n.store), iternames)
    Q.star         <- Q
    Qs             <- rep(Q, G)
    Q.store        <- provideDimnames(matrix(0, nr=G, nc=n.store), base=list(gnames, iternames))
    Q.large        <- Q.big <- Q.bigs <- FALSE
    err.z          <- z.err <- FALSE
    
    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    if(all(mu.zero == 0)) {
      mu.zero      <- matrix(0, nr=1, nc=G)
      cluster$l.switch[1]   <- FALSE
    }
    if(length(mu.zero)  == 1) {
      mu.zero      <- matrix(mu.zero,  nr=1, nc=G)
    }
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    pi.prop        <- cluster$pi.prop
    pi.alpha       <- cluster$pi.alpha
    .sim_psi.inv   <- switch(uni.type, unconstrained=.sim_psi.iu,  isotropic=.sim_psi.ii)
    .sim_psi.ip    <- switch(uni.type, unconstrained=.sim_psi.ipu, isotropic=.sim_psi.ipi)
    psi.beta       <- unique(round(psi.beta, min(nchar(psi.beta))))
    if(length(psi.beta) == 1) {
      psi.beta     <- matrix(psi.beta, nr=1, nc=G)
    }
    alpha.d1       <- cluster$alpha.d1
    alpha.d2       <- cluster$alpha.d2
    ad1.x          <- length(unique(alpha.d1)) == 1
    adk.x          <- length(unique(alpha.d2)) == 1
    mu0g           <- cluster$l.switch[1]
    psi0g          <- cluster$l.switch[2]
    delta0g        <- cluster$l.switch[3]
    qstar0g        <- cluster$l.switch[4]
    label.switch   <- any(cluster$l.switch)
    eta            <- .sim_eta.p(N=N, Q=Q)
    phi            <- lapply(Gseq, function(g) .sim_phi.p(Q=Q, P=P, nu=nu, plus1=nuplus1))
    delta          <- lapply(Gseq, function(g) c(.sim_delta.p(alpha=alpha.d1[g], beta=beta.d1), .sim_delta.p(Q=Q, alpha=alpha.d2[g], beta=beta.d2)))
    tau            <- lapply(delta, cumprod)
    lmat           <- lapply(Gseq, function(g) matrix(unlist(lapply(Pseq, function(j) .sim_load.ps(Q=Q, phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nr=P, byrow=TRUE))
    psi.inv        <- vapply(Gseq, function(g) .sim_psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    if(Q < .ledermann(N, P)) {
      fact.ind     <- nn    <= P
      fail.gs      <- which(fact.ind)
      for(g in which(!fact.ind)) {
        fact       <- try(stats::factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]      <- fact$scores
          lmat[[g]]         <- fact$loadings
          psi.inv[,g]       <- 1/fact$uniquenesses
        } else  {
          fail.gs  <- c(fail.gs, g)
        }               
      }
      fail.gs      <- fail.gs[Rfast::Order(fail.gs)]
      len.fail     <- length(fail.gs)
      if(len.fail   > 0)       message(paste0("Parameters of the following group", ifelse(len.fail > 2, "s ", " "), "were initialised by simulation from priors, not factanal: ", ifelse(len.fail > 1, paste0(paste0(fail.gs[-len.fail], sep="", collapse=", "), " and ", fail.gs[len.fail]), fail.gs), " - G=", G, ", Q=", Q))
    } else     {
      psi.tmp      <- psi.inv
      psi.inv      <- vapply(Gseq, function(g) if(nn[g] > 1) 1/Rfast::colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]      <- psi.tmp[inf.ind]
    }
    if(burnin       < 1)  {
      mu.store[,,1]         <- mu
      eta.store[,,1]        <- eta
      load.store[,,,1]      <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
      psi.store[,,1]        <- 1/psi.inv
      pi.store[,1]          <- pi.prop
      z.store[,1]           <- z
      ll.store[1]           <- .sim_z(data=data, mu=mu, Gseq=Gseq, N=N, pi.prop=pi.prop, sigma=lapply(Gseq, function(g) 
                               corpcor::make.positive.definite(tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g]))), Q0=Qs > 0)$log.like
      Q.store[,1]           <- Qs
    }
    init.time      <- proc.time() - start.time
    
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose   && iter   < burnin) utils::setTxtProgressBar(pb, iter)
      
    # Mixing Proportions
      pi.prop[]    <- .sim_pi(pi.alpha=pi.alpha, nn=nn)
      
    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      Q0           <- Qs  > 0
      Q1           <- Qs == 1
      z.log        <- utils::capture.output({ z.res <- try(.sim_z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0), silent=TRUE) })
      z.err        <- inherits(z.res, "try-error")
      if(z.err) {
        sigma      <- lapply(sigma, corpcor::make.positive.definite)
        z.res      <- .sim_z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0)
      }
      z            <- z.res$z
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn  > 0
      nn.ind       <- which(nn0)
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])
    
    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        eta        <- base::matrix(0, nr=N, nc=0)
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
        lmat       <- lapply(Gseq, base::matrix, 0, nr=P, nc=0)
      } else {
        eta.tmp    <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else base::matrix(0, nr=ifelse(Q0[g], 0, nn[g]), nc=Qs[g]))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) matrix(unlist(lapply(Pseq, function(j) .sim_load.s(Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]][,j], 
                             EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nr=P, byrow=TRUE) else 
                             matrix(unlist(lapply(Pseq, function(j) .sim_load.ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nr=P, byrow=FALSE))
        eta.tmp    <- if(length(unique(Qs)) != 1) lapply(Gseq, function(g) cbind(eta.tmp[[g]], matrix(0, nr=nn[g], nc=max(Qs) - Qs[g]))) else eta.tmp
        q0ng       <- !Q0 & nn0
        if(any(q0ng)) {
          eta.tmp[q0ng]     <- lapply(Gseq[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
        }
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      }
      
    # Means
      sum.data     <- vapply(dat.g, colSums, numeric(P))
      sum.eta      <- lapply(eta.tmp, colSums)
      mu           <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.zero=mu.zero[,g], sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                             sum.data=sum.data[,g], lmat=lmat[[g]], N=nn[g], P=P) else .sim_mu.p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero[,g]), numeric(P))
      
    # Uniquenesses
      psi.inv      <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta[,g], lmat=lmat[[g]],
                             P=P, eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    
    # Local Shrinkage
      load.2       <- lapply(lmat, .power2)
      phi          <- lapply(Gseq, function(g) if(nn0[g]) .sim_phi(Q=Qs[g], P=P, nu=nu, plus1=nuplus1,
                      tau=tau[[g]], load.2=load.2[[g]]) else .sim_phi.p(Q=Qs[g], P=P, nu=nu, plus1=nuplus1))
    
    # Global Shrinkage
      sum.terms    <- lapply(Gseq, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gseq)  {
        Qg         <- Qs[g]
        Q1g        <- Q1[g]
        if(nn0[g])    {
          for(k in seq_len(Qg)) { 
            delta[[g]][k]   <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2[g], beta.d2=beta.d2, delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P=P,
                                                    Q=Qg, k=k, sum.term.kq=sum.terms[[g]][k:Qg]) else .sim_delta1(Q=Qg, P=P, tau=tau[[g]], sum.term=sum.terms[[g]],
                                                                                                                 alpha.d1=ifelse(Q1g, alpha.d2[g], alpha.d1[g]), beta.d1=ifelse(Q1g, beta.d2, beta.d1), delta.1=delta[[g]][1])
            tau[[g]]        <- cumprod(delta[[g]])
          }
        } else {
          for(k in seq_len(Qg)) { 
            delta[[g]][k]   <- if(k > 1) .sim_delta.p(alpha=alpha.d2[g], beta=beta.d2) else .sim_delta.p(alpha=ifelse(Q1g, alpha.d2[g], alpha.d1[g]), beta=ifelse(Q1g, beta.d2, beta.d1))
            tau[[g]]        <- cumprod(delta[[g]])
          }
        }
      }
    
    # Adaptation  
      if(all(adapt, iter > adapt.at)) {      
        if(stats::runif(1)   < ifelse(iter < burnin, 0.5, 1/exp(b0 + b1 * (iter - adapt.at)))) {
          colvec   <- lapply(nn.ind, function(g) (if(Q0[g]) colSums(abs(lmat[[g]]) < epsilon)/P else 0) >= prop)
          nonred   <- lapply(colvec, .which0)
          numred   <- lengths(colvec) - lengths(nonred)
          notred   <- numred == 0
          ng.ind   <- seq_along(nn.ind)
          Qs.old   <- Qs[nn0]
          Qs[nn0]  <- vapply(ng.ind, function(h) if(notred[h]) Qs.old[h] + 1 else Qs.old[h] - numred[h], numeric(1))
          Q.big    <- Qs[nn0] > Q.star
          Q.bigs   <- any(Q.big)
          if(Q.bigs) {
            notred <- notred & !Q.big
            Qs[nn0][Q.big]  <- Q.star
          }
          phi[nn0]          <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(phi[[g]][,seq_len(Qs.old[h])], stats::rgamma(n=P, shape=nu + nuplus1, rate=nu)) else phi[[g]][,nonred[[h]], drop=FALSE])
          delta[nn0]        <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) c(delta[[g]][seq_len(Qs.old[h])], stats::rgamma(n=1, shape=alpha.d2, rate=beta.d2)) else delta[[g]][nonred[[h]]])  
          tau[nn0]          <- lapply(delta[nn.ind], cumprod)
          lmat[nn0]         <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(lmat[[g]][,seq_len(Qs.old[h])], stats::rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]])))) else lmat[[g]][,nonred[[h]], drop=FALSE])
          Qemp     <- Qs[!nn0]
          Fmax     <- max(Qs[nn0])
          Qmax     <- ifelse(all(Q.big), Fmax, max(Qs[nn0][!Q.big]))
          Qmaxseq  <- seq_len(Qmax)
          Qmaxold  <- max(Qs.old)
          eta      <- if(all(Fmax  > Qmaxold, !Q.bigs)) cbind(eta[,seq_len(Qmaxold)], stats::rnorm(N)) else eta[,seq_len(Fmax), drop=FALSE]
          if(Qmax   < max(Qemp, 0)) {
            Qs[Qmax < Qs & !nn0]  <- Qmax
            for(g  in Gseq[!nn0][Qemp > Qmax]) {  
              phi[[g]]      <- phi[[g]][,Qmaxseq,  drop=FALSE]
              delta[[g]]    <- delta[[g]][Qmaxseq]
              tau[[g]]      <- tau[[g]][Qmaxseq]
              lmat[[g]]     <- lmat[[g]][,Qmaxseq, drop=FALSE]
            }
          }
        }
      }
    
    # Label Switching
      if(label.switch)   {
        switch.lab <- .lab.switch(z.new=z, z.old=z.temp, Gs=Gseq)
        z          <- switch.lab$z
        z.perm     <- switch.lab$z.perm
        if(!identical(as.integer(z.perm), Gseq)) {
         if(length(unique(Qs)) != 1) {
          Qs       <- Qs[z.perm]  
         }
         mu        <- mu[,z.perm, drop=FALSE]
         lmat      <- lmat[z.perm]
         delta     <- delta[z.perm]
         phi       <- phi[z.perm]
         tau       <- tau[z.perm]
         psi.inv   <- psi.inv[,z.perm, drop=FALSE]
         pi.prop   <- pi.prop[z.perm]
         nn        <- nn[z.perm]
         if(mu0g)        {
          mu.zero  <- mu.zero[,z.perm, drop=FALSE]
         }
         if(psi0g)       {
          psi.beta <- psi.beta[,z.perm, drop=FALSE]
         }
         if(all(delta0g, 
                !ad1.x)) {
          alpha.d1 <- alpha.d1[z.perm]
         }
         if(all(delta0g, 
                !adk.x)) {
          alpha.d2 <- alpha.d2[z.perm]
         } 
        }
      }
      
      if(Q.bigs && !Q.large && iter > burnin) {        warning(paste0("Q has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star), call.=FALSE)
        Q.large    <- TRUE
      }
      if(z.err  && !err.z) {                           warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z      <- TRUE
      }
      if(is.element(iter, iters))   {
        if(verbose)   utils::setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])    mu.store[,,new.it]       <- mu  
        if(all(sw["s.sw"], 
           any(Q0)))  eta.store[,seq_len(max(Qs)),new.it]  <- eta
        if(sw["l.sw"])   {
          for(g in Gseq) {
            if(Q0[g]) load.store[,seq_len(Qs[g]),g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])   psi.store[,,new.it]      <- psi
        if(sw["pi.sw"])    pi.store[,new.it]        <- pi.prop
                           z.store[,new.it]         <- z 
                           ll.store[new.it]         <- z.res$log.like
                           Q.store[,new.it]         <- Qs
      }
    }
    close(pb)
    Qmax           <- seq_len(max(Q.store))
    eta.store      <- if(sw["s.sw"])  tryCatch(eta.store[,Qmax,, drop=FALSE],   error=function(e) eta.store)
    load.store     <- if(sw["l.sw"])  tryCatch(load.store[,Qmax,,, drop=FALSE], error=function(e) load.store)
    returns        <- list(mu       = if(sw["mu.sw"])  mu.store,
                           eta      = if(sw["s.sw"])   tryCatch(slam::as.simple_sparse_array(eta.store),  error=function(e) eta.store), 
                           load     = if(sw["l.sw"])   tryCatch(slam::as.simple_sparse_array(load.store), error=function(e) load.store),
                           psi      = if(sw["psi.sw"]) psi.store,
                           pi.prop  = if(sw["pi.sw"])  pi.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           Q.store  = Q.store,
                           time     = init.time)
    attr(returns, "Q.big")  <- Q.large
    return(returns)
  }