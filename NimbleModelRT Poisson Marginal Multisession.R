NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  lam0.fixed ~ dunif(0,15)
  sigma.fixed ~ dunif(0,10)
  #Expected density for marked + unmarked individuals
  D ~ dunif(0,10) #Expected density
  for(g in 1:N.session){
    #plug in shared df parameter for each session. Must use lam0[g] and sigma[g] here for custom update.
    #alternatively, can be estimated separately or with random effects.
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    lambda[g] <- D*area[g] #expected N
    N[g] ~ dpois(lambda[g]) #realized N
    #thinning prior, not shared across sessions here
    theta.thin[g] ~ dunif(0,1)
    #likelihoods (except for s priors)
    for(i in 1:M[g]) {
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      y.ID[g,i,1:J[g]] ~  dPoissonVector(lam[g,i,1:J[g]]*K1D[g,1:J[g]]*theta.thin[g],z=z[g,i]) #identified detections
    }
    #this speeds up lam0/sigma updates a bit by skipping z=0 inds
    bigLam[g,1:J[g]] <- GetbigLam(lam=lam[g,1:M[g],1:J[g]],z=z[g,1:M[g]])
    lam.noID[g,1:J[g]] <- bigLam[g,1:J[g]]*K1D[g,1:J[g]]*(1-theta.thin[g])
    y.noID[g,1:J[g]] ~ dPoissonVector(lam.noID[g,1:J[g]],z=1) #plug in z=1 to reuse dPoissonVector
  }
})# custom Metropolis-Hastings update for N/z