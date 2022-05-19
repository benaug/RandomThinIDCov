NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  lam0.fixed~dunif(0,15)
  sigma.fixed~dunif(0,10)
  theta.d~dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
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
    theta.thin[g]~dunif(0,1)
    #likelihoods (except for s priors)
    for(i in 1:M[g]) {
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      p[g,i,1:J[g]] <- theta.d/(theta.d+lam[g,i,1:J[g]])
      y.true[g,i,1:J[g]] ~ dNBVector(p=p[g,i,1:J[g]],theta.d=theta.d*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
      y.ID[g,i,1:J[g]] ~ dBinomialVector(theta.thin[g], y.true[g,i,1:J[g]],capcounts=capcounts[g,i])  # Model for ID process
    }
    capcounts[g,1:M[g]] <- Getcapcounts(y.true=y.true[g,1:M[g],1:J[g]]) #intermediate object to derive n
    n[g] <- Getncap(capcounts=capcounts[g,1:M[g]],ID=ID[g,1:n.samples[g]]) #number of captured individuals
  }
})# end model
