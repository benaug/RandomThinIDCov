NimModel <- nimbleCode({
  #detection function priors
  for(i in 1:n.levels[1]){ #make sure this is correct ID cov
    lam0[i] ~ dunif(0,2)
    sigma[i] ~ dunif(0,2)
  }
  #data augmentation prior
  psi ~ dunif(0,1)
  #thinning priors
  for(i in 1:n.levels[2]){ #make sure this is correct ID cov
    theta.thin[i] ~ dunif(0,1)
  }
  #categorical ID covariate priors
  for(m in 1:n.cat){
    alpha[m,1:n.levels[m]] <- 1 # prior parameters
    gammaMat[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]])
  }
  
  #likelihoods (except for s priors)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    for(m in 1:n.cat){
      G.true[i,m] ~ dcat(gammaMat[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma[G.true[i,1]], lam0=lam0[G.true[i,1]], z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J],z=z[i])  # Model for complete capture histories
    theta.i[i] <- theta.thin[G.true[i,2]]
    y.ID[i,1:J] ~ dBinomialVector(theta.i[i], y.true[i,1:J],capcounts=capcounts[i])  # Model for ID process
  }
  #calculate number of inds captured and abundance
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],G.latent=G.latent[1:M,1:n.cat])
  N <- sum(z[1:M])
})# end model
