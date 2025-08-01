NimModel <- nimbleCode({
  #detection function priors
  lam0 ~ dunif(0,20)
  sigma ~ dunif(0,20)
  #data augmentation prior
  psi ~ dunif(0,1)
  #thinning prior
  theta.thin ~ dunif(0,1)

  #likelihoods (except for s priors)
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J],z=z[i])  # Model for complete capture histories
    y.ID[i,1:J] ~ dBinomialVector(theta.thin, y.true[i,1:J],capcounts=capcounts[i])  # Model for ID process
  }
  #calculate number of inds captured
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M,capcounts.ID=capcounts.ID[1:M]) #intermediate object
  n <- Getncap(capcounts=capcounts[1:M])
  N <- sum(z[1:M])
})# end model
