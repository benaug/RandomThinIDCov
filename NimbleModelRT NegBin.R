NimModel <- nimbleCode({
  #detection function priors
  lam0 ~ dunif(0,20)
  sigma ~ dunif(0,20)
  theta.d ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  #data augmentation prior
  psi ~ dunif(0,1)
  #thinning prior
  theta.thin ~ dunif(0,1)

  #likelihoods (except for s priors)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    p[i,1:J] <- theta.d/(theta.d+lam[i,1:J])
    y.true[i,1:J] ~ dNBVector(p=p[i,1:J],theta.d=theta.d*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
    y.ID[i,1:J] ~ dBinomialVector(theta.thin, y.true[i,1:J],capcounts=capcounts[i])  # Model for ID process
  }
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J]) #intermediate object to derive n
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples]) #number of captured individuals
  N <- sum(z[1:M])
})# end model
