NimModel <- nimbleCode({
  #detection function priors
  lam0 ~ dunif(0,20)
  sigma ~ dunif(0,20)
  #thinning prior
  theta.thin ~ dunif(0,1)
  #abundance prior
  lambda.N ~ dunif(0,1000) #Expected N
  N ~ dpois(lambda.N) #realized N in state space
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.ID[i,1:J] ~ dPoissonVector(lam[i,1:J]*K1D[1:J]*theta.thin,z=z[i]) #identified detections
  }
  for(j in 1:J){
    bigLam[j] <- sum(lam[1:M,j]*z[1:M])
    lam.noID[j] <- bigLam[j]*K1D[j]*(1-theta.thin)
    y.noID[j] ~ dpois(lam.noID[j]) #unidentified detections, independent from identified with Poisson obsmod
  }
})# end model