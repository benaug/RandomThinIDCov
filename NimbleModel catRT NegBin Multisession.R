NimModel <- nimbleCode({
  #detection function priors as function of first ID cov
  #shared across sessions here
  for(df in 1:n.levels[1,1]){
    lam0.fixed[df] ~ dunif(0,15)
    sigma.fixed[df] ~ dunif(0,10)
  }
  theta.d ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  #Expected density for marked + unmarked individuals
  D ~ dunif(0,10) #Expected density
  for(g in 1:N.session){
    for(df in 1:n.levels[1,1]){ #assuming same levels across sessions so they can be shared!
      lam0[g,df] <- lam0.fixed[df]
      sigma[g,df] <- sigma.fixed[df]
    }
    lambda[g] <- D*area[g] #expected N
    N[g] ~ dpois(lambda[g]) #realized N
    #thinning priors
    for(i in 1:n.levels[g,2]){ #make sure this is correct ID cov
      theta.thin[g,i] ~ dunif(0,1)
    }
    #categorical ID covariate priors
    for(m in 1:n.cat[g]){ #skip first cat for marked/unmarked
      for(l in 1:n.levels[g,m]){
        alpha[g,m,l] <- 1 # prior parameters
      }
      gammaMat[g,m,1:n.levels[g,m]] ~ ddirch(alpha[g,m,1:n.levels[g,m]])
    }
    #likelihoods (except for s priors)
    for(i in 1:M[g]) {
      for(m in 1:n.cat[g]){
        G.true[g,i,m] ~ dcat(gammaMat[g,m,1:n.levels[g,m]])
      }
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g,G.true[g,i,1]], lam0=lam0[g,G.true[g,i,1]], z=z[g,i])
      p[g,i,1:J[g]] <- theta.d/(theta.d+lam[g,i,1:J[g]])
      y.true[g,i,1:J[g]] ~ dNBVector(p=p[g,i,1:J[g]],theta.d=theta.d*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
      theta.thin.i[g,i] <- theta.thin[g,G.true[g,i,2]] #make sure this is correct ID cov
      y.ID[g,i,1:J[g]] ~ dBinomialVector(theta.thin.i[g,i], y.true[g,i,1:J[g]],capcounts=capcounts[g,i])  # Model for ID process
    }
    capcounts[g,1:M[g]] <- Getcapcounts(y.true=y.true[g,1:M[g],1:J[g]]) #intermediate object to derive n
    n[g] <- Getncap(capcounts=capcounts[g,1:M[g]],ID=ID[g,1:n.samples[g]],G.latent=G.latent[g,1:M[g],1:n.cat[g]]) #number of captured individuals
  }
})# end model
