sim.RT.Dcov.multisession <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,xlim=NA,ylim=NA,
           p0=NA,lambda.d=NA,theta.thin=NA,K1D=NA,K2D=NA,obstype="poisson"){
    if(length(D.beta0)!=N.session)stop("D.beta0 must be of length N.session")
    if(length(D.beta1)!=N.session)stop("D.beta1 must be of length N.session")
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K)!=N.session)stop("K must be of length N.session")
    if(length(X)!=N.session)stop("X must be of length N.session")
    if(length(buff)!=N.session)stop("buff must be of length N.session")
    if(length(theta.thin)!=N.session)stop("theta.thin must be of length N.session")
    if(obstype%in%c("negbin","poisson")){
      if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    }
    if(obstype=="negbin"){
      if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
    }else{
      #make a dummy to pass to sim.RT
      theta.d <- rep(theta.d,N.session)
    }
    if(obstype=="negbinHurdle"){
      if(length(p0)!=N.session)stop("p0 must be of length N.session")
      if(length(lambda.d)!=N.session)stop("lambda.d must be of length N.session")
    }else{
      #make a dummy to pass to sim.RT
      lambda.d <- rep(lambda.d,N.session)
    }
    
    J <- rep(NA,N.session)
    for(g in 1:N.session){
      X[[g]] <- as.matrix(X[[g]])
      J[g] <- nrow(X[[g]])
    }
    
    #trap operation
    if(obstype!="negbinHurdle"){
      if(!any(is.na(K1D))){
        for(g in 1:N.session){
          if(any(K1D[[g]]>K[g])){
            stop("Some entries in K1D[[g]] are greater than K[g].")
          }
          if(is.null(dim(K1D[[g]]))){
            if(length(K1D[[g]])!=J[g]){
              stop("K1D[[g]] vector must be of length J[g].")
            }
          }
        }
      }else{
        K1D <- vector("list",N.session)
        for(g in 1:N.session){
          K1D[[g]] <- rep(K[g],J[g])
        }
      }
      K2D <- vector("list",N.session)
      for(g in 1:N.session){
        K2D[[g]] <- NA
      }
    }else{
      if(!any(is.na(K2D))){
        for(g in 1:N.session){
          if(nrow(K2D[[g]])!=J[g]){
            stop("K2D[[g]] must be a J[g] x K[g] matrix")
          }
          if(ncol(K2D[[g]])!=K[g]){
            stop("K2D[[g]] must be a J[g] x K[g] matrix")
          }
        }
      }else{
        K2D <- vector("list",N.session)
        for(g in 1:N.session){
          K2D[[g]] <- matrix(1,J[g],K[g])
        }
      }
      K1D <- vector("list",N.session)
      for(g in 1:N.session){
        K1D[[g]] <- NA
      }
    }
    
    #simulate sessions one at a time
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.RT.Dcov(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],InSS=InSS[[g]],res=res[g],
                               theta.thin=theta.thin[g],lam0=lam0[g],p0=p0[g],lambda.d=lambda.d[g],sigma=sigma[g],
                               K=K[g],X=X[[g]],obstype=obstype,theta.d=theta.d[g],xlim=xlim[g,],ylim=ylim[g,],
                               K1D=K1D[[g]],K2D=K2D[[g]])
    }
    return(data)
  }