e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

getArea <- function(X=X,buff=buff){
  N.session=length(X)
  area=rep(NA,N.session)
  for(a in 1:N.session){
    xlim=c(min(X[[a]][,1]),max(X[[a]][,1]))+c(-buff[[a]],buff[[a]])
    ylim=c(min(X[[a]][,2]),max(X[[a]][,2]))+c(-buff[[a]],buff[[a]])
    area[a]=diff(xlim)*diff(ylim)
  }
  return(area)
}

sim.RT.multisession <-
  function(N.session=NA,lambda=NA,lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,buff=NA,
           p0=NA,lambda.d=NA,theta.thin=NA,K1D=NA,K2D=NA,obstype="poisson"){
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
      theta.d=rep(theta.d,N.session)
    }
    if(obstype=="negbinHurdle"){
      if(length(p0)!=N.session)stop("p0 must be of length N.session")
      if(length(lambda.d)!=N.session)stop("lambda.d must be of length N.session")
    }else{
      #make a dummy to pass to sim.RT
      lambda.d=rep(lambda.d,N.session)
      lam0=rep(lam0,N.session)
    }
    
    #realized N
    N=rpois(N.session,lambda)
    
    library(abind)
    xlim=ylim=matrix(NA,N.session,2)
    s=D=vector("list",N.session)
    J=rep(NA,N.session)
    
    for(g in 1:N.session){
      X[[g]]=as.matrix(X[[g]])
      xlim[g,]=c(min(X[[g]][,1]),max(X[[g]][,1]))+c(-buff[g],buff[g])
      ylim[g,]=c(min(X[[g]][,2]),max(X[[g]][,2]))+c(-buff[g],buff[g])
      s[[g]]<- cbind(runif(N[g], xlim[g,1],xlim[g,2]), runif(N[g],ylim[g,1],ylim[g,2]))
      D[[g]]<- e2dist(s[[g]],X[[g]])
      J[g]=nrow(X[[g]])
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
        K1D=vector("list",N.session)
        for(g in 1:N.session){
          K1D[[g]]=rep(K[g],J[g])
        }
      }
      K2D=vector("list",N.session)
      for(g in 1:N.session){
        K2D[[g]]=NA
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
        K2D=vector("list",N.session)
        for(g in 1:N.session){
          K2D[[g]]=matrix(1,J[g],K[g])
        }
      }
      K1D=vector("list",N.session)
      for(g in 1:N.session){
        K1D[[g]]=NA
      }
    }
    
    #simulate sessions one at a time
    data=vector("list",N.session)
    for(g in 1:N.session){
      data[[g]]=sim.RT(N=N[g],theta.thin=theta.thin[g],
                        lam0=lam0[g],p0=p0[g],lambda.d=lambda.d[g],sigma=sigma[g],K=K[g],X=X[[g]],buff=buff[g],
                        obstype=obstype,theta.d=theta.d[g],K1D=K1D[[g]],K2D=K2D[[g]])
    }
    
    #combine session data
    n.samples=rep(NA,N.session)
    for(g in 1:N.session){
      n.samples[g]=length(data[[g]]$this.j)
    }
    n.samples.max=max(n.samples)
    this.j=this.k=ID=matrix(NA,N.session,n.samples.max)
    y.true=y.ID=y.true.full=s=vector("list",N.session)
    n.cap=n.ID=rep(NA,N.session)
    for(g in 1:N.session){
      this.j[g,1:n.samples[g]]=data[[g]]$this.j
      this.k[g,1:n.samples[g]]=data[[g]]$this.k
      ID[g,1:n.samples[g]]=data[[g]]$ID
      n.cap[g]=data[[g]]$n.cap
      n.ID[g]=data[[g]]$n.ID
      y.true[[g]]=data[[g]]$y.true
      y.true.full[[g]]=data[[g]]$y.true.full
      y.ID[[g]]=data[[g]]$y.ID
      s[[g]]=data[[g]]$s
    }
    
    out<-list(this.j=this.j,this.k=this.k,y.ID=y.ID, #observed data
              n.ID=n.ID,n.cap=n.cap,
              y.true.full=y.true.full,s=s,ID=ID,N=N,#true data
              X=X,K=K,K1D=K1D,K2D=K2D,buff=buff,xlim=xlim,ylim=ylim)
    return(out)
  }