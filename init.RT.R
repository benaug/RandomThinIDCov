e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.RT=function(data,inits=NA,M=NA,obstype="poisson"){
  library(abind)
  y.ID=data$y.ID
  this.j=data$this.j
  this.k=data$this.k
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- dim(y.ID)[3]
  buff<- data$buff
  n.ID=data$n.ID
  K1D=data$K1D
  K2D=data$K2D
  
  #data checks
  if(length(dim(y.ID))!=3){
    stop("dim(y.ID) must be 3. Reduced to 2 during initialization")
  }

  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0=inits$lam0
  p0=inits$p0
  sigma=inits$sigma
  lambda=inits$lambda
  theta.d=inits$theta.d
  
  #initialize IDs
  #This initialization algorithm matches samples if they are caught
  #at same trap. Requires larger data augmentation than necessary for sampling.
  n.samples=length(this.j)
  ID=rep(NA,n.samples)
  nextID=n.ID+1
  library(abind)
  y.ID=abind(y.ID,array(0,dim=c(M-n.ID,J,K)),along=1)
  y.true=y.ID
  for(l in 1:n.samples){
    sametrap=rowSums(y.true[,this.j[l],])>0
    if(any(sametrap)){
      ID[l]=which(sametrap)[1]
      y.true[ID[l],this.j[l],this.k[l]]=y.true[ID[l],this.j[l],this.k[l]]+1
    }else{#must be new ID
      ID[l]=nextID
      y.true[ID[l],this.j[l],this.k[l]]=y.true[ID[l],this.j[l],this.k[l]]+1
      nextID=nextID+1
    }
    if(nextID>M)stop("Need to raise M to initialize data.")
  }
  y.true2D=apply(y.true,c(1,2),sum)
  y.ID2D=apply(y.ID,c(1,2),sum)
  #initialize z
  z=1*(rowSums(y.true2D)>0)

  #intialize s
  s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y.true2D)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  D=e2dist(s, X)
  
  if(obstype=="poisson"){
    lamd <-lam0 * exp(-D * D/(2 * sigma * sigma))
    ll.y=y.true2D*0
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,]=dpois(y.true2D[i,],K1D*lamd[i,]*z[i],log=TRUE)
      }
    }
  }else if(obstype=="negbin"){
    lamd <-lam0 * exp(-D * D/(2 * sigma * sigma))
    ll.y=y.true2D*0
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,]=dnbinom(y.true2D[i,],mu=lamd[i,],size=theta.d*K1D,log=TRUE)
      }
    }
  }else if(obstype=="negbinHurdle"){
    pd <-p0 * exp(-D * D/(2 * sigma * sigma))
    ll.y=y.true*0
    if(is.null(data$K2D))stop("Must supply K2D for negbinHurdle model. J x K operation taking values 0 or 1.")
    K2D=data$K2D
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          for(k in 1:K){
            if(K2D[j,k]==1){
              if(y.true[i,j,k]==0){
                ll.y[i,j,k]=log(1-pd[i,j])
              }else{
                ll.y[i,j,k]=log(pd[i,j]) + log(dpois(y.true[i,j,k],lambda=lambda)/(1-exp(-lambda)))
              }
            }else{
              if(y.true[i,j,k]>0)stop("y.true[i,j,k]>0 but K2D[j,k]==0")
            }
          }
        }
      }
    }
  }else{
    stop("obstype not recognized")
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K1D (if supplied by user) or problem initializing data.")
  
  return(list(s=s,z=z,ID=ID,y.ID=y.ID2D,y.true=y.true2D,
              y.ID3D=y.ID,y.true3D=y.true,K1D=K1D,K2D=K2D,
         n.samples=n.samples,this.j=this.j,this.k=this.k,xlim=xlim,ylim=ylim))

}