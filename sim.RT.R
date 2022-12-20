e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.RT<-
  function(N=NA,lam0=NA,p0=NA,sigma=NA,theta.d=NA,lambda.d=NA,K=NA,X=NA,buff=NA,obstype="poisson",
           theta.thin=NA,K1D=NA,K2D=NA){
    library(abind)
    
    # simulate a population of activity centers
    X=as.matrix(X)
    xlim=c(min(X[,1]),max(X[,1]))+c(-buff,buff)
    ylim=c(min(X[,2]),max(X[,2]))+c(-buff,buff)
    s<- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D<- e2dist(s,X)
    lamd <-lam0 * exp(-D * D/(2 * sigma * sigma))
    J=nrow(X)
    
    #trap operation
    if(!any(is.na(K1D))){
      if(any(K1D>K)){
        stop("Some entries in K1D are greater than K.")
      }
      if(is.null(dim(K1D))){
        if(length(K1D)!=J){
          stop("K1D vector must be of length J.")
        }
      }
    }else{
      K1D=rep(K,J)
    }
    
    if(!any(is.na(K2D))){
      if(nrow(K2D)!=J){
        stop("K2D must be matrix of dimension J x K")
      }
      if(ncol(K2D)!=K){
        stop("K2D must be matrix of dimension J x K")
      }
    }else{
      K2D=matrix(1,J,K)
    }
    
    # Capture individuals
    y.true <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(lam0))stop("must provide lam0 for bernoulli obstype")
      pd=1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,1:K1D[j]]=rbinom(K1D[j],1,pd[i,j])
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,1:K1D[j]]=rpois(K1D[j],lamd[i,j])
        }
      }
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,1:K1D[j]]=rnbinom(K1D[j],mu=lamd[i,j],size=theta.d)
        }
      } 
    }else if(obstype=="negbinHurdle"){
      if(is.na(p0))stop("must provide p0 for negbinHurdle obstype")
      if(is.na(lambda.d))stop("must provide lambda.d for negbinHurdle obstype")
      if(is.na(theta.d))stop("must provide theta.d for negbinHurdle obstype")
      library(VGAM)
      pd<- p0*exp(-D*D/(2*sigma*sigma))
      y.det=array(0,dim=c(N,J,K))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.det[i,j,k]=rbinom(1,1,pd[i,j])
            if(y.det[i,j,k]==1){
              y.true[i,j,k]=rzanegbin(1,munb=lambda.d,size=theta.d,pobs0=0)
            }
          }
        }
      }
    }else{
      stop("obstype not recognized")
    }
    
    y.true.full=y.true
    n.cap=sum(rowSums(y.true.full)>0)
    
    captured=which(rowSums(y.true)>0)
    y.true=y.true[captured,,]
    n.cap=nrow(y.true)
   
    #split sightings into ID'd or not ID'd
    n.samples=sum(y.true)
    ID=rep(NA,n.samples)
    y.ID=y.true
    y.noID=array(0,dim=c(n.samples,J,K))
    idx=1
    for(i in 1:n.cap){
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y.true[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y.true[i,j,k]){ #then samples
              identify=rbinom(1,1,prob=theta.thin)
              if(identify==0){#no
                y.noID[idx,j,k]=1 #add to not identified data
                y.ID[i,j,k]=y.ID[i,j,k]-1 #subtract from identified data
                ID[idx]=i
                idx=idx+1
              }
            }
          }
        }
      }
    }
    n.samples.noID=sum(rowSums(y.noID))
    y.noID=y.noID[1:n.samples.noID,,]
    ID=ID[1:n.samples.noID]
    
    #reassemble to check for correctness
    y.test=y.ID
    for(i in 1:n.samples.noID){
      y.test[ID[i],,]=y.test[ID[i],,]+y.noID[i,,]
    }
    
   if(!all(y.test==y.true))stop("Error reassembling data. Bug in simulator.") #shouldn't happen 

    #reorder ID's if any y.ID inds end up with 0 samples after subsampling.
    #only needed to compare posterior ID match probs to correct ID #
    allIDs=1:nrow(y.ID) #all guys in y.ID with samples before subsampling
    ID.map=allIDs #what we'll use to map old IDs to new IDs
    fix=which(rowSums(y.ID)==0) #guys subsampled out
    if(length(fix)>0){
      ID.map[fix]=(max(allIDs)+1):(max(allIDs)+length(fix)) #give these guys new ID numbers not already used
      #then subtract the removed guys from ID numbers of other guy's samples
      for(i in 1:length(fix)){
        ID.map[ID.map>=fix[i]]=ID.map[ID.map>=fix[i]]-1
        fix=fix-1
      }
      ID=ID.map[ID]
    }
    
    #remove inds subsampled out of y.ID
    IDd=which(rowSums(y.ID)!=0)
    n.ID=length(IDd)
    y.ID=y.ID[IDd,,]
    
    #observed capture data can be represented by site and occasion of each count member
    this.j=this.k=rep(NA,n.samples.noID)
    for(i in 1:n.samples.noID){
      tmp=which(y.noID[i,,]==1,arr.ind=TRUE)
      this.j[i]=tmp[1]
      this.k[i]=tmp[2]
    }
    
    #reorder ID, this.j, this.k by ID (not required)
    ord=order(ID)
    ID=ID[ord]
    this.j=this.j[ord]
    this.k=this.k[ord]
    

    out<-list(y.true=y.true,y.ID=y.ID,this.j=this.j,this.k=this.k,
              X=X,K=K,K1D=K1D,K2D=K2D,buff=buff,s=s,xlim=xlim,ylim=ylim,
              ID=ID,n.ID=n.ID,y.true.full=y.true.full,n.cap=n.cap)
    return(out)
  }