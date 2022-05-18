e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.RT<-
  function(N=NA,lam0=NA,sigma=NA,theta.d=NA,K=10,X=NA,buff=NA,obstype="poisson",
           theta.thin=NA,tlocs=0){
    library(abind)
    
    # simulate a population of activity centers
    X=as.matrix(X)
    xlim=c(min(X[,1]),max(X[,1]))+c(-buff,buff)
    ylim=c(min(X[,2]),max(X[,2]))+c(-buff,buff)
    s<- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D<- e2dist(s,X)
    lamd <-lam0 * exp(-D * D/(2 * sigma * sigma))
    J=nrow(X)
    
    # Capture individuals
    y.true <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(lam0))stop("must provide lam0 for bernoulli obstype")
      pd=1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,]=rbinom(K,1,pd[i,j])
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,]=rpois(K,lamd[i,j])
        }
      }
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,]=rnbinom(K,mu=lamd[i,j],size=theta.d)
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

    #summarize captured ind data
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
    

    out<-list(y.true=y.true,y.ID=y.ID,this.j=this.j,this.k=this.k,
              X=X,K=K,buff=buff,s=s,xlim=xlim,ylim=ylim,
              ID=ID,n.ID=n.ID,
              y.true.full=y.true.full,n.cap=n.cap)
    return(out)
  }