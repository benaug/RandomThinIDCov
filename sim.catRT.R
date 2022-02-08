e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.catRT<-
  function(N=NA,lam0=NA,sigma=NA,theta.d=NA,K=10,X=NA,buff=NA,n.cat=NA,
           theta.cat=NA,IDcovs=NA,gamma=NA,theta.thin=NA,obstype="poisson",
           G.thin=NA){
    if(length(lam0)!=length(sigma))stop("Need the same number of values for each detection function parameter")
    if(length(IDcovs[[1]])!=length(lam0))stop("IDcovs[[1]] should have the same number of values as detection function parameter values")
    if(length(IDcovs[[G.thin]])!=length(theta.thin))stop("IDcovs element G.thin should have the same number of values as the number of thinning parameters")
    library(abind)
    
    #simulate IDcovs
    G.true=matrix(NA,nrow=N,ncol=n.cat) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:n.cat){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    
    # simulate a population of activity centers
    X=as.matrix(X)
    xlim=c(min(X[,1]),max(X[,1]))+c(-buff,buff)
    ylim=c(min(X[,2]),max(X[,2]))+c(-buff,buff)
    s<- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D<- e2dist(s,X)
    lamd <-lam0[1] * exp(-D * D/(2 * sigma[1] * sigma[1]))
    for(i in 2:length(sigma)){
      lamd[which(G.true[,1]==i), ]<-lam0[i]*exp(-D[which(G.true[,1]==i), ]^2/(2*sigma[i]^2))
    }
    J=nrow(X)
    
    # Capture individuals
    y.true <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,]=rbinom(K,1,pd[i,j])
        }
      }
    }else if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,]=rpois(K,lamd[i,j])
        }
      }
    }else if(obstype=="negbin"){
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
    G.full=G.true[captured,]
    n.cap=nrow(y.true)
   
    #split sightings into ID'd or not ID'd
    n.samples=sum(y.true)
    ID=rep(NA,n.samples)
    y.ID=y.true
    y.noID=array(0,dim=c(n.samples,J,K))
    G.noID=matrix(NA,nrow=n.samples,n.cat)
    idx=1
    for(i in 1:n.cap){
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y.true[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y.true[i,j,k]){ #then samples
              identify=rbinom(1,1,prob=theta.thin[G.full[i,G.thin]])
              if(identify==0){#no
                y.noID[idx,j,k]=1 #add to not identified data
                y.ID[i,j,k]=y.ID[i,j,k]-1 #subtract from identified data
                G.noID[idx,]=G.full[i,] #extract ID covs
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
    G.noID=G.noID[1:n.samples.noID,]
    ID=ID[1:n.samples.noID]
    
    #reassemble to check for correctness
    y.test=y.ID
    for(i in 1:n.samples.noID){
      y.test[ID[i],,]=y.test[ID[i],,]+y.noID[i,,]
      if(!all(G.noID[i,]==G.full[ID[i],]))stop("Error reassembling data. Bug in simulator") #shouldn't happen
    }
    
   if(!all(y.test==y.true))stop("Error reassembling data. Bug in simulator.") #shouldn't happen 

    #summarize captured ind data
    IDd=which(rowSums(y.ID)!=0)
    n.ID=length(IDd)
    y.ID=y.ID[IDd,,]
    G.ID=G.full[IDd,]
    
    #Observation failure for category levels in no ID. missing at random
    for(l in 1:n.cat){
      drop=which(rbinom(nrow(G.noID),1,theta.cat[l])==0)
      if(length(drop)>0){
        G.noID[drop,l]=0 #0 is dropout
      }
    }
    
    #observed capture data can be represented by site and occasion of each count member
    this.j=this.k=rep(NA,n.samples.noID)
    for(i in 1:n.samples.noID){
      tmp=which(y.noID[i,,]==1,arr.ind=TRUE)
      this.j[i]=tmp[1]
      this.k[i]=tmp[2]
    }
    

    out<-list(y.true=y.true,y.ID=y.ID,this.j=this.j,this.k=this.k,
              G.full=G.full,G.ID=G.ID,G.noID=G.noID,
              IDlist=list(n.cat=n.cat,IDcovs=IDcovs),
              X=X,K=K,buff=buff,s=s,xlim=xlim,ylim=ylim,
              ID=ID,n.ID=n.ID,
              y.true.full=y.true.full,n.cap=n.cap)
    return(out)
  }