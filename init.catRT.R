e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.catRT=function(data,inits=NA,M=NA,obstype="poisson"){
  library(abind)
  y.ID=data$y.ID
  this.j=data$this.j
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- dim(y.ID)[3]
  n.cat=data$IDlist$n.cat
  IDcovs=data$IDlist$IDcovs
  buff<- data$buff
  G.ID=data$G.ID
  G.noID=data$G.noID
  n.ID=data$n.ID
  
  n.levels=unlist(lapply(IDcovs,length))
  K1D=data$K1D
  if(!is.matrix(G.ID)){
    G.ID=matrix(G.ID)
  }
  if(!is.matrix(G.noID)){
    G.noID=matrix(G.noID)
  }
  if(!is.list(IDcovs)){
    stop("IDcovs must be a list")
  }
  if(ncol(G.ID)!=n.cat){
    stop("G.ID needs n.cat number of columns")
  }
  if(ncol(G.noID)!=n.cat){
    stop("G.noID needs n.cat number of columns")
  }
  
  #data checks
  if(length(dim(y.ID))!=3){
    stop("dim(y.ID) must be 3. Reduced to 2 during initialization")
  }
 
  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0=inits$lam0
  sigma<- inits$sigma
  gamma=inits$gamma
  if(!is.list(gamma)){
    stop("inits$gamma must be a list")
  }
  
  #initialize IDs
  #This initialization algorithm matches samples with consistent G.noID to same ID if they are caught
  #at same trap. Requires larger data augmentation than necessary for sampling. Improvement to implement:
  #combine if consistent G.noID and caught at traps within a distance consistent with initial sigma.
  G.true=matrix(0,nrow=M,ncol=n.cat)
  G.true[1:n.ID,]=G.ID
  n.samples=length(this.j)
  ID=rep(NA,n.samples)
  nextID=n.ID+1
  y.ID2D=apply(y.ID,c(1,2),sum)
  y.ID2D=rbind(y.ID2D,matrix(0,nrow=M-n.ID,ncol=J))
  y.true2D=y.ID2D
  
  for(l in 1:n.samples){
    # if(l==11)stop("A")
    #can you match an unmarked guy in same trap already assigned an ID?
    obsidx=which(G.noID[l,]!=0)
    matches=rep(FALSE,M)
    for(i2 in 1:M){#can match marked or unmarked
      obsidx2=which(G.true[i2,]!=0)
      sameobsidx=intersect(obsidx,obsidx2)
      matches[i2]=all(G.true[i2,sameobsidx]==G.noID[l,sameobsidx])
    }
    matches=which(matches)
    if(length(matches)==0){#must be new ID
      ID[l]=nextID
      G.true[ID[l],obsidx]=G.noID[l,obsidx]
      y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
      nextID=nextID+1
    }else if(length(matches)==1){
      if(y.true2D[matches,this.j[l]]>0){#caught at same trap?
        ID[l]=matches
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        #new sample for this ID might fill in some missing G.true indices
        notobsidx=which(G.true[ID[l],]==0)
        G.true[ID[l],notobsidx]=G.noID[l,notobsidx]
      }else{#must be new ID
        ID[l]=nextID
        G.true[ID[l],obsidx]=G.noID[l,obsidx]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        nextID=nextID+1
      }
    }else{
      sametrap=y.true2D[matches,this.j[l]]>0
      if(any(sametrap)){
        ID[l]=matches[which(sametrap)[1]]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        #new sample for this ID might fill in some missing G.true indices
        notobsidx=which(G.true[ID[l],]==0)
        G.true[ID[l],notobsidx]=G.noID[l,notobsidx]
      }else{#must be new ID
        ID[l]=nextID
        G.true[ID[l],obsidx]=G.noID[l,obsidx]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        nextID=nextID+1
      }
    }
    if(nextID>M)stop("Need to raise M to initialize data.")
  }
  
  
  #Check for correctness
  for(l in 1:n.samples){
    obsidx=which(G.noID[l,]!=0)
    if(!all(G.true[ID[l],obsidx]==G.noID[l,obsidx]))stop("Error in initialization")
  }
  
  #initialize remainder of G.true
  for(i in 1:M){
    for(c in 1:n.cat){
      if(G.true[i,c]==0){
        G.true[i,c]=sample(IDcovs[[c]],1,prob=inits$gamma[[c]])
      }
    }
  }
  
  #initialize match
  match=matrix(FALSE,nrow=n.samples,ncol=M)
  for(l in 1:n.samples){
    idx=which(G.noID[l,]!=0)
    if(length(idx)>1){#multiple observed
      match[l,]=apply(G.true[,idx],1,function(x){all(x==G.noID[l,idx])})
    }else if(length(idx)==1){#single observed
      match[l,]=G.true[,idx]==G.noID[l,idx]
    }else{#fully latent G.noID
      match[l,]=rep(TRUE,M)
    }
  }
  
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
  lamd <-lam0[1] * exp(-D * D/(2 * sigma[1] * sigma[1]))
  for(i in 1:length(sigma)){
    lamd[which(G.true[,1]==i), ]<-lam0[i]*exp(-D[which(G.true[,1]==i), ]^2/(2*sigma[i]^2))
  }
  
  if(obstype=="poisson"){
    ll.y=y.true2D*0
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,]=dpois(y.true2D[i,],K1D*lamd[i,]*z[i],log=TRUE)
      }
    }
  }else if(obstype=="negbin"){
    theta.d=inits$theta.d
    ll.y=y.true2D*0
    for(i in 1:M){
      if(z[i]==1){
        ll.y[i,]=dnbinom(y.true2D[i,],mu=lamd[i,],size=theta.d*K1D,log=TRUE)
      }
    }
  }else{
    stop("obstype not recognized")
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K1D (if supplied by user) or problem initializing data.")
  
  #Stuff category level probabilities into a ragged matrix. (may have different numbers of levels)
  gammaMat=matrix(0,nrow=n.cat,ncol=max(n.levels))
  for(l in 1:n.cat){
    gammaMat[l,1:n.levels[l]]=gamma[[l]]
  }
  
  #calculate G.latent, which indicies of G currently latent
  G.latent=matrix(1,nrow=M,ncol=n.cat)
  G.latent.ID=0*(data$G.ID>0)
  G.latent[1:n.ID,1:n.cat]=G.latent.ID
  for(l in 1:n.samples){
    for(m in 1:(n.cat)){
      if(G.noID[l,m]!=0){
        G.latent[ID[l],m]=0
      }
    }
  }
  return(list(s=s,z=z,G.true=G.true,ID=ID,G.noID=G.noID,y.ID=y.ID2D,y.true=y.true2D,K1D=K1D,xlim=xlim,ylim=ylim,
         n.samples=n.samples,gammaMat=gammaMat,this.j=this.j,G.latent=G.latent,G.latent.ID=G.latent.ID))

}