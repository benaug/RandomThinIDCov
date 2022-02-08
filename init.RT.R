e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.RT=function(data,inits=NA,M=NA){
  library(abind)
  y.ID=data$y.ID
  y.noID=data$y.noID
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- dim(y.ID)[3]
  buff<- data$buff
  n.ID=data$n.ID
  K1D=data$K1D
  
  #data checks
  if(length(dim(y.ID))!=3){
    stop("dim(y.ID) must be 3. Reduced to 2 during initialization")
  }
  if(length(dim(y.noID))!=3){
    stop("dim(y.noID) must be 3. Reduced to 2 during initialization")
  }
  
  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0=inits$lam0
  sigma<- inits$sigma
  
  #initialize IDs
  #This initialization algorithm matches samples with consistent G.noID to same ID if they are caught
  #at same trap. Requires larger data augmentation than necessary for sampling. Improvement to implement:
  #combine if consistent G.noID and caught at traps within a distance consistent with initial sigma.
  n.samples=nrow(y.noID)
  ID=rep(NA,n.samples)
  nextID=n.ID+1
  y.noID2D=apply(y.noID,c(1,2),sum)
  this.j=apply(y.noID2D,1,function(x){which(x>0)})
  y.ID2D=apply(y.ID,c(1,2),sum)
  y.ID2D=rbind(y.ID2D,matrix(0,nrow=M-n.ID,ncol=J))
  y.true2D=y.ID2D
  for(l in 1:n.samples){
    sametrap=y.true2D[,this.j[l]]>0
    if(any(sametrap)){
      ID[l]=which(sametrap)[1]
      y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
    }else{#must be new ID
      ID[l]=nextID
      y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
      nextID=nextID+1
    }
    if(nextID>M)stop("Need to raise M to initialize data.")
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
  lamd <-lam0 * exp(-D * D/(2 * sigma * sigma))
  ll.y=dpois(y.true2D,K1D*lamd*z,log=TRUE)
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K1D (if supplied by user) or problem initializing data.")
  
  return(list(s=s,z=z,ID=ID,y.ID=y.ID2D,y.true=y.true2D,K1D=K1D,
         n.samples=n.samples,this.j=this.j))

}