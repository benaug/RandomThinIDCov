e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.catRT.multisession=function(data,inits=NA,M=NA,obstype="poisson"){
  N.session=nrow(data$this.j)
  n.samples=rowSums(!is.na(data$this.j))
  init.session=vector("list",N.session)
  
  #split inits by session
  inits.use=vector("list",N.session)
  parms=names(inits)
  for(g in 1:N.session){
    inits.use[[g]]=vector("list",length(parms))
    names(inits.use[[g]])=parms
    inits.use[[g]]$lam0=inits$lam0[g,]
    inits.use[[g]]$sigma=inits$sigma[g,]
    if(obstype=="negbin"){
      inits.use[[g]]$theta.d=inits$theta.d[g]
    }
  }
  for(g in 1:N.session){
    #append gamma inits
    gamma=vector("list",data$IDlist[[g]]$n.cat) #population frequencies of each category level. Assume equal here.
    n.levels=unlist(lapply(data$IDlist[[g]]$IDcovs,length)) #number of levels per IDcat
    for(i in 1:data$IDlist[[g]]$n.cat){
      gamma[[i]]=rep(1/n.levels[i],n.levels[i])
    }
    inits.use[[g]]$gamma=gamma
  }
  
  #initialize sessions one by one
  for(g in 1:N.session){
    data.use=list(y.ID=data$y.ID[[g]],this.j=data$this.j[g,1:n.samples[g]],this.k=data$this.k[g,1:n.samples[g]],
                  n.ID=data$n.ID[g],X=data$X[[g]],buff=data$buff[g],K1D=data$K1D[[g]],
                  G.ID=data$G.ID[g,1:data$n.ID[g],],G.noID=data$G.noID[g,1:n.samples[g],],
                  IDlist=data$IDlist[[g]])
    init.session[[g]]=init.catRT(data.use,inits.use[[g]],M=M[g],obstype=obstype)
  }
  J=unlist(lapply(data$X,nrow))
  maxM=max(M)
  s=array(NA,dim=c(N.session,maxM,2))
  z=matrix(NA,N.session,maxM)
  ID=matrix(NA,N.session,max(n.samples))
  y.ID=array(NA,dim=c(N.session,maxM,max(J)))
  y.true=array(NA,dim=c(N.session,maxM,max(J)))
  K1D=matrix(NA,N.session,max(J))
  
  for(g in 1:N.session){
    s[g,1:M[g],]=init.session[[g]]$s
    z[g,1:M[g]]=init.session[[g]]$z
    ID[g,1:n.samples[g]]=init.session[[g]]$ID
    y.ID[g,1:M[g],1:J[g]]=init.session[[g]]$y.ID
    y.true[g,1:M[g],1:J[g]]=init.session[[g]]$y.true
    K1D[g,1:J[g]]=init.session[[g]]$K1D
  }
  
  #put X in ragged array
  X.new=array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],]=data$X[[g]]
  }
  
  #IDcov stuff
  n.cats=rep(NA,N.session)
  max.n.levels=rep(NA,N.session)
  for(g in 1:N.session){
    n.cats[g]=data$IDlist[[g]]$n.cat
    max.n.levels[g]=max(unlist(lapply(data$IDlist[[g]]$IDcovs,max)))
  }
  max.n.cat=max(n.cats)
  n.levels=matrix(NA,N.session,max.n.cat)
  for(g in 1:N.session){
    n.levels[g,]=unlist(lapply(data$IDlist[[g]]$IDcovs,max))
  }
  max.n.ID=max(data$n.ID)
  max.n.levels.all=max(n.levels)
  G.true=G.latent=array(NA,dim=c(N.session,maxM,max.n.cat))
  G.latent.ID=array(NA,dim=c(N.session,max.n.ID,max.n.cat))
  G.noID=array(NA,dim=c(N.session,max(n.samples),max.n.cat))
  gammaMat=array(NA,dim=c(N.session,max.n.cat,max.n.levels.all))
  for(g in 1:N.session){
    G.true[g,1:M[g],1:n.cats[g]]=init.session[[g]]$G.true
    G.latent[g,1:M[g],1:n.cats[g]]=init.session[[g]]$G.latent
    G.latent.ID[g,1:data$n.ID[g],1:n.cats[g]]=init.session[[g]]$G.latent.ID
    G.noID[g,1:n.samples[g],1:n.cats[g]]=init.session[[g]]$G.noID
    gammaMat[g,1:n.cats[g],1:max.n.levels[g]]=init.session[[g]]$gammaMat
  }
  
  return(list(s=s,z=z,ID=ID,y.ID=y.ID,y.true=y.true,K1D=K1D,J=J,X=X.new,
              n.samples=n.samples,this.j=data$this.j,xlim=data$xlim,ylim=data$ylim,
              G.true=G.true,G.noID=G.noID,G.latent=G.latent,G.latent.ID=G.latent.ID,
              gammaMat=gammaMat,n.levels=n.levels,n.cat=n.cats))
  
}