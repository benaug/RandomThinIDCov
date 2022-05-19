e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.RT.multisession=function(data,inits=NA,M=NA,obstype="poisson"){
  N.session=nrow(data$this.j)
  n.samples=rowSums(!is.na(data$this.j))
  init.session=vector("list",N.session)
  
  #split inits by session
  inits.use=vector("list",N.session)
  parms=names(inits)
  for(g in 1:N.session){
    inits.use[[g]]=vector("list",length(parms))
    names(inits.use[[g]])=parms
    for(i in 1:length(parms)){
      inits.use[[g]][[i]]=inits[[i]][g]
    }
  }
  
  #initialize sessions one by one
  for(g in 1:N.session){
    data.use=list(y.ID=data$y.ID[[g]],this.j=data$this.j[g,1:n.samples[g]],this.k=data$this.k[g,1:n.samples[g]],
                  n.ID=data$n.ID[g],X=data$X[[g]],buff=data$buff[g],
                  K1D=data$K1D[[g]])
    init.session[[g]]=init.RT(data.use,inits.use[[g]],M=M[g],obstype=obstype)
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
  
  return(list(s=s,z=z,ID=ID,y.ID=y.ID,y.true=y.true,K1D=K1D,J=J,X=X.new,
              n.samples=n.samples,this.j=data$this.j,
              xlim=data$xlim,ylim=data$ylim))
  
}