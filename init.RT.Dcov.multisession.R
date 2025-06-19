init.RT.Dcov.multisession <- function(data,inits=NA,M=NA,obstype="poisson"){
  N.session <- length(data)
  if(length(M)!=N.session)stop("Must supply an M for each session.")
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.RT.Dcov(data[[g]],inits.use,M=M[g],obstype="poisson")
  }
  
  n.cap <- unlist(lapply(data,function(x){x$n.cap}))
  n.samples <- unlist(lapply(data,function(x){length(x$this.j)}))
  J <- unlist(lapply(data,function(x){nrow(x$X)}))
  maxM <- max(M)
  s <- array(NA,dim=c(N.session,maxM,2))
  z <- matrix(NA,N.session,maxM)
  ID <- matrix(NA,N.session,max(n.samples))
  y.true <- array(NA,dim=c(N.session,maxM,max(J)))
  y.ID <- array(NA,dim=c(N.session,maxM,max(J)))
  K1D <- matrix(NA,N.session,max(J))
  
  n.cells <- unlist(lapply(data,function(x){x$n.cells}))
  n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
  n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
  n.cells.max <- max(n.cells)
  n.cells.x.max <- max(n.cells.x)
  n.cells.y.max <- max(n.cells.y)
  res <- unlist(lapply(data,function(x){x$res}))
  cellArea <- res^2
  xlim <- matrix(NA,N.session,2)
  ylim <- matrix(NA,N.session,2)
  x.vals <- matrix(NA,N.session,n.cells.x.max)
  y.vals <- matrix(NA,N.session,n.cells.y.max)
  dSS <- array(NA,dim=c(N.session,n.cells.max,2))
  InSS <- array(0,dim=c(N.session,n.cells.max))
  D.cov <- array(NA,dim=c(N.session,n.cells.max))
  cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
  
  for(g in 1:N.session){
    s[g,1:M[g],] <- init.session[[g]]$s
    z[g,1:M[g]] <- init.session[[g]]$z
    ID[g,1:n.samples[g]] <- init.session[[g]]$ID
    y.true[g,1:M[g],1:J[g]] <- init.session[[g]]$y.true
    y.ID[g,1:M[g],1:J[g]] <- init.session[[g]]$y.ID
    K1D[g,1:J[g]] <- init.session[[g]]$K1D
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
    K1D[g,1:J[g]] <- data[[g]]$K1D
    x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
    y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
    dSS[g,1:n.cells[g],] <- data[[g]]$dSS
    InSS[g,1:n.cells[g]] <- data[[g]]$InSS
    D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
    cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
  }
  
  #put X in ragged array
  X.new <- array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],] <- data[[g]]$X
  }
  dummy.data <- matrix(0,N.session,maxM) #dummy data not used, doesn't really matter what the values are
  
  z.data <- matrix(NA,N.session,maxM)
  for(g in 1:N.session){
    z.data[g,1:n.cap[g]] <- 1
  }
  
  #format data for marginalized sampler, too
  y.noID <- matrix(NA,N.session,max(J))
  for(g in 1:N.session){
    y.noID[g,1:J[g]] <- tabulate(init.session[[g]]$this.j,J[g]) #number of unidentified counts by trap
  }
  
  return(list(y.ID=y.ID,y.noID=y.noID,
              s.init=s,z.init=z,ID=ID,y.true=y.true,K1D=K1D,J=J,X=X.new,
              n.samples=n.samples,
              xlim=xlim,ylim=ylim,
              res=res,cellArea=cellArea,x.vals=x.vals,xlim=xlim,ylim=ylim,
              y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data,z.data=z.data))
  
}