#This version allows spatial density covariates and/or a habitat mask with multiple sessions.
#Poisson obsmod only.

library(nimble)
library(coda)
source("sim.RT.Dcov.multisession.R")
source("sim.RT.Dcov.R")
source("NimbleModelRT Poisson Dcov Marginal Multisession.R")
source("NimbleFunctionsRT Poisson Dcov Marginal Multisession.R")
source("init.RT.Dcov.multisession.R")
source("init.RT.Dcov.R")
source("sSampler Poisson Dcov Marginal Multisession.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

####Simulate some data####
N.session <- 3
#detection parameters
lam0 <- rep(0.5,N.session)
sigma <- rep(0.5,N.session)

K <- rep(10,N.session) #number of occasions
buff <- rep(2,N.session) #state space buffer
theta.thin <- rep(0.25,N.session) #sample thinning parameter

#make an SCR trapping array. Making the trapping array size vary by session
X <- vector("list",N.session)
X[[1]] <- as.matrix(expand.grid(1:10,1:10))
X[[2]] <- as.matrix(expand.grid(1:9,1:9))
X[[3]] <- as.matrix(expand.grid(1:11,1:11))

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- ylim <- matrix(NA,N.session,2)
for(g in 1:N.session){
  xlim[g,] <- range(X[[g]][,1]) + c(-buff[g],buff[g])
  ylim[g,] <- range(X[[g]][,2]) + c(-buff[g],buff[g])
}

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
for(g in 1:N.session){
  x.shift <- xlim[g,1]
  y.shift <- ylim[g,1]
  xlim[g,] <- xlim[g,] - x.shift
  ylim[g,] <- ylim[g,] - y.shift
  X[[g]][,1] <- X[[g]][,1]- x.shift
  X[[g]][,2] <- X[[g]][,2]- y.shift
}

res <- rep(0.25,N.session) #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- y.vals <- dSS <- cells <- vector("list",N.session)
n.cells <- n.cells.x <- n.cells.y <- rep(NA,N.session)
for(g in 1:N.session){
  x.vals[[g]] <- seq(xlim[g,1]+res[g]/2,xlim[g,2]-res[g]/2,res[g]) #x cell centroids
  y.vals[[g]] <- seq(ylim[g,1]+res[g]/2,ylim[g,2]-res[g]/2,res[g]) #y cell centroids
  dSS[[g]] <- as.matrix(cbind(expand.grid(x.vals[[g]],y.vals[[g]])))
  cells[[g]] <- matrix(1:nrow(dSS[[g]]),nrow=length(x.vals[[g]]),ncol=length(y.vals[[g]]))
  n.cells[g] <- nrow(dSS[[g]])
  n.cells.x[g] <- length(x.vals[[g]])
  n.cells.y[g] <- length(y.vals[[g]])
}

#create a density covariate - one for each session
library(geoR)
D.cov <- vector("list",N.session)
#need a simulated landscape with individuals living around traps to be captured
#these are pretty good
D.seeds <- c(13223,13216,13252)
for(g in 1:N.session){
  set.seed(D.seeds[g])
  D.cov.tmp <- grf(n.cells[g],grid=dSS[[g]],cov.pars=c(500,500),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
  D.cov.tmp <- as.numeric(scale(D.cov.tmp)) #scale
  par(mfrow=c(1,1),ask=FALSE)
  D.cov[[g]] <- D.cov.tmp
  image(x.vals[[g]],y.vals[[g]],matrix(D.cov[[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," D.cov"),xlab="X",ylab="Y",col=cols1)
}

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
InSS <- vector("list",N.session)
for(g in 1:N.session){
  dSS.tmp <- dSS[[g]] - res[g]/2 #convert back to grid locs
  InSS[[g]] <- rep(1,length(D.cov[[g]]))
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  image(x.vals[[g]],y.vals[[g]],matrix(InSS[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Habitat"))
}

#Density covariates
D.beta0 <- rep(-0.5,N.session)
D.beta1 <- rep(0.5,N.session)
#what is implied expected N in state space?
for(g in 1:N.session){
  lambda.cell <- InSS[[g]]*exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  print(sum(lambda.cell)) #expected N in state space
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," Expected Density"),col=cols1)
  points(X[[g]],pch=4,cex=0.75)
}


#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(14353244) #change seed for new data set
data <- sim.RT.Dcov.multisession(N.session=N.session,lam0=lam0,sigma=sigma,
                                 theta.thin=theta.thin,K=K,X=X,obstype="poisson",
                                 D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
                                 res=res,xlim=xlim,ylim=ylim)

#What is the observed data?
g <- 1 #session to look at
str(data[[g]]$y.ID) #the observed ID detections
head(data[[g]]$this.j) #the trap of detection for all unidentified detections
head(data[[g]]$this.k) #occasion of capture, but not used in this 2D data sampler

#Visualize activity centers
for(g in 1:N.session){
  lambda.cell <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Expected Density"))
  points(X[[g]],pch=4,cex=0.75)
  points(data[[g]]$s,pch=16)
}

for(g in 1:N.session){
  #function to test for errors in mask set up. 
  mask.check(dSS=data[[g]]$dSS,cells=data[[g]]$cells,n.cells=data[[g]]$n.cells,n.cells.x=data[[g]]$n.cells.x,
             n.cells.y=data[[g]]$n.cells.y,res=data[[g]]$res,xlim=data[[g]]$xlim,ylim=data[[g]]$ylim,
             x.vals=data[[g]]$x.vals,y.vals=data[[g]]$y.vals)
}

##Fit model in Nimble##
#data augmentation level
M <- rep(175,N.session)

J <- nrow(X) #number of detectors
K1D <- data$K1D #pull out trap operation

inits <- list(lam0=rep(1,N.session),sigma=rep(1,N.session)) #ballpark inits to build data, one per session

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.RT.Dcov.multisession(data,inits,M=M,obstype="poisson")

#inits for nimble
N.init <- rowSums(nimbuild$z.init,na.rm=TRUE) #N and z inits must be consistent
D0.init <- rowSums(nimbuild$z.init,na.rm=TRUE)/(rowSums(nimbuild$InSS)*nimbuild$cellArea)

Niminits <- list(N=N.init,z=nimbuild$z.init,s=nimbuild$s.init,theta.thin=rep(0.5,N.session),
                 lam0.fixed=0.5,sigma.fixed=0.5,D0=D0.init,D.beta1=rep(0,N.session))

#constants for Nimble
J <- unlist(lapply(data,function(x){nrow(x$X)}))
constants <- list(N.session=N.session,M=M,J=J,K1D=nimbuild$K1D,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,D.cov=nimbuild$D.cov,res=nimbuild$res,
                  cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells)

# Supply data to Nimble. marginalized data formatted in nimbuild object
Nimdata <- list(y.ID=nimbuild$y.ID, #ID detections
                y.noID=nimbuild$y.noID, #no ID detections
                X=nimbuild$X,cells=nimbuild$cells,
                dummy.data=nimbuild$dummy.data,InSS=nimbuild$InSS)

parameters <- c('D0','D.beta1','lambda.N','N','lam0.fixed',
                'sigma.fixed','theta.thin')
parameters2 <- c("lambda.cell",'D0') #record D0 here for plotting
nt <- 1 #thinning rate for parameters
nt2 <- 10 #thinning rate for parameters2

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#use block sampler below for 'D0.M','D.beta1'
config.nodes <- c('lam0.fixed','sigma.fixed','theta.thin')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2, thin2=nt2,
                      nodes=config.nodes)

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###*required* sampler replacement
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal.
for(g in 1:N.session){
  #nodes used for update
  y.ID.nodes <- Rmodel$expandNodeNames(paste("y.ID[",g,",1:",M[g],",1:",J[g],"]"))
  y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[",g,",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",1:",M[g],",1:",J[g],"]"))
  bigLam.nodes <- Rmodel$expandNodeNames(paste("bigLam[",g,",1:",J[g],"]")) #only need this in calcNodes
  lam.noID.nodes <- Rmodel$expandNodeNames(paste("lam.noID[",g,",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",1:",M[g],"]"))
  calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.noID.nodes,y.ID.nodes)
  inds.detected <- which(rowSums(nimbuild$y.ID[g,,],na.rm=TRUE)>0)
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(g=g,J=J[g],M=M[g],z.ups=z.ups[g],inds.detected=inds.detected,
                                                   lam.nodes=lam.nodes,lam.noID.nodes=lam.noID.nodes,
                                                   y.ID.nodes=y.ID.nodes,y.noID.nodes=y.noID.nodes,
                                                   N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),silent = TRUE)
}

# activity center sampler
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(g=g,i=i,J=J[g],n.cells=nimbuild$n.cells[g],
                                                   n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
                                                   scale=0.25),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#add lam0 and sigma samplers block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW.
conf$addSampler(target = c("lam0.fixed","sigma.fixed"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
#AF_slice pretty efficient for Dcov parameters. Block by session
for(g in 1:N.session){
  conf$addSampler(target = c(paste("D0[",g,"]"),paste("D.beta1[",g,"]")),
                  type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[-c(1:burnin),]))

exp(D.beta0)
unlist(lapply(data,function(x){x$lambda.N}))#expected Ns
unlist(lapply(data,function(x){x$N})) #realized Ns


cor(mvSamples[-c(1:burnin),])

#plot density surface, etc.
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
n.cells.max <- max(n.cells)
lambda.cell.idx <- matrix(lambda.cell.idx,N.session,n.cells.max)
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10 #consider nt2 thinning rate when setting burnin2

#compare expected D plot to truth
#image will show posterior means
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- array(NA,dim=c(N.session,n.cells.max,length(n.iter.use)))
lambda.cell <- lambda.cell.ests <- array(NA,dim=c(N.session,n.cells.max))
lambda.cell.HPDs <- array(NA,dim=c(N.session,n.cells.max,2))
for(g in 1:N.session){
  lambda.cell[g,1:n.cells[g]] <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  lambda.cell.post[g,1:n.cells[g],] <- t(cellArea[g]*mvSamples2[n.iter.use,D0.idx[g]]*
                                           mvSamples2[n.iter.use,lambda.cell.idx[g,1:n.cells[g]]])
  lambda.cell.ests[g,1:n.cells[g]] <- rowMeans(lambda.cell.post[g,1:n.cells[g],])
  lambda.cell.HPDs[g,1:n.cells[g],] <- HPDinterval(mcmc(t(lambda.cell.post[g,1:n.cells[g],])))
  #remove nonhabitat (or not, comment out)
  lambda.cell[g,InSS[[g]]==0] <- NA
  lambda.cell.ests[g,InSS[[g]]==0] <- NA
}

par(mfrow=c(1,1),ask=FALSE)
for(g in 1:N.session){
  zlim <- range(c(lambda.cell[g,],lambda.cell.ests[g,]),na.rm=TRUE) #use same zlim for plots below
  #truth
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"True Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
  #estimate, posterior means
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell.ests[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"Est Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
}
#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
for(g in 1:N.session){
  idx <- order(lambda.cell[g,1:n.cells[g]])
  plot(lambda.cell.ests[g,1:n.cells[g]][idx]~lambda.cell[g,1:n.cells[g]][idx],type="l",lwd=2,
       main=paste("Session",g,"True vs. Estimated Density"),ylim=range(lambda.cell.HPDs[g,1:n.cells[g],]))
  lines(lambda.cell.HPDs[g,1:n.cells[g],1][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  lines(lambda.cell.HPDs[g,1:n.cells[g],2][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  abline(0,1,col="darkred",lwd=2) #1:1 expectation
}
