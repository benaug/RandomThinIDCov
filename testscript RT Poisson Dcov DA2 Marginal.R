#this version uses an alternative data augmentation approach that runs faster and allows a poisson
#prior on N. Also uses observation model marginalized over individuals with results from Herliansyah et al.
#(2024) https://link.springer.com/article/10.1007/s13253-023-00598-3
#Marginal only possible with Poisson detections (among count models)
#This version allows spatial density covariates and/or a habitat mask.

library(nimble)
library(coda)
source("sim.RT.Dcov.R")
source("NimbleModelRT Poisson Dcov DA2 Marginal.R")
source("NimbleFunctionsRT Poisson Dcov DA2 Marginal.R")
source("init.RT.Dcov.R")
source("sSampler Poisson Dcov Marginal.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#detection parameters
lam0 <- 0.5
sigma <- 0.5

K <- 10 #number of occasions
buff <- 2 #state space buffer
X <- expand.grid(3:11,3:11) #make a trapping array
#sample thinning parameter
theta.thin <- 0.25

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
set.seed(152)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(500,500),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4)

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]>10] <- 0
InSS[dSS.tmp[,1]>10&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]>10&dSS.tmp[,2]>10] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -0.5 #data simulator uses intercept for marked + unmarked
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X,pch=4,cex=1,lwd=2)

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(143532) #change seed for new data set
data <- sim.RT.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,res=res,xlim=xlim,ylim=ylim,
                     lam0=lam0,sigma=sigma,theta.thin=theta.thin,K=K,X=X,obstype="poisson")
points(data$s,pch=16) #add activity centers


#What is the observed data?
str(data$y.ID) #the observed ID detections
head(data$this.j) #the trap of detection for all unidentified detections
head(data$this.k) #occasion of capture, but not used in this 2D data sampler

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

##Fit model in Nimble##
#data augmentation level
M <- 175

J <- nrow(X) #number of detectors
K1D <- data$K1D #pull out trap operation

inits <- list(lam0=1,sigma=1) #ballpark inits to build data

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.RT.Dcov(data,inits,M=M,obstype="poisson")

#inits for nimble
#init for D0
D0.init <- sum(nimbuild$z)/(sum(data$InSS)*data$res^2)
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #z and N inits must be consistent
                 D0=D0.init,D.beta1=0,
                 s=nimbuild$s,sigma=inits$sigma,lam0=inits$lam0,theta.thin=0.5)

#constants for Nimble
J <- nrow(data$X)
constants <- list(M=M,J=J,D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res)

# Supply data to Nimble. Note, y.true is completely latent.
y.noID <- tabulate(nimbuild$this.j,J) #number of unidentified counts by trap
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y.ID=nimbuild$y.ID,X=as.matrix(X),K1D=K1D,y.noID=y.noID,
                dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('lam0','sigma','theta.thin','N','lambda.N',"D0","D.beta1")
nt <- 1 #thinning rate
#record these for plotting density surface
parameters2 <- c("lambda.cell",'D0')
nt2 <- 10 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#use block sampler below for 'D0.M','D.beta1'
config.nodes <- c('lam0','sigma','theta.thin')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2, thin2=nt2,
                      nodes=config.nodes)

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###*required* sampler replacement
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal.
#nodes used for update
y.ID.nodes <- Rmodel$expandNodeNames(paste("y.ID[1:",M,",1:",J,"]"))
y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
bigLam.nodes <- Rmodel$expandNodeNames("bigLam") #only need this in calcNodes
lam.noID.nodes <- Rmodel$expandNodeNames("lam.noID")
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.noID.nodes,y.ID.nodes)
inds.detected <- which(rowSums(nimbuild$y.ID)>0)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,inds.detected=inds.detected,
                                                 lam.nodes=lam.nodes,lam.noID.nodes=lam.noID.nodes,
                                                 y.ID.nodes=y.ID.nodes,y.noID.nodes=y.noID.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)

#must use this activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,J=J,res=data$res,n.cells.x=data$n.cells.x,
                                                 n.cells.y=data$n.cells.y,xlim=nimbuild$xlim,
                                                 ylim=nimbuild$ylim,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#add lam0 and sigma samplers block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW.
conf$addSampler(target = c("lam0","sigma"),
                  type = 'RW_block',
                  control = list(adaptive=TRUE),
                  silent = TRUE)
#AF_slice block update pretty efficient for Dcovs. Using this only, no samplers assigned above.
conf$addSampler(target = c("D0","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
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
burnin <- 500
plot(mcmc(mvSamples[-c(1:burnin),]))

data$N #realized N
data$lambda.N #expected N
exp(D.beta0)

cor(mvSamples[-c(1:burnin),])

#plot density surface, etc.
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10 #should discard more burnin

#compare expected D plot to truth (for simulated data sets)
n.cells <- data$n.cells
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- t(cellArea*mvSamples2[n.iter.use,D0.idx]*mvSamples2[n.iter.use,lambda.cell.idx[1:n.cells]])
lambda.cell.ests <- rowMeans(lambda.cell.post[1:n.cells,])
lambda.cell.HPDs <- HPDinterval(mcmc(t(lambda.cell.post[1:n.cells,])))
#remove nonhabitat (or not, comment out)
lambda.cell[data$InSS==0] <- NA
lambda.cell.ests[data$InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)

#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
idx <- order(lambda.cell)
plot(lambda.cell.ests[1:n.cells][idx]~lambda.cell[1:n.cells][idx],type="l",lwd=2,
     main="True vs. Estimated Density",ylim=range(lambda.cell.HPDs[1:n.cells,]))
lines(lambda.cell.HPDs[1:n.cells,1][idx]~lambda.cell[1:n.cells][idx],lty=2)
lines(lambda.cell.HPDs[1:n.cells,2][idx]~lambda.cell[1:n.cells][idx],lty=2)
abline(0,1,col="darkred",lwd=2) #1:1 expectation
