library(nimble)
library(coda)
source("sim.RT.multisession.R")
source("sim.RT.R")
source("NimbleModelRT Poisson Marginal Multisession.R")
source("NimbleFunctionsRT Poisson Marginal Multisession.R")
source("init.RT.multisession.R")
source("init.RT.R")
source("sSampler Poisson Marginal Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session <- 3
D <- rep(0.4,N.session) #expected density in units of sigma and X
lam0 <- rep(0.5,N.session)
sigma <- rep(0.5,N.session)
K <- c(5,6,7) #number of occasions
buff <- rep(2,N.session) #state space buffer

#make trapping arrays
X1 <- expand.grid(3:11,3:11)
X2 <- expand.grid(3:12,3:12)
X3 <- expand.grid(3:13,3:13)
X <- list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area <- getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda <- D*area
lambda #expected N in each session

#sample thinning parameter
theta.thin <- rep(0.25,N.session)

data <- sim.RT.multisession(N.session=N.session,lambda=lambda,
                         lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,theta.thin=theta.thin,obstype="poisson")

#What is the observed data?
str(data$y.ID) #the observed ID detections. One list element per session
head(t(data$this.j)) #the trap of detection for all unidentified detections.
head(t(data$this.k)) #occasion of capture, but not used in this 2D data sampler


##Fit model in Nimble##

#data augmentation level
M <- c(155,165,175)
maxM <- max(M)

#trap operation matrix - plugging in perfect operation here
J <- unlist(lapply(X,nrow))
maxJ <- max(J)
K1D <- matrix(NA,nrow=N.session,ncol=maxJ)
for(g in 1:N.session){
  K1D[g,1:J[g]] <- rep(K[g],J[g])
}

inits <- list(lam0=rep(0.5,N.session),sigma=rep(1,N.session)) #ballpark inits to build data, one per session. Using simulated values here, don't do this in practice.

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.RT.multisession(data,inits,M=M,obstype="poisson")

#inits for nimble
N.init <- rowSums(nimbuild$z,na.rm=TRUE)

Niminits <- list(N=N.init,z=nimbuild$z,D=mean(N.init/area), #initialize D from N.init for faster convergence
                 s=nimbuild$s,lam0.fixed=0.5,sigma.fixed=0.5,theta.thin=rep(0.5,N.session))

#constants for Nimble
constants <- list(N.session=N.session,M=M,J=J,K1D=K1D,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,area=area)

y.noID <- matrix(NA,N.session,maxJ)
for(g in 1:N.session){
  y.noID[g,1:J[g]] <- tabulate(nimbuild$this.j[g,1:nimbuild$n.samples[g]],J[g]) #number of unidentified counts by trap
}

# Supply data to Nimble.
z.data <- matrix(NA,N.session,max(M))
for(g in 1:N.session){
  z.data[g,1:data$n.ID[g]] <- 1
}
Nimdata <- list(z=z.data,X=nimbuild$X,y.ID=nimbuild$y.ID,y.noID=y.noID)


# set parameters to monitor
parameters <- c('lam0.fixed','sigma.fixed','theta.thin','N','D','lambda')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)

config.nodes <- c('lam0.fixed','sigma.fixed','D','theta.thin')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy=FALSE,
                      nodes=config.nodes) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

## required z/N update
z.ups <- round(M*0.25) # how many z proposals per iteration per session?
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.ID.nodes <- Rmodel$expandNodeNames(paste("y.ID[",g,",","1:",M[g],",1:",J[g],"]"))
  y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[",g,",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  bigLam.nodes <- Rmodel$expandNodeNames(paste("bigLam[",g,",",1:J[g],"]"))#only need this in calcNodes
  lam.noID.nodes <- Rmodel$expandNodeNames(paste("lam.noID[",g,",",1:J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.ID.nodes,y.noID.nodes)
  inds.detected <- which(rowSums(nimbuild$y.ID[g,,],na.rm=TRUE)>0)
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups[g],J=J[g],M=M[g],inds.detected=inds.detected,
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.ID.nodes=y.ID.nodes,y.noID.nodes=y.noID.nodes,
                                                   lam.nodes=lam.nodes,lam.noID.nodes=lam.noID.nodes,
                                                   N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

#required s sampler
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,J=J[g],xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#add block update for lam0, sigma, and/or D if posteriors correlated
conf$addSampler(target = c("lam0.fixed",'sigma.fixed'),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(1500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time


mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[25:nrow(mvSamples),]))

data$N #realized abundance
data$lambda #expected abundance

