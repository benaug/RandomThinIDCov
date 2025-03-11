library(nimble)
library(coda)
source("sim.RT.R")
source("NimbleModelRT Poisson DA2.R")
source("NimbleFunctionsRT Poisson DA2.R")
source("init.RT.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
N <- 78
#detection parameters
lam0 <- 0.5
sigma <- 0.5

K <- 10 #number of occasions
buff <- 3 #state space buffer
X <- expand.grid(3:11,3:11) #make a trapping array
#sample thinning parameter
theta.thin <- 0.25

data <- sim.RT(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,theta.thin=theta.thin,obstype="poisson")

#What is the observed data?
str(data$y.ID) #the observed ID detections
head(data$this.j) #the trap of detection for all unidentified detections
head(data$this.k) #occasion of capture, but not used in this 2D data sampler


##Fit model in Nimble##

#data augmentation level
M <- 175

J <- nrow(X) #number of detectors
K1D <- data$K1D #pull out trap operation

inits <- list(lam0=1,sigma=1) #ballpark inits to build data

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.RT(data,inits,M=M,obstype="poisson")
capcounts.ID <- rowSums(nimbuild$y.ID)

#inits for nimble
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #z and N inits must be consistent
                 lambda.N=sum(nimbuild$z), #converges faster if you set N and lambda.N at similar values
                 s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,sigma=inits$sigma,lam0=inits$lam0,theta.thin=0.5)

#constants for Nimble
J <- nrow(data$X)
constants <- list(M=M,J=J,K1D=K1D,xlim=data$xlim,ylim=data$ylim,n.samples=nimbuild$n.samples)

# Supply data to Nimble. Note, y.true is completely latent.
z.data <- c(rep(1,data$n.ID),rep(NA,M-data$n.ID))

Nimdata <- list(y.true=matrix(NA,nrow=M,ncol=J),y.ID=nimbuild$y.ID,
              ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M),capcounts.ID=capcounts.ID)

# set parameters to monitor
parameters <- c('lam0','sigma','theta.thin','N','lambda.N','n')
nt <- 1 #thinning rate
#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2 <- 50 #thin more

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c('lam0','sigma','theta.thin','lambda.N')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,
                      nodes=config.nodes) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###One *required* sampler replacements

##Here, we remove the default samplers for y.true and y.event, which are not correct
#and replace it with the custom "IDSampler"
# conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K1D=K1D,n.samples=nimbuild$n.samples,
                                                  this.j=nimbuild$this.j),
                silent = TRUE)

###*required* sampler replacement
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal.
#nodes used for update
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 y.nodes=y.nodes,lam.nodes=lam.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  #                 type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW.
#Need to not use this update or modify it when using lam0 or sigma covariates.
conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c("lam0","sigma"),
                  type = 'RW_block',
                  control = list(adaptive=TRUE),
                  silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[-c(1:burnin),]))

data$n.cap #true number of captured individuals


#look at ID posteriors. Not removing any burnin here...
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)

check.sample <- 1
#posterior prob this sample belongs to each individual number
round(table(mvSamples2[,check.sample])/(nrow(mvSamples2)-1),2)
#truth (for simulated data sets)
data$ID[check.sample]
#These should match data$ID for all individuals with identified captures, 1, ... n.ID
#individuals with ID>n.ID will have different numbers than data$ID (but samples with same true ID should tend to have same posterior ID)
data$n.ID