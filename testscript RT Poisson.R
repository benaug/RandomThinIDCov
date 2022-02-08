library(nimble)
source("sim.RT.R")
source("NimbleModelRT Poisson.R")
source("NimbleFunctionsRT Poisson.R")
source("init.RT.R")

####Simulate some data####
N=78
#detection parameters
lam0=0.5
sigma=0.5

K=10 #number of occasions
buff=3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array
#sample thinning parameter
theta.thin=0.25

data=sim.RT(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,theta.thin=theta.thin,obstype="poisson")

#What is the observed data?
str(data$y.ID) #the observed ID detections
head(data$this.j) #the trap of detection for all unidentified detections
head(data$this.k) #occasion of capture, but not used in this 2D data sampler


##Fit model in Nimble##

#data augmentation level
M=175

#trap operation vector.
J=nrow(X)
K1D=rep(K,J)

inits=list(lam0=1,sigma=1) #ballpark inits to build data

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.RT(data,inits,M=M)

#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,sigma=inits$sigma,lam0=inits$lam0,theta.thin=0.5)

#constants for Nimble
J=nrow(data$X)
constants<-list(M=M,J=J,K=K,K1D=K1D,n.samples=nimbuild$n.samples,xlim=data$xlim,ylim=data$ylim)

# Supply data to Nimble. Note, y.true is completely latent.
z.data=c(rep(1,data$n.ID),rep(NA,M-data$n.ID))

Nimdata<-list(y.true=matrix(NA,nrow=M,ncol=J),y.ID=nimbuild$y.ID,
              ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M))

# set parameters to monitor
parameters<-c('psi','lam0','sigma','theta.thin','N','n')
nt=1 #thinning rate
#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2=50#thin more

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,nt2=nt2,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default samplers for y.true and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K1D=K1D,n.samples=nimbuild$n.samples,
                                                  n.ID=data$n.ID,this.j=nimbuild$this.j),
                silent = TRUE)

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#should be a little more efficient. Slice seems better than block random walk
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""), #do not adapt covariance bc s's not deterministically linked to unmarked individuals
  #                 type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=500),silent = TRUE)
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW.
#Need to not use this update or modify it when using lam0 or sigma covariates.
conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c("lam0","sigma"),
                  type = 'AF_slice',
                  control = list(adaptive=TRUE),
                  silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$n.cap #true number of captured individuals


mvSamples2 = as.matrix(Cmcmc$mvSamples2)
