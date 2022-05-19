library(nimble)
library(coda)
source("sim.RT.multisession.R")
source("sim.RT.R")
source("NimbleModelRT Poisson Multisession.R")
source("NimbleFunctionsRT Poisson Multisession.R")
source("init.RT.multisession.R")
source("init.RT.R")
source("sSampler Multisession.R")

#make sure to run this line!
nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches')

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session=3
D = rep(0.4,N.session) #expected density in units of sigma and X
lam0=rep(0.5,N.session)
sigma=rep(0.5,N.session)
K=c(5,6,7) #number of occasions
buff=rep(2,N.session) #state space buffer

#make trapping arrays
X1=expand.grid(3:11,3:11)
X2=expand.grid(3:12,3:12)
X3=expand.grid(3:13,3:13)
X=list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area=getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda=D*area
lambda #expected N in each session

#sample thinning parameter
theta.thin=rep(0.25,N.session)

data=sim.RT.multisession(N.session=N.session,lambda=lambda,
                         lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,theta.thin=theta.thin,obstype="poisson")

#What is the observed data?
str(data$y.ID) #the observed ID detections. One list element per session
head(t(data$this.j)) #the trap of detection for all unidentified detections.
head(t(data$this.k)) #occasion of capture, but not used in this 2D data sampler


##Fit model in Nimble##

#data augmentation level
M=c(155,165,175)

inits=list(lam0=lam0,sigma=sigma) #ballpark inits to build data, one per session. Using simulated values here, don't do this in practice.

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.RT.multisession(data,inits,M=M,obstype="poisson")

#inits for nimble
N.init=rowSums(nimbuild$z,na.rm=TRUE)

Niminits <- list(N=N.init,
                 z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=apply(nimbuild$y.true,c(1,2),sum,na.rm=TRUE),
                 y.true=nimbuild$y.true,
                 theta.thin=rep(0.5,N.session),lam0.fixed=0.5,sigma.fixed=0.5, #one init per lam0.fixed, sigma.fixed
                 D=0.5)

#constants for Nimble
J=unlist(lapply(data$X,nrow)) #number of detectors
constants<-list(N.session=N.session,M=M,J=J,K=K,K1D=nimbuild$K1D,n.samples=nimbuild$n.samples,
                xlim=data$xlim,ylim=data$ylim,area=area)

# Supply data to Nimble. Note, y.true is completely latent.
z.data=matrix(NA,N.session,max(M))
for(g in 1:N.session){
  z.data[g,1:data$n.ID[g]]=1
}

Nimdata<-list(y.ID=nimbuild$y.ID,y.true=array(NA,dim=c(N.session,max(M),max(J))),
              ID=matrix(NA,N.session,max(nimbuild$n.samples)),
              z=z.data,X=nimbuild$X,capcounts=matrix(NA,N.session,max(M)))

# set parameters to monitor
parameters<-c('lam0.fixed','sigma.fixed','theta.thin','N','n','D','lambda')
nt=1 #thinning rate
#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2=25#thin more

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default samplers for y.true and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.true")
for(g in 1:N.session){
  conf$addSampler(target = paste0("y.true[",g,",1:",M[g],",1:",J[g],"]"),
                  type = 'IDSampler',control = list(M=M[g],J=J[g],K1D=nimbuild$K1D[g,1:J[g]],
                                                    n.samples=nimbuild$n.samples[g],
                                                    this.j=nimbuild$this.j[g,1:nimbuild$n.samples[g]],
                                                    g=g),
                  silent = TRUE)
}

z.ups=c(25,25,25) # how many z proposals per iteration per session?
conf$removeSampler("N")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y.true[",g,",","1:",M[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,y.nodes,lam.nodes)
  
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(inds.detected=1:data$n.ID[g],z.ups=z.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}


###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler("s")
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW. 
#Need to not use this update or modify it when using lam0 or sigma covariates.
#This sampler is slower, so not worth it if data is not so sparse there is strong posterior correlation
#between lam0 and sigma.
conf$removeSampler(c("lam0.fixed","sigma.fixed"))
conf$addSampler(target = c("lam0.fixed","sigma.fixed"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)

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


mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$N #abundance
data$n.cap #true number of captured individuals

