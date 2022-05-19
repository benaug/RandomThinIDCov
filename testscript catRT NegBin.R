library(nimble)
library(coda)
source("NimbleModel catRT NegBin.R")
source("NimbleFunctions catRT NegBin.R")
source("init.catRT.R")
source("sim.catRT.R")
source("sSampler.R")

####Simulate some data####
N=78
#detection parameters
#vary by 1st ID cov (say, male-female), set up for 2 levels here
lam0=c(0.35,0.5)
sigma=c(0.75,0.5)
theta.d=0.1 #simulating 1 theta.d for all population groups (negbin overdispersion)

K=10 #number of occasions
buff=3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array

n.cat=2  #number of ID categories

gamma=IDcovs=vector("list",n.cat) #population frequencies of each category level. Assume equal here.
n.levels=rep(2,n.cat) #number of levels per IDcat
for(i in 1:n.cat){
  gamma[[i]]=rep(1/n.levels[i],n.levels[i])
  IDcovs[[i]]=1:n.levels[i]
}
theta.cat=rep(1,n.cat)#category level observation probabilities (probability you observe ID cov on capture)

#sample thinning
G.thin=2 #which ID cov is cov on thinning?
theta.thin=c(0.5,0.25) #must have n.levels[G.thin] thinning rates

##OK, what have we set up here? We have n.cat=2 ID covs. The data simulator will simulate
#from different detection parameters for the 1st cov. The first cov has 2 values, 1 and 2,
#and we have 2 lam0 and sigma parameters. We specify levels 1 and 2 are equal in the population
gamma[[1]]
#Then, since G.thin=2, the thinning rate varies by the 2nd ID cov, which also has 2 levels here and
#equal frequencies.
gamma[[2]]
#because there are 2 levels for this thinning covariate, we have 2 "theta.thin" parameters
theta.thin


data=sim.catRT(N=N,lam0=lam0,sigma=sigma,theta.d=theta.d,K=K,X=X,buff=buff,n.cat=n.cat,theta.cat=theta.cat,
               IDcovs=IDcovs,gamma=gamma,theta.thin=theta.thin,G.thin=G.thin,obstype="negbin")


#What are the observed data?
str(data$y.ID) #the identified samples
str(data$G.ID) #the observed ID covs for detected individuals (some may be NA)
head(data$this.j) #the trap of detection for all unidentified detections
head(data$this.k) #occasion of capture, but not used in this 2D data sampler
head(data$G.noID) #observed ID covs for every unobserved sample (some may be NA)
nrow(data$G.noID)==length(data$this.j) #equal

####Fit model in Nimble####

#Nimble code set up for G.thin=2 and ID covs a function of 1st ID covariate. Must modify if this is not the case.

#data augmentation level
M=175

J=nrow(X) #number of detectors
K1D=data$K1D #pull out trap operation

inits=list(lam0=lam0,sigma=sigma,theta.d=0.25,gamma=gamma) #ballpark inits to build data


#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.catRT(data,inits,M=M,obstype="negbin")

gammaMat=nimbuild$gammaMat
G.true.init=nimbuild$G.true
G.true.data=G.true.init*NA


#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,sigma=inits$sigma,lam0=inits$lam0,theta.d=inits$theta.d,
                 gammaMat=gammaMat,G.latent=nimbuild$G.latent,G.true=G.true.init)

#constants for Nimble
J=nrow(data$X)
constants<-list(M=M,J=J,K=K,K1D=data$K1D,n.samples=nimbuild$n.samples,xlim=data$xlim,ylim=data$ylim,
                n.cat=n.cat,n.levels=n.levels)

# Supply data to Nimble. Note, y.true is completely latent.
z.data=c(rep(1,data$n.ID),rep(NA,M-data$n.ID))

Nimdata<-list(y.true=matrix(NA,nrow=M,ncol=J),y.ID=nimbuild$y.ID,
              ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M),
              G.true=G.true.data)

# set parameters to monitor
parameters<-c('psi','lam0','sigma','theta.d','theta.thin','N','n','gammaMat')
nt=1 #thinning rate
#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2=25#thin more

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,nt2=nt2,useConjugacy = TRUE) 


###Two *required* sampler replacements

##Here, we remove the default samplers for y.true and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K1D=data$K1D,n.cat=n.cat,
                                          n.samples=nimbuild$n.samples,n.ID=data$n.ID,
                                      this.j=nimbuild$this.j,G.noID=nimbuild$G.noID,
                                      G.latent.ID=nimbuild$G.latent.ID),
                silent = TRUE)

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
conf$removeSampler("G.true")
for(i in 1:M){
  for(m in 1:n.cat){ #don't need to update first cat bc it is mark status
    conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
                    type = 'GSampler',
                    control = list(i = i,m=m,M=M,n.cat=n.cat,n.samples=nimbuild$n.samples,
                                   n.levels=n.levels,G.noID=nimbuild$G.noID), silent = TRUE)
  }
}

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#should be a little more efficient. Slice seems better than block random walk
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
for(i in 1:n.levels[1]){
  conf$addSampler(target = c(paste("lam0[",i,"]",sep=""),paste("sigma[",i,"]",sep="")),
                  type = 'RW_block',
                  control = list(adaptive=TRUE),
                  silent = TRUE)
}

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

data$n.cap #true number of captured individuals


#look at ID posteriors. Not removing any burnin here...
mvSamples2 = as.matrix(Cmcmc$mvSamples2)

check.sample=2
#posterior prob this sample belongs to each individual number
round(table(mvSamples2[,check.sample])/nrow(mvSamples2),2)
#truth
data$ID[check.sample]

#individual numbers larger than n.cap were not captured and identified
data$n.cap

