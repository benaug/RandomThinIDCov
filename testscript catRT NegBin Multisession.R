library(nimble)
library(coda)
source("sim.catRT.multisession.R")
source("sim.catRT.R")
source("NimbleModel catRT NegBin Multisession.R")
source("NimbleFunctions catRT NegBin Multisession.R")
source("init.catRT.multisession.R")
source("init.catRT.R")
source("sSampler Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session <- 3
D <- rep(0.4,N.session) #expected density in units of sigma and X

#here, we have N.session x 2 group detection functions. 
#Could be sex-specific parms shared across sessions.
lam0 <- cbind(rep(0.5,N.session),rep(0.25,N.session))
sigma <- cbind(rep(0.5,N.session),rep(0.75,N.session))
theta.d <- rep(0.1,N.session) # overdispersion parameter. Assuming shared across df groups here.

K <- c(8,9,10) #number of occasions
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

#categorical ID covariate stuff - not session-specific in data simulator
n.cat <- 2  #number of ID categories (not including marked status)
gamma <- IDcovs <- vector("list",n.cat) #population frequencies of each category level. Assume equal here.
n.levels <- rep(2,n.cat) #number of levels per IDcat
if(all(n.levels==1))stop("This specification has no categorial ID covariates. Use testscript for regular SMR.")
for(i in 1:n.cat){
  gamma[[i]] <- rep(1/n.levels[i],n.levels[i])
  IDcovs[[i]] <- 1:n.levels[i]
}
theta.cat <- rep(1,n.cat)#sample-level IDcov observation probabilities. Data missing at random if <1. 

#sample thinning
G.thin <- 2 #which ID cov is cov on thinning?
theta.thin <- cbind(rep(0.5,N.session),rep(0.25,N.session)) #must have n.levels[G.thin] thinning rates per session

##OK, what have we set up here? We have n.cat=2 ID covs. The data simulator will simulate
#from different detection parameters for the 1st cov. The first cov has 2 values, 1 and 2,
#and we have 2 lam0 and sigma parameters. We specify levels 1 and 2 are equal in the population
gamma[[1]]
#Then, since G.thin=2, the thinning rate varies by the 2nd ID cov, which also has 2 levels here and
#equal frequencies.
gamma[[2]]
#because there are 2 levels for this thinning covariate, we have 2 "theta.thin" parameters
theta.thin

data <- sim.catRT.multisession(N.session=N.session,lambda=lambda,
                         lam0=lam0,sigma=sigma,theta.d=theta.d,K=K,X=X,buff=buff,theta.thin=theta.thin,
                         obstype="negbin",
                         n.cat=n.cat,IDcovs=IDcovs,gamma=gamma,theta.cat=theta.cat,G.thin=G.thin)

#What is the observed data?
str(data$y.ID) #the observed ID detections. One list element per session
str(data$G.ID) #the observed ID covs for detected individuals (some may be NA)
head(t(data$this.j)) #the trap of detection for all unidentified detections.
head(t(data$this.k)) #occasion of capture, but not used in this 2D data sampler
#here are observed ID covariate data. Missing values coded with "0" instead of "NA"
dim(data$G.noID) #observed ID covs for every unobserved sample (some may be NA)

##Fit model in Nimble##

#data augmentation level
M <- c(155,165,175)

inits <- list(lam0=lam0,sigma=sigma,theta.d=theta.d) #ballpark inits to build data, one per session. Using simulated values here, don't do this in practice.

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.catRT.multisession(data,inits,M=M,obstype="negbin")
capcounts.ID <- apply(nimbuild$y.ID,c(1,2),sum)



#inits for nimble
N.init <- rowSums(nimbuild$z,na.rm=TRUE)
gammaMat.init <- nimbuild$gammaMat
theta.thin.init <- matrix(c(0.5,0.5),N.session,2,byrow=TRUE)
G.true.init <- nimbuild$G.true
G.true.data <- G.true.init*NA

Niminits <- list(N=N.init,
                 z=nimbuild$z,s=nimbuild$s,ID=nimbuild$ID,capcounts=apply(nimbuild$y.true,c(1,2),sum,na.rm=TRUE),
                 y.true=nimbuild$y.true,D=0.5,theta.d=c(0.2),
                 theta.thin=theta.thin.init,lam0.fixed=c(0.3,0.3),sigma.fixed=c(0.5,0.5), #one init per lam0.fixed, sigma.fixed
                 gammaMat=gammaMat.init,G.true=G.true.init,G.latent=nimbuild$G.latent)

#constants for Nimble
J <- unlist(lapply(data$X,nrow)) #number of detectors
constants<-list(N.session=N.session,M=M,J=J,K1D=nimbuild$K1D,n.samples=nimbuild$n.samples,
                xlim=data$xlim,ylim=data$ylim,area=area,n.cat=nimbuild$n.cat,n.levels=nimbuild$n.levels)

# Supply data to Nimble. Note, y.true is completely latent.
z.data <- matrix(NA,N.session,max(M))
for(g in 1:N.session){
  z.data[g,1:data$n.ID[g]] <- 1
}

Nimdata <- list(y.ID=nimbuild$y.ID,y.true=array(NA,dim=c(N.session,max(M),max(J))),
              ID=matrix(NA,N.session,max(nimbuild$n.samples)),
              z=z.data,X=nimbuild$X,capcounts=matrix(NA,N.session,max(M)),G.true=G.true.data,capcounts.ID=capcounts.ID)

# set parameters to monitor
parameters <- c('lam0.fixed','sigma.fixed','theta.d','theta.thin','N','n','D','lambda','gammaMat')
nt <- 1 #thinning rate
#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt2 <- 25#thin more

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#can use "nodes" argument in configureMCMC below to omit y.true and G.true that are replaced below for faster
#configuration
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Three *required* sampler replacements

##Here, we remove the default samplers for y.true and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.true")
for(g in 1:N.session){
  conf$addSampler(target = paste0("y.true[",g,",1:",M[g],",1:",J[g],"]"),
                  type = 'IDSampler',control = list(M=M[g],J=J[g],K1D=nimbuild$K1D[g,1:J[g]],
                                                    n.samples=nimbuild$n.samples[g],g=g,
                                                    this.j=nimbuild$this.j[g,1:nimbuild$n.samples[g]],
                                                    n.cat=nimbuild$n.cat[g],n.ID=data$n.ID[g],
                                                    G.noID=nimbuild$G.noID[g,1:nimbuild$n.samples[g],1:nimbuild$n.cat[g]],
                                                    G.latent.ID=nimbuild$G.latent.ID[g,1:data$n.ID[g],1:nimbuild$n.cat[g]]),
                  silent = TRUE)
}

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
conf$removeSampler("G.true")
for(g in 1:N.session){
  for(i in 1:M[g]){
    for(m in 1:nimbuild$n.cat[g]){
      conf$addSampler(target = paste("G.true[",g,",",i,",",m,"]", sep=""),
                      type = 'GSampler',
                      control = list(g=g,i=i,m=m,n.levels=nimbuild$n.levels[g,]), silent = TRUE) 
    }
  }
}

z.ups <- round(M*0.25) # how many z proposals per iteration per session?
conf$removeSampler("N")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y.true[",g,",","1:",M[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  p.nodes <- Rmodel$expandNodeNames(paste("p[",g,",","1:",M[g],",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,y.nodes,lam.nodes,p.nodes)
  
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(inds.detected=1:data$n.ID[g],z.ups=z.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,p.nodes=p.nodes,N.node=N.node,
                                                   z.nodes=z.nodes,calcNodes=calcNodes),
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
for(df in 1:nimbuild$n.levels[1,1]){
  conf$addSampler(target = c(paste("lam0.fixed[",df,"]"),paste("sigma.fixed[",df,"]")),type = 'RW_block',
                  control = list(adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time


mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$N #abundance
data$n.cap #true number of captured individuals

