###################################################################
# Custom nimbleFunctions to fit a categorical SPIM
###################################################################

#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)

#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dNBVector <- nimbleFunction(
  run = function(x = double(1), p = double(1), theta.d = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dnbinom(x, p = p, size = theta.d, log = TRUE))
      return(logProb)
    }
  }
)

#dummy random vector generator to make nimble happy
rNBVector <- nimbleFunction(
  run = function(n = integer(0),p = double(1),theta.d = double(1), z = double(0)) {
    returnType(double(1))
    J=nimDim(p)[1]
    out=numeric(J,value=0)
    return(out)
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), prob = double(0), size=double(1),capcounts = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(capcounts==0){
      return(0)
    }else{
      J=nimDim(size)[1]
      logProb <- 0
      for(j in 1:J){
        if(size[j]>0){
          logProb = logProb + dbinom(x[j], prob=prob, size=size[j], log = TRUE)
        }
      }
      return(logProb)
    }
  }
)
#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0),prob = double(0), size=double(1), capcounts = double(0)) {
    returnType(double(1))
    J=nimDim(size)[1]
    out=numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1),G.latent=double(2)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)


#------------------------------------------------------------------
# Custom sampler to update G.true, subject to constraints in G.latent
#------------------------------------------------------------------
#Can update from prior when detection function parameters do not depend on G.true
#Using MH here so that lam0 and sigma can vary by G.true values
GSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    n.cat <- control$n.cat
    n.levels <- control$n.levels
    n.samples <- control$n.samples
    G.noID <-control$G.noID
    calcNodes <- model$getDependencies(target)
    i <- control$i
    m <- control$m
  },
  run = function() {
    if(model$G.latent[i,m]==1){ #skip if not latent
      model.lp.initial <- model$getLogProb(calcNodes) #initial logProb
      prop.back <- model$gammaMat[m,model$G.true[i,m]] #backwards proposal prob
      G.prop <- rcat(1,model$gammaMat[m,1:n.levels[m]])
      prop.for <- model$gammaMat[m,G.prop] #forwards proposal prob
      model$G.true[i,m] <<- G.prop #store in model
      model.lp.proposed <- model$calculate(calcNodes)#proposed logProb
      log_MH_ratio <- (model.lp.proposed+log(prop.back)) - (model.lp.initial+log(prop.for))

      # log_MH_ratio
      accept <- decide(log_MH_ratio)
      if(accept) {
        copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      } else {
        copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      }
    }
  },
  methods = list( reset = function () {} )
)
#------------------------------------------------------------------
# Customer sampler to update latent IDs, and associated arrays
#------------------------------------------------------------------
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M<-control$M
    J <- control$J
    K1D <- control$K1D
    G.noID <- control$G.noID
    G.latent.ID <- control$G.latent.ID
    n.cat <- control$n.cat
    n.samples <- control$n.samples
    n.ID <- control$n.ID
    this.j <- control$this.j
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model$z
    G.true <- model$G.true
    y.true <- model$y.true
    y.ID <- model$y.ID
    
    #precalculate match. Does sample l match individual i?
    match=matrix(TRUE,nrow=n.samples,ncol=M) #start with all TRUE
    for(i in 1:M){
      if(z[i]==0){#no samples can match if z=0
        match[1:n.samples,i]=FALSE
      }else{
        for(l in 1:n.samples){ #each sample may conflict with certain z=1 individuals
          for(m in 1:n.cat){
            if(match[l,i]==TRUE){#so we abort when switched to false
              if(G.noID[l,m]!=0){#was sample l for cov m observed? If not, it matches
                if(G.noID[l,m]!=G.true[i,m]){#if so, does it match the value for individual i?
                  match[l,i]=FALSE
                }
              }
            }
          }
        }
      }
    }
    
    #precalculate log likelihoods. Can probably pull from NIMBLE somehow, but don't know how
    ll.y <- matrix(0,nrow=M,ncol=J)
    ll.y.ID <- matrix(0,nrow=M,ncol=J)
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          ll.y[i,j] <- dnbinom(y.true[i,j],size=model$theta.d[1]*model$K1D[j],prob=model$p[i,j],log=TRUE)
        }
      }
    }
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          if(y.true[i,j]>0){
            ll.y.ID[i,j] <- dbinom(y.ID[i,j],y.true[i,j],model$theta.i[i],log=TRUE)
          }
        }
      }
    }
    
    ll.y.cand <- ll.y
    ll.y.ID.cand <- ll.y.ID
    ID.curr <- model$ID

    ###update IDs
    for(l in 1:n.samples){#for all samples without known IDs
      ID.cand=ID.curr
      y.true.cand=y.true
      propprobs=model$lam[1:M,this.j[l]]
      for(i in 1:M){ #zero out nonmatches and z=0
        if(!match[l,i]){
          propprobs[i]=0
        }
      }
      denom=sum(propprobs) #abort if propprobs sum to 0. No matches anywhere nearby.
      if(denom>0){
        propprobs=propprobs/denom
        ID.cand[l]=rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped=c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob=propprobs[swapped[2]]
          backprob=propprobs[swapped[1]]
          #new y.true's - move sample from ID to ID.cand
          y.true.cand[ID.curr[l],this.j[l]]=y.true[ID.curr[l],this.j[l]]-1
          y.true.cand[ID.cand[l],this.j[l]]=y.true[ID.cand[l],this.j[l]]+1
          ll.y.cand[swapped[1],this.j[l]]=dnbinom(y.true.cand[swapped[1],this.j[l]],size=model$theta.d[1]*model$K1D[this.j[l]],prob=model$p[swapped[1],this.j[l]],log=TRUE)
          ll.y.cand[swapped[2],this.j[l]]=dnbinom(y.true.cand[swapped[2],this.j[l]],size=model$theta.d[1]*model$K1D[this.j[l]],prob=model$p[swapped[2],this.j[l]],log=TRUE)
          ll.y.ID.cand[swapped,this.j[l]] <- dbinom(y.ID[swapped,this.j[l]],y.true.cand[swapped,this.j[l]],model$theta.i[swapped],log=TRUE)
          
          #select sample to move proposal probabilities
          #P(select a sample of this type (not ID'd) for this ID)*P(select this j|sample of this type and this ID)
          #n.samples cancels out in MH ratio. Including for clarity
          focalprob=(sum(ID.curr==swapped[1])/n.samples)*
            (y.true[swapped[1],this.j[l]] - y.ID[swapped[1],this.j[l]])/sum(y.true[swapped[1],1:J] - y.ID[swapped[1],1:J])
          focalbackprob=(sum(ID.cand==swapped[2])/n.samples)*
            (y.true.cand[swapped[2],this.j[l]]- y.ID[swapped[2],this.j[l]])/sum(y.true.cand[swapped[2],1:J] - y.ID[swapped[2],1:J])

          #sum log likelihoods and do MH step
          lp_initial <- sum(ll.y[swapped,this.j[l]])+sum(ll.y.ID[swapped,this.j[l]])
          lp_proposed <- sum(ll.y.cand[swapped,this.j[l]])+sum(ll.y.ID.cand[swapped,this.j[l]])
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[swapped[1],this.j[l]]=y.true.cand[swapped[1],this.j[l]]
            y.true[swapped[2],this.j[l]]=y.true.cand[swapped[2],this.j[l]]
            ll.y[swapped[1],this.j[l]]=ll.y.cand[swapped[1],this.j[l]]
            ll.y[swapped[2],this.j[l]]=ll.y.cand[swapped[2],this.j[l]]
            ll.y.ID[swapped[1],this.j[l]]=ll.y.ID.cand[swapped[1],this.j[l]]
            ll.y.ID[swapped[2],this.j[l]]=ll.y.ID.cand[swapped[2],this.j[l]]
            ID.curr[l]=ID.cand[l]
          }
        }
      }
    }
    #update G.latent after ID changes
    G.latent=matrix(1,nrow=M,ncol=n.cat)
    G.latent[1:n.ID,1:n.cat]=G.latent.ID #known ID sample latent status contributions precalculated
    for(l in 1:n.samples){
      for(m in 1:n.cat){
        if(G.noID[l,m]!=0){ #if this sample is not latent
          G.latent[ID.curr[l],m]=0 #then its ID is not latent
        }
      }
    }
    
    #put everything back into the model$stuff after updating y.sight.true, y.sight.true.event
    model$y.true <<- y.true
    model$ID <<- ID.curr
    model$G.latent <<- G.latent
    model.lp.proposed <- model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)