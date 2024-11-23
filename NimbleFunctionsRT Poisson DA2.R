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

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lambda, log = TRUE))
      return(logProb)
    }
  }
)
#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lambda = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(lambda)[1]
    out <- numeric(J,value=0)
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
      J <- nimDim(size)[1]
      logProb <- 0
      for(j in 1:J){
        if(size[j]>0){
          logProb <- logProb + dbinom(x[j], prob=prob, size=size[j], log = TRUE)
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
    J <- nimDim(size)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts <- numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1)){ #don't need ID, but nimble requires is it used in a function 
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
# Customer sampler to update latent IDs, and associated arrays
#------------------------------------------------------------------
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    J <- control$J
    K1D <- control$K1D
    n.samples <- control$n.samples
    this.j <- control$this.j
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    z <- model$z
    y.true <- model$y.true
    y.ID <- model$y.ID
    
    #precalculate match. Does sample l match individual i?
    match <- matrix(TRUE,nrow=n.samples,ncol=M) #start with all TRUE
    for(i in 1:M){#can match any individual
      if(z[i]==0){#unless z is off
        match[1:n.samples,i] <- FALSE
      }
    }
    
    #precalculate log likelihoods. Can probably pull from NIMBLE somehow, but don't know how
    ll.y <- matrix(0,nrow=M,ncol=J)
    ll.y.ID <- matrix(0,nrow=M,ncol=J)
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          ll.y[i,j] <- dpois(y.true[i,j],K1D[j]*model$lam[i,j],log=TRUE)
        }
      }
    }
    for(i in 1:M){
      if(z[i]==1){
        for(j in 1:J){
          if(y.true[i,j]>0){
            ll.y.ID[i,j] <- dbinom(y.ID[i,j],y.true[i,j],model$theta.thin[1],log=TRUE)
          }
        }
      }
    }
    
    ll.y.cand <- ll.y
    ll.y.ID.cand <- ll.y.ID
    ID.curr <- model$ID

    ###update IDs
    for(l in 1:n.samples){#for all samples without known IDs
      ID.cand <- ID.curr
      y.true.cand <- y.true
      propprobs <- model$lam[1:M,this.j[l]]
      for(i in 1:M){ #zero out nonmatches and z=0
        if(!match[l,i]){
          propprobs[i] <- 0
        }
      }
      denom <- sum(propprobs) #abort if propprobs sum to 0. No matches anywhere nearby.
      if(denom>0){
        propprobs <- propprobs/denom
        ID.cand[l] <- rcat(1,prob=propprobs)
        if(ID.cand[l]!=ID.curr[l]){
          swapped <- c(ID.curr[l],ID.cand[l])
          #new sample proposal probabilities
          forprob <- propprobs[swapped[2]]
          backprob <- propprobs[swapped[1]]
          #new y.true's - move sample from ID to ID.cand
          y.true.cand[ID.curr[l],this.j[l]] <- y.true[ID.curr[l],this.j[l]]-1
          y.true.cand[ID.cand[l],this.j[l]] <- y.true[ID.cand[l],this.j[l]]+1
          ll.y.cand[swapped,this.j[l]] <- dpois(y.true.cand[swapped,this.j[l]],model$lam[swapped,this.j[l]]*K1D[this.j[l]],log=TRUE)
          ll.y.ID.cand[swapped,this.j[l]] <- dbinom(y.ID[swapped,this.j[l]],y.true.cand[swapped,this.j[l]],model$theta.thin,log=TRUE)
          
          #select sample to move proposal probabilities
          #P(select a sample of this type (not ID'd) for this ID)*P(select this j|sample of this type and this ID)
          #n.samples cancels out in MH ratio. Including for clarity
          focalprob <- (sum(ID.curr==swapped[1])/n.samples)*
            (y.true[swapped[1],this.j[l]] - y.ID[swapped[1],this.j[l]])/sum(y.true[swapped[1],1:J] - y.ID[swapped[1],1:J])
          focalbackprob <- (sum(ID.cand==swapped[2])/n.samples)*
            (y.true.cand[swapped[2],this.j[l]]- y.ID[swapped[2],this.j[l]])/sum(y.true.cand[swapped[2],1:J] - y.ID[swapped[2],1:J])

          #sum log likelihoods and do MH step
          lp_initial <- sum(ll.y[swapped,this.j[l]])+sum(ll.y.ID[swapped,this.j[l]])
          lp_proposed <- sum(ll.y.cand[swapped,this.j[l]])+sum(ll.y.ID.cand[swapped,this.j[l]])
          log_MH_ratio <- (lp_proposed+log(backprob)+log(focalbackprob)) - (lp_initial+log(forprob)+log(focalprob))
          accept <- decide(log_MH_ratio)
          if(accept){
            y.true[swapped[1],this.j[l]] <- y.true.cand[swapped[1],this.j[l]]
            y.true[swapped[2],this.j[l]] <- y.true.cand[swapped[2],this.j[l]]
            ll.y[swapped[1],this.j[l]] <- ll.y.cand[swapped[1],this.j[l]]
            ll.y[swapped[2],this.j[l]] <- ll.y.cand[swapped[2],this.j[l]]
            ll.y.ID[swapped[1],this.j[l]] <- ll.y.ID.cand[swapped[1],this.j[l]]
            ll.y.ID[swapped[2],this.j[l]] <- ll.y.ID.cand[swapped[2],this.j[l]]
            ID.curr[l] <- ID.cand[l]
          }
        }
      }
    }
    
    #put everything back into the model$stuff
    model$y.true <<- y.true
    model$ID <<- ID.curr
    model.lp.proposed <- model$calculate(calcNodes) #update logprob
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    lam.nodes <- control$lam.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(model$capcounts[pick]>0){#is this an individual with samples?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn pd off
          model$calculate(lam.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn pd on
          model$calculate(lam.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)