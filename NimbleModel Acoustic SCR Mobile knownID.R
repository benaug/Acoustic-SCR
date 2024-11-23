NimModel <- nimbleCode({
  ###priors###
  lambda.C ~ dunif(0,100) #mean calls/individual
  sigma.u ~ dunif(0,10) #call dispersion around individual
  #attenuation function parameters
  beta0 ~ dunif(50,100)
  beta1 ~ dunif(-100,0) #assuming must be negative (dB decays with distance)
  sigma.d ~ dunif(0,10)
  #data augmentation prior
  psi~dunif(0,1)
  
  #Model for individuals (parents)
  for(i in 1:M.s){
    z.s[i] ~ dbern(psi)
    #individual activity centers
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #call model
    calls[i] ~ dpois(lambda.C) #number of total calls per individual
  }
  #Model for calls (children)
  for(l in 1:M.u){
    #only update sigma.u from z.u=1 calls, requires custom u distribution and update
    u[l,1:2] ~ dnorm_censor(s[callID[l],1:2],sd=sigma.u,z.u = z.u[l])
    D[l,1:J] <- GetD(u=u[l,1:2], X=X[1:J,1:2], z.u=z.u[l])
    mu[l,1:J] <- GetMu(beta0=beta0, beta1=beta1, D=D[l,1:J], z.u=z.u[l])
    dB[l,1:J] ~ ddB(mu=mu[l,1:J], zeros=zeros[l,1:J], sigma.d=sigma.d, mindB=mindB,z.u=z.u[l])
  }
  N.ind <- sum(z.s[1:M.s]) #individual abundance
  N.call <- sum(z.u[1:M.u]) #call abundance
})# end model

#function to get distance from call to all ARUs (only if z.u=1)
GetD <- nimbleFunction(
  run = function(u = double(1), X=double(2),z.u=double(0)){
    returnType(double(1))
    J <- nimDim(X)[1]
    if(z.u==0){
      D <- rep(0,J)
    }else{
      D <- sqrt((u[1]-X[1:J,1])^2 + (u[2]-X[1:J,2])^2)
    }
    return(D)
  }
)

#function to get expected dB level for each call-ARU (only if z.u=1)
GetMu <- nimbleFunction(
  run = function(beta0=double(0),beta1=double(0), D = double(1),z.u=double(0)){
    returnType(double(1))
    J <- nimDim(D)[1]
    if(z.u==0){
      mu <- rep(0,J)
    }else{
      mu <- beta0 + beta1*D[1:J]
    }
    return(mu)
  }
)

#observation model (cannot be observed if z.u=0)
ddB <- nimbleFunction(
  run = function(x = double(1),mu = double(1), zeros = double(1),sigma.d = double(0),
                 mindB = double(0), z.u = double(0), log = integer(0)) {
    returnType(double(0))
    J <- nimDim(x)[1]
    logProb <- 0
    if(z.u==1){
      for(j in 1:J){
        if(zeros[j]){ #below dB threshold
          logProb <- logProb + pnorm(mindB, mu[j],sigma.d,log=TRUE)
        }else{ #above dB threshold
          logProb <- logProb + dnorm(x[j],mu[j],sigma.d,log=TRUE)
        }
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

#dummy RNG not used. makes nimble happy.
rdB <- nimbleFunction(
  run = function(n = integer(0),mu = double(1), zeros = double(1),sigma.d = double(0),
                 mindB = double(0), z.u = double(0)) {
    returnType(double(1))
    J <- nimDim(mu)[1]
    return(rep(0,J))
  }
)

#distrubution for dispersion of calls (children) around individuals (parents)
#bivariate normal, but only used if z.u=1
dnorm_censor <- nimbleFunction(
  run = function(x = double(1), s = double(1), sd= double(0),
                 z.u = double(0), log = integer(0)) {
    returnType(double(0))
    logProb=0
    if(z.u==1){
      logProb <- logProb + dnorm(x[1],s[1],sd=sd,log=log)
      logProb <- logProb + dnorm(x[2],s[2],sd=sd,log=log)
    }#otherwise, 0
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

#dummy function to make nimble happy
rnorm_censor <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sd= double(0),
                 z.u = double(0)) {
    returnType(double(1))
    return(rep(0,2))
  }
)


#Required custom update for z.s (parents)
#turn entire clusters on/off
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("z.s","dB","u"))
    calls.detected <- control$calls.detected
    J <- control$J
    M.s <- control$M.s
    M.u <- control$M.u
    mindB <- control$mindB
  },
  run = function() {
    #Nimble makes it hard to do this.... pulling out relevant likelihoods here, updating all individuals,
    #then handing it back to nimble
    ll.z <- model$logProb_z.s
    ll.dB <- model$logProb_dB[,1]
    ll.dB.cand <- ll.dB

    for(i in 1:M.s){
      if(model$z.s[i]==0){ #if z.s is off, turn on along with calls
        z.s.cand <- 1
        ll.z.cand <- dbinom(z.s.cand,1,model$psi[1],log=TRUE)
        empty.idx <- which(model$z.u==0) #find empty z.u's
        if(length(empty.idx)>=model$calls[i]){#skip if z.u is full. Need to raise M.u.  
          these.calls <- empty.idx[1:model$calls[i]]
          n.these.calls <- length(these.calls)
          
          #3 cases can occur that need to be dealt with slightly differently in nimble
          #1) There can be 0 calls to turn on, 2) there can be 1 call to turn on, 3) more than 1
          if(n.these.calls==0){#no calls, only z likelihoods to consider
            logProb.proposed <- ll.z.cand
            logProb.initial <- ll.z[i]
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              model$z.s[i] <<- z.s.cand
              ll.z[i] <-  ll.z.cand
            }
          }else if(n.these.calls==1){ #1 call
            u.cand <- c(rnorm(1,model$s[i,1],model$sigma.u[1]), rnorm(1,model$s[i,2],model$sigma.u[1]))
            #update D, mu from simulated u
            D.cand <- GetD(u=u.cand, X=model$X[1:J,1:2], z.u=1)
            mu.cand <- GetMu(beta0=model$beta0[1], beta1=model$beta1[1], D=D.cand,z.u=1)
            #get proposed logprobs
            ll.dB.cand[these.calls[1]] <- ddB(x=model$dB[these.calls[1],1:J],mu=mu.cand, 
                                              zeros=model$zeros[these.calls[1],1:J],sigma.d=model$sigma.d[1],
                                              mindB=mindB,z.u=1,log=TRUE)
            logProb.proposed <- ll.z.cand + ll.dB.cand[these.calls[1]]
            logProb.initial <- ll.z[i] + ll.dB[these.calls[1]]
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              model$z.s[i] <<- z.s.cand
              model$z.u[these.calls[1]] <<- 1
              #update log likelihoods
              ll.z[i] <-  ll.z.cand
              ll.dB[these.calls[1]] <-  ll.dB.cand[these.calls[1]]
              model$u[these.calls[1],] <<- u.cand
              model$D[these.calls[1],] <<- D.cand
              model$mu[these.calls[1],] <<- mu.cand
              model$callID[these.calls[1]] <<- i
            }
          }else{ #multiple calls to turn on
            #nimble doesn't like same object being defined more than once with
            #different dimensions, so using ".mat"
            u.cand.mat <- matrix(0,n.these.calls,2)
            D.cand.mat <- matrix(0,n.these.calls,J)
            mu.cand.mat <- matrix(0,n.these.calls,J)
            for(l in 1:n.these.calls){
              u.cand.mat[l,] <- c(rnorm(1,model$s[i,1],model$sigma.u[1]), rnorm(1,model$s[i,2],model$sigma.u[1]))
              #update D, mu from simulated u
              D.cand.mat[l,] <- GetD(u=u.cand.mat[l,], X=model$X[1:J,1:2], z.u=1)
              mu.cand.mat[l,] <- GetMu(beta0=model$beta0[1], beta1=model$beta1[1], D=D.cand.mat[l,],z.u=1)
              #get proposed logprobs
              ll.dB.cand[these.calls[l]] <- ddB(x=model$dB[these.calls[l],1:J],mu=mu.cand.mat[l,], zeros=model$zeros[these.calls[l],1:J], sigma.d=model$sigma.d[1],
                                                mindB=mindB,z.u=1,log=TRUE)
            }
            logProb.proposed <- ll.z.cand + sum(ll.dB.cand[these.calls])
            logProb.initial <- ll.z[i] + sum(ll.dB[these.calls])
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              ll.z[i] <-  ll.z.cand
              ll.dB[these.calls] <-  ll.dB.cand[these.calls]
              model$z.s[i] <<- z.s.cand
              model$z.u[these.calls] <<- rep(1,length(these.calls))
              model$u[these.calls,] <<- u.cand.mat
              model$D[these.calls,] <<- D.cand.mat
              model$mu[these.calls,] <<- mu.cand.mat
              model$callID[these.calls] <<- rep(i,n.these.calls)
            }
          }
        }else{
          print("Warning: Not enough z.u's left to turn on a z.s. Raise M.u if this happens past burn in.")
        }
      }else{ #if z is on, turn off with calls
        z.s.cand <- 0
        ll.z.cand <- dbinom(z.s.cand,1,model$psi[1],log=TRUE)
        these.calls <- which(model$callID==i&model$z.u==1)
        n.these.calls <- length(these.calls)
        #2 cases can occur that need to be dealt with slightly differently in nimble
        #1) There can be 0 calls to turn off, 2) there can be 1 call to turn off, 3) more than 1
        if(n.these.calls==0){#no calls, only z likelihoods to consider
          logProb.proposed <- ll.z.cand
          logProb.initial <- ll.z[i]
          log_MH_ratio <- (logProb.proposed) - (logProb.initial)
          accept <- decide(log_MH_ratio)
          if(accept) {
            model$z.s[i] <<- z.s.cand
            ll.z[i] <- ll.z.cand
          }
        }else if(n.these.calls==1){
          if(calls.detected[these.calls[1]]==0){#cannot turn z off if any calls detected
            ll.dB.cand[these.calls[1]] <- 0
            #add em up
            logProb.proposed <- ll.z.cand + ll.dB.cand[these.calls[1]]
            logProb.initial <- ll.z[i] + ll.dB[these.calls[1]]
            #MH
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              model$z.s[i] <<- z.s.cand
              model$z.u[these.calls[1]] <<- 0
              #update log likelihoods
              ll.z[i] <- ll.z.cand
              ll.dB[these.calls[1]] <- ll.dB.cand[these.calls[1]]
            }
          }
        }else{
          if(sum(calls.detected[these.calls])==0){#cannot turn z off if any calls detected
            ll.dB.cand[these.calls] <- rep(0,length(these.calls))
            #add em up
            logProb.proposed <- ll.z.cand + sum(ll.dB.cand[these.calls])
            logProb.initial <- ll.z[i] + sum(ll.dB[these.calls])
            #MH
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              model$z.s[i] <<- z.s.cand
              model$z.u[these.calls] <<- rep(0,length(these.calls))
              #update log likelihoods
              ll.z[i] <-  ll.z.cand
              ll.dB[these.calls] <- ll.dB.cand[these.calls]
            }
          }            # print("Warning: Not enough z.u's left to turn on a z.s. Raise M.u if this happens past burn in.")
        }
      }
    }
    #hand back to nimble
    model$N.ind[1] <<- sum(model$z.s[1:M.s])
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)


#Required custom update for number of calls (turn children on/off without changing parent status)
callSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(c("calls","dB","u"))
    calls.detected <- control$calls.detected
    call.ups <- control$call.ups
    J <- control$J
    M.s <- control$M.s
    M.u <- control$M.u
    mindB <- control$mindB
  },
  run = function() {
    #Nimble makes it hard to do this.... pulling out relevant likelihoods here, updating all individuals,
    #then handing it back to nimble
    ll.call <- model$logProb_calls
    ll.dB <- model$logProb_dB[,1]

    for(i in 1:M.s){
      if(model$z.s[i]==0){ #if z is off, propose from prior
        calls.cand <- rpois(1,model$lambda.C[1])
        model$calls[i] <<- calls.cand
      }else{ #if z is on, use Metropolis-Hastings
        #we'll propose to add or subtract 1 count with equal probability.
        #if adding, we can always choose a call to turn on
        #if subtracting, we need to consider that if we choose a detected call
        #it cannot be turned off. This could be incorporated into the proposal probs
        #but I'm just rejecting the update if you select a detected call to subtract.
        #You must account for this or it won't work.
        for(up in 1:call.ups){ #how many updates per iteration?
          #propose to add/subtract 1
          updown <- rbinom(1,1,0.5) #p=0.5 is symmetric
          reject <- FALSE #we reject if 1) proposed counts <0 or 2) you select a detected call
          if(updown==0){#subtract
            if(model$calls[i]==0){
              reject <- TRUE
            }else{
              #find all l indices associated with calls for individual i that are turned on
              these.calls <- which(model$callID==i&model$z.u==1)
              pick <- rcat(1,rep(1/model$calls[i],model$calls[i])) #select one of the calls
              l.idx <- these.calls[pick]
              if(calls.detected[l.idx]==1){ #is it one of the detected calls?
                reject <- TRUE #if so, we reject
              }
            }
            if(!reject){
              #propose new calls
              calls.cand <- model$calls[i]-1
              #get proposed logprobs
              ll.call.cand <- dpois(calls.cand,model$lambda.C[1],log=TRUE)
              ll.dB.cand <- 0 #always 0

              #add em up
              logProb.proposed <- ll.call.cand + ll.dB.cand
              logProb.initial <- ll.call[i] + ll.dB[l.idx]
              #MH
              log_MH_ratio <- (logProb.proposed) - (logProb.initial)
              accept <- decide(log_MH_ratio)
              if(accept) {
                model$calls[i] <<- calls.cand
                #turn these off
                model$z.u[l.idx] <<- 0
                model$D[l.idx,1:J] <<- rep(0,J)
                model$mu[l.idx,1:J] <<- rep(0,J)
                #update log likelihoods
                ll.call[i] <- ll.call.cand
                ll.dB[l.idx] <- ll.dB.cand
              }
            }
          }else{#add
            #find first empty u index
            if(sum(model$z.u)<M.u){ #cannot update if z.u maxed out. Need to raise M.u
              l.idx <- which(model$z.u==0)[1] #pick first z.u==0 to use

              #propose new calls
              calls.cand <- model$calls[i]+1

              #simulate a new u around this s
              u.cand <- c(rnorm(1,model$s[i,1],model$sigma.u[1]), rnorm(1,model$s[i,2],model$sigma.u[1]))

              #update z.u
              z.u.cand <- 1

              #update D, mu from simulated u
              D.cand <- GetD(u=u.cand, X=model$X[1:J,1:2], z.u=z.u.cand)
              mu.cand <- GetMu(beta0=model$beta0[1], beta1=model$beta1[1], D=D.cand,z.u=z.u.cand)

              #get proposed logprobs
              ll.call.cand <- dpois(calls.cand,model$lambda.C[1],log=TRUE)
              ll.dB.cand <- ddB(x=model$dB[l.idx,1:J],mu=mu.cand, zeros=model$zeros[l.idx,1:J], sigma.d=model$sigma.d[1],
                                mindB=mindB,z.u=z.u.cand,log=TRUE)

              #u logProb not used for MH update
              logProb.proposed <- ll.call.cand + ll.dB.cand
              logProb.initial <- ll.call[i] + ll.dB[l.idx]

              log_MH_ratio <- (logProb.proposed) - (logProb.initial)
              accept <- decide(log_MH_ratio)

              if(accept) {
                model$calls[i] <<- calls.cand
                #turn these on
                model$z.u[l.idx] <<- z.u.cand
                model$D[l.idx,1:J] <<- D.cand
                model$mu[l.idx,1:J] <<- mu.cand
                model$u[l.idx,1:2] <<- u.cand
                #update log likelihoods
                ll.call[i] <- ll.call.cand
                ll.dB[l.idx] <- ll.dB.cand
                # #update callID
                model$callID[l.idx] <<- i
              }
            }else{
              print("Warning: Not enough z.u's left to raise an individaul's call number. Raise M.u if this happens past burnin.")
            }
          }
        }
      }
    }
    #hand back to nimble
    model$N.call[1] <<- sum(model$z.u[1:M.u])
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#custom update for u so that they are only updated if z.u=1
uSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    scale <- control$scale
    M.u <- control$M.u
    mindB <- control$mindB
    J <- control$J
    calcNodes <- model$getDependencies("u","dB")
  },
  run = function() {
    ll.u <-  model$logProb_u[,1]
    ll.dB <-  model$logProb_dB[,1]
    callID <- model$callID
    s <- model$s
    for(l in 1:M.u){
      if(model$z.u[l]==1){ #only update z.u=1 calls. Let them leave state space
        u.cand <- c(rnorm(1,model$u[l,1],scale), rnorm(1,model$u[l,2],scale))
        #update D, mu
        D.cand <- GetD(u=u.cand, X=model$X[1:J,1:2], z.u=1)
        mu.cand <- GetMu(beta0=model$beta0[1], beta1=model$beta1[1], D=D.cand,z.u=1)
        
        #get proposed logprobs
        ll.u.cand <- dnorm_censor(u.cand,s[callID[l],1:2],sd=model$sigma.u[1],z.u = 1,log=TRUE)
        ll.dB.cand <- ddB(x=model$dB[l,1:J],mu=mu.cand, zeros=model$zeros[l,1:J], sigma.d=model$sigma.d[1],
                          mindB=mindB,z.u=1,log=TRUE)
      
        logProb.proposed <- ll.u.cand + ll.dB.cand
        logProb.initial <- ll.u[l] + ll.dB[l]
        
        log_MH_ratio <- (logProb.proposed) - (logProb.initial)
        accept <- decide(log_MH_ratio)
        
        if(accept) {
          model$u[l,1:2] <<- u.cand
          model$D[l,1:J] <<- D.cand
          model$mu[l,1:J] <<- mu.cand
          ll.u[l] <- ll.u.cand
          ll.dB[l] <- ll.dB.cand
        }
      }
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#custom update for s that only considers the z.u=1 calls
sSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    scale <- control$scale
    xlim <- control$xlim
    ylim <- control$ylim
    M.s <- control$M.s
    calcNodes <- model$getDependencies("u")
  },
  run = function() {
    ll.u <-  model$logProb_u[,1]
    ll.u.cand <- ll.u
    for(i in 1:M.s){
      #if z=0 or if z=1 and calls=0 we propose from prior.
      if(model$z.s[i]==0|model$calls[i]==0){#propose from uniform prior
        model$s[i, 1:2] <<- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
      }else{#MH
        s.cand <- c(rnorm(1,model$s[i,1],scale), rnorm(1,model$s[i,2],scale))
        inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
        if(inbox){
          #find z.u=1 calls
          these.calls <- which(model$callID==i&model$z.u==1)
          n.these.calls <- length(these.calls)
          if(n.these.calls>1){
            #calculate u likelihoods
            for(l in 1:n.these.calls){
              ll.u.cand[these.calls[l]] <- dnorm_censor(model$u[these.calls[l],1:2],s.cand[1:2],
                                                        sd=model$sigma.u[1],z.u = 1,log=TRUE)
            }
            logProb.proposed <- sum(ll.u.cand[these.calls])
            logProb.initial <- sum(ll.u[these.calls])
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              model$s[i,1:2] <<- s.cand
              ll.u[these.calls] <- ll.u.cand[these.calls] #dont really need to update, not used again
            }
          }else{
            #calculate u likelihoods
            ll.u.cand[these.calls[1]] <- dnorm_censor(model$u[these.calls[1],1:2],s.cand[1:2],
                                                      sd=model$sigma.u[1],z.u = 1,log=TRUE)
            logProb.proposed <- ll.u.cand[these.calls[1]]
            logProb.initial <- ll.u[these.calls[1]]
            log_MH_ratio <- (logProb.proposed) - (logProb.initial)
            accept <- decide(log_MH_ratio)
            if(accept) {
              model$s[i,1:2] <<- s.cand
              ll.u[these.calls[1]] <-  ll.u.cand[these.calls[1]] #dont really need to update, not used again
            }
          }
        }
      }
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)