NimModel <- nimbleCode({
  ###priors###
  lambda.C ~ dunif(0,100) #mean calls/individual
  #attenuation function parameters
  beta0 ~ dunif(50,100)
  beta1 ~ dunif(-100,0) #assuming must be negative (dB decays with distance)
  sigma.d ~ dunif(0,10)
  #data augmentation prior
  psi~dunif(0,1)
  
  ###likelihoods (except for s priors)###
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #call model
    calls[i] ~ dpois(lambda.C) #number of total calls per individual
    calls.zero[i] <- calls[i]-calls.detected[i] #number of undetected calls per individual
    #detection model
    D[i,1:J] <- GetD(s=s[i,1:2], X=X[1:J,1:2], z=z[i])
    mu[i,1:J] <- GetMu(beta0=beta0, beta1=beta1, D=D[i,1:J], z=z[i])
    #unobserved calls
    dB.unobs[i,1:J] ~ dNoDetect(z=z[i],calls.zero=calls.zero[i],
                                mu = mu[i,1:J], sigma.d=sigma.d, mindB=mindB)
  }
  #observed calls
  for(l in 1:n.calls){
    dB[l,1:J] ~ ddB(mu=mu[callID[l],1:J], zeros=zeros[l,1:J], sigma.d=sigma.d, mindB=mindB)
  }
  N.ind <- sum(z[1:M]) #individual abundance
  N.call <- sum(calls[1:M]*z[1:M]) #call abundance
})# end model

GetD <- nimbleFunction(
  run = function(s = double(1), X=double(2),z=double(0)){
    returnType(double(1))
    J <- nimDim(X)[1]
    if(z==0){
      D <- rep(0,J)
    }else{
      D <- sqrt((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
    }
    return(D)
  }
)

GetMu <- nimbleFunction(
  run = function(beta0=double(0),beta1=double(0), D = double(1),z=double(0)){
    returnType(double(1))
    J <- nimDim(D)[1]
    if(z==0){
      mu <- rep(0,J)
    }else{
      mu <- beta0 + beta1*D[1:J]
    }
    return(mu)
  }
)

ddB <- nimbleFunction(
  run = function(x = double(1),mu = double(1), zeros = double(1),sigma.d = double(0),
                 mindB = double(0), log = integer(0)) {
    returnType(double(0))
    J <- nimDim(x)[1]
    logProb <- 0
    for(j in 1:J){
      if(zeros[j]){ #below dB threshold
        logProb <- logProb + pnorm(mindB, mu[j],sigma.d,log=TRUE)
      }else{ #above dB threshold
        logProb <- logProb + dnorm(x[j],mu[j],sigma.d,log=TRUE)
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
                 mindB = double(0)) {
    returnType(double(1))
    J <- nimDim(mu)[1]
    return(rep(0,J))
  }
)

dNoDetect <- nimbleFunction(
  run = function(x = double(1), z = double(0), calls.zero = double(0),
                 mu = double(1),mindB = double(0), sigma.d = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      logProb <- 0
    }else{
      if(calls.zero<0){#keeps calls >= calls.detected
        logProb <- -Inf
      }else if(calls.zero==0){ #0 undetected calls for this individual
        logProb <- 0
      }else{ #calls.zero undetected calls for this individual
        J <- nimDim(x)[1]
        logProb.tmp <- 0
        for(j in 1:J){
          logProb.tmp <- logProb.tmp + pnorm(mindB, mu[j],sigma.d,log=TRUE)
        }
        logProb <- calls.zero*logProb.tmp
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
rNoDetect <- nimbleFunction(
  run = function(n=integer(0), z = double(0), calls.zero = double(0),
                 mu = double(1),mindB = double(0), sigma.d = double(0)) {
    returnType(double(1))
    J <- nimDim(mu)[1]
    return(rep(0,J))
  }
)

#Required custom update for number of calls
callSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    i <- control$i
    calls.detected <- control$calls.detected
    call.ups <- control$call.ups
  },
  run = function() {
    if(model$z[i]==0){ #if z is off, propose from prior
      calls.cand=rpois(1,model$lambda.C[1])
      model$calls[i] <<- calls.cand
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
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
        reject=FALSE #we reject if 1) proposed counts <0 or 2) you select a detected call
        if(updown==0){#subtract
          calls.cand <- model$calls[i]-1
          if(model$calls[i]==0){
            reject <- TRUE
          }else{
            tmp <- rcat(1,rep(1/model$calls[i],model$calls[i])) #select one of the calls
            if(tmp<=calls.detected){ #is it one of the detected calls?
              reject <- TRUE #if so, we reject
            }
          }
        }else{#add
          calls.cand <- model$calls[i]+1
        }
        if(!reject){
          model.lp.initial <- model$getLogProb(calcNodes)
          model$calls[i] <<- calls.cand
          model.lp.proposed <- model$calculate(calcNodes)#proposed logProb
          log_MH_ratio <- (model.lp.proposed) - (model.lp.initial)
          accept <- decide(log_MH_ratio)
          if(accept) {
            copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
          } else {
            copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
          }
        }
      }
    }
  },
  methods = list( reset = function () {} )
)