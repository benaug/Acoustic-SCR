library(nimble)
library(coda)
source("sim.Acoustic.ss.R")
source("build.cluster.R")
source("NimbleModel Acoustic SCR Stationary knownID.R")

#Simulation parameters

N <- 40 #Abundance
#Make an ARU array
X <- expand.grid(0:10,0:10) #regular grid
# X <- build.cluster(ntraps=64,clusterdim=4,spacingin=1,spacingout=4,plotit = TRUE) #cluster grid
buff <- 3 #ARU buffer around ARU array to define rectangular state space
lambda.C <- 10 #mean # calls/ind (Poisson call rate assumption, modifying requires modifying custom callSampler)
sigma.u <- 0 #BVN individual movement parameter. This sampler is for no movement, so sigma.u should be set to 0.
#data plot below doesn't illustrate movement if you simulate it.

#attenuation function parameters
beta0 <- 70 #dB at call source
beta1 <- -15 #slope for expected dB as function of distance from call source
sigma.d <- 2 #attenuation function error variance
mindB <- 60 #dB detection threshold (cannot detect call with received dB lower than this)

#simulate data
data <- sim.Acoustic.ss(N=N,beta0=beta0,beta1=beta1,sigma.d=sigma.d,sigma.u=sigma.u,mindB=mindB,X=X,
                     buff=buff,lambda.C=lambda.C)

#The observed data:
#1) the call by ARU detection history with NA for nondetections and a dB level for detections
data$dB.obs[1,] #call 1
#2) the individual IDs of each detected call
data$clusterID.obs
#Every call has an individual ID
nrow(data$dB.obs)==length(data$clusterID.obs)

#plot data
par(mfrow=c(1,1),ask=FALSE)
plot(NA,xlim=data$xlim,ylim=data$ylim,xlab="X",ylab="Y")
y2D <- apply(!is.na(data$dB.obs),c(1,2),sum)
cap <- which(rowSums(y2D)>0)
for(i in cap){
  trapcap <- which(y2D[i,]>0)
  for(j in trapcap){
    lines(x=c(data$s.obs[data$clusterID.obs[i],1],X[j,1]),y=c(data$s.obs[data$clusterID.obs[i],2],X[j,2]),lty=1,col="#999999",lwd=2)
  }
}
points(X,pch=4)
points(data$s,pch=16,col="#56B4E9",cex=2)
text(data$s[,1],data$s[,2],label=as.character(data$n.calls),cex=0.75)

#mean detections/detected call (how many ARUs registered call)
mean(rowSums(!is.na(data$dB.obs)))
#mean detected calls/parent
mean(table(data$clusterID))
#detected individuals
data$n.ind.obs  

#Process data for nimble
dummydB <- 99999 #dummy number for call-ARU nondetection events. NA's throw nimble warnings.
dB <- data$dB.obs
dB[is.na(dB)] <- dummydB
callID <- data$clusterID.obs
zeros <- dB==dummydB #zeros indicates the nondetection events
n.calls <- nrow(dB)
J <- nrow(X)

M <- 75 #data augmentation level. Raise this if N posterior ever hits M.

#build z.data. Every individual with a detected call gets a "1", all else NA to be estimated
z.data <- rep(NA,M)
z.data[unique(callID)] <- 1

#build counts.detected, detected counts per individual
calls.detected <- rep(M,0) 
for(i in 1:M){
  calls.detected[i] <- sum(callID==i)
}

#inits for nimble.
z.init <- z.data
z.init[is.na(z.init)] <- 0
#initialize the true number of calls per individual
#just adding 1 to the number of detected calls here. 
#Only important to make sure calls.init>=calls.detected
calls.init <- calls.detected+1

#initialize activity centers
s.init <- cbind(runif(M,data$xlim[1],data$xlim[2]),runif(M,data$ylim[1],data$ylim[2])) #start with uniform distribution
#update for detected individuals
y2D <- apply(!is.na(data$dB.obs),c(1,2),sum)
for(i in 1:data$n.ind.obs){
  these.calls <- which(callID==i)
  if(length(these.calls)>0){
    if(length(these.calls)>1){
      trapcaps <- which(colSums(dB[these.calls,]<dummydB)>0)
    }else{
      trapcaps <- which(dB[these.calls,]<dummydB)
    }
    if(length(trapcaps)>1){
      s.init[i,] <- colMeans(X[trapcaps,])
    }else{
      s.init[i,] <- as.numeric(X[trapcaps,])
    }
  }
}

#ballpark inits for attenuation function parameters - Check attenuation priors to make sure they are appropriate!
beta0.init <- max(data$dB.obs,na.rm=TRUE) #maximum observed received dB is a good starting value
beta1.init <- rnorm(1,beta1,1) #can let nimble draw these from prior, but may not converge if nowhere near truth.
sigma.d.init <- rnorm(1,sigma.d,0.1)
lambda.C.init <- mean(calls.init)

#supply stuff to nimble
Niminits <- list(z=z.init,calls=calls.init,s=s.init,beta0=beta0.init,beta1=beta1.init,
                 sigma.d=sigma.d.init,lambda.C=lambda.C.init,calls=calls.init)
constants <- list(M=M,xlim=data$xlim,ylim=data$ylim,J=J,n.calls=n.calls,zeros=zeros,callID=callID,mindB=mindB,
                calls.detected=calls.detected,X=as.matrix(X))
Nimdata <- list(dB=dB,dB.unobs=matrix(TRUE,M,J),z=z.data,calls=rep(NA,M))

# set parameters to monitor
parameters <- c('beta0','beta1','sigma.d','lambda.C','N.ind','N.call','psi')
parameters2 <- c("s") #monitor other things. with (possibly) different thinning rate...

start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1, monitors2=parameters2,thin2=1,useConjugacy = TRUE)

#One required sampler replacement. The nimble-selected slice sampler does not work correctly.
#Note, there is a tuning parameter there, "call.ups". The custom Metropolis-Hastings update
#only updates 1 call/individual at a time. "call.ups" determines how many times to do this per
#iteration. The only general suggestion I have is that "call.ups" should be larger when lambda.C is larger.
#if lambda.C mixing is poor, try raising "call.ups".
conf$removeSampler("calls")
for(i in 1:M){
  conf$addSampler(target = paste("calls[",i,"]", sep=""),
                  type = 'callSampler',
                  control=list(i=i,calls.detected=calls.detected[i],call.ups=5),silent = TRUE)
}


#Optional sampler replacements
#update x and y dimensions of s jointly.
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
}
#beta0 and beta1 posteriors are highly correlated
#RW_block is another option. AF slice mixes better, but is slower. Not sure which is more computationally efficient
conf$removeSampler(c("beta0","beta1"))
conf$addSampler(target = c("beta0","beta1"),
                type = 'AF_slice',control = list(adaptive=TRUE),silent = TRUE)
# conf$addSampler(target = c("beta0","beta1"),
#                 type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #short run for demonstration. Can run again to continue sampling where it stopped.
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

data$n.call #True number of calls


#if you monitor s (and don't thin), you can calculate the acceptance rate. Should remove some burnin
#after looking at some s posteriors
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
burnin <- 1000
1-rejectionRate(mcmc(mvSamples2[burnin:nrow(mvSamples2),]))

