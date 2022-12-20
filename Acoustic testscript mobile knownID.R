library(nimble)
library(coda)
source("sim.Acoustic.ss.R")
source("build.cluster.R")
source("NimbleModel Acoustic SCR Mobile knownID.R")
source("init.data.mobile.knownID.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE) 

#Simulation parameters

N=40#Abundance
#Make an ARU array
X=expand.grid(0:10,0:10) #regular grid
# X=build.cluster(ntraps=64,clusterdim=4,spacingin=1,spacingout=4,plotit = TRUE) #cluster grid
buff=3 #ARU buffer around ARU array to define rectangular state space
lambda.C=10 #mean # calls/ind (Poisson call rate assumption)
sigma.u=0.25 #BVN individual movement parameter. This sampler is for no movement, so sigma.u should be set to 0.
#data plot below doesn't illustrate movement if you simulate it.

#attenuation function parameters
beta0=70 #dB at call source
beta1=-15 #slope for expected dB as function of distance from call source
sigma.d=2 #attenuation function error variance
mindB=60 #dB detection threshold (cannot detect call with received dB lower than this)

#simulate data
data=sim.Acoustic.ss(N=N,beta0=beta0,beta1=beta1,sigma.d=sigma.d,sigma.u=sigma.u,mindB=mindB,X=X,
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
plot(data$s,xlim=data$xlim,ylim=data$ylim,pch=16,main=NA,xlab="X",ylab="Y")
points(data$u,col="#332288",pch=16)
points(data$u.obs,col="#6699CC",pch=16)
for(i in 1:data$n.call){
  lines(x=c(data$u[i,1],data$s[data$clusterID[i],1]),y=c(data$u[i,2],data$s[data$clusterID[i],2]))
}
points(X,pch=4)
points(data$s,col="#F0E442",pch=16)


#mean detections/detected call (how many ARUs registered call)
mean(rowSums(!is.na(data$dB.obs)))
#mean detected calls/parent
mean(table(data$clusterID))
#detected individuals
data$n.ind.obs  

#Process data for nimble
dummydB=99999 #dummy number for call-ARU nondetection events. NA's throw nimble warnings.
dB=data$dB.obs
dB[is.na(dB)]=dummydB
callID=data$clusterID.obs
zeros=dB==dummydB #zeros indicates the nondetection events
n.calls=nrow(dB)
J=nrow(X)

M.s=75 #individual data augmentation level. Raise this if N.ind posterior ever hits M.s.
M.u=750 #call data augmentation level. Raise this if N.call posterior ever hits M.u
#Really, it is more nuanced than that. You can never let N.call get close enough to M.u
#that a new individual cannot be turned on with it's expected number of calls. I have 
#nimble print warnings if you need to raise M.u.


#initialize some data structures and latent variables
inits=list(psi=0.5) #requires a psi init. Not really important to give this init to nimble.
nimbuild=init.data(data=data,inits=inits,M.s=M.s,M.u=M.u)

#build z.s.data. Every individual with a detected call gets a "1", all else NA to be estimated
z.s.data=rep(NA,M.s)
z.s.data[unique(callID)]=1

#build z.u.data. Every detected call gets a "1", all else NA to be estimated
z.u.data=rep(NA,M.u)
z.u.data[1:length(callID)]=1

#Augment zeros and dB
zeros.extra=matrix(rep(TRUE,J*(M.u-n.calls)),nrow=M.u-n.calls,byrow=TRUE)
zeros=rbind(zeros,zeros.extra)

dB.extra=matrix(rep(dummydB,J*(M.u-n.calls)),nrow=M.u-n.calls,byrow=TRUE)
dB=rbind(dB,dB.extra)

#which calls are detected vs. augmented
calls.detected=1*(rowSums(dB!=dummydB)>0)

#ballpark inits for attenuation function parameters - Check attenuation priors to make sure they are appropriate!
beta0.init=max(data$dB.obs,na.rm=TRUE) #maximum observed received dB is a good starting value
beta1.init=rnorm(1,beta1,1) #can let nimble draw this from prior, but may not converge if nowhere near truth.
sigma.d.init=rnorm(1,sigma.d,0.1)
lambda.C.init=mean(nimbuild$calls)
sigma.u.init=abs(rnorm(1,sigma.u,0.10))

#supply stuff to nimble
Niminits <- list(z.s=nimbuild$z.s,z.u=nimbuild$z.u,calls=nimbuild$calls,s=nimbuild$s,
                 u=nimbuild$u,calls=nimbuild$calls,beta0=beta0.init,beta1=beta1.init,
                 sigma.d=sigma.d.init,lambda.C=lambda.C.init,sigma.u=sigma.u.init)
constants<-list(M.s=M.s,M.u=M.u,xlim=data$xlim,ylim=data$ylim,J=J,n.calls=n.calls,zeros=zeros,
                mindB=mindB,X=as.matrix(X))
Nimdata<-list(dB=dB,z.s=z.s.data,z.u=z.u.data,calls=rep(NA,M.s),callID=nimbuild$callID)

# set parameters to monitor
parameters<-c('beta0','beta1','sigma.d','lambda.C','N.ind','N.call','psi','sigma.u')
#You must set the tuning parameters for s and u, so you should record them without thinning
parameters2=c("s","u","calls")

start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants,data=Nimdata,check=FALSE,inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1, monitors2=parameters2,thin2=1,useConjugacy = TRUE)

#Several required sampler replacements/additions
#call sampler has tuning parameter "call.ups", which determines how many times you propose to
#add/subtract calls to each detected individual per iteration. In general, if lambda.c is larger
#call.ups should be larger.
conf$removeSampler("calls")
conf$addSampler(target = paste("calls[",1,":",M.s,"]", sep=""),
                type = 'callSampler',
                control=list(calls.detected=calls.detected,J=J,call.ups=5,
                             M.s=M.s,M.u=M.u,mindB=mindB),silent = TRUE)

#Individual z.s sampler. This turns entire clusters on/off
conf$removeSampler("z.s")
conf$addSampler(target = paste("z.s[",1,":",M.s,"]", sep=""),
                type = 'zSampler',
                control=list(calls.detected=calls.detected,J=J,
                             M.s=M.s,M.u=M.u,mindB=mindB),silent = TRUE)

#u sampler. Must tune the scale parameter.
conf$removeSampler(paste("u[",1,":",M.u,"]", sep=""))
conf$addSampler(target = paste("u[",1,":",M.u,"]", sep=""),
                  type = 'uSampler',control=list(M.u=M.u,J=J,mindB=mindB,scale=0.2),silent = TRUE)

#s sampler. Must tune the scale parameter
conf$removeSampler(paste("s[",1,":",M.s,"]", sep=""))
conf$addSampler(target = paste("s[",1,":",M.s,"]", sep=""),
                type = 'sSampler',control=list(M.s=M.s,xlim=data$xlim,ylim=data$ylim,scale=0.15),silent = TRUE)

#Optional sampler replacement
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
start.time2<-Sys.time()
Cmcmc$run(1500,reset=FALSE) #short run for demonstration. Can run again to continue sampling where it stopped.
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples=as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$n.call #True number of calls


#if you monitor s (and don't thin), you can calculate the acceptance rate. Should remove some burnin
#after looking at some s posteriors
mvSamples2=as.matrix(Cmcmc$mvSamples2)

#check s and u acceptance rates
burnin=500 #discard some burnin

#we want to look at the minimum acceptance rates since there is only 1 tuning parameter for all u's
#pull out u indicies
idx=grep("u",colnames(mvSamples2))
#for u, only the detected calls are always updated. All other acceptance rates can be ignored.
idx=idx[1:data$n.call]
min(1-rejectionRate(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx]))) #shoot for 0.2 or so

#pull out s indices
idx=grep("\\bs\\b",colnames(mvSamples2)) #using exact match in grep here
#for s, only individuals with detected calls are always updated using MH. 
#All other acceptance rates can be ignored.
idx=idx[1:data$n.ind.obs]
min(1-rejectionRate(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx]))) #shoot for 0.2 or so

#posterior for each individual's true number of calls
idx=grep("calls",colnames(mvSamples2))
plot(mcmc(mvSamples2[5:nrow(mvSamples2),idx]))
