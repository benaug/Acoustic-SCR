init.data=function(data=NA,M.s=NA,M.u=NA,inits=NA){
  #pull out objects used to initialize data
  xlim <- data$xlim
  ylim <- data$ylim
  dB.obs <- data$dB.obs
  callID.obs <- data$clusterID.obs
  X <- as.matrix(data$X)
  #pull out inits
  psi <- inits$psi
  
  #initialize u and s
  u <- cbind(runif(M.u,xlim[1],xlim[2]),runif(M.u,ylim[1],ylim[2]))
  #Optimize starting locations given where they are heard
  idx <- which(rowSums(!is.na(dB.obs))>0) #switch for those actually caught
  n.u <- length(idx)
  
  for(i in 1:n.u){
    detections <- !is.na(dB.obs[i,])
    trps <-  matrix(X[detections,1:2],ncol=2,byrow=FALSE)
    weights <- dB.obs[i,detections]/sum(dB.obs[i,detections]) #weight by received dB
    if(nrow(trps)>1){
      u[i,] <- c(sum(trps[,1]*weights),sum(trps[,2]*weights))
    }else{
      u[i,] <- trps+rnorm(2,0,0.1) #put some noise in here so they are not directly on top of a trap
    }
  }
  s <- cbind(runif(M.s,xlim[1],xlim[2]),runif(M.s,ylim[1],ylim[2]))
  callID <- rep(NA,M.u)
  if(M.u<n.u)stop("More observed calls than M.u, raise M.u.")
  callID[1:n.u] <- callID.obs
  #then need to adjust s given known ID u inits
  capID <- unique(callID[1:n.u])
  if(max(capID)>M.s)stop("More detected inidivuals than M.s, raise M.s.")
  for(i in capID){
    u.tmp <- u[which(callID==i),]
    if(is.matrix(u.tmp)){
      s[i,] <- colMeans(u.tmp)
    }else{
      s[i,] <- u.tmp
    }
  }
  callID[(n.u+1):M.u] <- sample(1:M.s,M.u-n.u,replace=TRUE)
  z.u <- c(rep(1,n.u),rep(0,M.u-n.u))
  z.s <- rep(0,M.s)
  z.s[unique(callID[1:n.u])] <- 1
  Ntarget <- psi*M.s
  add <- Ntarget-sum(z.s)
  if(add>0){
    idx <- which(z.s==0)
    if(length(idx)==1){
      z.s[idx] <- 1
    }else{
      idx <- sample(idx,add)
      z.s[idx] <- 1
    }
  }
  calls <- rep(0,M.s)
  for(i in 1:M.s){
    if(z.s[i]==1){
      calls[i] <- sum(callID==i&z.u==1)
    }
  }
  out <- list(z.s=z.s,z.u=z.u,callID=callID,s=s,u=u,calls=calls)
  return(out)
}