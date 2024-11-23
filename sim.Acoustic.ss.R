e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.Acoustic.ss <-
  function(N=NA,beta0=NA,beta1=NA,sigma.d=NA,sigma.u=NA,mindB=NA,X=NA,buff=NA,lambda.C=NA){
    # simulate a population of activity centers
    xlim <- c(min(X[,1])-buff,max(X[,1])+buff)
    ylim <- c(min(X[,2])-buff,max(X[,2])+buff)
    s <- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    J <- nrow(X)
    # simulate number of calls. Must call at least 1 time.
    n.calls <- rpois(N,lambda.C)
    #simulate call locations
    u <- matrix(NA,nrow=sum(n.calls),ncol=2)
    clusterID <- matrix(NA,nrow=sum(n.calls),ncol=2)
    idx <- 1
    for(i in 1:N){
      if(n.calls[i]>0){
        for(j in 1:n.calls[i]){
          u[idx,] <- cbind(rnorm(1,s[i,1],sigma.u),rnorm(1,s[i,2],sigma.u))
          clusterID[idx,1] <- i
          clusterID[idx,2] <- j
          idx <- idx+1
        }
      }
    }
    #distances from all u to all ARUs
    D <- e2dist(u,X)
    #received decibel level
    dB <- beta0 + beta1*D + rnorm(prod(dim(D)), 0, sigma.d)
    dB.obs <- dB
    dB.obs[dB.obs<mindB] <- NA
    # How many detections for each call?
    n.loc.obs<- apply(!is.na(dB.obs),1,sum)
    # Matrix of observations
    dB.obs <- dB.obs[n.loc.obs>0,]
    u.obs <- u[n.loc.obs>0,]
    clusterID.obs <- clusterID[n.loc.obs>0,]
    n.obs <- sum(n.loc.obs>0)#how many detected calls/children
    
    #renumber clusterID.obs starting at 1
    clusterID.obs <- clusterID.obs[,1]
    these.guys <- unique(clusterID.obs)
    clusterID.obs2 <- rep(NA,n.obs)
    for(l in 1:n.obs){
      clusterID.obs2[l] <- which(these.guys==clusterID.obs[l])
    }
    s.obs <- s[these.guys,]
    n.calls.obs <- n.calls[these.guys]
    n.ind.obs <- length(these.guys) #how many detected parents/individuals
    
    out <- list(dB.obs=dB.obs,dB=dB,X=X,xlim=xlim,ylim=ylim,buff=buff,u=u,u.obs=u.obs,s=s,n.obs=n.obs,s.obs=s.obs,n.calls.obs=n.calls.obs,
              n.loc.obs=n.loc.obs,n.ind.obs=n.ind.obs,clusterID.obs=clusterID.obs2,clusterID=clusterID[,1],mindB=mindB,
              n.calls=n.calls,n.call=sum(n.calls))
    return(out)
  }
