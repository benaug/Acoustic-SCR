build.cluster=function(ntraps,clusterdim,spacingin,spacingout,plotit = F){
  size=ntraps/(clusterdim^2)
  nfullcluster=floor(size)
  coords=seq(0,spacingin*(clusterdim-1),by=spacingin)
  cluster=as.matrix(expand.grid(coords,coords),ncol=2)
  
  #Store all full clusters
  storecluster=list()
  for(l in 1:nfullcluster){
    storecluster[[l]]=cluster
  }
  #Are there traps that don't form a full cluster?
  remainder=ntraps-nfullcluster*clusterdim^2
  ncluster=nfullcluster
  if(remainder>0){
    #build remainder cluster
    x=y=0
    m=1
    X2=matrix(0,ncol=2,nrow=0)
    while(m<(remainder+1)){
      xy=getxy(x,y,m,spacingin)
      x=xy[1]
      y=xy[2]
      X2=rbind(X2,xy)
      m=m+1
    }
    #Put remainder cluster at the end
    X2=matrix(X2,ncol=2)
    storecluster=c(storecluster,list(X2))
    ncluster=ncluster+1
  }
  space=clusterdim*spacingin+spacingout
  X=matrix(0,nrow=0,ncol=2)
  x=y=0
  m=1
  
  #Build clusters
  while(m<(ncluster+1)){
    xy=getxy(x,y,m,space)
    x=xy[1]
    y=xy[2]
    newcluster=storecluster[[m]]
    newcluster[,1]=newcluster[,1]+x
    newcluster[,2]=newcluster[,2]+y
    X=rbind(X,newcluster)
    m=m+1
  }
  if(plotit){
    xlim=c(min(X[,1]-space),max(X[,1]+space))
    ylim=c(min(X[,2]-space),max(X[,2]+space))
    plot(X,xlim=xlim,ylim=ylim)
  }
  return(X)
}



getxy=function(x,y,m,space){
  #2x2
  if(m==2) y=y+space
  if(m==3) x=x+space
  if(m==4) y=y-space
  #3x3
  if(m==5){
    y=2*space
    x=0
  }
  if(m==6) x=x+space
  if(m==7) x=x+space
  if(m==8) y=y-space
  if(m==9) y=y-space
  #4x4
  if(m==10){
    y=3*space
    x=0
  }
  if(m==11) x=x+space
  if(m==12) x=x+space
  if(m==13) x=x+space
  if(m==14) y=y-space
  if(m==15) y=y-space
  if(m==16) y=y-space
  #5x5
  if(m==17){
    y=4*space
    x=0
  }
  if(m==18) x=x+space
  if(m==19) x=x+space
  if(m==20) x=x+space
  if(m==21) x=x+space
  if(m==22) y=y-space
  if(m==23) y=y-space
  if(m==24) y=y-space
  if(m==25) y=y-space
  #6x6
  if(m==26){
    y=5*space
    x=0
  }
  if(m==27) x=x+space
  if(m==28) x=x+space
  if(m==29) x=x+space
  if(m==30) x=x+space
  if(m==31) x=x+space
  if(m==32) y=y-space
  if(m==33) y=y-space
  if(m==34) y=y-space
  if(m==35) y=y-space
  if(m==36) y=y-space
  
  #7x7
  if(m==37){
    y=6*space
    x=0
  }
  if(m==38) x=x+space
  if(m==39) x=x+space
  if(m==40) x=x+space
  if(m==41) x=x+space
  if(m==42) x=x+space
  if(m==43) x=x+space
  if(m==44) y=y-space
  if(m==45) y=y-space
  if(m==46) y=y-space
  if(m==47) y=y-space
  if(m==48) y=y-space
  if(m==49) y=y-space
  
  #8x8
  if(m==50){
    y=7*space
    x=0
  }
  if(m==51) x=x+space
  if(m==52) x=x+space
  if(m==53) x=x+space
  if(m==54) x=x+space
  if(m==55) x=x+space
  if(m==56) x=x+space
  if(m==57) x=x+space
  if(m==58) y=y-space
  if(m==59) y=y-space
  if(m==60) y=y-space
  if(m==61) y=y-space
  if(m==62) y=y-space
  if(m==63) y=y-space
  if(m==64) y=y-space
  
  #9x9
  if(m==65){
    y=8*space
    x=0
  }
  if(m==66) x=x+space
  if(m==67) x=x+space
  if(m==68) x=x+space
  if(m==69) x=x+space
  if(m==70) x=x+space
  if(m==71) x=x+space
  if(m==72) x=x+space
  if(m==73) x=x+space
  if(m==74) y=y-space
  if(m==75) y=y-space
  if(m==76) y=y-space
  if(m==77) y=y-space
  if(m==78) y=y-space
  if(m==79) y=y-space
  if(m==80) y=y-space
  if(m==81) y=y-space
  
  #10x10
  if(m==82){
    y=9*space
    x=0
  }
  if(m==83) x=x+space
  if(m==84) x=x+space
  if(m==85) x=x+space
  if(m==86) x=x+space
  if(m==87) x=x+space
  if(m==88) x=x+space
  if(m==89) x=x+space
  if(m==90) x=x+space
  if(m==91) x=x+space
  if(m==92) y=y-space
  if(m==93) y=y-space
  if(m==94) y=y-space
  if(m==95) y=y-space
  if(m==96) y=y-space
  if(m==97) y=y-space
  if(m==98) y=y-space
  if(m==99) y=y-space
  if(m==100) y=y-space
  
  #11x11
  if(m==101){
    y=10*spacey
    x=0
  }
  if(m==102) x=x+space
  if(m==103) x=x+space
  if(m==104) x=x+space
  if(m==105) x=x+space
  if(m==106) x=x+space
  if(m==107) x=x+space
  if(m==108) x=x+space
  if(m==109) x=x+space
  if(m==110) x=x+space
  if(m==111) x=x+space
  if(m==112) y=y-space
  if(m==113) y=y-space
  if(m==114) y=y-space
  if(m==115) y=y-space
  if(m==116) y=y-space
  if(m==117) y=y-space
  if(m==118) y=y-space
  if(m==119) y=y-space
  if(m==120) y=y-space
  if(m==121) y=y-space
  
  #12x12
  if(m==122){
    y=11*space
    x=0
  }
  if(m==123) x=x+space
  if(m==124) x=x+space
  if(m==125) x=x+space
  if(m==126) x=x+space
  if(m==127) x=x+space
  if(m==128) x=x+space
  if(m==129) x=x+space
  if(m==130) x=x+space
  if(m==131) x=x+space
  if(m==132) x=x+space
  if(m==133) x=x+space
  if(m==134) y=y-space
  if(m==135) y=y-space
  if(m==136) y=y-space
  if(m==137) y=y-space
  if(m==138) y=y-space
  if(m==139) y=y-space
  if(m==140) y=y-space
  if(m==141) y=y-space
  if(m==142) y=y-space
  if(m==143) y=y-space
  if(m==144) y=y-space
  return(c(x,y))
}

build.cluster.alt=function(ntraps,ntrapsC,spacingin,spacingout,plotit = F){
  size=ntraps/ntrapsC
  nfullcluster=floor(size)
  #Build normal clusters
  #build remainder cluster
  x=y=0
  m=1
  cluster=matrix(0,ncol=2,nrow=0)
  while(m<(ntrapsC+1)){
    xy=getxy(x,y,m,spacingin)
    x=xy[1]
    y=xy[2]
    cluster=rbind(cluster,xy)
    m=m+1
  }
  
  #Store all full clusters
  storecluster=list()
  for(l in 1:nfullcluster){
    storecluster[[l]]=cluster
  }
  #Are there traps that don't form a full cluster?
  remainder=ntraps-nfullcluster*ntrapsC
  ncluster=nfullcluster
  if(remainder>0){
    #build remainder cluster
    x=y=0
    m=1
    X2=matrix(0,ncol=2,nrow=0)
    while(m<(remainder+1)){
      xy=getxy(x,y,m,spacingin)
      x=xy[1]
      y=xy[2]
      X2=rbind(X2,xy)
      m=m+1
    }
    #Put remainder cluster at the end
    X2=matrix(X2,ncol=2)
    storecluster=c(storecluster,list(X2))
    ncluster=ncluster+1
  }
  spacex=max(cluster[,1])+spacingout
  spacey=max(cluster[,2])+spacingout
  X=matrix(0,nrow=0,ncol=2)
  x=y=0
  m=1
  
  #Build clusters
  while(m<(ncluster+1)){
    xy=getxy.alt(x,y,m,spacex,spacey)
    x=xy[1]
    y=xy[2]
    newcluster=storecluster[[m]]
    newcluster[,1]=newcluster[,1]+x
    newcluster[,2]=newcluster[,2]+y
    X=rbind(X,newcluster)
    m=m+1
  }
  if(plotit){
    xlim=c(min(X[,1]-spacex),max(X[,1]+spacex))
    ylim=c(min(X[,2]-spacey),max(X[,2]+spacey))
    plot(X,xlim=xlim,ylim=ylim)
  }
  return(X)
}



getxy.alt=function(x,y,m,spacex,spacey){
  #2x2
  if(m==2) y=y+spacey
  if(m==3) x=x+spacex
  if(m==4) y=y-spacey
  #3x3
  if(m==5){
    y=2*spacey
    x=0
  }
  if(m==6) x=x+spacex
  if(m==7) x=x+spacex
  if(m==8) y=y-spacey
  if(m==9) y=y-spacey
  #4x4
  if(m==10){
    y=3*spacey
    x=0
  }
  if(m==11) x=x+spacex
  if(m==12) x=x+spacex
  if(m==13) x=x+spacex
  if(m==14) y=y-spacey
  if(m==15) y=y-spacey
  if(m==16) y=y-spacey
  #5x5
  if(m==17){
    y=4*spacey
    x=0
  }
  if(m==18) x=x+spacex
  if(m==19) x=x+spacex
  if(m==20) x=x+spacex
  if(m==21) x=x+spacex
  if(m==22) y=y-spacey
  if(m==23) y=y-spacey
  if(m==24) y=y-spacey
  if(m==25) y=y-spacey
  #6x6
  if(m==26){
    y=5*spacey
    x=0
  }
  if(m==27) x=x+spacex
  if(m==28) x=x+spacex
  if(m==29) x=x+spacex
  if(m==30) x=x+spacex
  if(m==31) x=x+spacex
  if(m==32) y=y-spacey
  if(m==33) y=y-spacey
  if(m==34) y=y-spacey
  if(m==35) y=y-spacey
  if(m==36) y=y-spacey
  
  #7x7
  if(m==37){
    y=6*spacey
    x=0
  }
  if(m==38) x=x+spacex
  if(m==39) x=x+spacex
  if(m==40) x=x+spacex
  if(m==41) x=x+spacex
  if(m==42) x=x+spacex
  if(m==43) x=x+spacex
  if(m==44) y=y-spacey
  if(m==45) y=y-spacey
  if(m==46) y=y-spacey
  if(m==47) y=y-spacey
  if(m==48) y=y-spacey
  if(m==49) y=y-spacey
  
  #8x8
  if(m==50){
    y=7*spacey
    x=0
  }
  if(m==51) x=x+spacex
  if(m==52) x=x+spacex
  if(m==53) x=x+spacex
  if(m==54) x=x+spacex
  if(m==55) x=x+spacex
  if(m==56) x=x+spacex
  if(m==57) x=x+spacex
  if(m==58) y=y-spacey
  if(m==59) y=y-spacey
  if(m==60) y=y-spacey
  if(m==61) y=y-spacey
  if(m==62) y=y-spacey
  if(m==63) y=y-spacey
  if(m==64) y=y-spacey
  
  #9x9
  if(m==65){
    y=8*spacey
    x=0
  }
  if(m==66) x=x+spacex
  if(m==67) x=x+spacex
  if(m==68) x=x+spacex
  if(m==69) x=x+spacex
  if(m==70) x=x+spacex
  if(m==71) x=x+spacex
  if(m==72) x=x+spacex
  if(m==73) x=x+spacex
  if(m==74) y=y-spacey
  if(m==75) y=y-spacey
  if(m==76) y=y-spacey
  if(m==77) y=y-spacey
  if(m==78) y=y-spacey
  if(m==79) y=y-spacey
  if(m==80) y=y-spacey
  if(m==81) y=y-spacey
  
  #10x10
  if(m==82){
    y=9*spacey
    x=0
  }
  if(m==83) x=x+spacex
  if(m==84) x=x+spacex
  if(m==85) x=x+spacex
  if(m==86) x=x+spacex
  if(m==87) x=x+spacex
  if(m==88) x=x+spacex
  if(m==89) x=x+spacex
  if(m==90) x=x+spacex
  if(m==91) x=x+spacex
  if(m==92) y=y-spacey
  if(m==93) y=y-spacey
  if(m==94) y=y-spacey
  if(m==95) y=y-spacey
  if(m==96) y=y-spacey
  if(m==97) y=y-spacey
  if(m==98) y=y-spacey
  if(m==99) y=y-spacey
  if(m==100) y=y-spacey
  
  #11x11
  if(m==101){
    y=10*spacey
    x=0
  }
  if(m==102) x=x+spacex
  if(m==103) x=x+spacex
  if(m==104) x=x+spacex
  if(m==105) x=x+spacex
  if(m==106) x=x+spacex
  if(m==107) x=x+spacex
  if(m==108) x=x+spacex
  if(m==109) x=x+spacex
  if(m==110) x=x+spacex
  if(m==111) x=x+spacex
  if(m==112) y=y-spacey
  if(m==113) y=y-spacey
  if(m==114) y=y-spacey
  if(m==115) y=y-spacey
  if(m==116) y=y-spacey
  if(m==117) y=y-spacey
  if(m==118) y=y-spacey
  if(m==119) y=y-spacey
  if(m==120) y=y-spacey
  if(m==121) y=y-spacey
  
  #12x12
  if(m==122){
    y=11*spacey
    x=0
  }
  if(m==123) x=x+spacex
  if(m==124) x=x+spacex
  if(m==125) x=x+spacex
  if(m==126) x=x+spacex
  if(m==127) x=x+spacex
  if(m==128) x=x+spacex
  if(m==129) x=x+spacex
  if(m==130) x=x+spacex
  if(m==131) x=x+spacex
  if(m==132) x=x+spacex
  if(m==133) x=x+spacex
  if(m==134) y=y-spacey
  if(m==135) y=y-spacey
  if(m==136) y=y-spacey
  if(m==137) y=y-spacey
  if(m==138) y=y-spacey
  if(m==139) y=y-spacey
  if(m==140) y=y-spacey
  if(m==141) y=y-spacey
  if(m==142) y=y-spacey
  if(m==143) y=y-spacey
  if(m==144) y=y-spacey
  return(c(x,y))
}