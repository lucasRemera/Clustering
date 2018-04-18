##Note : utilisation sur le wiki

###################################
# Potthoff-Whittinghill statistic #
###################################
PWstat=function(cas,temoins){
  P=sum(cas)/sum(cas+temoins)
  obs=cas+temoins
  expected=P*obs
  
  PW=sum(expected)*sum((cas*(cas-1)/expected))
  mPW=sum(cas)*(sum(cas)-1)
  vPW=2*(length(cas)-1)*mPW
  
  sPW=(PW-mPW)/sqrt(vPW)
  print(sPW)
  return(2*(1-pnorm(abs(sPW))))
}

replicatePW=function(ns=100,ra1=0.2,ra2=0.2,unif=FALSE){
  obs=sample(corine$obs,ns,TRUE)
  if(unif) cas=rbinom(ns,obs,runif(ns,ra1,ra2))
  else cas=rbinom(ns,obs,c(rep(ra1,ns/2),rep(ra2,ns/2)))
  return(PWstat(cas,obs-cas))
}



################################
# delta statistic              #
################################


kPPV=function(x0,y0,xref,yref, k=3,seuil=Inf,nmin=0,check=TRUE){
  d=sqrt((x0-xref)**2+(y0-yref)**2)
  i=order(d)[1:k]
  if(min(d)==0 & check) i=order(d)[2:(k+1)]
  
  if(min(d)>seuil){
    i=NULL
  } 
  else{
    i=i[which(d[i]<seuil)]
    
  }
  return(i)
}

deltaStat=function(data=corine,k=5,comparison=c("knn","random"),retour=c("pvalue","statistic","sum")){
  P0=sum(data$cas)/sum(data$obs)
  if(comparison[1]=="knn") deltaQ=sapply(1:nrow(data),function(i) sum((data$p[i] - data$p[kPPV(data$x[i],data$y[i],data$x,data$y,k = k)])**2/
                                                (P0*(1-P0)*(1/data$obs[i]+1/data$obs[kPPV(data$x[i],data$y[i],data$x,data$y,k=k)] ))   ))
  else{
    deltaQ=sapply(1:nrow(data),function(i){ 
      cc=sample(c(1:nrow(data))[-i],k)
      sum((data$p[i] - data$p[cc])**2/(P0*(1-P0)*(1/data$obs[i]+1/data$obs[cc] ))   )})

  }
  if(retour[1]=="sum") return(sum(deltaQ))
  else if(retour[1]=="statistic") return(deltaQ)
  else return(1-pchisq(sum(deltaQ),k*nrow(data)))
}

simDelta=function(p1=NA,p2=p1,pourc2=0, nsim=100,k=5,data=corine,retour=c("pvalue","statistic","sum")){
  P0=sum(data$cas)/(sum(data$temoins)+sum(data$cas))
  if(is.na(p1)) p1=P0
  pv=sapply(1:nsim,function(j){
    #pv=c()
    #for(j in 1:nsim){
    print(j)
    data0=data
    #print(data0)
    if((p1-p2)*pourc2==0){
      data0$cas=rbinom(nrow(data0),data0$obs,p1)
    } 
    else{
      ii=rep(0,nrow(data0))
      i=sample(1:nrow(data0),floor(pourc2*nrow(data0)))
      ii[i]=1
      data0$cas=rbinom(nrow(data0),data0$obs,p1+ii*(p2-p1))
    }
    data0$temoins=data0$obs-data0$cas
    data0$p=data0$cas/(data0$obs)
    P0=sum(data0$cas)/sum(data0$obs)
    deltaQ=sapply(1:nrow(data0),function(i) sum((data0$p[i] - data0$p[kPPV(data0$x[i],data0$y[i],data0$x,data0$y,k = k)])**2/
                                                  (P0*(1-P0)*(1/data0$obs[i]+1/data0$obs[kPPV(data0$x[i],data0$y[i],data0$x,data0$y,k=k)] ))   ))
    # hist(deltaQ,freq=F)
    # xx=seq(0,60,by=0.01)
    # lines(xx,dchisq(xx,df = k),col="green")
    #pv=c(pv,1-pchisq(sum(deltaQ),k*nrow(data0)))
    if(retour[1]=="sum") return(sum(deltaQ))
    else if(retour[1]=="statistic") return(deltaQ)
    else return(1-pchisq(sum(deltaQ),k*nrow(data0)))
    #}
  })
  return(pv)
}


