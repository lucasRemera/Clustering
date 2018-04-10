########################
# functions to return neighborhood as adjacent matrix
########################

#require deldir
adjFromTesselation=function(x,y){
  delaunay=deldir(x,y)$delsgs
  adj=matrix(0,nrow=length(x),ncol=length(x))
  
  for(k in 1:nrow(delaunay)){
    adj[delaunay[k,5],delaunay[k,6]]=1
    adj[delaunay[k,6],delaunay[k,5]]=1
  }
  return(adj)
  
}

adjFromDiscus=function(x,y,r){
  distM=as.matrix(dist(cbind(x,y),upper=T))<r
  distM=apply(distM,1,as.numeric)
  diag(distM)=0
  return(distM)
}

adjFromKNN=function(x,y,k=5){
  adj=matrix(0,nrow=length(x),ncol=length(x))
  for( i in 1:length(x)){
    vois=kPPV(x[i],y[i],x,y,k)
    
    adj[vois,i]=1
  }
  return(adj)
}

##############
# normalization functions
##############

getMod=function(x){
  xx=as.numeric(x==max(x))
  return(xx/sum(xx))
}

myNorm=function(M,mod=FALSE){
  if(mod) return(apply(M, 2, getMod))
  else return(apply(M,2,function(cc) cc/sum(cc)))
}

########################
# Likelihood functions
########################

vraisemblanceBinom=function(cas,obs,theta,llog=FALSE){
  if(!llog) ll=choose(obs,cas)*theta**cas*(1-theta)**(obs-cas)
  else ll=log(choose(obs,cas))+cas*log(theta)+(obs-cas)*log(1-theta)
  return(ll)
}

vraisemblanceVoisin=function(i,adj,etat,etatToTest,freq,correctFreq=T,theta=NULL,pow=1){
  #sum(etat[adj[,i]]==etatToTest)/sum(adj[,i])
  if(!is.null(theta)){
    w=(1-abs(theta[etatToTest]-theta))**pow
    l=sum(w%*%etat[,which(adj[,i]==1)])/freq[etatToTest]
  }
  else if(correctFreq) l=sum(etat[etatToTest,which(adj[,i]==1)])/freq[etatToTest]
  else l=sum(etat[etatToTest,which(adj[,i]==1)])
  return(l)
}


updateFreq=function(etat){
  freq=apply(etat,1,sum)
  freq=freq/sum(freq)
  return(freq)
}

# toMax=function(theta2,etat2=etat,cas=corine$cas,obs=corine$obs,freq2=freq,minus=-1){
#   yop=sapply(1:nrow(adj),function(iii){
#     sapply(1:K, function(kk){
#       vraisemblanceBinom(cas[iii],obs[iii],theta2[kk],TRUE)#+log(vraisemblanceVoisin(iii,adj,etat2,kk,freq2))
#     })
#   })
#   return(minus*sum(apply(yop*etat,1,sum))) #minus to transform maximisation into minimisation
# }


getInformationCriterion=function(theta2,cas=corine$cas,obs=corine$obs,etat2=etat,K=length(theta2),result=c("AIC","BIC","LIKELIHOOD")){
  yop=sapply(1:ncol(etat2),function(iii){
    sapply(1:K, function(kk){
      vraisemblanceBinom(cas[iii],obs[iii],theta2[kk],TRUE)
    })
  })
  vrais=(sum(apply(yop*etat2,1,sum))) 
  if(result[1]=="AIC") return(2*K-2*vrais)
  else if(result[1]=="BIC") return(log(ncol(etat2))*K-2*vrais)
  else return(vrais)
}


getIC=function(etat,cas=corine$cas,obs=corine$obs,precision=0.001,stat=TRUE,alpha=0.05){
  ic=function(k){
    vth=seq(precision,1-precision,by=precision)
    lth=rep(0,length(vth))
    for(i in 1:ncol(etat)){
      lth=lth+sapply(vth, function(th){
        etat[k,i]*vraisemblanceBinom(corine$cas[i],corine$obs[i],th,llog = TRUE)#+log(vraisemblanceVoisin(i,adj,etat,k,freq,correctFreq = T))
      })
    }
    if(!stat) return(10**lth/sum(10**lth))
    else{
      elth=cumsum(10**lth/sum(10**lth))
      p1=vth[which(elth>(alpha/2))[1]]-precision/2
      m=vth[which.max(lth)]
      p3=vth[which(elth>(1-(alpha/2)))[1]]-precision/2
      return(c(p1,m,p3))
    }
  }
  
  tor=sapply(1:nrow(etat), ic)
  if(!stat) return(cbind(seq(precision,1-precision,by=precision),tor))
  else return(tor)
}

EM=function(data=corine,theta=c(0.1,0.2),graphV=c("del","knn"),k.knn=5,mod=FALSE,pdsV=1,t=1:20,trace.iter=TRUE){
  
  toMax0=function(theta2,etat2=etat,cas=data$cas,obs=data$obs,freq2=freq,minus=-1){
    yop=sapply(1:nrow(adj),function(iii){
      sapply(1:K, function(kk){
        vraisemblanceBinom(cas[iii],obs[iii],theta2[kk],TRUE)#+log(vraisemblanceVoisin(iii,adj,etat2,kk,freq2))
      })
    })
    return(minus*sum(apply(yop*etat,1,sum))) #minus to transform maximisation into minimisation
  }
  
  ### INITIALISATION ###
  K=length(theta) #number of clusters
  
  if(graphV[1]=="knn") adj=adjFromKNN(data$x,data$y,k.knn)
  else adj=adjFromTesselation(data$x,data$y) #neighborhood definition
  
  posterior0=sapply(1:nrow(adj), function(i){
    sapply(1:length(theta), function(th){
      vraisemblanceBinom(data$cas[i],data$obs[i],theta[th])
    })
  }) #vector pi_i, init proba of being in each cluster
  
  if(is.null(dim(posterior0))) posterior0=t(as.matrix(posterior0))
  posterior0=myNorm(posterior0,mod)
  if(is.null(dim(posterior0))) posterior0=t(as.matrix(posterior0))
  freq=updateFreq(posterior0) #vector pi
  etat=posterior0
  
  for(lll in t){
    if(trace.iter) print(lll)
    etat=sapply(1:nrow(adj),function(i){
      sapply(1:K, function(k){
        freq[k]*vraisemblanceBinom(data$cas[i],data$obs[i],theta[k])*vraisemblanceVoisin(i,adj,etat,k,freq,correctFreq = T)**pdsV
      })
    }) #etape E
    if(is.null(dim(etat))) etat=t(as.matrix(etat))
    #etat=apply(etat, 2, function(cc) cc/sum(cc))
    etat=myNorm(etat,mod)
    if(is.null(dim(etat))) etat=t(as.matrix(etat))
    #print(etat)
    freq=updateFreq(etat)
    theta=optim(theta,toMax0)$par
    #print(freq)
    #print(theta)
    
  }
  
  return(list(theta=theta,pi_i=etat,pi_freq=freq))
}


EMArret=function(data=corine,theta=c(0.1,0.2),graphV=c("del","knn"),k.knn=5,mod=FALSE,pdsV=1,arret=0.01,trace.iter=TRUE){
  
  toMax0=function(theta2,etat2=etat,cas=data$cas,obs=data$obs,freq2=freq,minus=-1){
    yop=sapply(1:nrow(adj),function(iii){
      sapply(1:K, function(kk){
        vraisemblanceBinom(cas[iii],obs[iii],theta2[kk],TRUE)#+log(vraisemblanceVoisin(iii,adj,etat2,kk,freq2))
      })
    })
    return(minus*sum(apply(yop*etat,1,sum))) #minus to transform maximisation into minimisation
  }
  
  ### INITIALISATION ###
  K=length(theta) #number of clusters
  
  if(graphV[1]=="knn") adj=adjFromKNN(data$x,data$y,k.knn)
  else adj=adjFromTesselation(data$x,data$y) #neighborhood definition
  
  posterior0=sapply(1:nrow(adj), function(i){
    sapply(1:length(theta), function(th){
      vraisemblanceBinom(data$cas[i],data$obs[i],theta[th])
    })
  }) #vector pi_i, init proba of being in each cluster
  
  if(is.null(dim(posterior0))) posterior0=t(as.matrix(posterior0))
  posterior0=myNorm(posterior0,mod)
  if(is.null(dim(posterior0))) posterior0=t(as.matrix(posterior0))
  freq=updateFreq(posterior0) #vector pi
  etat=posterior0
  lll=1
  d=Inf
  while(d>arret){
    if(trace.iter) print(lll)
    print(theta)
    etat=sapply(1:nrow(adj),function(i){
      sapply(1:K, function(k){
        freq[k]*vraisemblanceBinom(data$cas[i],data$obs[i],theta[k])*vraisemblanceVoisin(i,adj,etat,k,freq,correctFreq = T)**pdsV
      })
    }) #etape E
    if(is.null(dim(etat))) etat=t(as.matrix(etat))
    #etat=apply(etat, 2, function(cc) cc/sum(cc))
    etat=myNorm(etat,mod)
    if(is.null(dim(etat))) etat=t(as.matrix(etat))
    #print(etat)
    freq=updateFreq(etat)
    theta2=optim(theta,toMax0)$par
    d=sqrt(sum((theta-theta2)**2))
    theta=theta2
    lll=lll+1
  }
  
  return(list(theta=theta,pi_i=etat,pi_freq=freq))
}

EMAuto=function(data=corine,K=2,nsim=5,graphV=c("del","knn"),k.knn=5,mod=FALSE,pdsV=1,arret=0.01,trace.iter=TRUE){
  best=-Inf
  l0=list()
  prange=data$cas/data$obs
  for(i in 1:nsim){
    l=EMArret(data,sort(runif(K,min(prange),max(prange))),graphV,k.knn,mod,pdsV,arret,trace.iter)
    print("END EM")
    likelihood=getInformationCriterion(l$theta,data$cas,data$obs,l$pi_i,K,"LIKELIHOOD" )
    if(likelihood>best){
      best=likelihood
      l0=l
    }
  }
  return(l0)
  
}

