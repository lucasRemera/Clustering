vraisemblance=function(idx,data=corine,coef=1,MC=FALSE){
  CAS=sum(data$cas)
  TEMOINS=sum(data$temoins)
  casIN=sum(data$cas[idx])
  temIN=sum(data$temoins[idx])
  if(MC){
    P=CAS/(CAS+TEMOINS)
    obsIN=casIN+temIN
    casIN=rbinom(1,obsIN,P)
    temIN=obsIN-casIN
  }
  
  
  
  casOUT=CAS-casIN
  temOUT=TEMOINS-temIN
  #print(paste(casIN,temIN))
  if((casIN+temIN)*(casOUT+temOUT)==0) return(Inf)
  #L0=(CAS/(CAS+TEMOINS))**CAS * (1-CAS/(CAS+TEMOINS))**(TEMOINS)
  LL0=CAS*(log(CAS)-log((CAS+TEMOINS))) + TEMOINS*(log(TEMOINS)-log((CAS+TEMOINS)))
  # 
  # if((casIN/(casIN+temIN))>(casOUT/(casOUT+temOUT))){
  #   L1=(casIN/(casIN+temIN))**casIN * (1-casIN/(casIN+temIN))**(temIN)
  #   LL1=casIN*(log(casIN)-log((casIN+temIN))) + temIN*(log(temIN)-log((casIN+temIN)))
  #   L2=(casOUT/(casOUT+temOUT))**casOUT * (1-casOUT/(casOUT+temOUT))**(temOUT)
  #   LL2=casOUT*(log(casOUT)-log((casOUT+temOUT))) + temOUT*(log(temOUT)-log((casOUT+temOUT)))
  # }
  # else{
  #   L1=L0
  #   L2=1
  #   LL1=LL0
  #   LL2=0
  # }
  if(TRUE){
    #L1=(casIN/(casIN+temIN))**casIN * (1-casIN/(casIN+temIN))**(temIN)
    if(casIN==0) L1=0
    else L1=casIN*(log(casIN)-log((casIN+temIN)))
    if(temIN==0) L2=0
    else L2=temIN*(log(temIN)-log((casIN+temIN)))
    LL1=L1+L2
    
    if(casOUT==0) L1=0
    else L1=casOUT*(log(casOUT)-log((casOUT+temOUT)))
    if(temOUT==0) L2=0
    else L2=temOUT*(log(temOUT)-log((casOUT+temOUT)))
    LL2=L1+L2
    
    #LL1=casIN*(log(casIN)-log((casIN+temIN))) + temIN*(log(temIN)-log((casIN+temIN)))
    #L2=(casOUT/(casOUT+temOUT))**casOUT * (1-casOUT/(casOUT+temOUT))**(temOUT)
    LL2=casOUT*(log(casOUT)-log((casOUT+temOUT))) + temOUT*(log(temOUT)-log((casOUT+temOUT)))
  }
  else{
    #L1=L0
    #L2=1
    LL1=LL0
    LL2=0
  }
  #print(paste(LL0,LL1,LL2))
  if((casIN/temIN)<(CAS/TEMOINS)) return(LL0)
  else return(LL1+LL2)
  #return(coef*(LL1+LL2-LL0))
  
}




vraisemblancePV=function(idx,data=corine,nsim=1000){
  L=vraisemblance(idx,data)
  Lsim=replicate(nsim,{vraisemblance(idx,data,MC=TRUE)})
  return(sum(Lsim>L)/nsim)
}


getProb=function(card,minn=3,maxx=Inf,opt=NULL,p=0.5){
  if(is.infinite(maxx)){
    if(card<=minn) return(0)
    else return(p)
  }
  else{
    if(card>=maxx) return(1)
    else{
      if(is.null(opt)) return((card-minn)/(maxx-minn))
      else{
        if(card>opt) return((1-p)*(card-opt)/(maxx-opt)+p)
        else return(p*(card-minn)/(opt-minn))
      }
    } 
    
  }
}


decreaseT=function(t,d=1,toLog=FALSE){
  if(toLog) return(1/(log(t+1)**d))
  else return(1/(t**d))
}

searchByConvex=function(n0=3,idx0=NULL,nsim=3000,d=1,toLog=FALSE,minn=3,maxx=Inf,opt=NULL,p=0.5,data=corine,trace=TRUE,pdel=1,llr=TRUE){
  if(is.null(idx0)){
    
    idx=sample(1:nrow(data),n0)
    
  }
  else idx=idx0
  X=data$x
  Y=data$y
  ECB=idx[chull(X[idx],Y[idx])]
  polyx=X[c(ECB,ECB[1])]
  polyy=Y[c(ECB,ECB[1])]
  EC=point.in.polygon(X,Y,polyx,polyy)
  EC=which(EC==1)
  l=list()
  tt=1:nsim
  for(t in tt){
    nE=length(EC)+length(ECB)
    l[[t]]=list(n=nE,EC=EC,ECB=ECB,L=L0)
    pu=runif(1,0,1)
    pp=getProb(nE,minn,maxx,opt,p)
    L0=vraisemblance(idx=c(ECB,EC),data =data)
    if(pu<pp){
      #if(TRUE){
      print("1") #ne pas toujours enlever un point
      if(length(ECB)>3 ){
        if(runif(1,0,1)<pdel){ #a verifier
          ECB1=ECB[-sample(1:length(ECB),1)]
          ECB1=c(ECB1,sample(EC,1)) 
          ECB1=ECB1[chull(X[ECB1],Y[ECB1])]
        }
        
        else ECB1=ECB[-sample(1:length(ECB),1)]
      }
      else{
        
        ECB1=ECB[-sample(1:length(ECB),1)]
        ECB1=c(ECB1,sample(EC,1))
      }
      
      polyx=X[c(ECB1,ECB1[1])]
      polyy=Y[c(ECB1,ECB1[1])]
      
      EC1=point.in.polygon(X,Y,polyx,polyy)
      EC1=which(EC1==1)
    }else{
      print("2")
      EC1=c(EC,ECB,sample(which(!1:nrow(data)%in%c(EC,ECB)),1))
      #print(EC1)
      ECB1=EC1[chull(X[EC1],Y[EC1])]
      #print(ECB1)
      polyx=X[c(ECB1,ECB1[1])]
      polyy=Y[c(ECB1,ECB1[1])]
      
      EC1=point.in.polygon(X,Y,polyx,polyy)
      EC1=which(EC1==1)
    }
    L1=vraisemblance(idx=c(ECB1,EC1),data =data)
    if(llr){
      if(L1>L0){
        EC=EC1
        ECB=ECB1
      }
      else{
        pp=exp(L1-L0)*decreaseT(t,d = d)
        pu=runif(1,0,1)
        if(pu<pp){
          ECB=ECB1
          EC=EC1
        }
      }
    }
    else{
      ECB=ECB1
      EC=EC1
    }
    #t=t+1
  }
  if(trace) return(l)
  else return(l[[nsim]])
  
}

