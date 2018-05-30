satscan1D=function(casIn,casOut,TIn,TOut,MC=FALSE){
  CAS=casIn+casOut
  TTot=TIn+TOut
  if(MC){
    casIn=rpois(1,TIn*CAS/TTot)
    casOut=CAS-casIn 
  } 
  if(casIn>1) lnk=sum(log(2:casIn))
  else lnk=0
  if(casOut>1) lnkout=sum(log(2:casOut))
  else lnkout=0
  LL0=casIn*log(CAS*TIn/TTot)-lnk-CAS*TIn/TTot
  L1=0
  L2=0
  if(casIn>0) L1=casIn*log(casIn)-lnk-casIn
  if(casOut>0) L2=casOut*log(casOut)-lnkout-casOut
  LL1=L1+L2
  #if((casIn/TIn)>(casOut/TOut)) return(LL1-LL0)
  #else return(0)
  return(LL1-LL0)
}


satScanPV=function(casIn,casOut,TIn,TOut,nsim=1000){
  scanTRUE=satscan1D(casIn,casOut,TIn,TOut)
  scanMC=replicate(nsim,{satscan1D(casIn,casOut,TIn,TOut,MC = TRUE)})
  return(1-sum(scanTRUE>scanMC)/nsim)
}

findScan1D=function(Tevent,t.begin=0,t.end=max(Tevent)){
  cas.tot=length(Tevent)
  t.tot=t.end-t.begin
  scan1D=function(intervalle){
    a=intervalle[1]
    b=intervalle[2]
    if(a<0) return(Inf)
    if(b>(max(Tevent)+1)) return(Inf)
    if(a>=b) return(Inf)
    cas.in=sum(Tevent>a&Tevent<b)
    cas.out=cas.tot-cas.in
    t.in=b-a
    t.out=t.tot-t.in
    return(-satscan1D(cas.in,cas.out,t.in,t.out))
  }
  opt=optim(t.begin+0.25*(t.end-t.begin)*c(1,3),scan1D)
  return(opt$par)
}

############
# examples #
############

dumb.t=seq(0,100)
in.cluster=as.numeric(dumb.t>40&dumb.t<60)
dumb.y=rpois(length(dumb.t),10+10*in.cluster)

satscan1D(sum(dumb.y*in.cluster),sum(dumb.y*(1-in.cluster)),sum(in.cluster),sum(1-in.cluster))
satscan1D(sum(dumb.y*in.cluster),sum(dumb.y*(1-in.cluster)),sum(in.cluster),sum(1-in.cluster),MC = TRUE)

satScanPV(sum(dumb.y*in.cluster),sum(dumb.y*(1-in.cluster)),sum(in.cluster),sum(1-in.cluster))

findScan1D(rep(dumb.t,dumb.y))
