##détails sur le wiki
#construit une matrice symétrique de similarité, basée sur la p-value du test de Fisher
#library(deldir)
getGNeighbourVoronoiMatrix=function(x,y,cas,tem,fisher=TRUE){
  m=matrix(0,ncol=length(x),nrow=length(x))
  ddir=deldir(x,y)
  #print(ddir)
  v1=ddir$dirsgs$ind1
  v2=ddir$dirsgs$ind2
  for(i in 1:nrow(ddir$dirsgs)){
    if(fisher) pp=fisher.test(matrix(c(sum(tem[v1[i]]),sum(cas[v1[i]]),sum(tem[v2[i]]),sum(cas[v2[i]])),ncol=2))$p.value #case and controls of adjacents vertices could be from the same distribution?
    else pp=1-abs(log((cas[v1[i]]/(cas[v1[i]]+tem[v1[i]]))/(cas[v2[i]]/(cas[v2[i]]+tem[v2[i]]))))
    #else pp=1-abs((cas[v1[i]]/(cas[v1[i]]+tem[v1[i]]))-(cas[v2[i]]/(cas[v2[i]]+tem[v2[i]])))
    m[v1[i],v2[i]]=pp
    m[v2[i],v1[i]]=pp
  }
  if(!fisher) m=m-min(m)
  return(m)
}

#effectue le partitionnement spectral
spectralPartitionnement=function(m,k=2,norme=FALSE,scaled=FALSE){
  N=dim(m)[1]
  Di=diag(apply(m,1,sum)) #matrice des degrés
  L=Di-m #laplacienne
  if(norme) L=solve(Di)%*%L
  #Di=diag(1/sqrt(apply(m,1,sum)))
  #L=(Di)%*%m%*%Di
  #vp=eigen(L)$values
  VP=eigen(L)$vectors #Vecteurs propres en colonne
  VP2=VP[,(N-k+1):N] #on récupère les k derniers VP
  if(scaled) VP2=scale(VP2)
  return(kmeans(VP2,k)$cluster)
}

###########
# exemple #
###########
x=runif(800)
y=runif(800)
inCluster=((x>.3)&(x<.6)&(y>.3)&(y<.6))
obs=rpois(800,500)
cas=rbinom(800,obs,0.2+as.numeric(inCluster)*0.15)
m=getGNeighbourVoronoiMatrix(x,y,cas,obs-cas)
cl=spectralPartitionnement(m)
plot(x,y,col=rainbow(2)[cl])#si tout marche bien, cluster au milieu
