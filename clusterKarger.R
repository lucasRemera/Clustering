##Clustering on Spatial data of coordinates @arg x,y to find cases aggregates in epidemiology (with @arg cas cases and @arg tem controls)
#Neighborhood defined by @arg voisin, voronoi tesseltation or by discus graph (with @arg radius r)
#probabilist algorithm if @arg random=TRUE, else determinist 
#stop criterion @arg alpha, the threshold in fisher test
## This algorithm makes cluster on a neighborhood graph, based on the Karger's min-cut algorithm
## Contiguous vertices/components are aggregated with a probability=fisher's exact test p-value (similarity score here),
##until each paired are too dissimilar
composanteConnexeKarger=function(x,y,cas,tem,r=1,alpha=0.05,voisin=c("radius","voronoi"),random=TRUE){
  ## INIT: neighborhood graph construction
  if(voisin[1]=="radius") liste=getGNeighbour(x,y,cas,tem,r)
  else liste=d=getGNeighbourVoronoi(x,y,cas,tem)
  
  # edge weighted by fisher's exact test p-value
  pv=sapply(liste, function(i) i$p) 
  while(length(which(pv>alpha))>0){ #if no adjacent vertices/clusters have equivalent prevalence, stop
    if(length(which(pv>alpha))>1){
      pv[which(pv<alpha)]=0 #don't group vertices if their p-value is under alpha risk

      weight=pv/sum(pv)
      if(random) idx=sample(1:length(pv),1,prob=weight) #select a random edge, weighted by the p-values
      else idx=which.max(weight)
    }
    else idx=which(pv>alpha)
    V1=liste[[idx]]$v1
    V2=liste[[idx]]$v2
    newcas=liste[[idx]]$c1+liste[[idx]]$c2
    newtem=liste[[idx]]$t1+liste[[idx]]$t2
    newpoint=c(V1,V2) #creation of a new vertex, union of two vertices, case and controls associated are summed
    liste[[idx]]=list(v1=NA,v2=NA,t1=NA,t2=NA,p=0) #delete the selected edge

    lv1=lapply(liste, function(i){ if(!is.na(i)){ i$v1}  else{NA}} )

    lv2=lapply(liste, function(i) { if(!is.na(i)){ i$v2}  else{NA}})
    toChange=which((sapply(lv1, function(j) prod(j==V1)| prod(j==V2) ) + sapply(lv2, function(j) prod(j==V1)| prod(j==V2) ))>0) #select vertices adjacent to the new vertices (one of the previous vertices)
    for(tc in toChange){
      if((prod(liste[[tc]]$v1==V1) + prod(liste[[tc]]$v1==V2)) & (prod(liste[[tc]]$v2==V1) + prod(liste[[tc]]$v2==V2))){
        liste[[tc]]=list(v1=NA,v2=NA,t1=NA,t2=NA,p=0) #delete the selected edge duplicates
      }else{
        if(prod(liste[[tc]]$v1==V1) + prod(liste[[tc]]$v1==V2)  ){ #replace vertices
          liste[[tc]]$v1=newpoint
          liste[[tc]]$c1=newcas
          liste[[tc]]$t1=newtem
          liste[[tc]]$p=fisher.test(matrix(c(newtem,newcas,liste[[tc]]$t2,liste[[tc]]$c2),ncol=2))$p.value #we have to renew edges and reset p-value
        }
        if(prod(liste[[tc]]$v2==V1) + prod(liste[[tc]]$v2==V2)){
          liste[[tc]]$v2=newpoint
          liste[[tc]]$c2=newcas
          liste[[tc]]$t2=newtem

          liste[[tc]]$p=fisher.test(matrix(c(newtem,newcas,liste[[tc]]$t1,liste[[tc]]$c1),ncol=2))$p.value
          #maybe don't do again fisher test, but keep sum/mean of last fisher's p-value
        }
        if(liste[[tc]]$v1==liste[[tc]]$v2){
          liste[[tc]]=list(v1=NA,v2=NA,t1=NA,t2=NA,p=0)
        }
      }
      
    }
    pv=sapply(liste, function(i) i$p)
  }
  #delete duplicates
  return(unique(c(sapply(liste[which(sapply(liste, function(i) prod(!is.na(i$v1)))==1)],function(j) j$v1),sapply(liste[which(sapply(liste, function(i) prod(!is.na(i$v2)))==1)],function(j) j$v2))))
  
  
}

#Construct a neighborhood graph (each pair of points with a distance< @arg r are linked)
getGNeighbour=function(x,y,cas,tem,r){
  liste=list()
  #frontiere=matrix(0,ncol=7,nrow=0)
  for(i in 1:(length(x)-1)){ #create an adjacence list
    for(j in (i+1):length(x)){
      #if(sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)<=r) frontiere=rbind(frontiere,matrix(c(i,j,cas[i],tem[i],cas[j],tem[j],fisher.test(matrix(c(sum(tem[i]),sum(cas[i]),sum(tem[j]),sum(cas[j])),ncol=7))$p.value),ncol=7))
      if(sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)<=r){
        if(cas[i]+tem[i]+cas[j]+tem[j]==0) pp=0
        else{ 
          pp=fisher.test(matrix(c(sum(tem[i]),sum(cas[i]),sum(tem[j]),sum(cas[j])),ncol=2))$p.value #case and controls of adjacents vertices could be from the same distribution?
        }
        distance=sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
        liste=c(liste,list(list(v1=i,v2=j,c1=cas[i],t1=tem[i],c2=cas[j],t2=tem[j],p=pp,d=distance)))
      } 
    }
  }
  return(liste)
}

#Construct a neighborhood graph with Voronoi tesselation
getGNeighbourVoronoi=function(x,y,cas,tem){
  liste=list()
  ddir=deldir(x,y)
  #print(ddir)
  v1=ddir$dirsgs$ind1
  v2=ddir$dirsgs$ind2
  for(i in 1:nrow(ddir$dirsgs)){
    pp=fisher.test(matrix(c(sum(tem[v1[i]]),sum(cas[v1[i]]),sum(tem[v2[i]]),sum(cas[v2[i]])),ncol=2))$p.value #case and controls of adjacents vertices could be from the same distribution?
    
    liste=c(liste,list(list(v1=v1[i],v2=v2[i],c1=cas[v1[i]],t1=tem[v1[i]],c2=cas[v2[i]],t2=tem[v2[i]],p=pp)))
    
  }
  return(liste)
}

#tranform main algorithm return into a vector of cluster numbers
listeToCLusters=function(liste,sorted=TRUE){
  cl=unlist(sapply(1:length(liste), function(i) rep(i, length(liste[[i]]))))
  names(cl)=unlist(liste)
  if(sorted) cl=cl[order(as.numeric(names(cl)))]
  return(cl)
}

#get prevalence of each clusters
clusterToPrev=function(cas,tem,clusters,idx=1:length(clusters), table=FALSE){
  toRet=clusters
  for(i in unique(clusters)){
    who=idx[which(clusters==i)]
    cc=sum(cas[who])
    tt=sum(tem[who])
    toRet[which(clusters==i)]=cc/(cc+tt)
  }
  
  if(table){
    ii=which(!duplicated(clusters))
    #toRet=unique(toRet)
    toRet=toRet[ii]
  } 
  return(toRet)
}


###############################
#   Examples                  #
###############################
karger=composanteConnexeKarger2(corine$x,corine$y,corine$cas,corine$temoins,voisin = "voronoi",alpha=0.02)
cluster.karger=listeToCLusters(karger2)
prev.karger=clusterToPrev(corine$cas,corine$temoins,cluster.karger)
g+geom_point(aes(x=corine$x,y=corine$y,col=prev.karger))+labs(col="")+scale_color_gradient2(low="green",mid="white",high="red",midpoint = 0.29)


