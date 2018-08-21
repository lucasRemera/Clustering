# data(iris)
# ##init
# nc=30
# carte=expand.grid(x=1:nc,y=1:nc)
# # carte$Sepal.Length=rnorm(nc**2,mean(iris$Sepal.Length),sd(iris$Sepal.Length))
# # carte$Sepal.Width=rnorm(nc**2,mean(iris$Sepal.Width),sd(iris$Sepal.Width))
# # carte$Petal.Length=rnorm(nc**2,mean(iris$Petal.Length),sd(iris$Petal.Length))
# # carte$Petal.Width=rnorm(nc**2,mean(iris$Petal.Width),sd(iris$Petal.Width))
# 
# carte$Sepal.Length=rnorm(nc**2)
# carte$Sepal.Width=rnorm(nc**2)
# carte$Petal.Length=rnorm(nc**2)
# carte$Petal.Width=rnorm(nc**2)
# 
# tf=nrow(iris)
# sigma0=2
# sigmaf=.1
# eps0=2
# epsf=.1
# classif=rep(0,tf)
# #set.seed(0)
# samp.iris=iris[sample(1:tf,tf),]
# #samp.iris[,1:4]=sapply(1:4,function(i) rnorm(tf,as.integer(samp.iris$Species),.1))
# samp.iris[,1:4]=scale(samp.iris[,1:4])
# 
# 
# for(tt in 1:tf){
#   eps=eps0*(epsf/eps0)**(tt/tf)
#   sigma=sigma0*(sigmaf/sigma0)**(tt/tf)
#   elemToCarte=apply(carte[,3:6], 1, function(i) sqrt(sum((i-samp.iris[tt,1:4])**2 )))
#   idx=which.min(elemToCarte)
#   nodesToNode=as.matrix(dist(carte[,1:2],upper = T))[idx,] #h
#   elemToNodes=t(t(-carte[,3:6])+as.numeric(samp.iris[tt,1:4]))
#   carte[,3:6]=carte[,3:6]+eps*exp(-nodesToNode**2/(2*sigma**2))*elemToNodes
#   classif[tt]=idx
# }

#@arg carte : les neurones : deux colonnes avec abscisse et ordonnées : connectés à ceux qui leur sont voisin à 1 u.a.
#@arg data les données (matrice ou df)
#@arg sigmaf et sigma0 : coefficient de voisinnage sigma(t)=sigma0*(sigmaf/sigma0)**t
#@arg epsf et eps0 : coefficient d'apprentissage epsilon(t)=eps0*(epsf/eps0)**t
#@arg tore : le plan est il un tore ?
#shuffle/sampling/turn : façon d'apprendre
SOM.learn=function(carte,data,sigma0=1,sigmaf=0.1,eps0=1,epsf=0.1,scaled=TRUE,shuffle=FALSE, sampling=FALSE,nsample=100,tore=FALSE,turn=FALSE, nturn=10,pca=FALSE){
  scaled=scaled | pca
  n.var=ncol(data)
  n.obs=nrow(data)
  n.neur=nrow(carte)
  if(scaled){
    data=scale(data)
    W=matrix(rnorm(n.var*n.neur),ncol = n.var,nrow=n.neur)
  }
  else W=sapply(1:n.var, function(i) rnorm(n.neur, mean(data[,i]),sd(data[,i])))
  if(shuffle & !turn){
    data.order=sample(1:n.obs,size = n.obs)
    data=data[data.order,]
  }
  classif=c()
  if(turn){
    for(kk in 1:nturn){
      data.order=sample(1:n.obs,size = n.obs)
      data2=data[data.order,]
      for(tt in 1:n.obs){
        eps=eps0*(epsf/eps0)**(((kk-1)*nobs+tt)/(n.obs*nturn))
        sigma=sigma0*(sigmaf/sigma0)**(((kk-1)*nobs+tt)/(n.obs*nturn))
        elemToCarte=apply(W, 1, function(i) sqrt(sum((i-data2[tt,])**2 )))
        idx=which.min(elemToCarte)
        #if(!tore) nodesToNode=as.matrix(dist(carte,upper = T))[idx,] #h
        #else nodesToNode=distTtore(carte,idx)
        nodesToNode=distTtore(carte,idx,!tore)#
        elemToNodes=t(t(-W)+as.numeric(data2[tt,]))
        #print(exp(-1**2/(2*sigma**2)))
        W=W+eps*exp(-nodesToNode**2/(2*sigma**2))*elemToNodes
        if(kk==nturn) classif=c(classif,idx)
      }
    }
    data=data2
  }
  else if(sampling){
    data.order=c()
    for(tt in 1:nsample){
      eps=eps0*(epsf/eps0)**(tt/nsample)
      sigma=sigma0*(sigmaf/sigma0)**(tt/nsample)
      idx.d=sample(1:n.obs,1)
      elemToCarte=apply(W, 1, function(i) sqrt(sum((i-data[idx.d,])**2 )))
      idx=which.min(elemToCarte)
      #if(!tore) nodesToNode=as.matrix(dist(carte,upper = T))[idx,] #h
      #else nodesToNode=distTtore(carte,idx)
      nodesToNode=distTtore(carte,idx,!tore)#
      elemToNodes=t(t(-W)+as.numeric(data[idx.d,]))
      W=W+eps*exp(-nodesToNode**2/(2*sigma**2))*elemToNodes
      classif=c(classif,idx)
      data.order=c(data.order,idx.d)
    }
  }
  else{
    for(tt in 1:n.obs){
      eps=eps0*(epsf/eps0)**(tt/n.obs)
      sigma=sigma0*(sigmaf/sigma0)**(tt/n.obs)
      elemToCarte=apply(W, 1, function(i) sqrt(sum((i-data[tt,])**2 )))
      idx=which.min(elemToCarte)
      #if(!tore) nodesToNode=as.matrix(dist(carte,upper = T))[idx,] #h
      #else nodesToNode=distTtore(carte,idx)
      nodesToNode=distTtore(carte,idx,!tore)#
      elemToNodes=t(t(-W)+as.numeric(data[tt,]))
      W=W+eps*exp(-nodesToNode**2/(2*sigma**2))*elemToNodes
      classif=c(classif,idx)
    }
    if(!shuffle) data.order=1:n.obs
  }
  return(list(learn=W,idx.data.order=data.order,class=classif,dataset=data,carte=carte))
}


#méthode pour calculer les distances entre neurones
distTtore=function(carte, idx,notore=FALSE){
  nx=max(carte[,1])
  ny=max(carte[,2])
  ii=carte[idx,1]
  jj=carte[idx,2]
  di=abs(carte[,1]-ii)
  dj=abs(carte[,2]-jj)
  if(notore) return(di+dj)
  dii=nx-pmax(ii,carte[,1])+pmin(ii,carte[,1])
  djj=ny-pmax(jj,carte[,2])+pmin(jj,carte[,2])
  return(pmin(di,dii)+pmin(dj,djj))
}

#calcule la U-Matrix
UMatrix=function(carte,W,tore=FALSE){
  nx=max(carte[,1])
  ny=max(carte[,2])
  vpx=seq(1.5,nx-.5,by=1)
  vpy=seq(1.5,nx-.5,by=1)
  vx=1:nx
  vy=1:ny
  xyy=expand.grid(x=vx,y=vpy)
  xxy=expand.grid(x=vpx,y=vy)
  deltax=sapply(1:nrow(xyy), function(i){
    idx=which(carte[,1]==xyy[i,1] & abs(xyy[i,2]-carte[,2])==.5)
    return(sum((W[idx[1],]-W[idx[2],])**2))
  })
  deltay=sapply(1:nrow(xxy), function(i){
    idx=which(carte[,2]==xxy[i,2] & abs(xxy[i,1]-carte[,1])==.5)
    return(sum((W[idx[1],]-W[idx[2],])**2))
  })
  return(cbind(rbind(xyy,xxy),delta=c(deltax,deltay) ))
}

# SOM.pred=function(learning,target,test,fuzzy=FALSE,burnin=0){
#   n.var=ncol(test)
#   n.obs=nrow(test)
#   n.neur=nrow(test)
#   sapply(1:n.obs, function(i){
#     idx=which.min(apply(learning, 1, function(i) sqrt(sum((i-target[tt,])**2 ))))
#   })
# }
# 
# nc=20
# xy=expand.grid(x=1:nc,y=1:(nc))
# prediction=SOM.learn(xy,iris[,1:4],sampling = T,nsample = 1000,sigma0 = 5,sigmaf = 1,eps0 = 1,epsf = .5,turn = TRUE)
# ggplot()+geom_point(aes(x=xy$x[prediction$class],y=xy$y[prediction$class],col=iris$Species[prediction$idx.data.order] ),alpha=.5)+labs(col="")
# 
# ggplot()+geom_point(aes(x=xy$x[prediction$class[800:1000]],y=xy$y[prediction$class[800:1000]],col=iris$Species[prediction$idx.data.order[800:1000]] ),alpha=.5)+labs(col="")
# 
# 
# cartecol=apply(prediction$learn,2,function(i) (i-min(i))/(max(i)-min(i)))
# ggplot()+geom_raster(aes(x=xy$x,y=xy$y),fill=rgb(cartecol[,1],cartecol[,2],cartecol[,3],1))
# 
# SOM.pred(prediction$learn,iris[,1:4],prediction$learn)
# 
# sapply(1:4,function(i) rnorm(tf,as.integer(samp.iris$Species),.1))
# 
# plot(samp.iris[,1],samp.iris[,2])

## EXAMPLES ##

nc=50 #le nombre de colonnes de neurones 
xy=expand.grid(x=1:nc,y=1:(nc)) # carte 50*50
#xy=expand.grid(x=1:50,y=1:(82))
freq.dyspatc=t(apply(count.dyspatc, 1, function(i) i/sum(i))) #calcul des fréquences des malformations pour chaque code ATC
#som.dyspatc=SOM.learn(xy,freq.dyspatc,sigma0 = 5,sigmaf = 1,eps0 = 1,epsf = .5,turn = TRUE,tore = FALSE)
som.dyspatc=SOM.learn(xy,freq.dyspatc,sigma0 = 5,sigmaf = 4,eps0 = 2,epsf = 1,turn = TRUE,tore = FALSE)
ggplot()+geom_text_repel(aes(x=xy$x[som.dyspatc$class],y=xy$y[som.dyspatc$class],label=rownames(freq.dyspatc)[som.dyspatc$idx.data.order],col=substr(rownames(freq.dyspatc)[som.dyspatc$idx.data.order],1,1) ))+labs(col="")

malfToAnalyze="Q60"
idxmalf=grep(malfToAnalyze, colnames(count.dyspatc))

meanmalf=som.dyspatc$learn[,idxmalf]
#meanmalf=apply(som.dyspatc$learn[,idxmalf],1,mean)
#sdmalf=apply(som.dyspatc$learn[,idxmalf],1,sd)
ggplot()+geom_raster(aes(x=xy$x,y=xy$y,fill=meanmalf))+geom_text_repel(aes(x=xy$x[som.dyspatc$class],y=xy$y[som.dyspatc$class],label=rownames(freq.dyspatc)[som.dyspatc$idx.data.order],col=substr(rownames(freq.dyspatc)[som.dyspatc$idx.data.order],1,1) ))+labs(col="")
#Quelles zones sont activées pour la malformation Q60 ?

umatrix=UMatrix(xy,som.dyspatc$learn)
ggplot()+geom_point(data=umatrix,aes(x=x,y=y,col=pmin(delta,50)))+scale_colour_gradient(low="white",high="black")
ggplot()+geom_raster(data=umatrix,aes(x=x,y=y,fill=pmin(delta,40)))+scale_fill_gradient(low="white",high="black")+
  geom_text_repel(aes(x=xy$x[som.dyspatc$class],y=xy$y[som.dyspatc$class],label=rownames(freq.dyspatc)[som.dyspatc$idx.data.order],col=substr(rownames(freq.dyspatc)[som.dyspatc$idx.data.order],1,1) ))+labs(col="")


#acp=dudi.pca((freq.dyspatc),scannf = F,nf = 2)
#acp$l1
