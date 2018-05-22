kernelAcp=function(data,scaled=TRUE,transformation=c("linear","polynomial","gaussian","laplacian"),sigma=1,gamma=1,d=1,alpha=1,centralize.K=FALSE,plot.value=FALSE){
  if(scaled) data=scale(data,scale = T,center = T)
  if(transformation[1]=="linear") varcovar=t(data)%*%(data)/(nrow(data)-1)
  else{
    varcovar=sapply(1:ncol(data), function(i){
      sapply(1:ncol(data), function(j){
        #if(transformation[1]=="polynomial") return((data[,i]%*%data[,j]/(nrow(data)-1)+gamma)**d) #polynomial
        if(transformation[1]=="polynomial") return((data[,i]%*%data[,j]+gamma)**d) #polynomial
        if(transformation[1]=="gaussian") return(exp(-sum((data[,i]-data[,j])**2)/(2*sigma**2)))
        if(transformation[1]=="laplacian") return(exp(-alpha*sqrt(sum((data[,i]-data[,j])**2))))
      })
    })
  }
  ##K=varcovar
  ##K'=K-1N%*%K-K%*%1N+1N%*%K%*%1N
  if(centralize.K){
    unN=matrix(1/nrow(varcovar),nrow = nrow(varcovar),ncol=nrow(varcovar))
    varcovar=varcovar-unN%*%varcovar-varcovar%*%unN+unN%*%varcovar%*%unN
  }
  Pas=eigen(varcovar)$vectors
  if(plot.value) plot(eigen(varcovar)$values,type='l')
  #dataprim=sapply(list, function)
  if(transformation[1]=="linear") dataprim=data%*%Pas
  else{
    dataprim=t(sapply(1:nrow(data), function(i){
      sapply(1:ncol(Pas), function(j){
        #if(transformation[1]=="polynomial") return((data[i,]%*%Pas[,j]/(nrow(data)-1)+gamma)**d) #polynomial
        if(transformation[1]=="polynomial") return((data[i,]%*%Pas[,j]+gamma)**d)
        if(transformation[1]=="gaussian") return(exp(-sum((data[i,]-Pas[,j])**2)/(2*sigma**2)))
        if(transformation[1]=="laplacian") return(exp(-alpha*sqrt(sum((data[i,]-Pas[,j])**2))))
      })
    }))
  } 
  return(dataprim)
}


###########
# Example #
###########

x=runif(500)
y=runif(500)
inCluster=((x>.3)&(x<.6)&(y>.3)&(y<.6))
obs=rpois(500,1000)
cas=rbinom(500,obs,0.2+as.numeric(inCluster)*0.1)

dumb=as.matrix(data.frame(x=x,y=y,p=cas/(cas+obs)))
dumb.acp=kernelAcp(dumb,transformation = "gaussian",sigma = 5)
ggplot()+geom_point(aes(dumb.acp[,1],dumb.acp[,2],col=factor(as.numeric(inCluster))))+labs(col="")

cl=cutree(hclust(dist(dumb.acp[,1:2]),method = "median"),2)
ggplot()+geom_point(aes(dumb[,1],dumb[,2],col=factor(cl)))+labs(col="")
