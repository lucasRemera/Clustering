kernelAcp=function(data,scaled=TRUE,transformation=c("linear","polynomial","gaussian","laplacian"),sigma=1,gamma=1,d=1,alpha=1,plot.value=FALSE){
  if(scaled) data=scale(data,scale = T,center = T)
  if(transformation[1]=="linear") varcovar=t(data)%*%(data)/(nrow(data)-1)
  else{
    varcovar=sapply(1:ncol(data), function(i){
      sapply(1:ncol(data), function(j){
        if(transformation[1]=="polynomial") return((data[,i]%*%data[,j]/(nrow(data)-1)+gamma)**d) #polynomial
        if(transformation[1]=="gaussian") return(exp(-sum((data[,i]-data[,j])**2)/(2*sigma**2)))
        if(transformation[1]=="laplacian") return(exp(-alpha*sqrt(sum((data[,i]-data[,j])**2))))
      })
    })
  }
  Pas=eigen(varcovar)$vectors
  if(plot.value) plot(eigen(varcovar)$values,type='l')
  dataprim=data%*%Pas
  return(dataprim)
}
