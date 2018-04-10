theta=seq(0,2*pi,length.out = 100)
r=theta

polaireToCartesien=function(r,theta){
  return(cbind(r*cos(theta),r*sin(theta)))
}

parameterCurve=function(xy,ordre=1:length(x)){
  xy=xy[ordre,]
  cumsum(c(0,sapply(1:(nrow(xy)-1), function(i) sqrt(sum((xy[i,]-xy[i+1])**2)))))
}


orderCurve=function(xy,method=c("nearest_insertion", "farthest_insertion", "cheapest_insertion", "arbitrary_insertion" )){
  library(TSP)
  tsp = TSP(dist(xy))
  
  tsp2=insert_dummy(tsp,1,label="cut")
  met=method[1]
  tour=solve_TSP(tsp2,method=met,two_opt=T)
  path=cut_tour(tour,"cut")
  return(path)
}


loess2D=function(xy,span=0.5,degree=1,method=c("nearest_insertion", "farthest_insertion", "cheapest_insertion", "arbitrary_insertion" )){
  path=orderCurve(xy)
  tt=parameterCurve(xy,path)
  xx=xy[path,1]
  yy=xy[path,2]
  xx.l=loess(xx~tt,span=span,degree = degree)
  xx.s=predict(xx.l)
  yy.l=loess(yy~tt,span = span,degree=degree)
  yy.s=predict(yy.l)
  cbind(xx.s,yy.s)
}

xysb=polaireToCartesien(r,theta)
et=0.3
xy=polaireToCartesien(r,theta)+cbind(rnorm(100,sd=et),rnorm(100,sd=et))
#path=orderCurve(xy)
#tt=parameterCurve(xy,path)

#plot(xy[,1],xy[,2])



gg1=ggplot(data.frame(xy))+geom_point(aes(x=X1,y=X2),col="black")
#gg1=gg1+geom_path(aes(x=X1[path],y=X2[path]),col="red")


# xx=xy[path,1]
# yy=xy[path,2]
# xx.l=loess(xx~tt)
# xx.s=predict(xx.l)
# yy.l=loess(yy~tt)
# yy.s=predict(yy.l)

xy.s=loess2D(xy,degree=2,span=0.4)

gg1+geom_path(aes(x=xy.s[,1],y=xy.s[,2]),col="green")+geom_path(aes(x=xysb[,1],y=xysb[,2]),lty=2,col="cyan")
