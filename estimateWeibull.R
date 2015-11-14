source('function.R')
# install.packages('FAdist')
library('FAdist')
library('R.matlab')
library('stats')
#load R-Data:
shen1=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
shen2=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")
logShen1=shen1
logShen1[,2]=log(shen1[,2])
bcshen1=shen1
bcshen1[,2]=shen1[,2]^(-0.5)
logShen2=shen2
logShen2[,2]=log(shen2[,2])

#estimating weibull3
#x <- rweibull3(200, shape = 3, scale = 1, thres = 100)
#suppressWarnings(fitdistr(x, function(x, shape, scale, thres) dweibull(x-thres, shape, scale), list(shape = 0.1, scale = 1, thres = 0)))


#estimating weibull2:
p=suppressWarnings(parameter(shen1,"weibull"))
p[,"est.mu.from.parameter"] <- p[,3]*gamma(1+1/p[,2])
p[,"est.sigma2.from.parameter"] <- p[,3]^2*(gamma(1+2/p[,2])-(gamma(1+1/p[,2]))^2)

#least square fit - for mu-function mu(s)=((g2-g1)/(s-g1))g3
difference<-function(g){#(stress,est.mu,g1,g2){
  g1=g[1]
  g2=g[2]
  g3=g[3]
  if((max(stress)>g2) || (min(stress)<g1)){
    returnValue=Inf
  }else{
    returnValue=0
    for(i in 1:length(stress)){
      #mu.function=log((g2-g1)/(stress[i]-g1))
      mu.function=((g2-g1)/(stress[i]-g1)-1)*g3
      returnValue=returnValue+((est.mu[i]-mu.function)/est.mu[i])^2
    }
  }
  returnValue
}

stress=p[,1]
est.mu=p[,4]
par=c(min(stress)*0.9,max(stress)*1.1,10000)
o.mu=optim(par,fn=difference)


#plot it
g1=o.mu$par[1];g2=o.mu$par[2];g3=o.mu$par[3]
par(mfrow=c(1,1))
xlim1=c(g1,max(stress)*1.1)
plot(stress,est.mu,xlim=xlim1)
xval=(0:100000)*(g2-g1)/100000+g1
#yval=log((g2-g1)/(xval-g1))
yval=((g2-g1)/(xval-g1)-1)*g3
polygon(c(xval,rev(xval)),c(yval,rev(yval)),border="red")


#same for variance and expectation:
difference2<-function(g){#(stress,est.mu,g1,g2){
  g1=g[1]
  g2=g[2]
  g3=g[3]
  g4=g[4]
  returnValue=0
  for(i in 1:length(stress)){
    #mu.function=log((g2-g1)/(stress[i]-g1))
    if(i*2<=length(stress)){
      par.function=((g2-g1)/(stress[i]-g1)-1)*g3
      returnValue=returnValue+((est.par[i]-par.function)/est.par[i])^2
    } else{
      par.function=((g2-g1)/(stress[i]-g1)-1)*g4
      returnValue=returnValue+(((est.par[i])^(1/1)-(par.function)^(1/1))/(est.par[i])^(1/1))^2
    }
  }
  returnValue
}

stress=c(p[,1],p[,1])
est.par=c(p[,4],p[,5])
par=c(min(stress)*0.9,max(stress)*1.1,10000,10000)
o.sigma2=optim(par,fn=difference2)

#plot it
g1=o.sigma2$par[1];g2=o.sigma2$par[2];g3=o.sigma2$par[3];g4=o.sigma2$par[4]
par(mfrow=c(1,2))
xlim1=c(g1,max(stress)*1.1)
xval=(0:10000)*(g2-g1)/10000+g1
yval1=((g2-g1)/(xval-g1)-1)*g3
yval2=((g2-g1)/(xval-g1)-1)*g4
plot(stress[1:length(stress)/2],est.par[1:length(stress)/2],xlim=xlim1,pch="x",main="mu")
polygon(c(xval,rev(xval)),c(yval1,rev(yval1)),border="red")
plot(stress[length(stress)/2+1:length(stress)],(est.par[length(stress)/2+1:length(stress)]),main="sigma^2",xlim=xlim1,pch="o")
polygon(c(xval,rev(xval)),(c(yval2,rev(yval2))),border="red")
