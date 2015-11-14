source('function.R')
library(R.matlab)
#load R-Data:
shen1=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
shen2=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")
logShen1=shen1
logShen1[,2]=log(shen1[,2])
bcshen1=shen1
bcshen1[,2]=shen1[,2]^(-0.5)
logShen2=shen2
logShen2[,2]=log(shen2[,2])

#returns 1. moment for each stress level
get.1moments<-function(data)
{
  stress=unique(data[,1])
  mu=c()
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    mu=c(mu,mean(dat))
  }
  r=data.frame(stress,mu)
}

#returns corrected 2. moment for each stress level
get.c2moments<-function(data)
{
  stress=unique(data[,1])
  cor.var=c()
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    cor.var=c(cor.var,sum((dat-mean(dat))^2)/(length(dat)-1))
  }
  r=data.frame(stress,cor.var)
}

#returns 2. moment for each stress level
get.2moments<-function(data)
{
  stress=unique(data[,1])
  var=c()
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    var=c(var,sum((dat-mean(dat))^2)/(length(dat)))
  }
  r=data.frame(stress,var)
}

#returns mean and corrected 2. moment for each stress level
get.moments<-function(data)
{
  d1=get.1moments(data)
  d2=get.c2moments(data)
  d3=get.2moments(data)
  r=data.frame(d1[1],d1[2],d2[2],d3[2])
}

f1=get.moments(shen1)
lf1=get.moments(logShen1)
f2=get.moments(shen2)
lf2=get.moments(logShen2)

#plot moments
plot(f1$stress,f1$mu,xlab="Belastung",ylab="Belastungszyklen",main="erstes empirisches Moment")
plot(f1$stress,f1$var,xlab="Belastung",ylab="Belastungszyklen",main="zweites zentrales empirisches Moment")

m.id=lm(f1$mu~f1$stress)
m.cl1=lm(1/f1$mu~f1$stress)
m.cl2=lm(1/f1$mu~f1$stress+I(f1$stress^2))

par(mfrow=c(1,2))
plot(f1$stress,f1$mu,xlab="stress",ylab="g(sample mean)",main="identity link")
plot(f1$stress,-1/f1$mu,xlab="stress",ylab="g(sample mean)",main="canoical link")
