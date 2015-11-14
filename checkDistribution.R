#set working directory before execution!!!!!
source('function.R')
source('basicMoments.R')
#load packages
library('stats')
library('stats4')
library('MASS')
library('evir')

#normal distribution boxcox transformation
boxcox(shen1[,2] ~ shen1[,1], lambda=seq(-1,0, len=18))
logShen1[,2]=log(shen1[,2])
bcshen1=shen1
bcshen1[,2]=shen1[,2]^(-0.5)
pval.normal=suppressWarnings(pfree(bcshen1,'normal','pnorm'))


#read dataset 1 -weibull
shen1=readRDS("shen1.rds")
pval1=round(p.val(shen1),2)
param1=param(shen1)
#Ann: shape parameter of wbl-dist constant
shape1=mean(param1$wbl.shape)
#p-values for given shape parameter (estimation done for both, shape replaced)
pnew1=suppressWarnings(pfree(shen1,'weibull','pweibull',shape1))


#read dataset 1
shen2=readRDS("shen2.rds")
pval2=round(p.val(shen2),2)
param2=param(shen2)
#Ann: shape parameter of wbl-dist constant
shape2=mean(param2$wbl.shape)
#p-values for given shape parameter (estimation done for both, shape replaced)
pnew2=suppressWarnings(pfree(shen2,'weibull','pweibull',shape2))


#compare variances in the normal distribution case
f.test<-function(data){
  set=unique(data[,1])
  df=data.frame(set)
  var=get.2moments(data)
  for(j in 1: length(set)){
    f=c()
    for(i in 1:length(set)){
      f=c(f,pf(var[j,2]/var[i,2],length(part(data[,2],data[,1]==set[j]))-1,length(part(data[,2],data[,1]==set[i]))-1))
    }
    df[,as.character(set[j])] <- f
  }
  round(df,3)
}
