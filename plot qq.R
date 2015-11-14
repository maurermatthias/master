source('function.R')
source('loadMATfileToR.R')
library(car)
library(nortest)

data=shen1




plot.qq2<-function(data,dist){
  stress=unique(data[,1])
  par(mfrow=c(3,4))
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    #qqnorm(dat,main=paste("Stress-level: ",as.character(stress[i])))
    #qqline(dat)
    qqPlot(dat,distribution=dist,ylab="theoretical quantiles",main=paste("Stress-level: ",as.character(stress[i])))
  }
  
}



plot.qq<-function(data,dist){
  plot.qq2(data,dist)
  a=1
  if(dist=="norm"){
    dat1=rnorm(20)
    dat2=rnorm(20)
  }else{
    a=0
  }
  
  if(a==1){
    qqPlot(dat1,distribution=dist,ylab="theoretical quantiles",main="generated sample 1")
    qqPlot(dat2,distribution=dist,ylab="theoretical quantiles",main="generated sample 2")
  }
}


plot.qq(data,"norm")
plot.qq(data,"gamma")

########################################
#lilliefors-test NV

lilliefors<-function(data){
  stress=unique(data[,1])
  p.val=c()
  p.val2=c()
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    p.val=c(p.val,lillie.test(dat)$p.value)
    p.val2=c(p.val2,suppressWarnings(ks.test(dat,'pnorm',mean(dat),sd(dat)))$p.val)
  }
  p.val
  a=list()
  a[['stress']]=stress
  a[['norm.lillie']]=p.val
  a[['norm.ks']]=p.val2
  a
}

l=lilliefors(shen1)
