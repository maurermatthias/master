source('function.R')
# install.packages('FAdist')
#library('FAdist')
#library('R.matlab')
#library('stats')

#####################################
#####################################



#calculate parameter for the fit of the structure: a...lower bound, b ... upper bound
#c1 ... multiplicator expectation, c2 ... multiplicator variance for EW=((b-a)/(x-a)-1)*c1
# and Var=((b-a)/(x-a)-1)*c2 with least square methode/maximum likelihood

#function for calculation sum of squared errors
error.function.ls<-function(stresslevel,moments, g){
  a=g[1]
  b=g[2]
  c1=g[3]
  c2=g[4]
  returnValue=0
  if(a >= min(stresslevel) || b <= max(stresslevel)){
    returnValue= Inf
  }else{
    mom1=moments[,1]
    mom2=moments[,2]
    
    for(i in 1:length(stresslevel)){
      par.function.ew=((b-a)/(stresslevel[i]-a)-1)*c1
      returnValue=returnValue+((mom1[i]-par.function.ew)/mom1[i])^2

      #par.function.var=((b-a)/(stresslevel[i]-a)-1)*c2
      par.function.var=((b-a)/(stresslevel[i]-a)-1)*c2*((b-a)/(stresslevel[i]-a))
      returnValue=returnValue+(((mom2[i])^(1/1)-(par.function.var)^(1/1))/(mom2[i])^(1/1))^2

    }
  }
  returnValue
}

#negative log likelihood function (https://stat.ethz.ch/R-manual/R-devel/library/stats4/html/mle.html)
error.function.ml.gamma<-function(data, g){
  a=g[1]
  b=g[2]
  c1=g[3]
  c2=g[4]
  
  x=data[,1]
  y=data[,2]
  
  #calculate gamma distribution parameter
  k=(1-(x-a)/(b-a))*(c1^2/c2)
  theta=(c2/c1)*((b-a)/(x-a))
  
  returnValue=0
  if( a<=0 || a >= min(x) || b <= max(x) || sum(k<=0)>0  || sum(theta<=0)>0){
    returnValue= -Inf
  }else{
    #calculate value of log-likelihoodfunction
    for(i in 1:length(x)){
      if(k[i]<=100)
        returnValue=returnValue+(k[i]-1)*log(x[i])-x[i]/theta[i]-k[i]*log(theta[i])-log(gamma(k[i]))
      else 
        returnValue=returnValue+(k[i]-1)*log(x[i])-x[i]/theta[i]-k[i]*log(theta[i])-(1/2)*log(2*pi*(k[i]-1))-(k[i]-1)*log((k[i]-1)/exp(1))
    }
  }
  returnValue=-returnValue
  returnValue
}


likelihood.par.est.gamma<-function(data){
  #par=c(min(data[,1])*0.9,max(data[,1])*1.1,10000,10000)
  par=moment.par.est(shen1)$estimates
  o.result=suppressWarnings(optim(par,fn=error.function.ml.gamma,data=data))
  
  a=list()
  a[["optimazation.result"]]=o.result
  a[["data"]]=data
  a[["estimates"]]=o.result$par
  a
}


moment.par.est<-function(data){
  stress=unique(data[,1])
  #calculate empirical moments for optimization:
  moments.for.optimization = get.moments(data)
  #2. empirical moment:4/corrected 2. emp. moment: 3
  used.moments=data.frame(mu=moments.for.optimization[,2], var=moments.for.optimization[,4])
  #initial values:
  par=c(min(stress)*0.9,max(stress)*1.1,10000,10000)
  o.result=suppressWarnings(optim(par,fn=error.function.ls,stresslevel=stress, moments = used.moments))
  a=list()
  a[["optimazation.result"]]=o.result
  a[["data"]]=data
  a[["emp.moments.used.for.fit"]]=used.moments
  a[["estimates"]]=o.result$par
  a
}


##helper functions:

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

get.3moments<-function(data){
  stress=unique(data[,1])
  m3=c()
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    m3=c(m3,sum((dat-mean(dat))^3)/(length(dat)))
  }
  r=data.frame(stress,m3)
}

get.skewness<-function(data){
  stress=unique(data[,1])
  var=get.2moments(data)[2]
  skew=c()
  for(i in 1:length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    skew=c(skew,sum(((dat-mean(dat))/(var[[1]][i])^(1/2))^3)/(length(dat)))
  }
  r=data.frame(stress,skew)
}

#returns mean and corrected 2. moment for each stress level
get.moments<-function(data)
{
  d1=get.1moments(data)
  d2=get.c2moments(data)
  d3=get.2moments(data)
  d4=get.3moments(data)
  d5=get.skewness(data)
  r=data.frame(d1[1],d1[2],d2[2],d3[2],d4[2],d5[2])
}

#plot moment estimation
plot.est.gamma=function(pa){
  g1=pa$opt$par[1];g2=pa$opt$par[2];g3=pa$opt$par[3];g4=pa$opt$par[4]
  stress=unique(shen1[,1])
  xlim1=c(g1,max(stress)*1.1)
  xval=(0:100000)*(g2-g1)/100000+g1
  yval1=((g2-g1)/(xval-g1)-1)*g3
  yval2=((g2-g1)/(xval-g1)-1)*g4*((g2-g1)/(xval-g1))
  
  par(mfrow=c(2,2))
  
  #1. moment
  x1=c(xval,rev(xval))
  y1=c(yval1,rev(yval1))
  y1.points=pa$emp.moments.used.for.fit$mu
  plot(stress,y1.points,xlim=xlim1,pch="o", xlab="stress", ylab="value")
  polygon(x1,y1,border="red")
  legend(x="topright", legend=c("mean est.       ","emp. mean       "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  #2. moment
  x2=c(xval,rev(xval))
  y2=c(yval2,rev(yval2))
  y2.points=pa$emp.moments.used.for.fit$var
  plot(stress,y2.points ,xlim=xlim1,pch="o",xlab="stress", ylab="value")
  polygon(x2,y2,border="red")
  legend(x="topright", legend=c("variance est.       ","emp. variance       "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))

  #quotient 2.moment/1.moment
  x3=c(xval,rev(xval))
  y3=y2/y1
  y3.points=y2.points/y1.points
  plot(stress,y3.points ,xlim=xlim1,pch="o",xlab="stress", ylab="value")
  polygon(x3,y3,border="red")
  legend(x="topright", legend=c("VAR/EV est.       ","emp. VAR/EV       "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  
  #1.moment^2/2.moment
  x4=c(xval,rev(xval))
  y4=y1^2/y2
  y4.points=y1.points^2/y2.points
  plot(stress,y4.points ,xlim=xlim1,pch="o",xlab="stress", ylab="value")
  polygon(x4,y4,border="red")
  legend(x="topright", legend=c("EV^2/Var est.       ","emp. EV^2/Var       "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
}









