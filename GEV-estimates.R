source('function.R')
library(evir)
library(stats4)

#negative log-likelihood function gev-distribution  -not used
nlh<-function(par,data){
  gamma=par[1]
  mu=par[2]
  sigma=par[3]
  values=data
  if(sigma <= 0){
    return(-Inf)
  }
  m<-length(values)
  val=-m*log(sigma)
  s1=0
  s2=0
  for(counter in 1:m){
    if(1+gamma*((values[counter]-mu)/sigma)<=0){
      return(-Inf)
    }
    if(gamma!=0){
      s1=s1+log(1+gamma*((values[counter]-mu)/sigma))
      s2=s2+(1+gamma*((values[counter]-mu)/sigma))^(-1/gamma)
    }
    if(gamma==0){
      s1=s1+((values[counter]-mu)/sigma)
      s2=s2+(exp(-(values[counter]-mu)/sigma))
    }
  }
  if(gamma!=0){
    s1=s1*(1+1/gamma)
  }
  val=val-s1-s2
  return(-val)
}

#ML-estimates for gev-distribution   - not used
likelihood.par.est.gev<-function(data){
  #par=c(min(data[,1])*0.9,max(data[,1])*1.1,10000,10000)
  #par=moment.par.est(shen1)$estimates
  par=gev(data,1)$par.est
  o.result=suppressWarnings(optim(par,fn=nlh,data=data))
  
  a=list()
  a[["optimazation.result"]]=o.result
  a[["data"]]=data
  a[["estimates"]]=o.result$par
  a
}


#gev-parameter estimation based on the evir package
gev.est.evir<-function(data){
  stress=unique(data[,1])
  xi=c()
  sigma=c()
  mu=c()
  estimate=list()
  for(i in 1: length(stress)){
    dat=part(data[,2],data[,1]==stress[i])
    est=gev(dat,1)
    xi=c(xi,est$par.est[[1]])
    sigma=c(sigma,est$par.est[[2]])
    mu=c(mu,est$par.est[[3]])
    estimate[[i]]=est
  }
  gev.evir.est=data.frame(stress,xi,sigma,mu)
  a=list()
  a[["param"]]=gev.evir.est
  a[["optim"]]=estimate
  a
}

plot.evir.est<-function(a){
  #par(mfrow=c(3,1))
  data=a$evir.est$param
  x=data$stress
  xi=data$xi
  sigma=data$sigma
  mu=data$mu
  plot.new()
  d=0.225
  
  est=a$est.par.struct
  k=est[1]
  a=est[2]
  b=est[3]
  c1=est[4]
  c2=est[5]
  xx=(0:100000)*(b-a)/100000+a
  xval=c(xx,rev(xx))
  yval.xi=rep(k,length(xval))
  yval.sigma=((b-a)/(xval-a)-1)*c1
  yval.mu=((b-a)/(xval-a)-1)*c2
  
  #xi plot
  par(fig=c(0,1,1-2*d,1), new=TRUE)
  plot(x,xi,xlab="",ylab=expression(xi),xaxt='n')
  polygon(xval,yval.xi,border="red")
  legend(x="topright", legend=c("estimate       ","observation       "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  #sigma plot
  par(fig=c(0,1,0.5-d,0.5+d), new=TRUE)
  plot(x,sigma,xlab="",ylab=expression(sigma),xaxt='n')
  polygon(xval,yval.sigma,border="red")
  
  #mu plot
  par(fig=c(0,1,0,2*d), new=TRUE)
  plot(x,mu,xlab="stress",ylab=expression(mu))
  polygon(xval,yval.mu,border="red")
  
}


#error function for fitting the parameter functions
error.gev.ls<-function(par,data){
  a=par[1]
  b=par[2]
  c1=par[3]
  c2=par[4]
  
  stress=data$stress
  
  returnValue=0
  if(a >= min(stress) || b <= max(stress)){
    returnValue= Inf
  }else{
    sigma=data$sigma
    mu=data$mu
    
    for(i in 1:length(stress)){
      par.function.sigma=((b-a)/(stress[i]-a)-1)*c1
      returnValue=returnValue+((sigma[i]-par.function.sigma)/sigma[i])^2
      
      par.function.mu=((b-a)/(stress[i]-a)-1)*c2
      #par.function.mu=((b-a)/(stresslevel[i]-a)-1)*c2*((b-a)/(stresslevel[i]-a))
      returnValue=returnValue+(((mu[i])^(1/1)-(par.function.mu)^(1/1))/(mu[i])^(1/1))^2
      
    }
  }
  returnValue
  
}

plot.evir.est<-function(a){
  #par(mfrow=c(3,1))
  data=a$evir.est$param
  x=data$stress
  xi=data$xi
  sigma=data$sigma
  mu=data$mu
  plot.new()
  d=0.225
  
  est=a$est.par.struct
  k=est[1]
  a=est[2]
  b=est[3]
  c1=est[4]
  c2=est[5]
  xx=(0:100000)*(b-a)/100000+a
  xval=c(xx,rev(xx))
  yval.xi=rep(k,length(xval))
  yval.sigma=((b-a)/(xval-a)-1)*c1
  yval.mu=((b-a)/(xval-a)-1)*c2
  
  #xi plot
  par(fig=c(0,1,1-2*d,1), new=TRUE)
  plot(x,xi,xlab="",ylab=expression(xi),xaxt='n')
  polygon(xval,yval.xi,border="red")
  legend(x="topright", legend=c("estimate       ","observation       "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  #sigma plot
  par(fig=c(0,1,0.5-d,0.5+d), new=TRUE)
  plot(x,sigma,xlab="",ylab=expression(sigma),xaxt='n')
  polygon(xval,yval.sigma,border="red")
  
  #mu plot
  par(fig=c(0,1,0,2*d), new=TRUE)
  plot(x,mu,xlab="stress",ylab=expression(mu))
  polygon(xval,yval.mu,border="red")
  
}

plot.est.moment<-function(A,type){
  #par(mfrow=c(3,1))
  data=A$evir.est$param
  xi=data$xi
  sigma=data$sigma
  mu=data$mu
  
  if(type == 'mle')
    est=A$mle.est.par.struct
  else
    est=A$est.par.struct
    
  k=est[1]
  a=est[2]
  b=est[3]
  c1=est[4]
  c2=est[5]
  xx=(0:100000)*(b-a)/100000+a
  xval=c(xx,rev(xx))
  yval.xi=rep(k,length(xval))
  yval.sigma=((b-a)/(xval-a)-1)*c1
  yval.mu=((b-a)/(xval-a)-1)*c2
  if(k>=0.5)
    stop('xi >= 0.5 ... expection/variance not finite')
  if(k>=1)
    stop('xi >= 1 ... expection not finite')
  if(k!=0)
    yval.EV=yval.mu+yval.sigma*(gamma(1-yval.xi)-1)/(yval.xi)
  else
    yval.EV=yval.mu+yval.sigma*0.5772156649015328606065120900824024310421 #eulers constant
  if(k!=0)
    yval.VAR=yval.sigma^2*(gamma(1-2*yval.xi)-(gamma(1-2*yval.xi))^2)/yval.xi^2
  else
    yval.VAR=yval.sigma^2*pi^2/6
  
  
  plot.new()
  d=0.59
  
  #plot expectation
  par(fig=c(0,1,1-d,1), new=TRUE)
  x=A$moments$stress
  y=A$moments$mu
  plot(x,y,xlab="",ylab="EV")
  polygon(xval,yval.EV,border="red")
  legend(x="topright", legend=c("estimate  EV         ","observation EV           "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  #plot variance
  par(fig=c(0,1,0,d), new=TRUE)
  x=A$moments$stress
  y=A$moments$var
  plot(x,y,xlab="stress",ylab="VAR")
  polygon(xval,yval.VAR,border="red")
  legend(x="topright", legend=c("estimate VAR          ","observation VAR          "), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
}

plot.est.moment.mle<-function(A){
  plot.est.moment(A,'mle')
}

plot.est.moment.ls<-function(A){
  plot.est.moment(A,'ls')
}

#plot expectation and quantile area containing q*100% of the data (theoretically)
plot.est.quantile<-function(A,q){
  stress=A$data[,1]
  
  est=A$est.par.struct
  k=est[1]
  a=est[2]
  b=est[3]
  c1=est[4]
  c2=est[5]
  xx=(0:10000)*(b-a)/10000+a
  xval=c(xx,rev(xx))
  yval.xi=rep(k,length(xval))
  yval.sigma=((b-a)/(xval-a)-1)*c1
  yval.mu=((b-a)/(xval-a)-1)*c2
  points.xi=rep(k,length(stress))
  points.sigma=((b-a)/(stress-a)-1)*c1
  points.mu=((b-a)/(stress-a)-1)*c2
  if(k>=1)
    stop('xi >= 1 ... expection not finite')
  if(k!=0){
    yval.EV=yval.mu+yval.sigma*(gamma(1-yval.xi)-1)/(yval.xi)
    points.EV=points.mu+points.sigma*(gamma(1-points.xi)-1)/(points.xi)
  }
  else{
    yval.EV=yval.mu+yval.sigma*0.5772156649015328606065120900824024310421 #eulers constant
    points.EV=points.mu+points.sigma*0.5772156649015328606065120900824024310421 #eulers constant
  }
  plot.new()
  par(fig=c(0,1,0,1), new=TRUE)
  x=A$data[,1]
  y=log(A$data[,2])
  boolv=(yval.EV!=0) & (yval.EV!=Inf)
  lyval.EV=part(log(yval.EV),boolv )
  xval2=part(xval,boolv)
  
  #plot points
  #plot(x,y,xlab="stress",ylab="log(N)")
  #quantile function exchange: gamma -> gev
  alpha=1-q
  lower=qgev(alpha/2,points.xi,points.mu,points.sigma)
  upper=qgev(1-alpha/2,points.xi,points.mu,points.sigma)
  #print(upper)
  plot(NULL, xlim=c(min(x)/1.01,max(x)*1.01), ylim=c(min(y/1.01),max(y)*1.01), ylab="log(N)", xlab="stress")
  count=0
  for(counter in 1:length(x)){
    if(log(lower[counter]) <= y[counter]  &&  log(upper[counter]) >= y[counter] ){
      points(x[counter],y[counter],pch="*",col="red", cex = 0.75)
      count=count+1;
    }
    else
      points(x[counter],y[counter],pch="o", col="blue", cex=0.5)
  }
  
  #plot expectation line
  polygon(xval2,lyval.EV,border="red")

  #plot quantile line 1
  q1=q+(1-q)/2
  yval.q=(yval.mu-(1-(-log(q1))^(-yval.xi))*(yval.sigma)/(yval.xi))
  boolv=(yval.q>0) & (yval.q<Inf) & !is.nan(yval.q)
  lyval.q=part(log(yval.q),boolv )
  xval3=part(xval,boolv )
  polygon(xval3,lyval.q,border="green",lty='dashed')
  
  #plot quantile line 2
  q2=(1-q)/2
  yval.q2=(yval.mu-(1-(-log(q2))^(-yval.xi))*(yval.sigma)/(yval.xi))
  boolv2=(yval.q2>0) & (yval.q2<Inf) & !is.nan(yval.q2)
  lyval.q2=part(log(yval.q2),boolv2 )
  xval4=part(xval,boolv2 )
  polygon(xval4,lyval.q2,border="green",lty='dashed')
  

  
  legend(x="topright", legend=c("estimate  EV         ","observations           ","observations           ","quantile line"), col = c("red","red","blue","green"), lty = c("solid",NA,NA,"dashed"), pch=c(NA,"*",'o',NA))
  
  inside=count/length(x)
  inside

}

#nll-function for gev (xi != 0)
nll.gev<-function(data, xiv, parameter){
  #get function parameters
  k=parameter[1]
  a=parameter[2]
  b=parameter[3]
  c1=parameter[4]
  c2=parameter[5]
  
  #get stress-levels
  x.vec=unique(data[,1])
  
  #if(a>=min(x.vec) || b<=max(x.vec) || c1<=0 || c2<=0)
  #  return(Inf)
  
  #get distribution parameter
  xi=k
  sigma.vec=((b-a)/(x.vec-a)-1)*c1
  mu.vec=((b-a)/(x.vec-a)-1)*c2
  
  if(xi==0 || xi<=min(xiv) || xi >=max(xiv))
    return(Inf)
  
  ml.val=0;
  
  
  for(i in 1:length(x.vec)){
    x=part(data[,2],data[,1]==x.vec[i])
    sigma=sigma.vec[i]
    mu=mu.vec[i]
    s1=0
    s2=0
    for(j in 1: length(x)){
      if(sigma<=0 || (1+xi*(x[j]-mu)/sigma)<=0){
        return(Inf)
      }
      s1=s1+log(1+xi*(x[j]-mu)/sigma)
      s2=s2+(1+xi*(x[j]-mu)/sigma)^(-1/xi)
    }
    ml.val=ml.val+length(x)*sigma+(1+1/xi)*s1+s2
  }
  ml.val
}


#try to fit the parameters with functions
#xi-constant(k) .... no need for a fit - take mean
#mu/sigma -1/x*c_1/c_2
gev.est.par.structure<-function(a, control){
  observations=a$evir.est$param
  stress=observations$stress
  par=c(min(stress)*0.9,max(stress)*1.1,10000,10000)
  o.result=suppressWarnings(optim(par,fn=error.gev.ls,data=observations, control=control))
  xi=mean(observations$xi)
  
  a[["optimazation.result.par.struct"]]=o.result
  a[["est.par.struct"]]=c(xi,o.result$par)
  a
  
}


mle.gev<-function(a,control){
  par=c(-0.3,48,367,13800,5566)
  xiv=a$evir.est$param$xi
  o.result=suppressWarnings(optim(par,fn=nll.gev, method = "SANN",data=a$data, xiv=xiv, control=control))
  #start=list(parameter=par)
  #fixed=list(data=a$data, xiv=xiv)
  #o.result=suppressWarnings(mle(minuslogl=nll.gev, start, method = "BFGS", fixed))
  a[["mle.optimazation.result.par.struct"]]=o.result
  a[["mle.est.par.struct"]]=o.result$par
  a
}


#####over all function
gev.approach<-function(data){
  a=list()
  a[["data"]]=data
  a[["moments"]]=get.moments(data)
  #calculate ml-estimates for gev distribution for each stress-level
  a[["evir.est"]]=suppressWarnings(gev.est.evir(data))
  #calculate ls-estimates for parameter structure
  control = list(maxit = 20000)
  #control=list()
  a=gev.est.par.structure(a, control)
  a=mle.gev(a, control)
  
  #plot structure and save plot data in a
  #plot.evir.est(a)
  #plot.est.moment(a)
  a
}

