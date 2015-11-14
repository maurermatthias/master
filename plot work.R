source('checkDistribution.R')


#plot data, y-log transformiert
data=shen1
data[,2]=log(data[,2])
plot(data[,1],log(data[,2]), xlab="stress", ylab = "log(N)")

#plot data, y-log transformiert / 1/X transformed
data=shen1
data[,2]=log(data[,2])
plot(1/data[,1],log(data[,2]), xlab="1/stress", ylab = "log(N)")

#linear model
data=shen1
y=log(data[,2])
x=1/data[,1]
mod =lm(y ~ x)
intercept=mod$coef[[1]]
slope=mod$coef[[2]]
plot(x,y,xlab="1/stress", ylab = "log(N)")
xval=c(0:100)*(max(x)-min(x))/100+min(x)
yval=intercept+slope*xval
#polygon(xval,yval,border="green",lty="dotdash")

#squared model
data=shen1
y=log(data[,2])
x=1/data[,1]
x2=x^2
mod =lm(y ~ x+x2)
intercept=mod$coef[[1]]
slope.x=mod$coef[[2]]
slope.x2=mod$coef[[3]]
#plot(x,y)
xval=c(0:100)*(max(x)-min(x))/100+min(x)
yval=intercept+slope.x*xval+slope.x2*xval^2
polygon(c(xval,rev(xval)),c(yval,rev(yval)),border="blue",lty="dashed")

#squared model mit interaktion
data=shen1
y=log(data[,2])
x=1/data[,1]
x2=x^2
mod =lm(y ~ x*x2)
intercept=mod$coef[[1]]
slope.x=mod$coef[[2]]
slope.x2=mod$coef[[3]]
slope.x.x2=mod$coef[[4]]
#plot(x,y)
xval=c(0:100)*(max(x)-min(x))/100+min(x)
yval=intercept+slope.x*xval+slope.x2*xval^2+slope.x.x2*xval*xval^2
polygon(c(xval,rev(xval)),c(yval,rev(yval)),border="red",lty="dotted")

legend(min(x),max(y), c("x+x^2","x+x^2+x^3"),bty="n", lty=c("dashed","dotted"), lwd=c(1.5,1.5),col=c("blue","red"))


m2 =lm(y ~ x+x2)
par(mfrow=c(2,2))
plot(m2)
m3 =lm(y ~ x+x2+x3)
par(mfrow=c(2,2))
plot(m3)


############################################
#plot of confidence, prediction interval:
#model 2: x+x^2
data=shen1
a=1
if(a==1){
  y=log(data[,2])
  x=1/data[,1]
}else{
  #grouped data
  y=log(data[,2])
  x=1/data[,1]
  uni=unique(x)
  x.new=c()
  y.new=c()
  for(counter in 1 :(length(uni))){
    dat=part(y,x==uni[counter])
    x.new=c(x.new,uni[counter])
    y.new=c(y.new,mean(dat))
  }
  x.old=x
  y.old=y
  x=x.new
  y=y.new
}
x2=x^2
m2 =lm(y ~ x+x2)
intercept=m2$coef[[1]]
slope.x=m2$coef[[2]]
slope.x2=m2$coef[[3]]

#calculate mean
xval=c(0:100)*(max(x)-min(x))/100+min(x)
xval2=xval^2
yvalMean=intercept+slope.x*xval+slope.x2*xval2

#calculate model parameter/quantities
one=rep(1,length(x))
X=matrix(c(one,x,x2), nrow = length(x))
Y=matrix(y,nrow=length(y))
beta=matrix(c(m2$coefficients[[1]],m2$coefficients[[2]],m2$coefficients[[3]]),nrow=3)
XtXn=solve(t(X)%*%X)
b=XtXn%*%t(X)%*%Y
if(a==1){
  S2=(t(Y-X %*% beta) %*% (Y-X %*% beta))/(length(x)-length(beta)-1)
}else{
  S2=20*(t(Y-X %*% beta) %*% (Y-X %*% beta))/(length(x)-length(beta)-1)
}

#calculate confidence interval for the mean
alphaConf=0.05
one2=rep(1,length(xval))
XVAL=matrix(c(one2,xval,xval2),nrow=length(one2))

conf.interval<-function(XVAL,XtXn,S2,q)
{
  sol=c()
  for(i in 1:length(XVAL[,1]))
    sol=c(sol,q*sqrt(S2*t(XVAL[i,])%*%XtXn%*%XVAL[i,]))

  sol
}
q=qt(alphaConf/2,length(x)-length(beta))
delta.con=conf.interval(XVAL,XtXn,S2,q)

#calculate prediction interval for new observation
alphaPred=0.05
one2=rep(1,length(xval))
XVAL=matrix(c(one2,xval,xval2),nrow=length(one2))

pred.interval<-function(XVAL,XtXn,S2,q)
{
  sol=c()
  for(i in 1:length(XVAL[,1]))
    sol=c(sol,q*sqrt(S2*(1+t(XVAL[i,])%*%XtXn%*%XVAL[i,])))
  
  sol
}
q=qt(alphaPred/2,length(x)-length(beta))
delta.pred=pred.interval(XVAL,XtXn,S2,q)

#plot everything
par(mfrow=c(1,1))
if(a==0){
  plot(x,y, xlab="1/stress", ylab = "log(N)" , pch="-", cex=3)
  points(x.old,y.old, cex=0.6)
}else{
  plot(x,y, xlab="1/stress", ylab = "log(N)" )
}
polygon(c(xval,rev(xval)),c(yvalMean,rev(yvalMean)),border="blue",lty="solid")
polygon(c(xval,rev(xval)),c(yvalMean+delta.con,rev(yvalMean+delta.con)),border="green",lty="dashed")
polygon(c(xval,rev(xval)),c(yvalMean-delta.con,rev(yvalMean-delta.con)),border="green",lty="dashed")
polygon(c(xval,rev(xval)),c(yvalMean+delta.pred,rev(yvalMean+delta.pred)),border="red",lty="dotted")
polygon(c(xval,rev(xval)),c(yvalMean-delta.pred,rev(yvalMean-delta.pred)),border="red",lty="dotted")
legend(min(x),max(y), c("mean","confidence","prediction"),bty="n", lty=c("solid","dashed","dotted"), lwd=c(1.5,1.5),col=c("blue","green","red"))


#model3: x+x^2+x^3
x3=x^3
m3 =lm(y ~ x+x2+x3)
xval3=xval^3


#calculate model parameter/quantities
one=rep(1,length(x))
X=matrix(c(one,x,x2,x3), nrow = length(x))
Y=matrix(y,nrow=length(y))
beta=matrix(c(m3$coefficients[[1]],m3$coefficients[[2]],m3$coefficients[[3]],m3$coefficients[[4]]),nrow=4)
XtXn=solve(t(X)%*%X)
b=XtXn%*%t(X)%*%Y
if(a==1){
  S2=(t(Y-X %*% beta) %*% (Y-X %*% beta))/(length(x)-length(beta))
S2}else{
  S2=20*(t(Y-X %*% beta) %*% (Y-X %*% beta))/(length(x)-length(beta))
}
yvalMean=beta[1]+beta[2]*xval+beta[3]*xval2+beta[4]*xval3

#calculate confidence interval for the mean
alphaConf=0.05
one2=rep(1,length(xval))
XVAL=matrix(c(one2,xval,xval2,xval3),nrow=length(one2))
q=qt(alphaConf/2,length(x)-length(beta))
delta.con=conf.interval(XVAL,XtXn,S2,q)

#calculate prediction interval for new observation
alphaPred=0.05
one2=rep(1,length(xval))
XVAL=matrix(c(one2,xval,xval2,xval3),nrow=length(one2))
q=qt(alphaPred/2,length(x)-length(beta))
delta.pred=pred.interval(XVAL,XtXn,S2,q)

#plot everything
#par(mfrow=c(1,1))
if(a==0){
  plot(x,y, xlab="1/stress", ylab = "log(N)" , pch="-", cex=3)
  points(x.old,y.old, cex=0.6)
}else{
  plot(x,y, xlab="1/stress", ylab = "log(N)" )
}
polygon(c(xval,rev(xval)),c(yvalMean,rev(yvalMean)),border="blue",lty="solid")
polygon(c(xval,rev(xval)),c(yvalMean+delta.con,rev(yvalMean+delta.con)),border="green",lty="dashed")
polygon(c(xval,rev(xval)),c(yvalMean-delta.con,rev(yvalMean-delta.con)),border="green",lty="dashed")
polygon(c(xval,rev(xval)),c(yvalMean+delta.pred,rev(yvalMean+delta.pred)),border="red",lty="dotted")
polygon(c(xval,rev(xval)),c(yvalMean-delta.pred,rev(yvalMean-delta.pred)),border="red",lty="dotted")
legend(min(x),max(y), c("mean","confidence","prediction"),bty="n", lty=c("solid","dashed","dotted"), lwd=c(1.5,1.5),col=c("blue","green","red"))



#calculates the p-values for given stresslevels for a lm fit
p.val.norm<-function(data,lm){
  set=unique(data[,1])
  pval=c()
  means=c()
  pval=c()
  s=summary(lm)$sigma
  for(i in 1:length(set)){
    m=lm$coefficients[[1]]
    for(j in 2: length(lm$coefficients)){
      m=m+lm$coefficients[[j]]*(1/set[i])^(j-1)
    }
    means=c(means,m)
    pp=suppressWarnings(ks.test(log(part(data[,2],data[,1]==set[i])),'pnorm',m,s))
    pval=c(pval,pp$p.val)
  }
  df=data.frame(stress=set,mean=means,sigma2=rep(s^2,length(means)),pval)
  df
}

q=p.val.norm(shen1,m3);q

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
#GLM

data=shen1
#y=log(data[,2])
y.all=data[,2]
x.all=data[,1]

#grouped data
y=data[,2]
x=data[,1]
uni=unique(x)
x.new=c()
y.new=c()
for(counter in 1 :(length(uni))){
  dat=part(y,x==uni[counter])
  x.new=c(x.new,uni[counter])
  y.new=c(y.new,mean(dat))
}
x.old=x
y.old=y
x.mean=x.new
y.mean=y.new

mod.mean<-glm(y.mean~x.mean+I(x.mean^2),family=Gamma)
mod.all<-glm(y.all~x.all+I(x.all^2)+I(x.all^3)+I(x.all^4)+I(x.all^5),family=Gamma)

xVx<-function(xVal,V){
  nr.of.pred=(length(V))^0.5
  ret=c()
  for(counter in 1 : length(xVal)){
    vec=c(1)
    for(c.nr in 2: nr.of.pred)
      vec=c(vec,xVal[counter]^(c.nr-1))
    x.mat=matrix(vec)
    ret=c(ret,t(x.mat)%*%V%*%(x.mat))
  }
  ret
}


plot.glm<-function(mod){
  y=mod$model[[1]];
  x=mod$model[[2]];
  #inclusives intercept:
  nr.of.pred = length(mod$coefficients);
  coef = mod$coefficients;
  
  xval=c(0:100)*(max(x)-min(x))/100+min(x)
  yvalMeanLinear=rep(coef[1],length(xval))
  for(counter in 2 : length(coef)){
    yvalMeanLinear = yvalMeanLinear + coef[counter]*(xval)^(counter-1)
  }
  yvalMean=1/yvalMeanLinear;
  ymean=c(yvalMean,rev(yvalMean));
  
  #calculate X
  Xvec=rep(1,length(mod$model[[1]]))
  for(counter in 2: length(mod$model))
    Xvec=c(Xvec,mod$model[[counter]])
  X<-matrix(Xvec,nrow=length(mod$model[[1]]))
  #calculate V
  V<-matrix(rep(0,length(mod$model[[1]])^2),nrow=length(mod$model[[1]]))
  k=1/(summary(mod)$dispersion);
  #lambda = k*unique(mod$linear.predictors);
  lambda = k*(mod$linear.predictors);
  lambda.all=k*mod$linear.predictors;
  mu=k/lambda
  for(counter in 1: length(mod$model[[1]]))
    V[counter,counter]=(k*(mu[counter])^2);
  #Var Beta
  varB=solve(t(X)%*%V%*%X)
  alpha=0.05
  delta.con=qnorm(alpha/2)*xVx(xval,varB)
  conf.up=1/c(yvalMeanLinear+delta.con,rev(yvalMeanLinear+delta.con))
  conf.low=1/c(yvalMeanLinear-delta.con,rev(yvalMeanLinear-delta.con))
  
  lambda.val = k*yvalMeanLinear;
  theta.val=1/lambda.val;
  
  
  ######PLOTTING
  alpha=0.05;
  y1=c(k/lambda.val,rev(k/lambda.val));
  plot(NULL, xlim=c(min(x)/1.01,max(x)*1.01), ylim=c(min(log(y)/1.01),max(log(y))*1.01), ylab="log(N)", xlab="stress")
  y.95=c(qgamma(1-alpha/2,k,lambda.val),qgamma(1-alpha/2,k,rev(lambda.val)))
  y.05=c(qgamma(alpha/2,k,lambda.val),qgamma(alpha/2,k,rev(lambda.val)))
  polygon(c(xval,rev(xval)),log(y1),border="blue",lty="solid")
  polygon(c(xval,rev(xval)),log(y.95),border="green",lty="dashed")
  polygon(c(xval,rev(xval)),log(y.05),border="green",lty="dashed")
  #plot(x,log(y), xlab="stress", ylab = "log(N)")
  for(counter in 1:length(x)){
    if(qgamma(alpha/2,k,lambda.all[counter]) <= y[counter]  &&  qgamma(1-alpha/2,k,lambda.all[counter]) >= y[counter] )
      points(x[counter],log(y[counter]),pch="*",col="red", cex = 0.75)
    else
      points(x[counter],log(y[counter]),pch="o", col="blue", cex=0.5)
  }
  #polygon(c(xval,rev(xval)),ymean,border="blue",lty="solid")
  #polygon(c(xval,rev(xval)),conf.up,border="green",lty="dashed")
  #polygon(c(xval,rev(xval)),conf.low,border="red",lty="dashed")
  
  legend(160,14.5, c("mean","quantile lines (0.025/0.975)"),bty="n", lty=c("solid","dashed"), lwd=c(1.5,1.5),col=c("blue","green"))
  
  l=length(conf.low)/2
  df=data.frame(conf.low[1:l],ymean[1:l],conf.up[1:l],yvalMeanLinear,delta.con)
  a=list()
  a[["lines"]]=df
  a[["xval"]]=xval
  a[["varB"]]=varB
  a[["X"]]=X
  a[["V"]]=V
  a
  
}


plot.glm(mod.all);
plot2.glm(mod.mean);

plot2.glm<-function(mod){
  y=mod$model[[1]];
  x=mod$model[[2]];
  #inclusives intercept:
  nr.of.pred = length(mod$coefficients);
  coef = mod$coefficients;
  
  xval=c(0:100)*(max(x)-min(x))/100+min(x)
  yvalMeanLinear=rep(coef[1],length(xval))
  for(counter in 2 : length(coef)){
    yvalMeanLinear = yvalMeanLinear + coef[counter]*(xval)^(counter-1)
  }
  yvalMean=1/yvalMeanLinear;
  ymean=c(yvalMean,rev(yvalMean));
  
  #calculate X
  Xvec=rep(1,length(mod$model[[1]]))
  for(counter in 2: length(mod$model))
    Xvec=c(Xvec,mod$model[[counter]])
  X<-matrix(Xvec,nrow=length(mod$model[[1]]))
  #calculate V
  V<-matrix(rep(0,length(mod$model[[1]])^2),nrow=length(mod$model[[1]]))
  k=1/(summary(mod)$dispersion);
  lambda = k*unique(mod$linear.predictors);
  lambda.all=k*mod$linear.predictors;
  mu=k/lambda
  for(counter in 1: length(mod$model[[1]]))
    V[counter,counter]=(k*(mu[counter])^2);
  #Var Beta
  varB=solve(t(X)%*%V%*%X)
  alpha=0.05
  delta.con=qnorm(alpha/2)*xVx(xval,varB)
  conf.up=1/c(yvalMeanLinear+delta.con,rev(yvalMeanLinear+delta.con))
  conf.low=1/c(yvalMeanLinear-delta.con,rev(yvalMeanLinear-delta.con))
  
  lambda.val = k*yvalMeanLinear;
  theta.val=1/lambda.val;
  
  
  ######PLOTTING
  alpha=0.05;
  y1=c(k/lambda.val,rev(k/lambda.val));
  
  y.95=c(qgamma(1-alpha/2,k,lambda.val),qgamma(1-alpha/2,k,rev(lambda.val)))
  y.05=c(qgamma(alpha/2,k,lambda.val),qgamma(alpha/2,k,rev(lambda.val)))
  polygon(c(xval,rev(xval)),log(y.95),border="black",lty="dotted")
  polygon(c(xval,rev(xval)),log(y.05),border="black",lty="dotted")
  #plot(x,log(y), xlab="stress", ylab = "log(N)")
  #polygon(c(xval,rev(xval)),ymean,border="blue",lty="solid")
  #polygon(c(xval,rev(xval)),conf.up,border="green",lty="dashed")
  #polygon(c(xval,rev(xval)),conf.low,border="red",lty="dashed")
  
 legend(160,14.5, c("mean","quantile lines (0.025/0.975)","quantile lines mean (0.025/0.975)"),bty="n", lty=c("solid","dashed","dotted"), lwd=c(1.5,1.5,1.5),col=c("blue","green","black"))
  
  
  l=length(conf.low)/2
  df=data.frame(conf.low[1:l],ymean[1:l],conf.up[1:l],yvalMeanLinear,delta.con)
  
}


