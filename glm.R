source('function.R')
source('plot work.R')
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


#poisson model
data=shen1
y=log(data[,2])
x=data[,1]
m.poi1<-glm(y~x,family=poisson)
m.poi2<-glm(y~x+I(x^2),family=poisson)
m.poi3<-glm(y~x*I(x^2),family=poisson)
summary(m.poi2)
anova(m.poi1,m.poi2,m.poi3, test="Chisq")

#gamma model
data=shen1
#y=log(data[,2])
y=data[,2]
x=data[,1]
a=0
if(a==1){
  y=data[,2]
  x=data[,1]
}else{
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
  x=x.new
  y=y.new
}
m.gam1<-glm(y~x,family=Gamma)
m.gam2<-glm(y~x+I(x^2),family=Gamma)
m.gam3<-glm(y~x+I(x^2)+I(x^3),family=Gamma)
m.gam4<-glm(y~x+I(x^2)+I(x^3)+I(x^4),family=Gamma)
m.gam5<-glm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5),family=Gamma)
m.gam6<-glm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6),family=Gamma)
m.gam7<-glm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6)+I(x^7),family=Gamma)
m.gam8<-glm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6)+I(x^7)+I(x^8),family=Gamma)
anova(m.gam1,m.gam2,m.gam3,m.gam4, test="Chisq")

v.aic=c(m.gam1$aic,m.gam2$aic,m.gam3$aic,m.gam4$aic,m.gam5$aic,m.gam6$aic,m.gam7$aic,m.gam8$aic)

m.par1=m.gam1$coefficients
m.par2=m.gam2$coefficients
m.par3=m.gam3$coefficients
m.par4=m.gam4$coefficients

summary(m.gam1)
summary(m.gam2)

fit = m.gam1$fitted.values

beta=m.gam1$coefficients
#linear predictor nu
nu.1=beta[[1]]+beta[[2]]*f1$stress
#faster
nu.2 =  m.gam1$linear.predictors

#dispersion parameter
fi.1 = sum((f1$mu-m.gam1$fitted.values)^2/m.gam1$fitted.values^2)/8
#faster
fi.2 = summary(m.gam1)$dispersion


#param = gamma.mom(dat)


#calculate parameter for the gamma distribution 
g.par<-function(data,mod1,mod2){
  k.ind=c();
  lambda.ind=c()
  k1=1/(summary(mod1)$dispersion);
  k2=1/(summary(mod2)$dispersion);
  lambda.m1 = k1*unique(mod1$linear.predictors);
  lambda.m2 = k2*unique(mod2$linear.predictors);
  k.m1=rep(k1,length(lambda.m1));
  k.m2=rep(k2,length(lambda.m2));
  
  set=unique(data[,1])
  for(i in 1 : length(set)){
    dat=part(data[,2],data[,1]==set[i])
    param = gamma.mom(dat)
    k.ind=c(k.ind,param[1])
    lambda.ind=c(lambda.ind,param[2])
  }
  df=data.frame(set, k.ind,1/lambda.ind,k.ind/lambda.ind, k.m1, 1/lambda.m1,k.m1/lambda.m1, k.m2, 1/lambda.m2,k.m2/lambda.m2)
}

#calculate Gamma Parameter for GLM Model
getPar<-function(mod,data){
  k1=1/(summary(mod)$dispersion);
  lambda = k1*unique(mod$linear.predictors);
  k=rep(k1,length(lambda))
  
  set=unique(data[,1])
  pval=c()
  for(i in 1:length(set)){
    dat=part(data[,2],data[,1]==set[i])
    param = gamma.mom(dat)
    test=suppressWarnings(ks.test(dat,"pgamma",k[i],lambda[i]))
    pval[length(pval)+1]=round(test$p.value,3)
  }
  
  df = data.frame(stress=set,k,lambda,theta=1/lambda,pval)
}
r=getPar(m.gam1,shen1); r
m1=getPar(m.gam1,shen1)$pval
m2=getPar(m.gam2,shen1)$pval
m3=getPar(m.gam3,shen1)$pval
m4=getPar(m.gam4,shen1)$pval
m5=getPar(m.gam5,shen1)$pval
m6=getPar(m.gam6,shen1)$pval
m7=getPar(m.gam7,shen1)$pval
m8=getPar(m.gam8,shen1)$pval
pvalues=data.frame(stress=r$stress,m1,m2,m3,m4,m5,m6,m7,m8)
df2latex(pvalues)


#plot model xx
model = m.gam2
names(model)
w1=model$fitted.values-model$y
w2=model$residuals
xval=c(0:100)*(max(x)-min(x))/100+min(x)
yvalMeanLinear=model$coefficients[[1]]+model$coefficients[[2]]*xval
if(length(model$coefficients)>2){
  for(int in 3:(length(model$coefficients))){
    yvalMeanLinear = yvalMeanLinear + model$coefficients[[int]]*xval^(int-1)
  }
}
yvalMean=1/yvalMeanLinear

#for gamma distribution:
#a=1/k for X_i ~ Beta(k,l)
#theta=-l/k
a<-function(){return 2;}
#b(t)=-log(-t) => ddb(x)=-(1/x^2)
ddb<-function(x){ret = -(1/x)^(2); }
#g(t)=(b')^-1=(1/t)^-1=1/t => dg(x)=-1/x^2
db<-function(x){ret=-1/x^2;}

one=rep(1,length(x))
X=matrix(c(one,x,x2), nrow = length(x))
Y=matrix(y,nrow=length(y))
XtXn=solve(t(X)%*%X)
b=XtXn%*%t(X)%*%Y


if(a==0){
  plot(x,y, xlab="1/stress", ylab = "log(N)" , pch="-", cex=3)
  points(x.old,y.old, cex=0.6)
}else{
  plot(x,y, xlab="1/stress", ylab = "log(N)" )
}
polygon(c(xval,rev(xval)),c(yvalMean,rev(yvalMean)),border="blue",lty="solid")
#polygon(c(xval,rev(xval)),c(yvalMean+delta.con,rev(yvalMean+delta.con)),border="green",lty="dashed")
#polygon(c(xval,rev(xval)),c(yvalMean-delta.con,rev(yvalMean-delta.con)),border="green",lty="dashed")
#polygon(c(xval,rev(xval)),c(yvalMean+delta.pred,rev(yvalMean+delta.pred)),border="red",lty="dotted")
#polygon(c(xval,rev(xval)),c(yvalMean-delta.pred,rev(yvalMean-delta.pred)),border="red",lty="dotted")
#legend(min(x),max(y), c("mean","confidence","prediction"),bty="n", lty=c("solid","dashed","dotted"), lwd=c(1.5,1.5),col=c("blue","green","red"))
legend(min(x),max(y), c("mean"),bty="n", lty=c("solid"), lwd=c(1.5,1.5),col=c("blue"))

