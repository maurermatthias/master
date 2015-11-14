source('basicMoments.R')

#overall normality test
param=fitdistr(shen1[,2],"normal")
shape=param$estimate[1]
scale=param$estimate[2]
test=ks.test(shen1[,2],"pnorm",shape,scale)


alpha=0.05
n=length(x)
data=shen1
data[,1]=1/shen1[,1]
data[,2]=log(shen1[,2])
moments=get.moments(data)
x=data[,1]
y=data[,2]
#b1=sum((x-mean(x))*(y-mean(y)))/sum((x-mean(x))^2)
#b0=mean(y)-b1*mean(x)
#the same as above:
fm <- lm(y~x)
intercept=fm$coef[[1]]
slope=fm$coef[[2]]
plot(x, y, xlab="1/Stress", ylab = "log(N)", main="Dataset shen1 - normal distribution")
#abline(intercept,slope, col="red")
polygon(c(min(x),max(x)),c(intercept+slope*min(x),intercept+slope*max(x)),border="red")   #only valid for linear approach
mse=sum((y-mean(y))^2)/(n-2)
pos=(0:100)*((max(x)-min(x))/100)+min(x)
f=qf(1-alpha/2,2,n-2)
sx=sum((x-mean(x))^2)
dist=sqrt(2*f*mse*(1/n+(pos-mean(x))^2/sx))
means=intercept+slope*pos       #only valid for linear approach
polygon(c(pos,rev(pos)),c(means+dist,rev(means+dist)),border="green")
polygon(c(pos,rev(pos)),c(means-dist,rev(means-dist)),border="green")
legend(min(x),max(y), c("mean","Working-Hotelling"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red","green"))
