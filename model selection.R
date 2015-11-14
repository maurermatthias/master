source('checkDistribution.R')


#plot data, y-log transformiert
data=shen1
data[,2]=log(data[,2])
plot(data[,1],log(data[,2]), xlab="stress", ylab = "log(N)")

#plot data, y-log transformiert / 1/X transformed
data=shen1

plot(1/data[,1],log(data[,2]), xlab="1/stress", ylab = "log(N)")

#linear model
if(1==0){
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
  x=x.new
  y=y.new
}


x2=x^2
x3=x^3
x4=x^4
x5=x^5
x6=x^6
x7=x^7
x8=x^8
x9=x^9
x10=x^10

m0=lm(y~1)
m1 =lm(y ~ x)
m2 =lm(y ~ x+x2)
m3 =lm(y ~ x+x2+x3)
m4 = lm(y ~ x+x2+x3+x4)
m5 = lm(y ~ x+x2+x3+x4+x5)
m6 = lm(y ~ x+x2+x3+x4+x5+x6)
m7 = lm(y ~ x+x2+x3+x4+x5+x6+x7)
m8 = lm(y ~ x+x2+x3+x4+x6+x7+x8)
m9 = lm(y ~ x+x2+x3+x4+x6+x7+x8+x9)
m10 = lm(y ~ x+x2+x3+x4+x6+x7+x8+x9+x10)
s0=summary(m0)
s1=summary(m1)
s2=summary(m2)
s3=summary(m3)
s4=summary(m4)
s5=summary(m5)
s6=summary(m6)
s7=summary(m7)
s8=summary(m8)
s9=summary(m9)
s10=summary(m10)
r2=c(
  s0$r.squared,
  s1$r.squared,
  s2$r.squared,
  s3$r.squared,
  s4$r.squared,
  s5$r.squared,
  s6$r.squared,
  s7$r.squared,
  s8$r.squared,
  s9$r.squared,
  s10$r.squared)
r2.adj=c(
  s0$adj.r.squared,
  s1$adj.r.squared,
  s2$adj.r.squared,
  s3$adj.r.squared,
  s4$adj.r.squared,
  s5$adj.r.squared,
  s6$adj.r.squared,
  s7$adj.r.squared,
  s8$adj.r.squared,
  s9$adj.r.squared,
  s10$adj.r.squared)

dff=data.frame(r2,r2.adj)
df2latex(dff,3)

#histogram:
#hist(m2$residuals)