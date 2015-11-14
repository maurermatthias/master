source('function.R')
library(evir)
library(stats4)

#load R-Data:
shen1=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
shen2=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")


x=shen1[,1]
y=shen1[,2]
lx=log(x)
ly=log(y)


#mode: ly=a+b*lx+epsilon

m1=lm(ly~lx)
a=m1$coefficients[[1]]
b=m1$coefficients[[2]]

plot(lx,ly)
abline(a,b,col="red")







#vers.mathematische abschätzen  y->y^(-b) ; y reißt nach unten aus -> yy reißt nach oben aus

#b=1
x=shen1[,1]
y=1/shen1[,2]
plot(x,y)

#gev-distribution
gev.ind.est<-function(x,y){
  xi=c()
  sigma=c()
  mu=c()
  stress=c()
  for(i in 1:length(unique(x))){
    parTmp=c(NA,NA,NA)
    try({
      ytmp=part(y,x==unique(x)[i]);
      parTmp=suppressWarnings(gev(ytmp)$par.ests)
      xi=c(xi,parTmp[[1]])
      sigma=c(sigma,parTmp[[2]])
      mu=c(mu,parTmp[[3]])
      stress=c(stress,unique(x)[i])
    },TRUE)
  }
  par=data.frame(stress,xi,mu,sigma)
  return(par)
}

##parameter-estimate function vor gev-input, when having same xi
#gev.xi.same()$y.std contains standaduzed values (~iid)
gev.xi.same<-function (data, control=NA) 
{
  ######calculate start value for optim
  get.start.val<-function(data){
    x=data[,1]
    stress=unique(x)
    y=data[,2]
    sigma=c()
    mu=c()
    for(j in 1:length(stress)){
      y.part=part(y,x==stress[j])
      sigma.tmp=sqrt(6 * var(y.part))/pi
      mu.tmp <- mean(y.part) - 0.57722 * sigma.tmp
      sigma=c(sigma,sigma.tmp)
      mu=c(mu,mu.tmp)
    }
    xi = 0.1
    return(c(xi,sigma,mu))
  }
  x=data[,1]
  y=data[,2]
  theta <- get.start.val(data)
  #####negativeloglikelihod function vor gev-input, when having same xi
  #tmp contains x and y values
  negloglik <- function(theta, tmp) {
    x=tmp[,1]
    stress=unique(x)
    y=tmp[,2]
    
    xi=theta[1]
    sigma.v=theta[2:(length(stress)+1)]
    mu.v=theta[(length(stress)+2):(2*length(stress)+1)]
    
    out=0
    for(j in 1: length(stress)){
      mu=mu.v[j]
      sigma=sigma.v[j]
      y.part=part(y,x==stress[j])
      v <- 1 + (xi * (y.part - mu))/sigma
      #print(v)
      #print(sigma)
      if ((sigma < 0) || (min(v) < 0)) 
        out = out + 1e+06
      else {
        term1 <- length(y.part) * logb(sigma)
        term2 <- sum((1 + 1/xi) * logb(v))
        term3 <- sum(v^(-1/xi))
        out = out + term1 + term2 + term3
      }
    }
    out
  }
  if (is.na(control)) {
    fit = optim(theta, negloglik, tmp = data)
  } else {
    fit = optim(theta, negloglik, tmp = data,  control=control) 
  }
  if (fit$convergence) 
    warning("optimization may not have succeeded")
  par.ests.all <- fit$par
  par.ests = list(xi=par.ests.all[1], sigma=par.ests.all[1:length(unique(x))+1], mu = par.ests.all[(length(unique(x))+2):(2*length(unique(x))+1)])
  standardize<-function(x,y,xi,sigma,mu){
    y.new=c()
    stress=unique(x)
    for(j in 1:length(x)){
      pos=which(x[j]==stress)
      y.new=c(y.new,(y[j]-mu[pos])/sigma[pos])
    }
    return(y.new)
  }
  y.new=standardize(x,y,par.ests.all[1],par.ests.all[1:length(unique(x))+1],par.ests.all[(length(unique(x))+2):(2*length(unique(x))+1)])
  out <- list(x=x, y=y, y.std=y.new,
              par.ests = par.ests, counts=fit$counts,  
              converged = fit$convergence, nllh.final = fit$value)
  out
}

ks.test.all=function(x,y){
  stress=unique(x)
  m=matrix(rep(NA,length(stress)^2),nrow=length(stress))
  for(j in 1:length(stress)){
    d1=part(y,x==stress[j])
    for(i in 1:length(stress)){
      d2=part(y,x==stress[i])
      m[i,j]=suppressWarnings(ks.test(d1,d2)$p.value)
    }
  }
  return(m)
}

#stdandardize values
control = list(maxit=12000)
z=gev.xi.same(shen1, control)

#ks-test
ks=ks.test.all(shen1[,1],z$y.std)
ks

#fit std values
ind.all=gev.ind.est(shen1[,1],z$y.std)
ind=gev(z$y.std)$par.ests

#plot std values
plot(shen1[,1],z$y.std)
abline(0,0,col="red")



