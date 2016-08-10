# install.packages()

#setwd("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R")

#functions:

library('stats')
library('stats4')
library('MASS')
library('evir')
library('optimx')
library('evir')
library('ismev')
library('psych')   #df2latex  convert dataframe to latex table (http://artax.karlin.mff.cuni.cz/r-help/library/psych/html/df2latex.html)

############################################################################################################
###common functions

#convert array to hard-coded text array  
cr<-function(x){
  lines=dim(x)[1]
  col=dim(x)[2]
  txt=""
  for(i in 1:lines){
    for(j in 1:col){
      txt=paste(txt,"data[",as.character(i),",",as.character(j),"]=",as.character(x[i,j]),"; ",sep="")
    }
  }
  return(txt)
}

#takes two vektors of the same size, returns those elements of vector
#one where vector two contains true
part<-function(data,bool)
{
  if(length(bool)!=length(data))
    stop("vectors need to be of the same size!") 
  ret=c()
  for(i in 1:length(bool))
    if(bool[i])
      ret[length(ret)+1]=data[i]
  ret
}

#removes one observation from data matrix
rem.obs<-function(data,obs){
  if(obs>1 && obs<dim(data)[1]){
    d2=matrix(c(data[,1][1:(obs-1)],data[,1][(obs+1):dim(data)[1]],data[,2][1:(obs-1)],data[,2][(obs+1):dim(data)[1]]),ncol=2)
    return(d2)
  }else if(obs==1){
    d2=matrix(c(data[,1][(obs+1):dim(data)[1]],data[,2][(obs+1):dim(data)[1]]),ncol=2)
    return(d2)
  }else if(obs==dim(data)[1]){
    d2=matrix(c(data[,1][1:(obs-1)],data[,2][1:(obs-1)]),ncol=2)
    return(d2)
  }
}



############################################################################################################
###p-value calculation

#returns p-values of kolmogorov smirnof test if distribution is "dist"
#p>alpha => assumption of wbl distribution cannot be thrown away
#usage: pfree(datavector,'weibull','pweibull')
pfree<-function(data,dist,pdist){
  set=unique(data[,1])
  pval=c()
  for(i in 1:length(set)){
    dat=part(data[,2],data[,1]==set[i])
    if(pdist=='pgev')
    {
      #ge = gev.fit(dat);
      ge = gev(dat,1)  #blocksize=1
      par = ge$par.est
      xi=par[[1]]
      sigma=par[[2]]
      mu=par[[3]]
      test=ks.test(dat,pdist,xi,mu,sigma)
      #http://www.inside-r.org/packages/cran/dgof/docs/ks.test (gev)
      pval[length(pval)+1]=test$p.value 
    }else if(pdist=='pgamma'){
      param = gamma.mom(dat)
      test=ks.test(dat,pdist,param[1],param[2])
      pval[length(pval)+1]=test$p.value
    } else {
      param=fitdistr(dat,dist)
      if(length(param$estimate)==2)
      {
        shape=param$estimate[[1]]
        scale=param$estimate[[2]]
        test=ks.test(dat,pdist,shape,scale)
        pval[length(pval)+1]=test$p.value
      }
      else if (length(param$estimate)==1)
      {
        test=ks.test(dat,pdist,param$estimate[[1]])
        pval[length(pval)+1]=test$p.value
      }
      else
      {
        print("number of parameter not supportet")
      }
    }
  }
  pval
}

#calls pfree for weibull and logn distribution
p.val<-function(data){
  weibull=suppressWarnings(pfree(data,'weibull','pweibull'))
  logn=suppressWarnings(pfree(data,'lognormal','plnorm'))
  normal=suppressWarnings(pfree(data,'normal','pnorm'))
  poisson=suppressWarnings(pfree(data,'poisson','ppois'))
  gam.mom=suppressWarnings(pfree(data,'gamma','pgamma'))
  gev=suppressWarnings(pfree(data,'gev','pgev'))
  data2=data
  data2[,2]=log(data2[,2])
  #lweibull = suppressWarnings(pfree(data2,'weibull','pweibull'))
  stress = unique(data[,1])
  #r=format(round(data.frame(stress,weibull,logn,normal,gev),3),nsmall = 3)
  r=data.frame(stress,weibull,logn,normal,gam.mom,poisson,gev)
}

#eins=fitdistr(part(shen1[,2],shen1[,1]==unique(shen1[,1])[1]),'gamma')
#t.eins=ks.test(part(shen1[,2],shen1[,1]==unique(shen1[,1])[1]),'pgamma',eins$estimate[1],eins$estimate[2])
#

############################################################################################################
###Parameter estimation

#ans=gamma.mom(part(shen1[,2],shen1[,1]==unique(shen1[,1])[1]))
#t.ans=ks.test(part(shen1[,2],shen1[,1]==unique(shen1[,1])[1]),'pgamma',ans[1],1/ans[2])
gamma.mom<-function(data){
  m1=sum(data)/length(data)
  m2=sum(data^2)/length(data)
  l=m1^2/(m2-m1^2)
  a=m1/(m2-m1^2)
  ret=c(l,a)
  ret
}

##estimating GEV parameter for one level of stress
gev.est<-function(data){
  u = unique(data[,1])
  xi=c()
  sigma=c()
  mu=c()
  for(i in 1:length(u))
  {
    dat=part(data[,2],data[,1]==u[i])
    ge = gev(dat,1)  #blocksize=1
    par = ge$par.ests
    xi=c(xi,par[1])
    sigma=c(sigma,par[2])
    mu=c(mu,par[3])
  }
  p.est=data.frame(u,xi,sigma,mu)
}

#estimate params for distribution for each stress level
#for 2 param-distribution
parameter<-function(data,dist){
  stresslevel=unique(data[,1])
  param1=c()
  param2=c()
  for(i in 1:length(stresslevel)){
    dat=part(data[,2],data[,1]==stresslevel[i])
    param=fitdistr(dat,dist)
    param1[length(param1)+1]=param$estimate[1]
    param2[length(param2)+1]=param$estimate[2]
  }
  pval=data.frame(stresslevel,param1,param2)
}
#calls fkt parameter for weibull and logn distribution
param<-function(data){
  stress=unique(data[,1])
  weibull=suppressWarnings(parameter(data,'Weibull'))
  wbl.shape=weibull$param1
  wbl.scale=weibull$param2
  logn=parameter(data,'lognormal')
  logn.ew=logn$param1
  logn.var=logn$param2
  #gev=suppressWarnings(gev(data,1))
  #gev.xi=gev$par.ests[[1]]
  #gev.sigma=gev$par.ests[[2]]
  #gev.mu=gev$par.ests[[3]]
  #d=data.frame(stress,wbl.shape,wbl.scale,logn.ew,logn.var,gev.xi,gev.sigma,gev.mu)
  d=data.frame(stress,wbl.shape,wbl.scale,logn.ew,logn.var)
}

#negative log-log likelihood weibull function (n sample size)
nllwbl<-function(shape, scale){
  data=dataForFit
  n=length(data)
  ret=0
  for(i in 1: n){
    ret=ret-log(log(scale*shape))-log((shape-1)*log((scale*data[i])))+shape*log((scale*data[i]))
  }
  ret
}

#calculates mle estimators of weibull distribution
mlewbl<-function(data){
  param=suppressWarnings(fitdistr(dat,'Weibull'))
  p1=param$estimate[1]
  p2=param$estimate[2]
  dataForFit=data
  est=mle(nllwbl,list(shape=p1,scale=p2))
}


############################################################################################################
###visualisation functions

plot.latex.table<-function(data){
  stresslevel=unique(data[,1])
  reflen=length(part(data[,2],data[,1]==stresslevel[2]))
  nr=1:reflen
  df=data.frame(nr)
  for(i in 1:length(stresslevel)){
    dat=part(data[,2],data[,1]==stresslevel[i])/100
    if(length(dat) != reflen){
      dif=reflen-length(dat)
      for(j in 1:dif){
        dat=c(dat,99999)
      }
    }
    df[,as.character(stresslevel[i])] <- dat
  }
  df$nr <- NULL
  df2latex(df)
  #df
}

plot.par<-function(data){

}
