
#########GEV
#defining the input list
gev=list();
gev[["distr"]]="gev"
data=load.data()
gev[["xval"]]=data[,1]
gev[["yval"]]=data[,2]
gev[["control"]]=list(maxit=20000);

#defining additional input-list fields for fitting
gev[["type"]]="fit"
gev[["error.type"]]="rel"
gev[["validity.fun"]]="val.gev"
gev[["struct.fun"]]=c("gev1","gev2","gev3")
gev[["struct.start.parameter"]]=c(-0.031,49.67,20809,18.8,69.41)
#gev[["struct.start.parameter"]]=c(-0.2157395,49.6430314,22760.4432796,17.2437222,63.4111482)
gev[["quantiles"]]=1:9/10

#validity function
val.gev<-function(stress,parameter){
  k=parameter[1];
  a=parameter[2];
  b=parameter[3];
  c1=parameter[4];
  c2=parameter[5];
  
  if(b<=max(stress) || a>=min(stress) || c1<=0 || c2<=0 || k==0 || min(gev2(stress,parameter))<=0 ){
    return(FALSE);
  }else{
    return(TRUE);
  }
}


#define parameter-regression function, dependendt on structure parameter 
#xi gev-distibution
gev1<-function(stress,parameter){
  k=parameter[1];
  return(k);
}

#sigma^2 gev-distribution
gev2<-function(stress,parameter){
  a=parameter[2];
  b=parameter[3];
  c1=parameter[4];
  return(((b-a)/(stress-a)-1)*c1)
}

#mu gev-distribution
gev3<-function(stress,parameter){
  a=parameter[2];
  b=parameter[3];
  c2=parameter[5];
  return(((b-a)/(stress-a)-1)*c2)
}

gev[["ML"]]=TRUE

#perform simulation (observe)
simGev=list()
simGev[["times"]]=10000
simGev[["ratio"]]=0.9
simGev[["plot"]]=FALSE
#v=pr.sim(gev, sim)

#perform simulation (generate)
simGev[["type"]]="generate"
simGev[["par"]]=c(-0.22,49.64,22760,17.25,63.4)
simGev[["xval"]]=unique(data[,1])
simGev[["n"]]=rep(20,length(simGev[["xval"]]))


#b=pr.sim(gev,sim)
#saveRDS(b, file="simGev10000ml")




#####################MORMAL
#normal
norm=list();
norm[["distr"]]="norm"
norm[["control"]]=list(maxit=1000);
norm[["xval"]]=data[,1]
norm[["yval"]]=data[,2]
norm[["ML"]]=TRUE
#norm[["validity.fun"]]="val.norm"
norm[["struct.fun"]]=c("norm1","norm2")
norm[["struct.start.parameter"]]=c(0.9*min(data[,1]),1.1*max(data[,1]),1000,1000)

#validity function
val.norm<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  c2=parameter[4];
  sigma.min=min(norm2(stress,parameter))
  if(b<=max(stress) || a>=min(stress) || c1<=0 || c2<=0 || sigma.min <= 0){
    return(FALSE);
  }else{
    return(TRUE);
  }
}

#mu normal
norm1<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  return(((b-a)/(stress-a)-1)*c1)
}

#sigma^2 normal
norm2<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c2=parameter[4];
  return(((b-a)/(stress-a)-1)*c2)
}

#norm.result=pr(norm)

#perform simulation (observe)
simNorm=list()
simNorm[["times"]]=10000
simNorm[["ratio"]]=0.9
simNorm[["plot"]]=FALSE
#v=pr.sim(gev, sim)

#perform simulation (generate)
simNorm[["type"]]="generate"
simNorm[["par"]]=c(49.6,346,7261,1401)
simNorm[["xval"]]=unique(data[,1])
simNorm[["n"]]=rep(20,length(simNorm[["xval"]]))

norm[["error.type"]]="rel"
#aa=pr.sim(norm,sim)
#saveRDS(aa, file="simNorm10000ml")



###################LOGN

logn=list();
logn[["distr"]]="logn"
logn[["error.type"]]="wei"
logn[["xval"]]=data[,1]
logn[["yval"]]=data[,2]
logn[["control"]]=list(maxit=20000);

logn[["ML"]]=TRUE
logn[["struct.fun"]]=c("logn1","logn2")
logn[["validity.fun"]]="val.logn"
logn[["quantiles"]]=1:9/10
logn[["struct.start.parameter"]]=c(0.9*min(data[,1]),1.1*max(data[,1]),10,0.3)

#validity function
val.logn<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  c2=parameter[4];
  if(b<=max(stress) || a>=min(stress) || c1<=0 || c2<=0){
    return(FALSE);
  }else{
    return(TRUE);
  }
}

#mu logn
logn1<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  return(((b-a)/(stress-a)-1)*c1+9)
}

#sigma^2 logn
logn2<-function(stress,parameter){
  c2=parameter[4];
  return(c2)
}

#logn.result=pr(logn)


#perform simulation (observe)
simLogn=list()
simLogn[["type"]]="observe"
simLogn[["times"]]=10000
simLogn[["ratio"]]=0.9
simLogn[["plot"]]=FALSE

#perform simulation (generate)
simLogn[["type"]]="generate"
simLogn[["par"]]=c(47.3,323,0.092,0.293)
simLogn[["xval"]]=unique(data[,1])
simLogn[["n"]]=rep(20,length(simLogn[["xval"]]))

#logn[["error.type"]]="rel"
#aaa=pr.sim(logn,sim)
#saveRDS(aaa, file="simLogn10000ml")

###################do simulations
n=100
times=10000

#NORM
norm[["error.type"]]="rel"
norm[["ML"]]=FALSE
norm[["times"]]=times
simNorm[["n"]]=rep(n,length(simNorm[["xval"]]))
normRel=pr.sim(norm,simNorm)
saveRDS(normRel,"NormRel")

norm[["error.type"]]="abs"
norm[["ML"]]=FALSE
norm[["times"]]=times
simNorm[["n"]]=rep(n,length(simNorm[["xval"]]))
normRel=pr.sim(norm,simNorm)
saveRDS(normRel,"NormAbs")

norm[["error.type"]]="wei"
norm[["ML"]]=FALSE
norm[["times"]]=times
simNorm[["n"]]=rep(n,length(simNorm[["xval"]]))
normRel=pr.sim(norm,simNorm)
saveRDS(normRel,"NormWei")

norm[["ML"]]=TRUE
norm[["times"]]=times
simNorm[["n"]]=rep(n,length(simNorm[["xval"]]))
normRel=pr.sim(norm,simNorm)
saveRDS(normRel,"NormMl")

#LNORM
logn[["error.type"]]="abs"
logn[["ML"]]=FALSE
logn[["times"]]=times
simLogn[["n"]]=rep(n,length(simLogn[["xval"]]))
val=pr.sim(logn,simLogn)
saveRDS(val,"LognAbs")

logn[["error.type"]]="rel"
logn[["ML"]]=FALSE
logn[["times"]]=times
simLogn[["n"]]=rep(n,length(simLogn[["xval"]]))
val=pr.sim(logn,simLogn)
saveRDS(val,"LognRel")

logn[["error.type"]]="wei"
logn[["ML"]]=FALSE
logn[["times"]]=times
simLogn[["n"]]=rep(n,length(simLogn[["xval"]]))
val=pr.sim(logn,simLogn)
saveRDS(val,"LognWei")

logn[["ML"]]=TRUE
logn[["times"]]=times
simLogn[["n"]]=rep(n,length(simLogn[["xval"]]))
val=pr.sim(logn,simLogn)
saveRDS(val,"LognMl")

#GEV
gev[["error.type"]]="abs"
gev[["ML"]]=FALSE
gev[["times"]]=times
simGev[["n"]]=rep(n,length(simGev[["xval"]]))
val=pr.sim(gev,simGev)
saveRDS(val,"GevAbs")

gev[["error.type"]]="rel"
gev[["ML"]]=FALSE
gev[["times"]]=times
simGev[["n"]]=rep(n,length(simGev[["xval"]]))
val=pr.sim(gev,simGev)
saveRDS(val,"GevRel")

gev[["error.type"]]="wei"
gev[["ML"]]=FALSE
gev[["times"]]=times
simGev[["n"]]=rep(n,length(simGev[["xval"]]))
val=pr.sim(gev,simGev)
saveRDS(val,"GevWei")

gev[["ML"]]=TRUE
gev[["times"]]=times
simGev[["n"]]=rep(n,length(simGev[["xval"]]))
val=pr.sim(gev,simGev)
saveRDS(val,"GevMl")

####################GAMMA

#gam=list();
#gam[["distr"]]="gamma"
#gam[["error.type"]]="rel"
#gam[["error.type"]]="wei"
#gam[["xval"]]=data[,1]
#gam[["yval"]]=data[,2]
#gam[["validity.fun"]]="val.gam"
#gam[["struct.fun"]]=c("gam1","gam2")
#gam[["struct.start.parameter"]]=c(49.67,20809,18.8,20,0)
#gam[["quantiles"]]=1:9/10
#gam[["ML"]]=TRUE
#gam[["control"]]=list(maxit=20000)

#validity function
#val.gam<-function(stress,parameter){
#  a=parameter[1];
#  b=parameter[2];
#  c1=parameter[3];
#  k=parameter[4];
#  d=parameter[5];
#  if(b<=max(stress) || a>=min(stress) || c1<=0 ){
#    return(FALSE);
#  }else{
#    return(TRUE);
#  }
#}


#k
#gam1<-function(stress,parameter){
#  k=parameter[4];
#  d=parameter[5];
#  return(k*stress+d);
#}

#theta
#gam2<-function(stress,parameter){
#  a=parameter[1];
#  b=parameter[2];
#  c1=parameter[3];
#  return(((b-a)/(stress-a)-1)*c1)
#}



#gam.result=pr(gam)

#perform simulation (observe)
#sim=list()
#sim[["type"]]="observe"
#sim[["times"]]=10000
#sim[["ratio"]]=0.9
#sim[["plot"]]=FALSE

#perform simulation (generate)
#sim[["type"]]="generate"
#sim[["par"]]=c(42.9,2283,14.14,-0.01048,21.32)
#sim[["xval"]]=unique(data[,1])
#sim[["n"]]=rep(20,length(sim[["xval"]]))

#gam[["error.type"]]="rel"
#aaaa=pr.sim(gam,sim)
#saveRDS(aaaa, file="simGam10000rel")


##############################################
#plot changes

plot.changes<-function(mat,par=NULL){
  for(j in 1: length(mat[1,])){
    m=c()
    s=c()
    for(i in 2:length(mat[,1])){
      v=mat[1:i,j]
      m=c(m,mean(v))
      s=c(s,sd(v))
    }
    plot.new()
       
    
    #plot first parameter
    par(fig=c(0,1,0.4,1),new=TRUE)
    if(!is.null(par)){
      ymin=round(min(m,par[j])*100)/100
      ymax=round(max(m,par[j])*100)/100
      plot(m,ylab="mean",xlab="",main=paste("Parameter ",as.character(j),sep=""),type="l",ylim=c(ymin,ymax))
      abline(h=par[j],col="red")
    }else{
      plot(m,ylab="mean",xlab="",main=paste("Parameter ",as.character(j),sep=""),type="l")
    }
    par(fig=c(0,1,0,0.6), new=TRUE)
    plot(s,ylab="sd",xlab="",type="l")
    cat ("Press [enter] to continue")
    line <- readline()
  }
}

#plot.changes(a,sim[["par"]])

#plot result
#gamma
#parGam=c(42.9,2283,14.14,-0.01048,21.32)
#simGam10000wei=readRDS("simGam10000wei")
#plot.changes(simGam10000wei,parGam);
#simGam10000rel=readRDS("simGam10000rel")
#plot.changes(simGam10000rel,parGam);
#simGam10000abs=readRDS("simGam10000abs")
#plot.changes(simGam10000abs,parGam);
#gev
#parGev=c(-0.31,49.67,20809,18.8,69.41)
#simGev10000wei=readRDS("n50times10000/GevMl")
#plot.changes(simGev10000wei,parGev);
#simGev10000rel=readRDS("simGam10000rel")
#plot.changes(simGev10000rel,parGev);
#simGev10000abs=readRDS("simGam10000abs")
#plot.changes(simGev10000abs,parGev);
#norm
#data=load.data()
#parNorm=c(0.9*min(data[,1]),1.1*max(data[,1]),1000,1000)
#simNorm10000wei=readRDS("simNorm10000wei")
#plot.changes(simNorm10000wei,parNorm);
#simNorm10000rel=readRDS("simNorm10000rel")
#plot.changes(simNorm10000rel,parNorm);
#simNorm10000abs=readRDS("simNorm10000abs")
#plot.changes(simNorm10000abs,parNorm);
#lnorm
#parLogn=c(0.9*min(data[,1]),1.1*max(data[,1]),10,0.3)
#simLogn10000wei=readRDS("simLogn10000wei")
#plot.changes(simLogn10000wei,parLogn);
#simLogn10000rel=readRDS("simLogn10000rel")
#plot.changes(simLogn10000rel,parLogn);
#simLogn10000abs=readRDS("simLogn10000abs")
#plot.changes(simLogn10000abs,parLogn);

#simLogn10000ml=readRDS("simLogn10000ml")
#plot.changes(simLogn10000ml,parLogn);
#simNorm10000ml=readRDS("simNorm10000ml")
#plot.changes(simNorm10000ml,parNorm);
#simGev10000ml=readRDS("simGev10000ml")
#plot.changes(simGev10000ml,parGev);


#############################################################################
####GEV mit \xi = 0
#Sinnlos: wollen nicht fitten sondern schaun ob fit ok ist

gev=list();
gev[["distr"]]="gev"
data=load.data()
gev[["xval"]]=data[,1]
gev[["yval"]]=data[,2]
gev[["control"]]=list(maxit=20000);
gev[["type"]]="fit"
gev[["validity.fun"]]="val.gev"
gev[["struct.fun"]]=c("gev1","gev2","gev3")
gev[["struct.start.parameter"]]=c(49.67,20809,18.8,69.41)
#gev[["struct.start.parameter"]]=c(-0.2157395,49.6430314,22760.4432796,17.2437222,63.4111482)
gev[["quantiles"]]=1:9/10

#validity function
val.gev<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  c2=parameter[4];
  
  if(b<=max(stress) || a>=min(stress) || c1<=0 || c2<=0  || min(gev2(stress,parameter))<=0 ){
    return(FALSE);
  }else{
    return(TRUE);
  }
}


#define parameter-regression function, dependendt on structure parameter 
#xi gev-distibution
gev1<-function(stress,parameter){
  k=0;
  return(k);
}

#sigma^2 gev-distribution
gev2<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  return(((b-a)/(stress-a)-1)*c1)
}

#mu gev-distribution
gev3<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c2=parameter[4];
  return(((b-a)/(stress-a)-1)*c2)
}



#perform simulation (observe)
simGev=list()
simGev[["times"]]=10000
simGev[["ratio"]]=0.9
simGev[["plot"]]=FALSE
#v=pr.sim(gev, sim)

#perform simulation (generate)
simGev[["type"]]="generate"
simGev[["par"]]=c(-0,49.64,22760,17.25,63.4)
simGev[["xval"]]=unique(data[,1])



gev[["ML"]]=TRUE
gev[["error.type"]]="rel"
simGev[["n"]]=rep(20,length(simGev[["xval"]]))

#b=pr.sim(gev,sim)
#saveRDS(b, file="simGev10000ml")
