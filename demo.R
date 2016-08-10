install.packages("par.reg_1.0.tar.gz", repos=NULL,type="source")
library(par.reg)


data=load.data()
plot(data[,1],log(data[,2]))


gev=list();
gev[["distr"]]="gev"
gev[["xval"]]=data[,1]
gev[["yval"]]=data[,2]


#print individual parameter estimations - decide on the parameter functions
gev[["type"]]="diag"
gev.result=pr(gev)


#defining additional input-list fields for fitting
gev[["type"]]="fit"
gev[["xval"]]=data[,1]
gev[["yval"]]=data[,2]
gev[["error.type"]]="rel"
gev[["validity.fun"]]="val.gev"
gev[["struct.fun"]]=c("gev1","gev2","gev3")
gev[["struct.start.parameter"]]=c(-0.31,49.67,20809,18.8,69.41)
gev[["quantiles"]]=1:10/20

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


#xi gev-distribution
gev1<-function(stress,parameter){
  k=parameter[1];
  return(k);
}

#sigma gev-distribution
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

gev[["control"]]=list(maxit=1000)

#perform fit
gev.result=pr(gev)

#perform simulation with observed data
sim=list()
sim[["type"]]="observe"
sim[["times"]]=50
sim[["ratio"]]=0.9
v=pr.sim(gev, sim)

#perform simulation with generated data
sim[["type"]]="generate"
sim[["par"]]=gev.result$struct.par.est
sim[["xval"]]=unique(data[,1])
sim[["n"]]=rep(20,length(sim[["xval"]]))

a=pr.sim(gev,sim)
