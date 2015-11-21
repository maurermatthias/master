source('ParameterRegression.R')

###################################################################
#TEST:






#normal
norm=list();
norm[["distr"]]="norm"
norm[["control"]]=list(maxit=20000);
norm[["xval"]]=shen1[,1]
norm[["yval"]]=shen1[,2]
#norm[["validity.fun"]]="val.norm"
norm[["struct.fun"]]=c("norm1","norm2")
norm[["struct.start.parameter"]]=c(0.9*min(shen1[,1]),1.1*max(shen1[,1]),1000,1000)

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

norm.result=pr(norm)
pr.parplot(norm.result)




######################################################################################
#lognormal
logn=list();
logn[["distr"]]="logn"
logn[["error.type"]]="wei"
logn[["xval"]]=shen1[,1]
logn[["yval"]]=shen1[,2]
logn[["struct.fun"]]=c("logn1","logn2")
logn[["validity.fun"]]="val.logn"
logn[["quantiles"]]=1:9/10
logn[["struct.start.parameter"]]=c(0.9*min(shen1[,1]),1.1*max(shen1[,1]),10,0.3)

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

logn.result=pr(logn)
pr.parplot(logn.result)


#################################################################################
#gev
gev=list();
gev[["distr"]]="gev"
gev[["error.type"]]="rel"
gev[["xval"]]=shen1[,1]
gev[["yval"]]=shen1[,2]
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
  if(b<=max(stress) || a>=min(stress) || c1<=0 || c2<=0 || k==0){
    return(FALSE);
  }else{
    return(TRUE);
  }
}


#xi gev
gev1<-function(stress,parameter){
  k=parameter[1];
  return(k);
}

#sigma^2 gev
gev2<-function(stress,parameter){
  a=parameter[2];
  b=parameter[3];
  c1=parameter[4];
  return(((b-a)/(stress-a)-1)*c1)
}

#mu gev
gev3<-function(stress,parameter){
  a=parameter[2];
  b=parameter[3];
  c2=parameter[5];
  return(((b-a)/(stress-a)-1)*c2)
}

gev.result=pr(gev)
pr.parplot(gev.result)

gev[["quantiles"]]=1:9/10
v=pr.sim(gev, 0.5, 500)



###########################################################
#Gamma

gam=list();
gam[["distr"]]="gamma"
#gam[["error.type"]]="rel"
#gam[["error.type"]]="wei"
gam[["xval"]]=shen1[,1]
gam[["yval"]]=shen1[,2]
gam[["validity.fun"]]="val.gam"
gam[["struct.fun"]]=c("gam1","gam2")
gam[["struct.start.parameter"]]=c(49.67,20809,18.8,20,0)
gam[["quantiles"]]=1:9/10

#validity function
val.gam<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  k=parameter[4];
  d=parameter[5];
  if(b<=max(stress) || a>=min(stress) || c1<=0 ){
    return(FALSE);
  }else{
    return(TRUE);
  }
}


#k
gam1<-function(stress,parameter){
  k=parameter[4];
  d=parameter[5];
  return(k*stress+d);
}

#theta
gam2<-function(stress,parameter){
  a=parameter[1];
  b=parameter[2];
  c1=parameter[3];
  return(((b-a)/(stress-a)-1)*c1)
}



gam.result=pr(gam)
pr.parplot(gam.result)
