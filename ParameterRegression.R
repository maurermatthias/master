#source('function.R')
library('stats')
library('evir')

###################################################################
#helper functions

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

#evaluates function with name "fun" at each point in vector xval and parameter par 
myapply<-function(fun,xval,par){
  func=match.fun(fun);
  val=c()
  for(i in 1:length(xval)){
    val=c(val,func(xval[i],par));
  }
  return(val)
}



###################################################################
#checks input for the function pr 
checkInput<-function(x,allowedDistributions){
  #check distribution input
  if(is.null(x[["distr"]]))
    stop("Distribution (distr) of observations not specified!");
  if(typeof(x[["distr"]]) != "character")
    stop("Distribution (distr) needs to be a character-type!")
  if(sum(allowedDistributions==x[["distr"]])==0)
    stop(paste("Distribution (distr) ",x[["distr"]]," not supported!",sep=""))
  #check data (observations/predictors) input
  if(is.null(x[["xval"]]))
    stop("Predictor values (xval) of observations not specified!");
  if(is.null(x[["yval"]]))
    stop("observation values (yval) not specified!");
  if(length(x[["xval"]]) != length(x[["yval"]]))
    stop("Predictor values (xval) and observation values (yval) vectorlength differ!"); 
  if(typeof(x[["xval"]]) != "double" || typeof(x[["yval"]]) != "double")
    stop("Predictor values (xval) and observation values (yval) needs to be a double-type!");
  
}

###################################################################
#Parameter estimators for each stress level

#estimate normal-distribution parameter
#adjusted ML-estimators (adjusted sample variance instead of sample variance)
estimateParameters.norm<-function(x){
  predictor.levels=unique(x[["xval"]])
  p1=c();
  p2=c();
  for(i in 1:length(predictor.levels)){
    data=part(x[["yval"]],x[["xval"]]==predictor.levels[i]);
    if(length(data)<=2)
      stop("For each predictor-value (xval) there need to be at least two observations for the estimation!");
    p1=c(p1,mean(data));
    p2=c(p2,sd(data));
  }
  ret=list();
  ret[["distribution"]]="norm";
  ret[["numberOfParameters"]]=2;
  ret[["predictorlevels"]]=predictor.levels;
  ret[["parameter.name"]]=c("mu","sigma")
  ret[["p1"]]=p1;
  ret[["p2"]]=p2;
  ret[["p"]]=t(array(c(p1,p2),dim=c(length(p1),2)));
  return(ret);
}

#estimate log normal-distribution parameter
#adjusted ML-estimators (adjusted sample variance instead of sample variance)
estimateParameters.logn<-function(x){
  predictor.levels=unique(x[["xval"]])
  p1=c();
  p2=c();
  for(i in 1:length(predictor.levels)){
    data=part(x[["yval"]],x[["xval"]]==predictor.levels[i]);
    if(length(data)<=2)
      stop("For each predictor-value (xval) there need to be at least two observations for the estimation!");
    p1=c(p1,mean(log(data)));
    p2=c(p2,sd(log(data)));
  }
  ret=list();
  ret[["distribution"]]="logn";
  ret[["numberOfParameters"]]=2;
  ret[["predictorlevels"]]=predictor.levels;
  ret[["parameter.name"]]=c("mu","sigma")
  ret[["p1"]]=p1;
  ret[["p2"]]=p2;
  ret[["p"]]=t(array(c(p1,p2),dim=c(length(p1),2)));
  return(ret);
}

#estimate gev distribution parameter
estimateParameters.gev<-function(x){
  predictor.levels=unique(x[["xval"]])
  p1=c();
  p2=c();
  p3=c();
  for(i in 1:length(predictor.levels)){
    data=part(x[["yval"]],x[["xval"]]==predictor.levels[i]);
    if(length(data)<=3)
      stop("For each predictor-value (xval) there need to be at least three observations for the estimation!");
    p.est=gev(data,1)$par.est
    p1=c(p1,p.est[[1]]);
    p2=c(p2,p.est[[2]]);
    p3=c(p3,p.est[[3]]);
  }
  ret=list();
  ret[["distribution"]]="gev";
  ret[["numberOfParameters"]]=3;
  ret[["predictorlevels"]]=predictor.levels;
  ret[["parameter.name"]]=c("xi","sigma","mu")
  ret[["p1"]]=p1;
  ret[["p2"]]=p2;
  ret[["p3"]]=p3;
  ret[["p"]]=t(array(c(p1,p2,p3),dim=c(length(p1),3)));
  return(ret);
}

#estimate gamma distribution parameter
estimateParameters.gamma<-function(x){
  predictor.levels=unique(x[["xval"]])
  p1=c();
  p2=c();
  for(i in 1:length(predictor.levels)){
    data=part(x[["yval"]],x[["xval"]]==predictor.levels[i]);
    if(length(data)<=3)
      stop("For each predictor-value (xval) there need to be at least three observations for the estimation!");
    #moment estimators:
    m1=sum(data)/length(data)
    m2=sum(data^2)/length(data)
    k=m1^2/(m2-m1^2)
    beta=m1/(m2-m1^2)
    theta=1/beta  #f(x,k,theta)=(x^(k-1)e^(-x/theta))/(Gamma(k)theta^k)
    #moment est. end
    p1=c(p1,k);
    p2=c(p2,theta);
  }
  ret=list();
  ret[["distribution"]]="gamma";
  ret[["numberOfParameters"]]=2;
  ret[["predictorlevels"]]=predictor.levels;
  ret[["parameter.name"]]=c("k","theta")
  ret[["p1"]]=p1;
  ret[["p2"]]=p2;
  ret[["p"]]=t(array(c(p1,p2),dim=c(length(p1),2)));
  return(ret);
}

#estimates parameter for the observations for each predictor level
estimateParameters<-function(x){
  if(x[["distr"]]=="norm"){
    return(estimateParameters.norm(x));
  }else if(x[["distr"]]=="logn"){
    return(estimateParameters.logn(x));
  }else if(x[["distr"]]=="gev"){
    return(estimateParameters.gev(x));
  }else if(x[["distr"]]=="gamma"){
    return(estimateParameters.gamma(x));
  }else{
    stop(paste("Distribution (distr) ",x[["distr"]]," - no parameter estimation implemented so far!",sep=""));
  }
}


###################################################################
#estimation of structure parameter with LS methode

#least square error function for optimization
error.ls<-function(parameter,data){
  x=data[[1]];
  if(!is.null(x[["validity.fun"]])){
    val.fun=match.fun(x[["validity.fun"]]);
    if(!val.fun(x[["xval"]],parameter)){
      return(Inf);
    }
  }
  ind.par.est=data[[2]];
  parameter.fun.name = x[["struct.fun"]]
  error.value=0;
  p=ind.par.est[["p"]]
  predictorlevels=ind.par.est[["predictorlevels"]]
  if(length(parameter.fun.name) != ind.par.est[["numberOfParameters"]])
    stop("Number of Parameter do not fit number of Parameter-structure-function supplied!");
  for(i in 1:length(parameter.fun.name)){
    fun=match.fun(parameter.fun.name[i]);
    for(j in 1:length(predictorlevels)){
      structure.value=fun(predictorlevels[j],parameter);
      #if(structure.value < 0)
      #  return(Inf);
      if(!is.null(x[["error.type"]]) && x[["error.type"]]=="rel"){
        #relative error
        error.value=error.value+((p[i,j]-structure.value)/p[i,j])^2;
      }else if(!is.null(x[["error.type"]]) && x[["error.type"]]=="wei"){
        #weigthed error - weigth are mean of values
        error.value=error.value+((p[i,j]-structure.value)/sum(p[i,]))^2;
      }else{
        #absolute error
        error.value=error.value+((p[i,j]-structure.value))^2;
      }
    }   
  }
  return(error.value);
}

#uses grid search to find error.ls-function parameter for finite function value 
grid.search<-function(data){
  max.it=25;
  start.delta=10;
  x=data[[1]];
  ind.par.est=data[[2]];
  stop("Grid search not implemented yet.");
}

find.start.par<-function(data){
  if(is.null(data[[1]][["struct.start.parameter"]])){
    #grid.search(data);
    stop("No initial structure-parameter (struct.start.parameter) for optimization provided!");
  }else{
    return(data[[1]][["struct.start.parameter"]]);
  }
}

estimateStructureParameter<-function(x,ind.par.est){
  if(is.null(x[["control"]])){
    control = list(maxit = 1000);
  }else{
    control = x[["control"]];
  }
  data=list(x,ind.par.est);
  start.par=find.start.par(data);
  o.result=suppressWarnings(optim(start.par,fn=error.ls,data=data, control=control))
  return(o.result);
}


###################################################################
#function gets a named list
#NEEDED:
#distr .... distribution of the observations
#xval ..... predictor values as vector
#yval ..... observation values as vector
#struct.fun...list of function names for modelling parameter structure (one name for each parameter)
#struct.start.parameter....start parameter vector for ls-optimization (see par - optim)
#OPTIONAL:
#error.type ... if rel then relative errors are used for the ls-optimization
#control.......  named-list for optim (see control - optim)
#validity.fun .. name of function returning true, if parameter fulfill requirements, false otherwise
pr<-function(x){
  
  #check the input
  allowedDistributions=c("norm","logn","gev","gamma");
  checkInput(x,allowedDistributions);
  
  #estimate parameters for each predictor level
  ind.par.est = estimateParameters(x);
  
  #define variable for holding information
  val = list();
  val[["input"]]=x;
  val[["ind.par.est"]]=ind.par.est;
  
  
  
  #estimate structure parameter
  struct.par.est = estimateStructureParameter(x,ind.par.est);
  val[["struct.par.opt.result"]]=struct.par.est;
  val[["struct.par.est"]]=struct.par.est$par;
  
  #calculate chi2 test p-value for goodness of fit
  if(is.null(x[["quantiles"]])){
    no.observations=length(val[["input"]][["xval"]])
    no.iqr.areas=round(no.observations/20)
    val[["input"]][["quantiles"]]=1:(no.iqr.areas-1)/no.iqr.areas
  }
  chi2=chi2.test(val);
  val[["chi2.test"]]=chi2

  
  return(val);
}


#####################PLOT
#plot the 2 parameter calculatet - x is the return value from pr 
pr.parplot2<-function(x){
  input=x[["input"]];
  ind.par.est=x[["ind.par.est"]];
  struct.par.est=x[["struct.par.est"]];

  x.points=ind.par.est[["predictorlevels"]]
  y.points1=ind.par.est[["p"]][1,]
  y.points2=ind.par.est[["p"]][2,]
  #xmin=max(min(input[["xval"]])*0.9,struct.par.est[1]+0.0000001);
  xmin=min(input[["xval"]])*0.9;
  xmax=max(input[["xval"]])*1.1
  x.lineTMP1=(0:1000)*(xmax-xmin)/1000+xmin
  x.lineTMP=c(x.lineTMP1,rev(x.lineTMP1));
  
  parameter.fun.name = input[["struct.fun"]];
  #fun1=match.fun(parameter.fun.name[1]);
  #fun2=match.fun(parameter.fun.name[2]);
  y.line1TMP=myapply(parameter.fun.name[1], x.lineTMP, struct.par.est);
  y.line2TMP=myapply(parameter.fun.name[2], x.lineTMP, struct.par.est);
  #y.line1TMP=fun1(x.lineTMP,struct.par.est);
  #y.line2TMP=fun2(x.lineTMP,struct.par.est);
  
  boolv=(!is.na(y.line1TMP)) & (!is.na(y.line2TMP))
  if(input[["distr"]]=="norm"  || input[["distr"]]=="logn")
    boolv= boolv & ((y.line1TMP>=0) & ( y.line2TMP>=0));
  x.line=part(x.lineTMP,boolv);
  y.line1=part(y.line1TMP,boolv);
  y.line2=part(y.line2TMP,boolv);
  
  if(!is.null(ind.par.est[["parameter.name"]])){
    name1=ind.par.est[["parameter.name"]][1]
    name2=ind.par.est[["parameter.name"]][2]
  }else{
    name1="parameter1";
    name2="parameter2";
  }
  
  #x.line=x.lineTMP
  #y.line1=y.line1TMP
  #y.line2=y.line2TMP
  
  #create new plot
  plot.new()
  
  #plot first parameter
  par(fig=c(0,1,0.4,1), new=TRUE)
  plot(x.points,y.points1,xlab="",ylab=name1)
  polygon(x.line,y.line1,border="red")
  legend(x="topright", legend=c(paste("estimate ",name1, sep=""),paste("observation ",name1, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  
  #plot second parameter
  par(fig=c(0,1,0,0.6), new=TRUE)
  plot(x.points,y.points2,xlab="stress",ylab=name2)
  polygon(x.line,y.line2,border="red")
  legend(x="topright", legend=c(paste("estimate ",name2, sep=""),paste("observation ",name2, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
}

#plot the 3 parameter calculatet - x is the return value from pr 
pr.parplot3<-function(x){
  input=x[["input"]];
  ind.par.est=x[["ind.par.est"]];
  struct.par.est=x[["struct.par.est"]];
  
  x.points=ind.par.est[["predictorlevels"]]
  y.points1=ind.par.est[["p"]][1,]
  y.points2=ind.par.est[["p"]][2,]
  y.points3=ind.par.est[["p"]][3,]
  xmin=min(input[["xval"]])*0.9;
  xmax=max(input[["xval"]])*1.1
  x.lineTMP1=(0:1000)*(xmax-xmin)/1000+xmin
  x.lineTMP=c(x.lineTMP1,rev(x.lineTMP1));
  
  parameter.fun.name = input[["struct.fun"]];
  #fun1=match.fun(parameter.fun.name[1]);
  #fun2=match.fun(parameter.fun.name[2]);
  y.line1TMP=myapply(parameter.fun.name[1], x.lineTMP, struct.par.est);
  y.line2TMP=myapply(parameter.fun.name[2], x.lineTMP, struct.par.est);
  y.line3TMP=myapply(parameter.fun.name[3], x.lineTMP, struct.par.est);
  #y.line1TMP=fun1(x.lineTMP,struct.par.est);
  #y.line2TMP=fun2(x.lineTMP,struct.par.est);
  

  boolv = (!is.na(y.line1TMP)) & (!is.na(y.line2TMP)) & (!is.na(y.line3TMP))
  if(input[["distr"]]=="gev")
    boolv = boolv & ((y.line2TMP>=0) & ( y.line3TMP>=0))

  x.line=part(x.lineTMP,boolv);
  y.line1=part(y.line1TMP,boolv);
  y.line2=part(y.line2TMP,boolv);
  y.line3=part(y.line3TMP,boolv);
  
  #x.line=x.lineTMP
  #y.line1=y.line1TMP
  #y.line2=y.line2TMP
  #y.line3=y.line3TMP
  
  if(!is.null(ind.par.est[["parameter.name"]])){
    name1=ind.par.est[["parameter.name"]][1]
    name2=ind.par.est[["parameter.name"]][2]
    name3=ind.par.est[["parameter.name"]][3]
  }else{
    name1="parameter1";
    name2="parameter2";
    name3="parameter3";
  }
  
  #create new plot
  plot.new()
  d=0.225
  
  #plot first parameter
  par(fig=c(0,1,1-2*d,1), new=TRUE)
  plot(x.points,y.points1,xlab="",ylab=name1)
  polygon(x.line,y.line1,border="red")
  legend(x="topright", legend=c(paste("estimate ", name1, sep=""),paste("observation ", name1, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  
  #plot second parameter
  par(fig=c(0,1,0.5-d,0.5+d), new=TRUE)
  plot(x.points,y.points2,xlab="",ylab=name2)
  polygon(x.line,y.line2,border="red")
  legend(x="topright", legend=c(paste("estimate ", name2, sep=""),paste("observation ", name2, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
  #plot third parameter
  par(fig=c(0,1,0,2*d), new=TRUE)
  plot(x.points,y.points3,xlab="stress",ylab=name3)
  polygon(x.line,y.line3,border="red")
  legend(x="topright", legend=c(paste("estimate ", name3, sep=""),paste("observation ", name3, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
  
}

#plot the parameter calculatet - x is the return value from pr
pr.parplot<-function(x){
  if(x[["ind.par.est"]][["numberOfParameters"]]==2){
    pr.parplot2(x);
  }else if(x[["ind.par.est"]][["numberOfParameters"]]==3){
    pr.parplot3(x);
  }else{
    stop("Printing of parameter only implemented for two and three parameters!");
  }
}



##########################Chi2 test

chi2.test<-function(val){
  a=list();
  a[["quantiles"]]=val[["input"]][["quantile"]]
  a[["quantile.values"]]=getQuantileValues(val);
  a[["quantile.numbers"]]=getQuantileNumbers(val,a[["quantile.values"]]);
  a[["chi2.value"]]=getChi2Value(val[["input"]][["quantiles"]],a[["quantile.numbers"]])
  a[["df"]]=length(val[["input"]][["xval"]])
  warning("Degrees of freedom for chi2 test not corrected yet!");
  a[["p.value"]]=pchisq(a[["chi2.value"]],a[["df"]]);
  
  return(a)
}

getChi2Value = function(quantiles,quantile.numbers){
  no.stress=dim(quantile.numbers)[1]-1
  no.quantile.areas=dim(quantile.numbers)[2]-1
  no.observations=sum(quantile.numbers)-sum(quantile.numbers[1])
  observed=c()
  expected=c()
  for(i in 1:no.quantile.areas){
    observed=c(observed,sum(quantile.numbers[[i+1]]))
    if(i==1){
      expected=c(expected,quantiles[1])
    }else if(i==no.quantile.areas){
      expected=c(expected,1-quantiles[i-1])
    }else{
      expected=c(expected,quantiles[i]-quantiles[i-1])
    }
  }
  expected=expected*no.observations
  val=sum((observed-expected)^2/expected)
  return(val)
}

getQuantileNumbers<-function(val,quantileValues){
  stress=unique(val[["input"]][["xval"]]);
  quantiles=val[["input"]][["quantiles"]];  
  A=matrix(rep(0,(length(stress)+1)*(length(quantiles)+1)),nrow=(length(stress)+1),ncol=(length(quantiles)+1));
  for(i in 1:length(stress)){
    data=part(val[["input"]][["yval"]],val[["input"]][["xval"]]==stress[i])
    for(j in 1:(length(quantiles)+1)){
      if(j==1){
        A[i,j]=sum(data<quantileValues[i,j+1])
      }else if(j==(length(quantiles)+1)){
        A[i,j]=sum(data >= quantileValues[i,j])
      }else{
        A[i,j]=sum(data<quantileValues[i,j+1] & data >= quantileValues[i,j])
      }
      
    }
  }
  for(j in 1:(length(quantiles)+1)){
    A[length(stress)+1,j]=sum(A[,j])
  }
  
  df=data.frame(c(stress,0))
  for(i in 1:(length(quantiles)+1)){
    if(i==1){
      df[[paste("q_0-q_",as.character(quantiles[i]), sep= "")]]=A[,i];
    }else if(i==(length(quantiles)+1)){
      df[[paste("q_",as.character(quantiles[i-1]),"-q_1", sep= "")]]=A[,i];
    }else{
      df[[paste("q_",as.character(quantiles[i-1]),"-q_",as.character(quantiles[i]), sep= "")]]=A[,i];
    }
  }
  return(df);

}


getQuantileValues<-function(val){
  stress=unique(val[["input"]][["xval"]]);
  quantiles=val[["input"]][["quantiles"]];
  struct.parameter=val[["struct.par.est"]]
  parameter.function.name=val[["input"]][["struct.fun"]]
  nr.dist.parameter=val[["ind.par.est"]][["numberOfParameters"]]
  
  #calculate distribution parameters for each stress level
  dist.parameter=matrix(rep(0,nr.dist.parameter*length(stress)),nrow=length(stress), ncol=nr.dist.parameter)
  for(j in 1:nr.dist.parameter){
    fun=match.fun(parameter.function.name[j]);
    for(i in 1:length(stress)){
      dist.parameter[i,j]=fun(stress[i],struct.parameter);
    }
  }
  
  #calculate quantile values for each stress level 
  A=matrix(rep(0,length(stress)*length(quantiles)),nrow=length(stress),ncol=length(quantiles));
  for(i in 1:length(stress)){
    parameter=dist.parameter[i,];
    for(j in 1:length(quantiles)){
      if(val[["input"]][["distr"]]=="norm"){
        A[i,j]=qnorm(quantiles[j],parameter[1],parameter[2]);
      }else if(val[["input"]][["distr"]]=="logn"){
        A[i,j]=qlnorm(quantiles[j],parameter[1],parameter[2]);
      }else if(val[["input"]][["distr"]]=="gev"){
        A[i,j]=qgev(quantiles[j],parameter[1],parameter[3],parameter[2]);
      }else if(val[["input"]][["distr"]]=="gamma"){
        A[i,j]=qgamma(quantiles[j],parameter[1],scale=parameter[2]);
      }else{
        stop("Chi2 test not implemented with this distribution!");
      }
    }
  }
  df=data.frame(stress)
  for(i in 1:length(quantiles)){
    df[[paste("q_",as.character(quantiles[i]), sep= "")]]=A[,i];
  }
  return(df);
}


########################################################################################
########################################################################################
