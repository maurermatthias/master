#source('function.R')
library('stats')
library('evir')

###################################################################
#function gets a named list
#NEEDED:
#distr .... distribution of the observations
#xval ..... predictor values as vector
#yval ..... observation values as vector
#struct.fun...list of function names for modelling parameter structure (one name for each parameter)
#struct.start.parameter....start parameter vector for ls-optimization (see par - optim)
#OPTIONAL:
#type  .........diag/fit; std=fit; if type is diag no fit is performed - only the individual estimates are ploted
#error.type ... if rel then relative errors are used for the ls-optimization
#control.......  named-list for optim (see control - optim)
#validity.fun .. name of function returning true, if parameter fulfill requirements, false otherwise
pr<-function(x){
  
  ###################################################################
  ###################################################################
  #Functions needed for pr
  
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
  
  
  ##########################Chi2 test
  
  chi2.test<-function(val){
    a=list();
    a[["quantiles"]]=val[["input"]][["quantile"]]
    a[["quantile.values"]]=getQuantileValues(val);
    a[["quantile.numbers"]]=getQuantileNumbers(val[["input"]][["xval"]],val[["input"]][["yval"]],val[["input"]][["quantiles"]],a[["quantile.values"]]);
    a[["chi2.value"]]=getChi2Value(val[["input"]][["quantiles"]],a[["quantile.numbers"]])
    a[["df"]]=length(val[["input"]][["xval"]])
    if(is.null(val[["input"]][["sim"]]) || ( !is.null(val[["input"]][["sim"]]) && !val[["input"]][["sim"]]))
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
  
  #x...xval to be tested
  #y...yval to be tested
  #quantiles....vector of qunatile area e.g.: c(0.25,0.5,0.75)
  #qunatileValues....values of the corresponding quantile (per stresslevel)
  getQuantileNumbers<-function(x,y,quantiles,quantileValues){
    stress=unique(x);
    A=matrix(rep(0,(length(stress)+1)*(length(quantiles)+1)),nrow=(length(stress)+1),ncol=(length(quantiles)+1));
    for(i in 1:length(stress)){
      data=part(y,x==stress[i])
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
  #functionality pr start:
  
  #check the input
  allowedDistributions=c("norm","logn","gev","gamma");
  checkInput(x,allowedDistributions);
  
  if(is.null(x[["type"]])){
    x[["type"]]="fit";
  }else{
    if(x[["type"]]!="diag" && x[["type"]] !="fit"){
      x[["type"]]="fit";
    }
  }
  
  #estimate parameters for each predictor level
  ind.par.est = estimateParameters(x);
  
  
  #define variable for holding information
  val = list();
  val[["input"]]=x;
  val[["ind.par.est"]]=ind.par.est;
  
  #plot and return in case of diagnosis plot
  if(x[["type"]]=="diag"){
    pr.parplot(val)
    return(NULL);
  }
  
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


#plot the parameter calculatet - x is the return value from pr
pr.parplot<-function(x){
  ############################################################################
  #pr.parplot needed functions
  
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
  
  #####################PLOT
  #plot the 2 parameter calculatet - x is the return value from pr 
  pr.parplot2<-function(x){
    input=x[["input"]];
    type=input[["type"]]
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
    
    if(type=="fit"){
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
    }
    
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
    if(type=="fit"){
      polygon(x.line,y.line1,border="red")
      legend(x="topright", legend=c(paste("estimate ",name1, sep=""),paste("observation ",name1, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
    }
    
    #plot second parameter
    par(fig=c(0,1,0,0.6), new=TRUE)
    plot(x.points,y.points2,xlab="stress",ylab=name2)
    if(type=="fit"){
      polygon(x.line,y.line2,border="red")
      legend(x="topright", legend=c(paste("estimate ",name2, sep=""),paste("observation ",name2, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
    }
  }
  
  #plot the 3 parameter calculatet - x is the return value from pr 
  pr.parplot3<-function(x){
    input=x[["input"]];
    type=input[["type"]]
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
    
    if(type=="fit"){
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
    }
    
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
    if(type=="fit"){
      polygon(x.line,y.line1,border="red")
      legend(x="topright", legend=c(paste("estimate ", name1, sep=""),paste("observation ", name1, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
    }
    
    #plot second parameter
    par(fig=c(0,1,0.5-d,0.5+d), new=TRUE)
    plot(x.points,y.points2,xlab="",ylab=name2)
    if(type=="fit"){
      polygon(x.line,y.line2,border="red")
      legend(x="topright", legend=c(paste("estimate ", name2, sep=""),paste("observation ", name2, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
    }
    
    #plot third parameter
    par(fig=c(0,1,0,2*d), new=TRUE)
    plot(x.points,y.points3,xlab="stress",ylab=name3)
    if(type=="fit"){
      polygon(x.line,y.line3,border="red")
      legend(x="topright", legend=c(paste("estimate ", name3, sep=""),paste("observation ", name3, sep="")), col = c("red","black"), lty = c("solid",NA), pch=c(NA,"o"))
    }
  }
  
  
  ############################################################################
  #functionality pr.parplot start
  if(x[["ind.par.est"]][["numberOfParameters"]]==2){
    pr.parplot2(x);
  }else if(x[["ind.par.est"]][["numberOfParameters"]]==3){
    pr.parplot3(x);
  }else{
    stop("Printing of parameter only implemented for two and three parameters!");
  }
}




########################################################################################
########################################################################################
#SIMULATION
#input....input list like for pr()
#ratio....p.u. (0<ratio<1) amount of observations per stress level used for fitting
#times...number of times evaluation is done
pr.sim<-function(input, ratio, times){
  ############################################################################
  #needed function for pr.sim
  
  
  ##############
  #chi2-test functions
  
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
  
  #x...xval to be tested
  #y...yval to be tested
  #quantiles....vector of qunatile area e.g.: c(0.25,0.5,0.75)
  #qunatileValues....values of the corresponding quantile (per stresslevel)
  getQuantileNumbers<-function(x,y,quantiles,quantileValues){
    stress=unique(x);
    A=matrix(rep(0,(length(stress)+1)*(length(quantiles)+1)),nrow=(length(stress)+1),ncol=(length(quantiles)+1));
    for(i in 1:length(stress)){
      data=part(y,x==stress[i])
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
  
  
  #supress only special warning
  suppressWarnings2<-function(expr, regex=character()){
    withCallingHandlers(expr, warning=function(w) {
      if (length(regex) == 1 && length(grep(regex, conditionMessage(w)))) {
        invokeRestart("muffleWarning")
      }
    })
  } 
  
  
  
  #creates random bool matrix
  create.rand.boolmatrix<-function(input, ratio, times){
    x=input[["xval"]]
    stress=unique(x)
    no.lines=length(x)
    ratio.new=round(ratio*length(x))
    #warning("Random observation generation not implemented properly just yet");
    #return(t(replicate(length(input[["xval"]]),rnorm(times))<=0.85))
    max.times=1;
    for(gr in 1:length(stress)){
      gr.len=sum(stress[gr]==x)
      gr.ratio=round((1-ratio)*gr.len) 
      max.times=max.times*choose(gr.len,gr.ratio)
    }
    if(max.times<times)
      stop(paste("number of simulation too big for this datasample"));
    vectors=list()
    old.percentage=0;
    for(col in 1:times){
      if(round(100*col/times) != old.percentage){
        old.percentage=round(100*col/times)
        msg=paste("creating simulation matrix - ",as.character(old.percentage),"% done.", sep="")
        status.update(msg);
      }
      #do-while loop
      repeat{
        bool.vec=c()
        for(gr in 1:length(stress)){
          gr.len=sum(stress[gr]==x)
          gr.ratio=round((1-ratio)*gr.len) #number of FALSE
          if(gr.ratio<=0 || gr.ratio>=gr.len)
            stop("ratio too small or too large - need to be between 1 and 0 and be adjusted to the group size!");
          tmp.vec=rep(TRUE,gr.len)
          pos.to.change=sample(1:gr.len,gr.ratio)
          for(i in 1:length(pos.to.change)){
            tmp.vec[pos.to.change]=FALSE
          }
          bool.vec=c(bool.vec,tmp.vec)
        }
        stop.rep=TRUE
        if(length(vectors)!=0){
          for(i in 1:length(vectors)){
            if(sum(bool.vec==vectors[[i]])==length(bool.vec)){
              stop.rep=FALSE
            }
          }
        }
        if(stop.rep){
          break
        }
      }
      vectors[[col]]=bool.vec
    }
    #warning("no check for multiple appereance of the same vector");
    mat=c()
    for(i in 1:times){
      mat=c(mat,vectors[[i]])
    }
    return(matrix(mat,nrow=no.lines,byrow=FALSE))
    
  }
  
  #one simulation step 
  #input ... input list like for pr()
  #bool ... boolvector of size length(input[["xval"]]) - specifies which values are takern for fitting (TRUE)
  #         and which ones are taking for evaluation (FALSE)
  sim<-function(bool,input){
    x.fit=part(input[["xval"]],bool)
    y.fit=part(input[["yval"]],bool)
    x.eval=part(input[["xval"]],bool==FALSE)
    y.eval=part(input[["yval"]],bool==FALSE)
    input.new=input
    input.new[["xval"]]=x.fit
    input.new[["yval"]]=y.fit
    result=suppressWarnings2(pr(input.new),"NaNs produced")
    #evaluation
    quantile.numbers=getQuantileNumbers(x.eval,y.eval,input[["quantiles"]],result$chi2.test[["quantile.values"]]);
    chi2.value=getChi2Value(input[["quantiles"]],quantile.numbers)
    df=length(x.eval)
    #warning("Degrees of freedom for chi2 test not corrected yet!");
    p.value=pchisq(chi2.value,df);
    return(p.value)
  }

  status.update<-function(string){
    cat("                          ", " \r")
    flush.console();
    cat(string, " \r")
    flush.console();
  }
  
  #print gained EDF for p-values
  print.edf<-function(x, points=FALSE){
    x.new=unique(sort(x))
    #points=FALSE
    no.ob=length(x.new)
    y.new=c()
    for(j in 1:no.ob){
      y.tmp=sum(x.new<=x.new[j])/no.ob
      y.new=c(y.new,y.tmp)
    }
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
    #abline(h=0,col="grey")
    #abline(v=0,col="grey")
    if(x.new[1]>0)
      segments(0,0,x.new[1],0,col="red")
    if(x.new[no.ob]<1){
      segments(x.new[no.ob],1,1,1,col="red")
    }else if(points){
      points(1,1,col="red",pch=20)
    }
    for(j in 1:(no.ob-1)){
      segments(x.new[j],y.new[j],x.new[j+1],y.new[j],col="red")
    }
    if(points)
      points(x.new,y.new,col="red",pch=20)
    
  }
  
  
  get.time<-function(sec){
    if(sec<60){
      return(paste(as.character(sec),"[s]",sep=""))
    }else{
      n.min=sec%/%60
      n.sec=as.character(sec-n.min*60)
      if(nchar(n.sec)==1)
        n.sec=paste("0",n.sec,sep="")
      return(paste(as.character(n.min),"[m] ",n.sec,"[s]",sep=""))
    }
  }
  
  
  ###########################################################################
  #pr.sim functionality start
  if(is.null(input[["quantiles"]]))
    stop("Quantiles are needed for this approach!")
  input[["sim"]]=TRUE
  p.val=c()
  display.steps=0.1
  status.update("creating simulation matrix");
  bool.m=create.rand.boolmatrix(input, ratio, times)
  status.update("start calculation");
  start.time = proc.time()[["elapsed"]]
  old.val=-1;
  for(c in 1:times){
    if(TRUE || (c/times)%%0.1==0){
      val=round(100*(c/times))
      if(old.val != val){
        old.val=val;
        time.string=get.time(round(proc.time()[["elapsed"]]-start.time))
        msg = paste("Calculated ",as.character(val),"% in ",time.string ,". Abort with ESC.",sep="");
        status.update(msg);
      }
    }
    bool=bool.m[,c]
    try({
      p.val.tmp = sim(bool,input);
      p.val=c(p.val,p.val.tmp)
      }, silent=TRUE)
  }
  msg=paste("Simulation done in ",get.time(round(proc.time()[["elapsed"]]-start.time))," with ",as.character(round(100*length(p.val)/times)),"% success rate.",sep="");
  status.update(msg);
  par(mfrow=c(2,1))
  print.edf(p.val)
  hist(p.val,main="")#
  par(mfrow=c(1,1))
  return(p.val)
}



