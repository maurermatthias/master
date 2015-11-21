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
#error.type ... if rel then relative errors are used, if wei the error is weighted with the mean, abs is standard for the ls-optimization
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
  
  #supress only special warning
  suppressWarnings2<-function(expr, regex=character()){
    withCallingHandlers(expr, warning=function(w) {
      if (length(regex) == 1 && length(grep(regex, conditionMessage(w)))) {
        invokeRestart("muffleWarning")
      }
    })
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
  #functionality pr start:
  
  #check the input
  allowedDistributions=c("norm","logn","gev","gamma");
  checkInput(x,allowedDistributions);
  
  if(is.null(x[["sim"]]))
    x[["sim"]]=FALSE
  
  if(is.null(x[["type"]])){
    x[["type"]]="fit";
  }else{
    if(x[["type"]]!="diag" && x[["type"]] !="fit"){
      x[["type"]]="fit";
    }
  }
  
  #estimate parameters for each predictor level
  ind.par.est = suppressWarnings2(estimateParameters(x),"NaNs produced");
  
  
  #define variable for holding information
  val = list();
  val[["input"]]=x;
  val[["ind.par.est"]]=ind.par.est;
  
  #plot and return in case of diagnosis plot
  if(x[["type"]]=="diag"){
    pr.parplot(val)
    return(val);
  }
  
  #estimate structure parameter
  struct.par.est = estimateStructureParameter(x,ind.par.est)
  val[["struct.par.opt.result"]]=struct.par.est;
  val[["struct.par.est"]]=struct.par.est$par;
  

  #calculate chi2 test p-value for goodness of fit
  if(is.null(x[["quantiles"]])){
    no.observations=length(val[["input"]][["xval"]])
    no.iqr.areas=round(no.observations/20)
    val[["input"]][["quantiles"]]=1:(no.iqr.areas-1)/no.iqr.areas
  }
  chi2=chi2.test(val)
  val[["chi2.test"]]=chi2
  
  if(x[["sim"]]!=TRUE)
    pr.parplot(val)
  
  return(val);
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
        msg=paste("creating simulation matrix - ",as.character(old.percentage),"% done. Abort with ESC.", sep="")
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
    result=pr(input.new)
    #evaluation
    quantile.numbers=getQuantileNumbers(x.eval,y.eval,input[["quantiles"]],result$chi2.test[["quantile.values"]]);
    chi2.value=getChi2Value(input[["quantiles"]],quantile.numbers)
    df=length(x.eval)
    #warning("Degrees of freedom for chi2 test not corrected yet!");
    p.value=pchisq(chi2.value,df);
    return(p.value)
  }

  status.update<-function(string){
    cat("                                                                           ", " \r")
    flush.console();
    cat(paste(string,"                                                "), " \r")
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
  
  #supress only special warning
  suppressWarnings2<-function(expr, regex=character()){
    withCallingHandlers(expr, warning=function(w) {
      if (length(regex) == 1 && length(grep(regex, conditionMessage(w)))) {
        invokeRestart("muffleWarning")
      }
    })
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
  input[["type"]]="fit"
  
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
      p.val.tmp = suppressWarnings2(sim(bool,input),"NaNs produced");
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

#load example dataset
load.data<-function(){
  data=matrix(rep(NA,600),nrow=200)
  data[1,1]=294.3; data[1,2]=5300; data[1,3]=0; data[2,1]=294.3; data[2,2]=6200; 
  data[2,3]=0; data[3,1]=294.3; data[3,2]=6500; data[3,3]=0; data[4,1]=294.3; 
  data[4,2]=6600; data[4,3]=0; data[5,1]=294.3; data[5,2]=7000; data[5,3]=0; 
  data[6,1]=294.3; data[6,2]=7500; data[6,3]=0; data[7,1]=294.3; data[7,2]=8000; 
  data[7,3]=0; data[8,1]=294.3; data[8,2]=8400; data[8,3]=0; data[9,1]=294.3; 
  data[9,2]=8700; data[9,3]=0; data[10,1]=294.3; data[10,2]=8800; data[10,3]=0; 
  data[11,1]=294.3; data[11,2]=9000; data[11,3]=0; data[12,1]=294.3; data[12,2]=9200; 
  data[12,3]=0; data[13,1]=294.3; data[13,2]=9200; data[13,3]=0; data[14,1]=294.3; 
  data[14,2]=9400; data[14,3]=0; data[15,1]=294.3; data[15,2]=9500; data[15,3]=0; 
  data[16,1]=294.3; data[16,2]=9500; data[16,3]=0; data[17,1]=294.3; data[17,2]=9800; 
  data[17,3]=0; data[18,1]=294.3; data[18,2]=10000; data[18,3]=0; data[19,1]=294.3; 
  data[19,2]=10500; data[19,3]=0; data[20,1]=294.3; data[20,2]=11800; data[20,3]=0; 
  data[21,1]=220.7; data[21,2]=5100; data[21,3]=0; data[22,1]=220.7; data[22,2]=6100; 
  data[22,3]=0; data[23,1]=220.7; data[23,2]=7000; data[23,3]=0; data[24,1]=220.7; 
  data[24,2]=7700; data[24,3]=0; data[25,1]=220.7; data[25,2]=8600; data[25,3]=0; 
  data[26,1]=220.7; data[26,2]=9000; data[26,3]=0; data[27,1]=220.7; data[27,2]=9100; 
  data[27,3]=0; data[28,1]=220.7; data[28,2]=9300; data[28,3]=0; data[29,1]=220.7;
  data[29,2]=9600; data[29,3]=0; data[30,1]=220.7; data[30,2]=9700; data[30,3]=0; 
  data[31,1]=220.7; data[31,2]=9700; data[31,3]=0; data[32,1]=220.7; data[32,2]=10100; 
  data[32,3]=0; data[33,1]=220.7; data[33,2]=10300; data[33,3]=0; data[34,1]=220.7; 
  data[34,2]=11200; data[34,3]=0; data[35,1]=220.7; data[35,2]=11500; data[35,3]=0; 
  data[36,1]=220.7; data[36,2]=11600; data[36,3]=0; data[37,1]=220.7; data[37,2]=12300; 
  data[37,3]=0; data[38,1]=220.7; data[38,2]=12500; data[38,3]=0; data[39,1]=220.7; 
  data[39,2]=13400; data[39,3]=0; data[40,1]=220.7; data[40,2]=15900; data[40,3]=0; 
  data[41,1]=176.6; data[41,2]=6200; data[41,3]=0; data[42,1]=176.6; data[42,2]=9400; 
  data[42,3]=0; data[43,1]=176.6; data[43,2]=10000; data[43,3]=0; data[44,1]=176.6; 
  data[44,2]=10000; data[44,3]=0; data[45,1]=176.6; data[45,2]=10200; data[45,3]=0; 
  data[46,1]=176.6; data[46,2]=10800; data[46,3]=0; data[47,1]=176.6; data[47,2]=11300; 
  data[47,3]=0; data[48,1]=176.6; data[48,2]=12600; data[48,3]=0; data[49,1]=176.6; 
  data[49,2]=12800; data[49,3]=0; data[50,1]=176.6; data[50,2]=13900; data[50,3]=0; 
  data[51,1]=176.6; data[51,2]=14000; data[51,3]=0; data[52,1]=176.6; data[52,2]=14200; 
  data[52,3]=0; data[53,1]=176.6; data[53,2]=14300; data[53,3]=0; data[54,1]=176.6; 
  data[54,2]=14700; data[54,3]=0; data[55,1]=176.6; data[55,2]=15100; data[55,3]=0; 
  data[56,1]=176.6; data[56,2]=15200; data[56,3]=0; data[57,1]=176.6; data[57,2]=16600; 
  data[57,3]=0; data[58,1]=176.6; data[58,2]=16900; data[58,3]=0; data[59,1]=176.6; 
  data[59,2]=17000; data[59,3]=0; data[60,1]=176.6; data[60,2]=18200; data[60,3]=0; 
  data[61,1]=134.9; data[61,2]=9100; data[61,3]=0; data[62,1]=134.9; data[62,2]=9300; 
  data[62,3]=0; data[63,1]=134.9; data[63,2]=9400; data[63,3]=0; data[64,1]=134.9; 
  data[64,2]=9700; data[64,3]=0; data[65,1]=134.9; data[65,2]=14500; data[65,3]=0; 
  data[66,1]=134.9; data[66,2]=15900; data[66,3]=0; data[67,1]=134.9; data[67,2]=16000; 
  data[67,3]=0; data[68,1]=134.9; data[68,2]=16200; data[68,3]=0; data[69,1]=134.9; 
  data[69,2]=17900; data[69,3]=0; data[70,1]=134.9; data[70,2]=18500; data[70,3]=0; 
  data[71,1]=134.9; data[71,2]=19800; data[71,3]=0; data[72,1]=134.9; data[72,2]=20800; 
  data[72,3]=0; data[73,1]=134.9; data[73,2]=21000; data[73,3]=0; data[74,1]=134.9; 
  data[74,2]=21800; data[74,3]=0; data[75,1]=134.9; data[75,2]=22100; data[75,3]=0; 
  data[76,1]=134.9; data[76,2]=22400; data[76,3]=0; data[77,1]=134.9; data[77,2]=22400; 
  data[77,3]=0; data[78,1]=134.9; data[78,2]=25700; data[78,3]=0; data[79,1]=134.9; 
  data[79,2]=25800; data[79,3]=0; data[80,1]=134.9; data[80,2]=27800; data[80,3]=0; 
  data[81,1]=105.4; data[81,2]=12800; data[81,3]=0; data[82,1]=105.4; data[82,2]=15600;
  data[82,3]=0; data[83,1]=105.4; data[83,2]=17400; data[83,3]=0; data[84,1]=105.4; 
  data[84,2]=19000; data[84,3]=0; data[85,1]=105.4; data[85,2]=19000; data[85,3]=0; 
  data[86,1]=105.4; data[86,2]=19700; data[86,3]=0; data[87,1]=105.4; data[87,2]=20000; 
  data[87,3]=0; data[88,1]=105.4; data[88,2]=21000; data[88,3]=0; data[89,1]=105.4;
  data[89,2]=21300; data[89,3]=0; data[90,1]=105.4; data[90,2]=24400; data[90,3]=0; 
  data[91,1]=105.4; data[91,2]=25100; data[91,3]=0; data[92,1]=105.4; data[92,2]=25400;
  data[92,3]=0; data[93,1]=105.4; data[93,2]=26700; data[93,3]=0; data[94,1]=105.4; 
  data[94,2]=26800; data[94,3]=0; data[95,1]=105.4; data[95,2]=26900; data[95,3]=0; 
  data[96,1]=105.4; data[96,2]=28300; data[96,3]=0; data[97,1]=105.4; data[97,2]=28500; 
  data[97,3]=0; data[98,1]=105.4; data[98,2]=29500; data[98,3]=0; data[99,1]=105.4; 
  data[99,2]=30900; data[99,3]=0; data[100,1]=105.4; data[100,2]=38200; data[100,3]=0;
  data[101,1]=83.4; data[101,2]=18200; data[101,3]=0; data[102,1]=83.4; data[102,2]=25000; 
  data[102,3]=0; data[103,1]=83.4; data[103,2]=25700; data[103,3]=0; data[104,1]=83.4;
  data[104,2]=28600; data[104,3]=0; data[105,1]=83.4; data[105,2]=29000; data[105,3]=0; 
  data[106,1]=83.4; data[106,2]=33700; data[106,3]=0; data[107,1]=83.4; data[107,2]=35000; 
  data[107,3]=0; data[108,1]=83.4; data[108,2]=36400; data[108,3]=0; data[109,1]=83.4; 
  data[109,2]=39900; data[109,3]=0; data[110,1]=83.4; data[110,2]=40000; data[110,3]=0; 
  data[111,1]=83.4; data[111,2]=40700; data[111,3]=0; data[112,1]=83.4; data[112,2]=44000; 
  data[112,3]=0; data[113,1]=83.4; data[113,2]=45100; data[113,3]=0; data[114,1]=83.4; 
  data[114,2]=46000; data[114,3]=0; data[115,1]=83.4; data[115,2]=46100; data[115,3]=0; 
  data[116,1]=83.4; data[116,2]=46800; data[116,3]=0; data[117,1]=83.4; data[117,2]=48700; 
  data[117,3]=0; data[118,1]=83.4; data[118,2]=50000; data[118,3]=0; data[119,1]=83.4; 
  data[119,2]=54300; data[119,3]=0; data[120,1]=83.4; data[120,2]=55600; data[120,3]=0; 
  data[121,1]=73.6; data[121,2]=12000; data[121,3]=0; data[122,1]=73.6; data[122,2]=40000; 
  data[122,3]=0; data[123,1]=73.6; data[123,2]=45000; data[123,3]=0; data[124,1]=73.6; 
  data[124,2]=48000; data[124,3]=0; data[125,1]=73.6; data[125,2]=62000; data[125,3]=0; 
  data[126,1]=73.6; data[126,2]=65000; data[126,3]=0; data[127,1]=73.6; data[127,2]=65000;
  data[127,3]=0; data[128,1]=73.6; data[128,2]=67000; data[128,3]=0; data[129,1]=73.6; 
  data[129,2]=70000; data[129,3]=0; data[130,1]=73.6; data[130,2]=80000; data[130,3]=0;
  data[131,1]=73.6; data[131,2]=81000; data[131,3]=0; data[132,1]=73.6; data[132,2]=83000;
  data[132,3]=0; data[133,1]=73.6; data[133,2]=88000; data[133,3]=0; data[134,1]=73.6; 
  data[134,2]=91000; data[134,3]=0; data[135,1]=73.6; data[135,2]=92000; data[135,3]=0; 
  data[136,1]=73.6; data[136,2]=94000; data[136,3]=0; data[137,1]=73.6; data[137,2]=95000; 
  data[137,3]=0; data[138,1]=73.6; data[138,2]=104000; data[138,3]=0; data[139,1]=73.6;
  data[139,2]=108000; data[139,3]=0; data[140,1]=73.6; data[140,2]=112000; data[140,3]=0;
  data[141,1]=56.4; data[141,2]=114000; data[141,3]=0; data[142,1]=56.4; data[142,2]=130000; 
  data[142,3]=0; data[143,1]=56.4; data[143,2]=157000; data[143,3]=0; data[144,1]=56.4;
  data[144,2]=157000; data[144,3]=0; data[145,1]=56.4; data[145,2]=159000; data[145,3]=0;
  data[146,1]=56.4; data[146,2]=170000; data[146,3]=0; data[147,1]=56.4; data[147,2]=180000;
  data[147,3]=0; data[148,1]=56.4; data[148,2]=201000; data[148,3]=0; data[149,1]=56.4; 
  data[149,2]=205000; data[149,3]=0; data[150,1]=56.4; data[150,2]=210000; data[150,3]=0;
  data[151,1]=56.4; data[151,2]=230000; data[151,3]=0; data[152,1]=56.4; data[152,2]=244000; 
  data[152,3]=0; data[153,1]=56.4; data[153,2]=250000; data[153,3]=0; data[154,1]=56.4; 
  data[154,2]=251000; data[154,3]=0; data[155,1]=56.4; data[155,2]=257000; data[155,3]=0; 
  data[156,1]=56.4; data[156,2]=266000; data[156,3]=0; data[157,1]=56.4; data[157,2]=273000; 
  data[157,3]=0; data[158,1]=56.4; data[158,2]=287000; data[158,3]=0; data[159,1]=56.4; 
  data[159,2]=296000; data[159,3]=0; data[160,1]=56.4; data[160,2]=309000; data[160,3]=0;
  data[161,1]=54; data[161,2]=285000; data[161,3]=0; data[162,1]=54; data[162,2]=308000; 
  data[162,3]=0; data[163,1]=54; data[163,2]=336000; data[163,3]=0; data[164,1]=54; 
  data[164,2]=377000; data[164,3]=0; data[165,1]=54; data[165,2]=380000; data[165,3]=0; 
  data[166,1]=54; data[166,2]=396000; data[166,3]=0; data[167,1]=54; data[167,2]=427000; 
  data[167,3]=0; data[168,1]=54; data[168,2]=497000; data[168,3]=0; data[169,1]=54; 
  data[169,2]=510000; data[169,3]=0; data[170,1]=54; data[170,2]=551000; data[170,3]=0; 
  data[171,1]=54; data[171,2]=560000; data[171,3]=0; data[172,1]=54; data[172,2]=595000; 
  data[172,3]=0; data[173,1]=54; data[173,2]=617000; data[173,3]=0; data[174,1]=54; 
  data[174,2]=660000; data[174,3]=0; data[175,1]=54; data[175,2]=668000; data[175,3]=0;
  data[176,1]=54; data[176,2]=685000; data[176,3]=0; data[177,1]=54; data[177,2]=714000;
  data[177,3]=0; data[178,1]=54; data[178,2]=733000; data[178,3]=0; data[179,1]=54; 
  data[179,2]=849000; data[179,3]=0; data[180,1]=54; data[180,2]=895000; data[180,3]=0; 
  data[181,1]=51.5; data[181,2]=820000; data[181,3]=0; data[182,1]=51.5; data[182,2]=839000; 
  data[182,3]=0; data[183,1]=51.5; data[183,2]=938000; data[183,3]=0; data[184,1]=51.5; 
  data[184,2]=1024000; data[184,3]=0; data[185,1]=51.5; data[185,2]=1040000; data[185,3]=0; 
  data[186,1]=51.5; data[186,2]=1048000; data[186,3]=0; data[187,1]=51.5; data[187,2]=1100000;
  data[187,3]=0; data[188,1]=51.5; data[188,2]=1103000; data[188,3]=0; data[189,1]=51.5; 
  data[189,2]=1136000; data[189,3]=0; data[190,1]=51.5; data[190,2]=1145000; data[190,3]=0;
  data[191,1]=51.5; data[191,2]=1147000; data[191,3]=0; data[192,1]=51.5; data[192,2]=1150000; 
  data[192,3]=0; data[193,1]=51.5; data[193,2]=1151000; data[193,3]=0; data[194,1]=51.5; 
  data[194,2]=1163000; data[194,3]=0; data[195,1]=51.5; data[195,2]=1200000; data[195,3]=0;
  data[196,1]=51.5; data[196,2]=1210000; data[196,3]=0; data[197,1]=51.5; data[197,2]=1319000; 
  data[197,3]=0; data[198,1]=51.5; data[198,2]=1320000; data[198,3]=0; data[199,1]=51.5;
  data[199,2]=1321000; data[199,3]=0; data[200,1]=51.5; data[200,2]=1630000; data[200,3]=0;
  return(data)
}



