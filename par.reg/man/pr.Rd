\name{pr}
\alias{pr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Parameter Regression for fatigue data
}
\description{ 
Does either a LS- or ML-fit for user-spezified parameter-function.
}
\usage{
pr(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 
The function takes one argument "x", a named list with the following members:
  \item{x$distr}{Distribution of the observations, element out of \{norm, lnorm, gamma, gev\} for normal-,lognormal-, gamma- and extreme-value-distribution.}
  \item{x$xval}{Predictor values (stress-level) as vector.}
  \item{x$yval}{Observation values (circles to failure) as vector.}
  \item{x$ML}{If TRUE a ML-fit instead of a LS-fit is done.}
  \item{x$struct.fun}{List of names for parameter functions (one name for each parameter).}
  \item{x$struct.start.parameter}{Start parameter vector for Least-Square optimization (see par - optim).}
  
  
  \item{x$type (optional)}{Has values in \{fit (default),diag\}. If type has value diag no fit is performed - only the individual estimates are ploted to help deciding on the parameter regression functions. }
  \item{x$error.type (optional)}{Has values in \{abs (default),rel,wei\} and specifies the error function for the Least-Square error. If "abs" is selected the error can be calculated as the sum over all parameters of sum((observed.Parameter-calculated.Parameter)^2). If "rel" is selected the error can be calculated as the sum over all parameters of sum(((observed.Parameter-calculated.Parameter)/observed.Parameter)^2). If "wei" is selected the error can be calculated as the sum over all parameters of sum(((observed.Parameter-calculated.Parameter)/mean(observed.Parameter))^2) }
  \item{x$control (optional)}{Named-list for optim (see argument control - optim)}
  \item{x$validity.fun (optional)}{
  Name of function returning TRUE, if parameters fulfill requirements, FALSE otherwise.   
  }
  \item{x$quantiles (optional)}{Vector with values v: 0<v<1, v[i]<v[i+1]. describing the quantile areas for the chi-square goodness of fit test.}
}



%\details{ details test
%%  ~~ If necessary, more details than the description above ~~
%}
\value{ 
This function returns a list with the following fields:
\item{input }{List with input fields to pr() supplemented with default values.}
\item{ind.par.est}{List containing the individual parameter estimations for each stress-level - if LS-fit is done.}
\item{struct.par.opt.result }{Complete Optimization-result gained from optim.}
\item{struct.par.est}{Parameter values gained from the Least-Square estimation of the structural parameter used in the parameter functions.}
\item{chi2.test}{List containing values gained from the chi-squared test performed based on the Least-Square fit.}
}
%\references{ references test
%% ~put references to the literature/web site here ~
%}
\author{ 
Matthias Maurer
}
%\note{ note test
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ seealso test
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{ 
#defining the input list
gev=list();
gev[["distr"]]="gev"
data=load.data()
gev[["xval"]]=data[,1]
gev[["yval"]]=data[,2]

#print individual parameter estimations - decide on the parameter functions
gev[["type"]]="diag"
gev.result=pr(gev)

#defining additional input-list fields for fitting
gev[["type"]]="fit"
gev[["error.type"]]="rel"
gev[["validity.fun"]]="val.gev"
gev[["struct.fun"]]=c("gev1","gev2","gev3")
gev[["struct.start.parameter"]]=c(-0.31,49.67,20809,18.8,69.41)
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

#perform fit
gev.result=pr(gev)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
