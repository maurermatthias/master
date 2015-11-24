\name{pr.sim}
\alias{pr.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Perform Parameter-Regression Quality Check by simulation
}
\description{ 
Takes repeatedly a part of the data as input for the function pr() and the other part to perform a chi-squared goodness of fit test. It can also be used to generate sample data with given structural parameters and estimate those parameters afterwards. 
}
\usage{
pr.sim(input, ratio, times, sim=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ arguments test
  \item{input}{ Input List like used with the function pr().
}
  \item{ratio}{ Number for which 0<ratio<1 holds. Describes how many observations are used for the fitting (per stress-level: rounding is done).
}
  \item{times}{ Describes how often fitting and evaluation is done.
}
  \item{sim}{ A named list describing simulation behaviour in case data sample is generated. It must contain the variables "dist" specifiing the distribution for generating the sample, a vector "xval" describing the predictor values of the sample, a vector n specifiing the number of observations at each predictor-value and par describing the structur-parameter for generating the data sample. If this parameter is not used a quality check with the passed observations in the "input"-list is done. 
}
}
%\details{ details test
%%  ~~ If necessary, more details than the description above ~~
%}
\value{ 
A list containing the p-values for all performed chi-squared goodness of fit test is returned, in case a quality check for an observation is done (sim is not passed to the function).
}
%\references{ references test
%% ~put references to the literature/web site here ~
%}
\author{ 
Matthias Maurer
}
\note{ 
If a fit does not succeed, the p-value is skipped. The success-rate prompted at the end of the function-execution shows how many fits succeeded.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ seealso test
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{ 
#define input-list
gev=list();
gev[["distr"]]="gev"
gev[["error.type"]]="rel"
data=load.data()
gev[["xval"]]=data[,1]
gev[["yval"]]=data[,2]
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

gev[["control"]]=list(maxit=1000)

#perform simulation
v=pr.sim(gev, 0.9, 50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
