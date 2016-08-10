\name{pr.sim}
\alias{pr.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Perform Parameter-Regression simulation
}
\description{ 
Takes repeatedly a part of the data as input for the function pr() and the other part to perform a chi-squared goodness of fit test. It can also be used to generate sample data with given structural parameters and estimate those parameters afterwards. 
}
\usage{
pr.sim(input, sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ The following arguments are needed:
  \item{input}{ Input List like used with the function pr().
}
  \item{sim}{ A named list describing simulation behaviour in case data sample is generated.
}
\item{sim$type}{Specifying the basic simulation behaviour: in case it has value "observe" a part of the data is taken for a fit and the other part of the data is taken to evaluate the fit for several times. In case it has the value "generate" observations following the assumed model are generated multiple times.}
  \item{sim$ratio}{ Number for which 0<ratio<1 holds. Describes how many observations are used for the fitting in case "observe"-simulation (see sim$type) is done (per stress-level: rounding is done).
}
  \item{sim$times}{ Describes how often the fitting is done.
}
\item{sim$xval}{Only needed in case "generate"-simulation (see sim$type) is done. A vector describing for which predictor-values the values should be generated.}
\item{sim$n}{Only needed in case "generate"-simulation (see sim$type) is done. A vector describing how many values should be generated at each predictor-value.}
\item{sim$par}{Only needed in case "generate"-simulation (see sim$type) is done. A vector describing the structur-parameter vector used for the assumption (Parameters used for the function spezified in the input-list describing the parameter behaviour).}
\item{sim$plot (optional)}{Only needed in case "generate"-simulation (see sim$type) is done. Default-value is TRUE, spezifies if a plot of the mean an sd of the fitted data over iterations should be done at the end of the simulation. }
}
%\details{ details test
%%  ~~ If necessary, more details than the description above ~~
%}
\value{ There are two possible outputs, depending on the passed sim[["type"]]-value:
\item{sim[["type"]]=="observe"}{A list containing the p-values for all performed chi-squared goodness of fit test is returned.}
\item{sim[["type"]]=="generate"}{A matrix containing all structur-parameter estimations. Each line describes the parametervector for one fit.}
}
%\references{ references test
%% ~put references to the literature/web site here ~
%}
\author{ 
Matthias Maurer
}
\note{ 
If a fit does not succeed, the p-value is skipped in case sim[["type"]]=="observe" - the success-rate prompted at the end of the function-execution shows how many fits succeeded.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
