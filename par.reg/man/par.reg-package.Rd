\name{par.reg-package}
\alias{par.reg-package}
\alias{par.reg}
\docType{package}
\title{
Least-Square Parameter Regression
}
\description{
Package for performing a Least-Square fit for an user defined distribution of fatigue data.
}
\details{
\tabular{ll}{
Package: \tab par.reg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-11-21\cr
%License: \tab What license is it under?\cr
}
Within this package the function pr(...) is used to estimate parameters for fatigue-data. Therefore the data needs to contain multiple observations (cyrcles to failure) for each stress-level. First the parameters for the user-defined distributed observations are calculated for each stress-level and then a Least-Square fit is done. Also a chi-squared goodness of fit test is done. The function pr.sim(...) performes multiple time the function pr(...) with a subset of the observations and uses the other subset to perform a chi-squared goodness of fit test.
}
\author{
Matthias Maurer
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%~~ directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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

#perform simulation
v=pr.sim(gev, 0.9, 50)
}
