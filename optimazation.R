source('function.R')
shen1=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
shen2=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")

#########################################################################################################################
################################################STD######################################################################
#########################################################################################################################

#data...struct of observations
#output: first scale parameter, then shape parameter
para.est<-function(data,joinedDensity){
  FUN<-match.fun(joinedDensity)
  #estimate parameter of weibull distribution for independent parameter for each stress level
  #as starting parameter for the optimazation
  pa=suppressWarnings(parameter(data,"Weibull"))
  optim(c(pa[,2],pa[,1]),FUN)
}

#negative likelihood density without anny assumptions (only iid for same stress, independent for different stress)
nlogWeibull<-function(params){
  observations=data[,2]
  grouping=data[,1]
  univ=unique(grouping)
  nrOfGroups=length(univ)
  val=0
  for(counter in 1:length(observations)){
    scale=params[which(univ == grouping[counter])]
    shape=params[which(univ == grouping[counter])+nrOfGroups]
    val=val+log(scale)+log(shape)+(shape-1)*log(scale*observations[counter])-(scale*observations[counter])^(shape)
  }
  val=-val
  val
}


##std
data=shen1
p=suppressWarnings(para.est(shen1,nlogWeibull))
p$par
p.ind=suppressWarnings(parameter(data,"Weibull"))
pind=c(p.ind[,2],p.ind[,1])
val=pind==p$par
dif=abs(pind-p$par)

#########################################################################################################################
################################################CONST#SHAPE##############################################################
#########################################################################################################################

#negative likelihood density without anny assumptions (only iid for same stress, independent for different stress)
nlogWeibullConstShape<-function(params){
  observations=data[,2]
  grouping=data[,1]
  univ=unique(grouping)
  nrOfGroups=length(univ)
  val=0
  for(counter in 1:length(observations)){
    scale=params[which(univ == grouping[counter])]
    shape=params[nrOfGroups+1]
    val=val+log(scale)+log(shape)+(shape-1)*log(scale*observations[counter])-(scale*observations[counter])^(shape)
  }
  val=-val
  val
}

#data...struct of observations
#output: first scale parameter then one shape parameter
para.est.constShape<-function(data,joinedDensity){
  FUN<-match.fun(joinedDensity)
  #estimate parameter of weibull distribution for independent parameter for each stress level
  #as starting parameter for the optimazation
  pa=suppressWarnings(parameter(data,"Weibull"))
  optim(c(pa[,2],mean(pa[,1])),FUN)
}

##const shape
data=shen1
p.const=suppressWarnings(para.est.constShape(shen1,nlogWeibullConstShape))
p.const$par
p.ind.const=suppressWarnings(parameter(data,"Weibull"))
pind.const=c(p.ind.const[,2],mean(p.ind.const[,1]))
val.const=pind.const==p.const$par
dif.const=abs(pind.const-p.const$par)

