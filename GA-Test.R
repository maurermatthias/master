source('GA-EstimateMomentStructure.R')
source('GA-DistributionParameterEstimation.R')

#load R-Data:
shen1=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
shen2=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")

ls.est.gamma = moment.par.est(shen1)
plot.est.gamma(ls.est.gamma)
pb=distribution.parameter.est(ls.est.gamma,unique(shen1[,1]))


ml.est.gamma=likelihood.par.est.gamma(shen1)
plot.est.gamma(ml.est.gamma)

error.function.ml.gamma(shen1,ml.est.gamma$estimates)
error.function.ml.gamma(shen1,ls.est.gamma$estimates)
