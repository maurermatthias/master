source('function.R')
library(R.matlab)
#annealed aluminium wire (AAW)
shen1=readMat("C:/Users/moja/Dropbox/12 sem/Masterarbeit/MATLAB/shen1.mat")$dat
#2024-T4 aluminum alloy fatigue data
shen2=readMat("C:/Users/moja/Dropbox/12 sem/Masterarbeit/MATLAB/shen2.mat")$dat


data=shen1
#Scatterplot 
plot(data[,2],data[,1], main="Scatterplot",  xlab="N ", ylab="Belastung [MPa]",)
#logN - Scatterplot
plot(log(data[,2]),data[,1], main="Scatterplot",  xlab="log(N) ", ylab="Belastung [MPa]",)

par.normal=parameter(data,"normal")
par.lognormal=parameter(data,"lognormal")
par.weibull=parameter(data,"weibull")
par.gev = gev.est(data)
