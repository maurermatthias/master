#install.packages("R.matlab")
#reading .mat file in:
library(R.matlab)
shen1=readMat("C:/Users/moja/Dropbox/12 sem/Masterarbeit/MATLAB/shen1.mat")
shen2=readMat("C:/Users/moja/Dropbox/12 sem/Masterarbeit/MATLAB/shen2.mat")
#save R-Data:
saveRDS(shen1$dat,"C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
saveRDS(shen2$dat,"C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")
#load R-Data:
shen1=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen1.rds")
shen2=readRDS("C:/Users/moja/Dropbox/12 sem/Masterarbeit/R/shen2.rds")