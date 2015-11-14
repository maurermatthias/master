source('function.R')


distribution.parameter.est=function(pa,xval){
  gamma.parameter=list()
  a=pa$opt$par[1]
  b=pa$opt$par[2]
  c1=pa$opt$par[3]
  c2=pa$opt$par[4]
  gamma.parameter[["theta"]]=(c2/c1)*((b-a)/(xval-a))
  gamma.parameter[["k"]]=(1-(xval-a)/(b-a))*(c1^2/c2)
  
  stress=unique(pa$data[,1])
  k.stress=(1-(stress-a)/(b-a))*(c1^2/c2)
  theta.stress=(c2/c1)*((b-a)/(stress-a))
  pval=c()
  for(counter in 1 :(length(stress))){
    dat=part(pa$data[,2],pa$data[,1]==stress[counter])
    test=suppressWarnings(ks.test(dat,'pgamma',k.stress[counter],1/theta.stress[counter]))
    #test=suppressWarnings(ks.test(dat,'pgamma',29.44,1/290.28))
    pval[length(pval)+1]=test$p.value
  }
  gamma.parameter[["p.val.stresslevel"]]=pval
  
  gamma.parameter
}
