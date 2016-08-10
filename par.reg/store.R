library(psych)


info.load<-function(){
  return(readRDS("info"))
}

info.clear<-function(){
  a=list()
  saveRDS(a,"info")
}

info.store<-function(n,distr,type,data,p.vec=NULL){
  
  getValuesSim<-function(val){
    a=list()
    a[["mean"]]=list()
    a[["sd"]]=list()
    a[["mean"]][["set"]]=c()
    a[["sd"]][["set"]]=c()
    for(j in (1+round(length(val)*0.99)):length(val)){
      a[["mean"]][["set"]]=c(a[["mean"]][["set"]],mean(val[1:j]))
      a[["sd"]][["set"]]=c(a[["sd"]][["set"]],sd(val[1:j]))
    }
    a[["mean"]][["intercept"]]  = lm(a[["mean"]][["set"]]~1)[["coefficients"]][[1]]
    a[["sd"]][["intercept"]]  = lm(a[["sd"]][["set"]]~1)[["coefficients"]][[1]]
    a[["mean"]][["slope"]]  = lm(a[["mean"]][["set"]]~I(1:length(a[["mean"]][["set"]])))[["coefficients"]][[2]]
    a[["sd"]][["slope"]]  = lm(a[["sd"]][["set"]]~I(1:length(a[["sd"]][["set"]])))[["coefficients"]][[2]]
    a[["mean"]][["range"]]  = max(a[["mean"]][["set"]])-min(a[["mean"]][["set"]])
    a[["sd"]][["range"]]  =max(a[["sd"]][["set"]])-min(a[["sd"]][["set"]])
    #a[["bigger"]]  = sum()
    return(a);
  }
  
  info=readRDS("info")
  if(is.null(info[[as.character(n)]])){
    info[[as.character(n)]]=list()
    info[[as.character(n)]][[distr]]=list()
    info[[as.character(n)]][[distr]][[type]]=list()
  }else if(is.null(info[[as.character(n)]][[distr]])){
    info[[as.character(n)]][[distr]]=list()
    info[[as.character(n)]][[distr]][[type]]=list()
  }else if(is.null(info[[as.character(n)]][[distr]][[type]])){
    info[[as.character(n)]][[distr]][[type]]=list()
  }
  no.par=length(data[1,])
  for(i in 1:no.par){
    info[[as.character(n)]][[distr]][[type]][[paste("p",as.character(i),sep="")]]=getValuesSim(data[,i])
  }
  
  saveRDS(info,"info")
}


info.fill<-function(){
  info.clear()
  d=info.load()
  d[["parameter"]]=list()
  d[["parameter"]][["Gev"]]=c(-0.22,49.64,22760,17.25,63.4)
  d[["parameter"]][["Norm"]]=c(49.6,346,7261,1401)
  d[["parameter"]][["Logn"]]=c(47.3,323,0.092,0.293)
  saveRDS(d,"info")
  
  folder=c("n5times10000","n10times10000","n20times10000","n30times10000","n40times10000","n50times10000","n100times10000")
  folderNr=c(5,10,20,30,40,50,100)
  dist=c("Gev","Norm","Logn")
  type=c("Abs","Ml","Rel","Wei");
  for(j in 1:length(folder)){
    for(i in 1:length(dist)){
      for(k in 1:length(type)){
        dat=readRDS(paste(folder[j],"/",dist[i],type[k],sep=""));
        info.store(folderNr[j],dist[i],type[k],dat,p.vec=d[["parameter"]][[dist]])
      }
    }
  }
  


}

info.table.n<-function(n){
  data=info.load()
  if(is.null(data[[as.character(n)]])){
    return(NULL)
  }
}

info.table.dist<-function(dist,pnr){
  data=info.load()
  original=data[["parameter"]][[dist]][pnr]
  df=data.frame(type=names(a[[1]][[1]]))
  for(i in 1:(length(data)-1)){
    m.intercept=c()
    s.intersept=c()
    m.slope=c()
    s.slope=c()
    m.range=c()
    s.range=c()
    for(j in 1:length(data[[i]][[dist]])){
      d=data[[i]][[dist]][[j]][[paste("p",as.character(pnr),sep="")]]
      m.intercept=c(m.intercept,d[["mean"]][["intercept"]])
      s.intersept=c(s.intersept,d[["sd"]][["intercept"]])
      m.slope=c(m.slope,d[["mean"]][["slope"]])
      s.slope=c(s.slope,d[["sd"]][["slope"]])
      m.range=c(m.range,d[["mean"]][["range"]])
      s.range=c(s.range,d[["sd"]][["range"]])
    }
    df[[paste("n=",names(data)[i]," m.i",sep="")]]=m.intercept
  }
  return(df)
}

info.get<-function(n,dist,type,pnr,quantity,value){
  return(info.load()[[as.character(n)]][[dist]][[type]][[paste("p",as.character(pnr),sep="")]][[quantity]][[value]])
}

info.table<-function(dist,type,pnr){
  data=info.load()
  original=data[["parameter"]][[dist]][pnr]
  name=names(data)
  a=list()
  
  for(i in 1:(length(data)-1)){

    d=data[[i+1]][[dist]][[type]][[paste("p",as.character(pnr),sep="")]]
    if(i!=1){
      a[[1]]=c(a[[1]],d[["mean"]][["intercept"]])
      a[[2]]=c(a[[2]],d[["sd"]][["intercept"]])
      a[[3]]=c(a[[3]],d[["mean"]][["slope"]])
      a[[4]]=c(a[[4]],d[["sd"]][["slope"]])
      a[[5]]=c(a[[5]],d[["mean"]][["range"]])
      a[[6]]=c(a[[6]],d[["sd"]][["range"]])
    }else{
      a[[1]]=c(d[["mean"]][["intercept"]])
      a[[2]]=c(d[["sd"]][["intercept"]])
      a[[3]]=c(d[["mean"]][["slope"]])
      a[[4]]=c(d[["sd"]][["slope"]])
      a[[5]]=c(d[["mean"]][["range"]])
      a[[6]]=c(d[["sd"]][["range"]])
    }
    
  }

  
  df=data.frame(m.i=a[[1]])
  df[[" s.i"]]=a[[2]]
  df[[" m.s"]]=a[[3]]
  df[[" s.s"]]=a[[4]]
  df[[" m.r"]]=a[[5]]
  df[[" s.r"]]=a[[6]]

  
  mat=t(data.matrix(df))
  
  df2=data.frame(nr=1:length(mat[,1]))
  
  for(i in 1:(length(data)-1)){
    df2[[paste("n=",as.character(name[i+1]),sep="")]]=mat[,i]
  }
  print("Zeilen(jeweils abwechselnd  mean,sd): intercept,slope,range")
  return(df2)
}


info.table.change<-function(dist,type){
  data=info.load()
  anz.par=length(data[[1]][[dist]][[1]])
  real.par=data[["parameter"]][[dist]]
  names=names(data)
  
  a=data.frame(real=real.par)
  #n
  for(i in 1:(length(data)-1)){
    par=c()
    #pnr
    for(j in 1:length(data[[2]][[dist]][[type]])){
      d=data[[i+1]][[dist]][[type]][[paste("p",as.character(j),sep="")]][["mean"]][["intercept"]]
      par=c(par,d)
    }
    a[[paste("r=",names[i+1],sep="")]]=par
  }
  
  b=data.frame(real=real.par)
  for(i in 2:length(a)){
    b[[names(a)[i]]]=abs((a[[i]]-a[[1]])/a[[1]])*100
  }
  
  
  #a...absolute werte intercept
  #b...relative werte intercept
  return(b)
}


h=info.table.change("Norm","Abs")




info.table.change.sd<-function(dist,type){
  data=info.load()
  anz.par=length(data[[1]][[dist]][[1]])
  real.par=data[["parameter"]][[dist]]
  names=names(data)
  
  a=data.frame(real=real.par)
  #n
  for(i in 1:(length(data)-1)){
    par=c()
    #pnr
    for(j in 1:length(data[[2]][[dist]][[type]])){
      d=data[[i+1]][[dist]][[type]][[paste("p",as.character(j),sep="")]][["sd"]][["intercept"]]
      par=c(par,d)
    }
    a[[paste("r=",names[i+1],sep="")]]=par
  }
  
  b=data.frame(real=real.par)
  for(i in 2:length(a)){
    b[[names(a)[i]]]=abs((a[[i]])/a[[1]])*100
  }
  
  
  #a...absolute werte sd
  #b...relative werte sd
  return(a)
}


info.table.change.sddiff<-function(dist,type){
  data=info.load()
  anz.par=length(data[[1]][[dist]][[1]])
  real.par=data[["parameter"]][[dist]]
  names=names(data)
  
  a=data.frame(real=real.par)
  s=data.frame(real=real.par)
  #n
  for(i in 1:(length(data)-1)){
    par=c()
    par2=c()
    #pnr
    for(j in 1:length(data[[2]][[dist]][[type]])){
      d=data[[i+1]][[dist]][[type]][[paste("p",as.character(j),sep="")]][["mean"]][["intercept"]]
      par=c(par,d)
      d2=data[[i+1]][[dist]][[type]][[paste("p",as.character(j),sep="")]][["sd"]][["intercept"]]
      par2=c(par2,d2)
    }
    a[[paste("r=",names[i+1],sep="")]]=par
    s[[paste("r=",names[i+1],sep="")]]=par2
  }
  
  b=data.frame(real=real.par)
  for(i in 2:length(a)){
    b[[names(a)[i]]]=abs((a[[i]]-a[[1]]))
  }
  
  c=data.frame(real=real.par)
  for(i in 2:length(a)){
    c[[names(a)[i]]]=b[[i]]/s[[i]]
  }
  ns=data.frame(real=real.par)
  for(i in 2:length(a)){
    ns[[names(a)[i]]]=abs(s[[i]]/a[[i]])
  }
  
  #a...absolute werte intercept
  #b...difference real-estimate abs
  #s...standard derivation
  #c...b/s
  #ns..s/a
  return(c)
}

info.table.change.sddiff("Norm","Ml")



