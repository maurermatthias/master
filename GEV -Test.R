source('GEV-estimates.R')

#does a parameter estimation for each stress-level and a
#least square fit for the parameter
g=gev.approach(shen1)

#plot least square fit for parameter
plot.evir.est(g)

#plot resulting 1. and 2. moment
plot.est.moment.ls(g)

#plot iq-area
plot.est.quantile(g,0.9)


#####################
a=mle.gev(g)

#par=c(-0.3,48,367,13800,5566)
#nll.gev(par,a$data)

plot.est.moment.mle(g)


g$mle.est.par.struct
g$est.par.struct

g$est.par.struct
g$mle.est.par.struct

g$optimazation.result.par.struct
g$mle.optimazation.result.par.struct

#########################################
a0=rgev(100,-0.3,300,20)
p0=gev(a0)$par.est;p0
a1=(a0-p0[3])/p0[2]
p1=gev(a1)$par.est;p1

