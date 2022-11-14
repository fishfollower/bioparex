library(TMB)

sw<-stockassessment::read.ices("sw.dat")

data<-list(Y=log(sw))

for(i in 1:10)data$Y<-rbind(data$Y, NA) # for predictions

param<-list()
param$logPhi<-c(0,0)
param$mu<-numeric(ncol(data$Y)) # one mean for each age group
param$logSdProc<-0
param$logSdObs<-0
param$P<-matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))

compile("gmrf.cpp")
dyn.load(dynlib("gmrf"))

obj<-MakeADFun(data,param,random="P", DLL="gmrf")
opt<- nlminb(obj$par,obj$fn, obj$gr)
sdr<-sdreport(obj)
pl<-as.list(sdr,"Est")
plsd<-as.list(sdr,"Std")

matplot(exp(pl$P), type="l")
matplot(exp(pl$P-2*plsd$P), type="l", add=T)
matplot(exp(pl$P+2*plsd$P), type="l", add=T)

matplot(exp(data$Y), add=TRUE)

