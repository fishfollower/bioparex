library(TMB)

sw<-stockassessment::read.ices("sw.dat")

data<-list(Y=log(sw))

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

matplot((pl$P), type="l", col="gray")
matplot((pl$P-2*plsd$P), type="l", add=T, col="gray")
matplot((pl$P+2*plsd$P), type="l", add=T, col="gray")
matplot((data$Y), add=TRUE)

runwith<-function(include){
  thisdat<-data
  thisdat$Y[-include,]<-NA    
  obj<-MakeADFun(thisdat,param,random="P", DLL="gmrf")
  opt<- nlminb(obj$par,obj$fn, obj$gr)
  sdr<-sdreport(obj)
  pl<-as.list(sdr,"Est")
  plsd<-as.list(sdr,"Std")
  matplot((pl$P), type="l", add=TRUE, lwd=2)
  matplot((pl$P-2*plsd$P), type="l", add=T, lwd=2)
  matplot((pl$P+2*plsd$P), type="l", add=T, lwd=2)
}


runwith(1:10)


