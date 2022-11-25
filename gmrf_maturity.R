library(TMB)
## library(glmmTMB)
library(abind)

## e <- new.env()
## load("mat_data.Rdata",e)
## e$mat_data

## dat <- with(subset(e$mat_data, SubArea=="NW"),
##             data.frame(y = Maturity,
##                        n = NoAtAlk,
##                        fage = factor(Age),
##                        fyear = factor(Year),
##                        wt = wi)
##             )

## m4 <- glmmTMB(cbind(y,n-y) ~ fage - 1+ (1|fyear) + (1|(fage:fyear)), data = dat,
## weight=wt,family = "betabinomial") 

## ndat <- expand.grid(fage = factor(levels(dat$fage),levels(dat$fage)),
##                                   fyear = factor(levels(dat$fyear),levels(dat$fyear)),
##                                   wt = 1)
## pr <- predict(m4, newdata = ndat,
##         se.fit = TRUE
##         )


## Y <- xtabs(m4$modelInfo$family$linkinv(pr$fit) ~ ndat$fyear + ndat$fage)
## w <- xtabs(pr$se.fit ~ ndat$fyear + ndat$fage)

## write.table(Y,"maturity.tab")
## write.table(w,"maturity_w.tab")

Y <- read.table("maturity.tab")
w <- read.table("maturity_w.tab")

compile("gmrf_maturity.cpp")
dyn.load(dynlib("gmrf_maturity"))

tmb_data <- list(Y = rbind(Y,array(NA,dim=c(10,ncol(Y)))),
                 w = w)
tmb_pars <- list(logPhi = numeric(3),
                 mu = numeric(ncol(Y)),
                 logSdProc = 0,
                 logSdObs = 0,
                 P = array(0, dim = dim(tmb_data$Y)))

obj <- MakeADFun(tmb_data, tmb_pars, random = "P", DLL = "gmrf_maturity")

obj$gr()

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$he <- optimHess(opt$par, obj$n, obj$gr)

round(cov2cor(solve(opt$he)),3)

sdr<-sdreport(obj, opt$par, opt$he)
pl<-as.list(sdr,"Est")
plsd<-as.list(sdr,"Std")

matplot(plogis(pl$P), type="l")
matplot(plogis(pl$P-2*plsd$P), type="l", add=T)
matplot(plogis(pl$P+2*plsd$P), type="l", add=T)
matplot((tmb_data$Y), add=TRUE)

opt$par
