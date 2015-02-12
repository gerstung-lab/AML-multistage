load("sim2Data.RData")
source("../../mg14/R/mg14.R")
library(CoxHD)
library(parallel)
nData <- c(100, 1000, 10000)
nJobs <- as.numeric(Sys.getenv("NJOBS"))

jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))


set.seed(jobIndex)
simCoef <- CoxHD:::SimCoef(coxRFXFitOsTDGGc, groups = simGroups)
simRisk <- as.matrix(simDataFrame[names(whichRFXOsTDGG)]) %*% simCoef[names(whichRFXOsTDGG)]
simSurv <- SimSurvNonp(simRisk, os)

w100 <- sample(1:nrow(simDataFrame), 100)
fit100 <- CoxRFX(simDataFrame[w100,names(whichRFXOsTDGG)], simSurv[w100], simGroups[names(whichRFXOsTDGG)], nu=1,which.mu=mainGroups)
fit100$X <- NULL

w1000 <- sample(1:nrow(simDataFrame), 1000)
fit1000 <- CoxRFX(simDataFrame[w1000,names(whichRFXOsTDGG)], simSurv[w1000], simGroups[names(whichRFXOsTDGG)],which.mu=mainGroups)
fit1000$X <- NULL

fit10000 <- CoxRFX(simDataFrame[names(whichRFXOsTDGG)], simSurv, simGroups[names(whichRFXOsTDGG)],which.mu=mainGroups)
fit10000$X <- NULL

rm(simDataFrame)

save.image(file=paste("simRFX/",Sys.getenv("LSB_JOBNAME"),"_",Sys.getenv("LSB_JOBINDEX"),".RData", sep=""))