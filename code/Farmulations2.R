load("sim2Data.RData")
library(mg14)
library(CoxHD)
library(parallel)
nData <- c(100, 200, 500, 1000, 2000, 5000, 10000)
nJobs <- as.numeric(Sys.getenv("NJOBS"))

jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))


set.seed(jobIndex)
simCoef <- CoxHD:::SimCoef(coxRFXFitOsTDGGc, groups = simGroups)
simRisk <- as.matrix(simDataFrame[names(whichRFXOsTDGG)]) %*% simCoef[names(whichRFXOsTDGG)]
simSurv <- SimSurvNonp(simRisk, os)

for(n in nData){
	s <- if(n < 10000) sample(1:nrow(simDataFrame), n) else 1:10000
	f <- CoxRFX(simDataFrame[s,names(whichRFXOsTDGG)], simSurv[s], simGroups[names(whichRFXOsTDGG)], nu=1,which.mu=mainGroups)
	f$X <- NULL
	assign(paste0("w",n), s)
	assign(paste0("fit",n),f)	
}

rm(simDataFrame)

save.image(file=paste("simRFX/",Sys.getenv("LSB_JOBNAME"),"_",Sys.getenv("LSB_JOBINDEX"),".RData", sep=""))