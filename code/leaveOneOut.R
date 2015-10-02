load("loo.RData")
library(mg14)
library(CoxHD)
library(Rcpp)

#save(dataFrame, nrdData, crGroups, mainGroups, prdData, relData, prdData, osData, cr, dataFrameOsTD, osTD, tplSplitOs, groups, data, whichRFXOsTDGG, clinicalData, MultiRFX5,  file="../../code/loo.RData")

jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

cvIdx <- 1:nrow(dataFrame)
whichTrain <- which(cvIdx != jobIndex)
rfxNrs <- CoxRFX(nrdData[nrdData$index %in% whichTrain, names(crGroups)], Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
rfxNrs$coefficients["transplantRel"] <- 0
#prsData$time1[!is.na(prsData$time1)] <- 0
rfxPrs <-  CoxRFX(prdData[prdData$index %in% whichTrain, names(crGroups)], Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% whichTrain], groups=crGroups, nu=1, which.mu = intersect(mainGroups, unique(crGroups)))
rfxRel <-  CoxRFX(relData[relData$index %in% whichTrain, names(crGroups)], Surv(relData$time1, relData$time2, relData$status)[relData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
rfxRel$coefficients["transplantRel"] <- 0
rfxCr <- CoxRFX(osData[whichTrain, names(crGroups)], Surv(cr[,1], cr[,2]==2)[whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
rfxEs <- CoxRFX(osData[whichTrain, names(crGroups)], Surv(cr[,1], cr[,2]==1)[whichTrain], groups=crGroups, which.mu = NULL)

ix <- tplSplitOs %in% whichTrain
rfxOs <- CoxRFX(dataFrameOsTD[ix,whichRFXOsTDGG], osTD[ix], groups[whichRFXOsTDGG], which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 

xx <- 0:2000
coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])[prdData$index %in% whichTrain,]) 
tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						

coxphOs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))[osData$index %in% whichTrain,]) 
exp(pmin(predict(coxphOs, newdata=data.frame(time0=500)),predict(coxphOs, newdata=data.frame(time0=xx[-1])))) ## cap predictions at induction length 500 days.
multiRfx5 <- MultiRFX5(rfxEs, rfxCr, rfxNrs, rfxRel, rfxPrs, data[cvIdx == jobIndex,,drop=FALSE], tdPrmBaseline = tdPrmBaseline, tdOsBaseline = tdOsBaseline, x=2000)

save(rfxEs, rfxCr, rfxEs, rfxNrs, rfxPrs, rfxRel, rfxOs, multiRfx5, file=paste0("loo/",jobIndex,".RData"))

