#' AML simulations
#' =================

#+ Preliminaries, echo=FALSE
options(width=120)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('png','pdf'), fig.ext=c('png','pdf'), fig.width=4, fig.height=4, smallMar=TRUE)

#' Libraries
library(RColorBrewer)
library(survival)
library(parallel)
library(glmnet)
set1 <- brewer.pal(8, "Set1")

#+ load
load("2014-05-12-AML.RData", envir = globalenv())
library(CoxHD)
source("../../CoxHD/CoxHD/R/functions.R")
ls()

#' 1. Subsets
#' ---------------------------
#' We subset and compute the predictive power of the model
#+ subsetConcordance, cache=TRUE, warning=FALSE, collapse=TRUE
library(survivalROC)
set.seed(42)
subsets <- c(300,600,900,1200,1500)
subsetConcordance <- lapply(subsets, function(s){
			lapply(1:5, function(i){
						trn <- 1:nrow(dataFrame) %in% sample(nrow(dataFrame), s)
						tst <-  !trn #1:nrow(dataFrame) %in% sample(setdiff(1:nrow(dataFrame), trn), 300)
						stbCx <- CoxCPSS(dataFrame[!is.na(survival) & trn,mainIdx], survival[!is.na(survival) & trn], bootstrap.samples = 50, control="BH") ## Only main effects
						slctd <- which(stbCx$Pi > .8)
						cxFt <- coxph(as.formula(paste("survival[!is.na(survival) & trn] ~ ", paste(names(slctd), collapse="+"))), data=dataFrame[!is.na(survival) & trn,mainIdx])
						prdRsk <- predict(cxFt, newdata = dataFrame[!is.na(survival) & tst,mainIdx])
						C <- survConcordance(survival[!is.na(survival) & tst] ~ prdRsk)
						ROC <- survivalROC(Stime=survival[!is.na(survival) & tst,1], status=survival[!is.na(survival) & tst,2], marker = prdRsk, predict.time = 278, method="KM", cut.values=seq(-5,5,0.1))
						list(C, ROC, trn, tst, slctd)})
		})

#+ subsetConcordancePlot, fig.width=4, fig.height=4
#pdf("subsetConcordance.pdf", 2.5,2.5, pointsize=8)
par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
plot(NA,NA, xlim=c(0,1),ylim=c(0,1), xlab="FPR",ylab="TPR")
abline(0,1, lty=3)
for(i in seq_along(subsets)){
	x <- sapply(subsetConcordance[[i]], function(x) x[[2]]$FP)
	y <- sapply(subsetConcordance[[i]], function(x) x[[2]]$TP)
	lines(rowMeans(x),rowMeans(y), col=col1[i], type="l")
}
legend("bottomright", legend=rev(subsets), lty=1, col=col1[5:1], bty="n")

rangeplot <- function(x, y, col = 1, pch = 19, lty = 1, ...){
	plot(x, colMeans(y), col = col, pch=pch, ylim = range(y), ..., xaxt="n")
	axis(at = x, labels=x, side=1)
	segments(x,apply(y,2,min),x,apply(y,2,max), col=col, lty = lty)
}

rangeplot(x=subsets, y = sapply(subsetConcordance, function(x) sapply(1:5, function(i) x[[i]][[2]]$AUC)) , col=col1, xlab="Cohort", ylab="AUC", lty=1)
rangeplot(x=subsets, y = sapply(subsetConcordance, function(x) sapply(1:5, function(i) x[[i]][[1]]$concordance)) , col=col1, xlab="Cohort", ylab="Concordance", lty=1)
rangeplot(x=subsets, y = sapply(subsetConcordance, function(x) sapply(1:5, function(i) length(x[[i]][[5]]))) , col=col1, xlab="Cohort", ylab="Selected variables", lty=1)

#boxplot(as.numeric(sapply(subsetConcordance, function(x) sapply(1:5, function(i) x[[i]][[2]]$AUC))) ~ rep(subsets, each=5), border=col1, xlab="Cohort", ylab="AUC", boxwex=.66, staplewex=0, lty=1, pch=16)
#boxplot(as.numeric(sapply(subsetConcordance, function(x) sapply(1:5, function(i) x[[i]][[1]]$concordance))) ~ rep(subsets, each=5), border=col1, xlab="Cohort", ylab="Concordance", boxwex=.66, staplewex=0, lty=1, pch=16)
#boxplot(as.numeric(sapply(subsetConcordance, function(x) sapply(1:5, function(i) length(x[[i]][[5]])))) ~ rep(subsets, each=5), border=col1, xlab="Cohort", ylab="Selected variables", boxwex=.66, staplewex=0, lty=1, pch=16)
#dev.off()

#' Assuming an oracle told us the variables selected on the entire cohort..
#+ subsetOracle, cache=TRUE
subsetOracle <- lapply(subsetConcordance, function(x) lapply(x, function(y){
						trn <- y[[3]]
						tst <- y[[4]]
						cxFt <- coxph(as.formula(paste("survival[!is.na(survival) & trn] ~ ", paste(selectedInt, collapse="+"))), data=dataFrame[!is.na(survival) & trn,mainIdx])
						prdRsk <- predict(cxFt, newdata = dataFrame[!is.na(survival) & tst,mainIdx])
						C <- survConcordance(survival[!is.na(survival) & tst] ~ prdRsk)
						ROC <- survivalROC(Stime=survival[!is.na(survival) & tst,1], status=survival[!is.na(survival) & tst,2], marker = prdRsk, predict.time = 278, method="KM", cut.values=seq(-5,5,0.1))
						list(C, ROC, trn, tst, selectedInt)})
)

#+ subsetOraclePlot, fig.width=4, fig.height=4
#pdf("subsetOracle.pdf", 2.5,2.5, pointsize=8)
par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
plot(NA,NA, xlim=c(0,1),ylim=c(0,1), xlab="FPR",ylab="TPR")
abline(0,1, lty=3)
for(i in seq_along(subsets)){
	x <- sapply(subsetOracle[[i]], function(x) x[[2]]$FP)
	y <- sapply(subsetOracle[[i]], function(x) x[[2]]$TP)
	lines(rowMeans(x),rowMeans(y), col=col1[i], type="l")
}
legend("bottomright", legend=rev(subsets), lty=1, col=col1[5:1], bty="n")

rangeplot(x=subsets, y = sapply(subsetOracle, function(x) sapply(1:5, function(i) x[[i]][[2]]$AUC)) , col=col1, xlab="Cohort", ylab="AUC", lty=1)
rangeplot(x=subsets, y = sapply(subsetOracle, function(x) sapply(1:5, function(i) x[[i]][[1]]$concordance)) , col=col1, xlab="Cohort", ylab="Concordance", lty=1)
rangeplot(x=subsets, y = sapply(subsetOracle, function(x) sapply(1:5, function(i) length(x[[i]][[5]]))) , col=col1, xlab="Cohort", ylab="Selected variables", lty=1)

#boxplot(as.numeric(sapply(subsetOracle, function(x) sapply(1:5, function(i) x[[i]][[2]]$AUC))) ~ rep(subsets, each=5), border=col1, xlab="Cohort", ylab="AUC", boxwex=.66, staplewex=0, lty=1, pch=16)
#boxplot(as.numeric(sapply(subsetOracle, function(x) sapply(1:5, function(i) x[[i]][[1]]$concordance))) ~ rep(subsets, each=5), border=col1, xlab="Cohort", ylab="Concordance", boxwex=.66, staplewex=0, lty=1, pch=16)
#boxplot(as.numeric(sapply(subsetOracle, function(x) sapply(1:5, function(i) length(x[[i]][[5]])))) ~ rep(subsets, each=5), border=col1, xlab="Cohort", ylab="Selected variables", boxwex=.66, staplewex=0, lty=1, pch=16)
#dev.off()

#' 2. Survival only
#' ----------------
#' Non-parametric function to simulate survival based on an observed survival curve
SimSurvNonp

#' Simulate s
set.seed(42)
partRisk <- PartialRisk(coxRFXFit)
totalRisk <- rowSums(partRisk)
simSurvival <- SimSurvNonp(risk=totalRisk, H0=basehaz(coxph(survival ~ totalRisk), centered=FALSE), surv=survival)
plot(survfit(simSurvival ~ 1))
lines(survfit(survival ~ 1), col="red")
survConcordance(simSurvival[trainIdx] ~ totalRisk[trainIdx])
survConcordance(simSurvival[trainIdx][-na.action(coxFit)] ~ predict(coxFit))


#+ simCoxRFXFit, cache=TRUE
simCoxRFXFit <- CoxRFX(dataFrame[,whichRFX], simSurvival, groups=groups[whichRFX], sigma0 = 0.1)

#+ simCoxRFXFitBox
par(mar=c(5,7,1,1))
boxplot(coef(simCoxRFXFit)~groups[whichRFX], border=col1, horizontal=TRUE, las=2)
abline(v=0)

#+ simCoxRFXFitScatter
plot(coef(simCoxRFXFit), coef(coxRFXFit), col = col1[groups[whichRFX]])

#+ simPredictedRisk
simPredictedRisk <- PredictRiskMissing(simCoxRFXFit)
plot(totalRisk, simPredictedRisk[,1], xlab="True risk (simulated)",ylab="Estimated risk")
segments(totalRisk, simPredictedRisk[,1] + sqrt(simPredictedRisk[,2]), totalRisk, simPredictedRisk[,1] - sqrt(simPredictedRisk[,2]), col="#00000022")
abline(0,1, col="red")

#+ simPartialRisk, fig.width=8, fig.height=8
par(mfrow=c(3,3))
p <- PartialRisk(simCoxRFXFit)
v <- PartialRiskVar(simCoxRFXFit)
for(i in 1:8){
	plot(partRisk[,i], p[,i], main=colnames(partRisk)[i], xlab="True risk component (simulated)",ylab="Estimated risk component")
	segments(partRisk[,i], p[,i] + sqrt(v[,i]), partRisk[,i], p[,i] - sqrt(v[,i]), col="#00000022")
	abline(0,1, col="red")
}

#' #### Stability selection with simulated survival
#+ simCoxCPSS, cache=TRUE
set.seed(42)
s <- simSurvival
s[s[,1]==0,1] <- .5 #Machine$double.eps
simCoxCPSS <- CoxCPSS(dataFrame[!is.na(survival),], s[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

#' Plots
#+ simCoxCPSSPlot
plot(simCoxCPSS)
plot(simCoxCPSS$Pi, coxCPSS$Pi)

#' Selected variables
simSelected <- which(simCoxCPSS$Pi > .8)
simSelected
simCoxFit <- coxph(simSurvival[trainIdx] ~ ., data=dataFrame[trainIdx,][,names(simSelected)])
summary(simCoxFit)

#' #### Stability selection with interactions
#+ simCoxCPSSInt, cache=TRUE
set.seed(42)
simCoxCPSSInt <- CoxCPSSInteractions(dataFrame[groups %in% mainGroups], s, bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))

#+ simCoxCPSSIntPlot
plot(simCoxCPSSInt)
plot(coxCPSSInt$Pi, simCoxCPSSInt$Pi[names(coxCPSSInt$Pi)])


#' Selected variables
simSelectedInt <- which(simCoxCPSSInt$Pi > .8)
simSelectedInt
simCoxFitInt <- coxph(simSurvival[trainIdx] ~ ., data=dataFrame[trainIdx,][,names(simSelectedInt)])
summary(simCoxFitInt)


#' Coefficients
set.seed(42)
v <- coxRFXFit$sigma2 * coxRFXFit$df[-8] / as.numeric(table(groups[whichRFXTD])-1)
v["GeneGene"] <- v["GeneTreat"] <- v["CytoTreat"] <- v["GeneCyto"] <-  0.1
simCoef <- numeric(ncol(simDataFrame))
names(simCoef) <- colnames(simDataFrame)
simCoef[whichRFXTD] <- rnorm(length(groups[whichRFXTD]), coxRFXFit$mu[groups[whichRFXTD]], sqrt(v[groups[whichRFXTD]]))

#' Survival
set.seed(42)
simDataRisk <- (as.matrix(simDataFrame) %*% simCoef)[,1] 
simDataRisk <- simDataRisk #- mean(simDataRisk)
simDataSurvival <- SimSurvNonp(simDataRisk, H0=basehaz(coxph(survival ~ totalRisk), centered=FALSE), surv=survival)
f <- as.formula(paste("simDataSurvival ~", names(which.max(simCoef))))
f
plot(survfit(f, data=simDataFrame))

#' ### 1. RFX model
#+ simDataCoxRFXFit, cache=TRUE
simDataCoxRFXFit <- CoxRFX(simDataFrame[whichRFXTD], simDataSurvival, groups = groups[whichRFXTD], sigma0 = 0.1)
plot(simCoef[whichRFXTD],coef(simDataCoxRFXFit), col=col1[groups[whichRFXTD]])
plot(sapply(split(simCoef[whichRFXTD], groups[whichRFXTD]), var),simDataCoxRFXFit$sigma2, col=col1, pch=19)
plot(sapply(split(simCoef[whichRFXTD], groups[whichRFXTD]), mean),simDataCoxRFXFit$mu, col=col1, pch=19)


#' ### 2. Stability selection
#+ simDataCoxCPSS, cache=TRUE
set.seed(42)
s <- simDataSurvival
s[s[,1]==0,1] <- .Machine$double.eps
simDataCoxCPSS <- CoxCPSS(simDataFrame[!is.na(survival),], s[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)
#save(simDataCoxCPSS, file="simDataCoxCPSS")

#' Plots
plot(simDataCoxCPSS)
plot(simDataCoxCPSS$Pi, coxCPSS$Pi)
plot(simCoef, simDataCoxCPSS$Pi)
plot(colMeans(simDataFrame), simDataCoxCPSS$Pi, cex=sqrt(abs(simCoef))*2+.1, log="x")


#' Selected variables
simSelectedInt <- which(simDataCoxCPSS$Pi > .8)
simSelectedInt
simDataCoxFit <- coxph(simDataSurvival[trainIdx] ~ ., data=simDataFrame[trainIdx,][,names(simSelectedInt)])
summary(simDataCoxFit)


#' ### Test interactions without main effects
#+ simDataPIntNoMain, cache=TRUE
set.seed(42)
intIdx <- groups %in% c("GeneGene","GeneTreat","GeneCyto")
scope <- strsplit(colnames(dataFrame)[intIdx],":")

simCoefInt <- simCoef
simCoefInt[sample(which(!mainIdx), 500)] <- 0
size <- 10000#nrow(dataFrame)
s <- sample(10000, size)
simDataRiskInt <- (as.matrix(simDataFrame[s,]) %*% simCoefInt)
simDataRiskInt <- simDataRiskInt - mean(simDataRiskInt)
simDataSurvivalInt <- SimSurvNonp(simDataRiskInt, H0=basehaz(coxph(survival ~ totalRisk), centered=FALSE), surv=survival)
simDataPIntNoMain <- TestInteractions(simDataFrame[s,mainIdx], simDataSurvivalInt, scope, whichMain=NULL)

idx <- simCoefInt[intIdx] == 0;

#+ simCoefIntWald
plot( simCoefInt[intIdx], simDataPIntNoMain$pWald, log='y', col=ifelse(idx,2,1))
abline(h=max(simDataPIntNoMain$pWald[p.adjust(simDataPIntNoMain$pWald, "BH")<0.1], na.rm=TRUE))
qqplot(simDataPIntNoMain$pWald[idx], simDataPIntNoMain$pWald[!idx], log='xy', xlab="Coef == 0", ylab="Coef != 0")
abline(0,1)
title("QQ-plot")

#+ redHerring, fig.width=6, fig.height=8
redHerring <- rownames(simDataPIntNoMain)[idx][head(order(simDataPIntNoMain$pWald[idx]))]
redHerring
par(mfrow=c(3,2))
for(h in redHerring)
	plot(survfit(as.formula(paste("simDataSurvivalInt ~", sub(":","+", h))), data=simDataFrame[s,] ), col=col1)

#+ simCoefIntLR
par(mfrow=c(1,1))
plot( simCoefInt[intIdx], simDataPIntNoMain$pLR, log='y', col=ifelse(simCoefInt[intIdx]==0,2,1))
qqplot(simDataPIntNoMain$pLR[idx], simDataPIntNoMain$pLR[!idx], log='xy', xlab="Coef == 0", ylab="Coef != 0")
abline(0,1)
title("QQ-plot")

#+ simCoefIntRisk
meanRisk <- sapply(scope, function(p) mean(simDataRiskInt[simDataFrame[s,p[1]] * simDataFrame[s,p[2]]==1]) - simCoefInt[p[1]] - simCoefInt[p[2]])
plot(simCoefInt[intIdx][simDataPIntNoMain$warn==0], simDataPIntNoMain$coef[simDataPIntNoMain$warn==0], xlab="True coefficient", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(simCoefInt[intIx][simCoefInt[intIx]!=0], simDataPIntNoMain$coef[simCoefInt[intIx]!=0], use="complete", method="spearman"),2)))
plot(meanRisk[simDataPIntNoMain$warn==0], simDataPIntNoMain$coef[simDataPIntNoMain$warn==0], xlab="Residual risk", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(meanRisk, simDataPIntNoMain$coef, use="complete", method="spearman"),2)))

#' ### Test interactions with main effects
#+ simDataPIntMain, cache=TRUE
simDataPIntMain <- TestInteractions(simDataFrame[s,mainIdx], simDataSurvivalInt, scope, whichMain=colnames(dataFrame)[mainIdx], mc.cores=8)

#load("temp.RData")
plot(simCoefInt[intIx], simDataPIntMain$coef, xlab="True coefficient", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(simCoefInt[intIx][simCoefInt[intIx]!=0], simDataPIntMain$coef[simCoefInt[intIx]!=0], use="complete", method="spearman"),2)))
plot(meanRisk, simDataPIntMain$coef, xlab="Residual risk", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(meanRisk, simDataPIntMain$coef, use="complete", method="spearman"),2)))

plot(simCoefInt[intIx],simDataPIntNoMain$pWald, log="y", main="Marginal only", col=ifelse(simCoefInt[intIx]==0,2,1), ylab="P-value")
abline(h=max(simDataPIntNoMain$pWald[p.adjust(simDataPIntNoMain$pWald, "BH")<0.1], na.rm=TRUE), lty=3)
abline(h=max(simDataPIntNoMain$pWald[p.adjust(simDataPIntNoMain$pWald, "holm")<0.05], na.rm=TRUE), lty=2)

plot(simCoefInt[intIx],simDataPIntMain$pWald, log="y", main="Incl. all mains", col=ifelse(simCoefInt[intIx]==0,2,1), ylab="P-value")
abline(h=max(simDataPIntMain$pWald[p.adjust(simDataPIntMain$pWald, "BH")<0.1], na.rm=TRUE), lty=3)
abline(h=max(simDataPIntMain$pWald[p.adjust(simDataPIntMain$pWald, "holm")<0.05], na.rm=TRUE), lty=2)
#
#
#simResults <- list()
#for(size in c(100,1000,10000)){
#	method <- "BH"
#	simResults[[as.character(size)]]$tpr.bh <- sapply(simulations[[as.character(size)]], function(sim){
#				c <- cut(colSums(simDataFrame[sim$s,!mainIdx]), breaks = c(0,1,10,50,100,500,1000, 5000), include.lowest=TRUE)
#				idx <- sim$simCoef[!mainIdx] != 0
#				sapply(split((p.adjust(sim$simDataPInt$pLR,method) < 0.1)[idx], c[idx]), mean)
#			})
#	simResults[[as.character(size)]]$fpr.bh <- sapply(simulations[[as.character(size)]], function(sim){
#				c <- cut(colSums(simDataFrame[sim$s,!mainIdx]), breaks = c(0,1,10,50,100,500,1000, 5000), include.lowest=TRUE)
#				idx <- sim$simCoef[!mainIdx] == 0
#				sapply(split((p.adjust(sim$simDataPInt$pLR,method) < 0.1)[idx], c[idx]), mean)
#			})
#	method <- "bonf"
#	simResults[[as.character(size)]]$tpr.bonf <- sapply(simulations[[as.character(size)]], function(sim){
#				c <- cut(colSums(simDataFrame[sim$s,!mainIdx]), breaks = c(0,1,10,50,100,500,1000, 5000), include.lowest=TRUE)
#				idx <- sim$simCoef[!mainIdx] != 0
#				sapply(split((p.adjust(sim$simDataPInt$pLR,method) < 0.1)[idx], c[idx]), mean)
#			})
#	simResults[[as.character(size)]]$fpr.bonf <- sapply(simulations[[as.character(size)]], function(sim){
#				c <- cut(colSums(simDataFrame[sim$s,!mainIdx]), breaks = c(0,1,10,50,100,500,1000, 5000), include.lowest=TRUE)
#				idx <- sim$simCoef[!mainIdx] == 0
#				sapply(split((p.adjust(sim$simDataPInt$pLR,method) < 0.1)[idx], c[idx]), mean)
#			})
#}
#
#plot(sapply(split(p.adjust(simDataPInt$pWald,"BH") < 0.1 & simCoef[-(1:ncol(X))]!=0, c), mean), col="red", ylab="TPR", type="l", lty=3, xaxt="n", xlab="")
#axis(side=1, levels(c), at=1:nlevels(c))
#lines(sapply(split(p.adjust(simDataPInt$pWald,"BH") < 0.1 & simCoef[-(1:ncol(X))]==0, c), mean), lty=3)
#lines(sapply(split(p.adjust(simDataPInt$pWald) < 0.1 &  c), mean), col="red", ylab="TPR", type="l")
#lines(sapply(split(p.adjust(simDataPInt$pWald) < 0.1 & simCoef[-(1:ncol(X))]==0, c), mean))

#' ### Newer simulations: Keep all others fixed
SimTreatmentSurvival <- function(X, coxRFXFit, treatEffect, whichGene, geneTreatEffect) {
	treatment <- rbinom(nrow(X), 1, 0.5) ## assume randomize treatment
	Y <- cbind(X, Test = treatment)
	risk <- as.matrix(X) %*% coxRFXFit$coef + treatment * treatEffect +  X[,whichGene] * treatment * geneTreatEffect 
	simSurvival <-  SimSurvNonp(risk, H0=basehaz(coxph(coxRFXFit$surv ~ predict(coxRFXFit, newdata = as.data.frame(coxRFXFit$X))), centered=FALSE), surv=coxRFXFit$surv)
	return(list(Y=Y, simSurvival=simSurvival))
}
SimAndTest <- function(X, coxRFXFit, whichGene, treatEffect = -0.5, geneTreatEffect = .5, nSim = 1, whichMain = colnames(X)){
	sim <- SimTreatmentSurvival(X = X, coxRFXFit = coxRFXFit, treatEffect = treatEffect, whichGene = whichGene, geneTreatEffect = geneTreatEffect)
	TestInteractions(sim$Y, sim$simSurvival, list(c("Test", whichGene)), whichMain=whichMain)
}

#+ simGenes, cache=TRUE, fig.width=6
set.seed(42)
effectSizes <- seq(-1,1,0.1)
whichMain <- colnames(simDataFrame)[groups %in% c("Clinical","Genetics","Cytogenetics","Treatment")]
simASXL1 <- lapply(c(1000,10000), function(nSim) sapply(effectSizes, function(e) sapply(1:10, function(i) SimAndTest(simDataFrame[1:nSim,whichRFX], coxRFXFit, "ASXL1", geneTreatEffect = e, whichMain = whichMain)[1,1])))
simNPM1 <- lapply(c(1000,10000), function(nSim) sapply(effectSizes, function(e) sapply(1:10, function(i) SimAndTest(simDataFrame[1:nSim,whichRFX], coxRFXFit, "NPM1", geneTreatEffect = e, whichMain = whichMain)[1,1])))

par(mfrow=c(1,2))
for(x in c("simASXL1","simNPM1")){
	tmp <- get(x)
	boxplot(tmp[[1]], log="y", names=effectSizes, border="grey", ylim=range(tmp[[2]]))
	boxplot(tmp[[2]], log="y", names=effectSizes, add=TRUE, col=NA)
	abline(h=0.05)
	abline(h=0.05 / length(pairs))
	title(xlab='Interaction effect size', ylab="P-value", main=x)
}



#' 4. Comprehensive simulations
#' ------------------
#' In this section we will use comprehensive simulations over a range of parameter values to assses the power for identifying gene-treatment interactions.
#' The treatment will be assumed to be randomized and given in on of two arm of equal size. Parameters are
#' * Cohort size 100,1000,10000
#' * Mutation frequencies as observed
#' * Baseline treatment effect -2 ... 0
#' * Gene-treatment interactions -2 ... +2
#' We run 100 simulations for each case
#' ### Prepare data
save("coxRFXFit","simDataFrame","whichRFX", "groups","TestInteractions", "SimSurvNonp","SimAndTest","SimTreatmentSurvival", file="../simulations/simulationData.RData")

#' ### Simulations
#' The following R code is executed on the farm
#+ farmulations, eval=FALSE, cache=FALSE
read_chunk('Farmulations.R')

#' ### Load data 
#' Now load the data
#+ simResults, cache=TRUE
simResults <- do.call("rbind",lapply(dir("../simulations/output-2014-05-20", pattern = "*.RData", full.names = TRUE), function(file){
			load(file)
			return(simResults)
		}))
simResults <- simResults[order(simResults$indeces),]
simResults$pWald[simResults$warn==1] <- NA

interactionEffects <- seq(-2,2,0.25)
nData <- c(100, 1000, 10000)
treatmentEffects <- seq(0,-2,-0.25)
nSim <- 100

simPwaldArray <-  array(simResults$pWald, dim = (c(length(nData), sum(groups=="Genetics"),length(treatmentEffects), length(interactionEffects), nSim)))
simPwaldArray <- aperm(simPwaldArray, 5:1)
x <- colMeans(simDataFrame[,groups=="Genetics"])
o <- order(x)
l <- (length(interactionEffects) +1)/2

#' The following plot show the P-values as a function of mutations frequency
#+ frequencyPval
plot(NA,NA, xlim=range(x), ylim=c(1e-40,1), log="y", xlab="Mutation frequency", ylab="P-value")
for(n in seq_along(nData))
	for(e in 1:l){
		lines(x=x[o],y=colMeans(simPwaldArray[,2*e-1,,o,n], dims=2, na.rm=TRUE), col=brewer.pal(l,"Spectral")[e], lty=4-n, type="b", pch=n)
	}
legend("bottomleft", legend=interactionEffects[seq(1,17,2)], col = brewer.pal(l,"Spectral"), lty=1, bty='n')
legend("bottom", pch=1:3, legend=nData, bty="n")

#' There are some regularities, which can be collapsed into the invariant quantity mutation frequency x cohort size x effect^2
#+ effectPval
plot(NA,NA, xlim=c(1,2500), ylim=c(1e-21,1), log="xy", xlab="Mutation number x effect^2", ylab="P-value")
for(n in seq_along(nData))
	for(e in 1:l)
		points(x=x[o]*nData[n]*abs(interactionEffects[2*e-1])^2,y=(colMeans(simPwaldArray[,2*e-1,,o,n], dims=2, na.rm=TRUE)), col=brewer.pal(9,"Spectral")[e], lty=4-n, pch=n)
xx <- x[o]*nData[n]*abs(interactionEffects[2*e-1])^2
lines(xx, 10^-(xx/100))
legend("bottomleft", legend=interactionEffects[seq(1,17,2)], col = brewer.pal(l,"Spectral"), lty=1, bty='n')
legend("bottom", pch=1:3, legend=nData, bty="n")

#' Power as a function of relative mutation frequency
#+ frequencyPower
plot(NA,NA, xlim=range(x)+ 1e-3, ylim=c(0,1), xlab="Mutation frequency", ylab="Power", log="x")
for(n in seq_along(nData))
	for(e in 1:l){
		y <- colMeans(simPwaldArray[,2*e-1,,o,n] < 5e-2 / (length(x)), dims=2, na.rm=TRUE)
		points(x[o],y, col=brewer.pal(l,"Spectral")[e], lty=4-n)
		s <- smooth.spline(log(x[o] + .5/max(nData))[!is.na(y)],na.omit(y), spar=.5)
		lines(exp(s$x),s$y, col=brewer.pal(11,"Spectral")[e], lty=4-n)
	}
legend("topleft", legend=interactionEffects[seq(1,length(interactionEffects),2)], col = brewer.pal(l,"Spectral"), lty=1, bty='n')

#' Power as a function of the invariant effect size
#+ effectPower
plot(NA,NA, xlim=c(1,2500), ylim=c(0,1), xlab="Mutation number x effect^2", ylab="Power", log="x")
for(n in seq_along(nData))
	for(e in 1:l){
		y <- colMeans(simPwaldArray[,2*e-1,,o,n] < 5e-2 / (length(x)), dims=2, na.rm=TRUE)
		points(x[o] *nData[n]* abs(interactionEffects[2*e-1])^2, y, col=brewer.pal(l,"Spectral")[e], lty=4-n, pch=n)
		s <- smooth.spline(log(x[o]+ .5/max(nData))[!is.na(y)], na.omit(y), spar=.66)
		lines(exp(s$x)*nData[n]*abs(interactionEffects[2*e-1])^2,s$y, col=brewer.pal(11,"Spectral")[e], lty=4-n)
	}
legend("topleft", legend=interactionEffects[seq(1,length(interactionEffects),2)], col = brewer.pal(11,"Spectral"), lty=1, bty='n')
legend("bottomleft", pch=1:3, legend=nData, bty="n")


#' Effect size vs. mutation freq
lev <-  c(1,2,5) * rep(10^(-3:-1), each=3)
c <- cut(x,lev)
tmp <- lapply( seq_along(nData), function(n){
			sapply(levels(c), function(lev){
						sapply(1:l, function(e){
									mean(simPwaldArray[,2*e-1,,which(c==lev),n, drop=FALSE] < 5e-2 / (length(x)), na.rm=TRUE)
								})
					})})

#+ hazardFreq, fig.width=6, fig.height=6
par(mfrow=c(2,2))
for(n in 1:3){
	plot(NA,NA, xlim=exp(c(-1,1)), ylim=range(log10(lev)), xlab="Interaction hazard", ylab="Mutation frequency",  yaxt="n", log="x")
	.filled.contour(z=tmp[[n]], x = exp(interactionEffects[seq(1,length(interactionEffects),2)]), y=log10(lev[-1]),levels=seq(0.1,1.1,0.1), col=colorRampPalette(brewer.pal(9,"Reds"))(11))
	contour(z=tmp[[n]], x = exp(interactionEffects[seq(1,length(interactionEffects),2)]), y=log10(lev[-1]), levels=seq(0.1,1,0.1), col="black", add=TRUE, axes=FALSE)
	axis(2, at=log10(lev[-9]), labels=lev[-9])
	title(paste('n =' ,nData[n]))
	abline(v=c(1/seq(1,2.5,.5),seq(1,2.5,.5)), col="grey")
}

#' Interaction effect vs main effect
lev <-  c(0,c(1,2,5) * rep(10^(0:3), each=3),10000, 20000)
X <- expand.grid(nData=nData, gene=colnames(simDataFrame)[which(groups=="Genetics")], treatmentEffects=treatmentEffects, interactionEffects=interactionEffects, nSim=1:nSim)
numTimesEffect <- cut(X$nData * x[X$gene] * X$interactionEffects^2, breaks=lev)
tmp <- 	sapply(levels(numTimesEffect), function(l){
						w <- numTimesEffect == l
						sapply(rev(treatmentEffects), function(t){
									mean(simResults$pWald[w & X$treatmentEffects==t] < 5e-2 / (length(x)), na.rm=TRUE)
								})
					})

#+ interactionMain, fig.width=3, fig.height=3
par(mfrow=c(1,1))
plot(NA,NA, xlim=range(treatmentEffects), ylim=range((lev[-1])), xlab="Treatment risk", ylab="Mutation frequency", log="y", xaxs="i", yaxs="i")
image(z=tmp, x = rev(treatmentEffects), y=lev,breaks=c(seq(0,1,0.01)), col=colorRampPalette(brewer.pal(9,"Reds")[-9])(100), add=TRUE, right=FALSE) 
contour(z=tmp, x = rev(treatmentEffects), y=lev[-1] - diff(lev), levels=seq(0.1,1,0.1), col="black", add=TRUE, axes=FALSE)
abline(v=c(1/seq(1,2.5,.5),seq(1,2.5,.5)), col="grey")


set.seed(42)
sigma2 <- coxRFXFitTD$sigma2
sigmaX <- apply(simDataFrame,2,var)

c <- sapply(1:20, function(i){
			#simRisk <- as.matrix(simDataFrame) %*% rnorm(length(groups), mean = coxRFXFitTD$mu[groups], sd = sqrt(sigma2)[groups])
			simRisk <- as.matrix(simDataFrame[whichRFX]) %*% rnorm(length(groups[whichRFX]), mean = coxRFXFitTD$mu[groups[whichRFX]], sd = sqrt(sigma2)[groups[whichRFX]])
			s <- SimSurvNonp(simRisk, coxRFXFit$surv, H0 = basehaz(coxRFXFit, centered = FALSE))
			c(var(simRisk),survConcordance(s ~ simRisk)$concordance)
		})

v <- sapply(1:100, function(i){
			simRisk <- as.matrix(simDataFrame) %*% rnorm(length(groups), mean = coxRFXFitTD$mu[groups], sd = sqrt(sigma2)[groups])
			#simRisk <- as.matrix(simDataFrame[whichRFX]) %*% rnorm(length(groups[whichRFX]), mean = coxRFXFitTD$mu[groups[whichRFX]], sd = sqrt(sigma2)[groups[whichRFX]])
			var(simRisk)
		})

simRisk <- as.matrix(simDataFrame) %*% rnorm(length(groups), mean = coxRFXFitTD$mu[groups], sd = sqrt(sigma2)[groups])
c <- sapply(seq(0,5,0.1), function(i){
			#simRisk <- as.matrix(simDataFrame[whichRFX]) %*% rnorm(length(groups[whichRFX]), mean = coxRFXFitTD$mu[groups[whichRFX]], sd = sqrt(sigma2)[groups[whichRFX]])
			s <- SimSurvNonp(simRisk*i, coxRFXFit$surv, H0 = basehaz(coxRFXFit, centered = FALSE))
			c(var(simRisk * i),survConcordance(s ~ simRisk)$concordance)
		})