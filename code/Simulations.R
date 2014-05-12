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
simSurvival <- SimSurvNonp(totalRisk, basehaz(coxph(survival ~ totalRisk), centered=FALSE), survival, span=.1)
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

#' 3. Extrapolations
#' -----------------
#' Simulate data 
#+ simData, cache=TRUE
set.seed(42)
SimDataNonp
data <- data.frame(Clinical=dataList$Clinical, Genetics=dataList$Genetics, Cytogenetics=dataList$Cytogenetics, Treatment=dataList$Treatment)
simData <- SimDataNonp(data, nData = 10000, m=10)
names(simData) <- names(data)
save(file="simData.RData", simData)
#load(file="simData.RData")

g <- sub("\\..+","", colnames(data))
colnames(simData) <- gsub("(Clinical|Genetics|Treatment|Cytogenetics).","", colnames(simData))
simDataFrame <- data.frame(simData,
		MakeInteractions(simData[,g=="Genetics"], simData[,g=="Genetics"])[,as.vector(upper.tri(matrix(0,ncol=sum(g=="Genetics"), nrow=sum(g=="Genetics"))))],
		MakeInteractions(simData[,g=="Genetics"], simData[,g=="Cytogenetics"]),
		MakeInteractions(simData[,g=="Genetics"], simData[,g=="Treatment"]),
		MakeInteractions(simData[,g=="Cytogenetics"], simData[,g=="Treatment"]), check.names = FALSE
)
for(n in unique(which(is.na(simDataFrame), arr.ind = TRUE)[,2]))
	simDataFrame[[n]] <- poorMansImpute(simDataFrame[[n]])
simDataFrame <- StandardizeMagnitude(simDataFrame)
simDataFrame <- simDataFrame[,colnames(dataFrame)]

#simDataFrameTD <- simDataFrame[tplSplit,]
#simDataFrameTD[which(tplIndex), grep("TPL", colnames(simDataFrameTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

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
simDataRisk <- simDataRisk - mean(simDataRisk)
simDataSurvival <- SimSurvNonp(simDataRisk, basehaz(coxRFXFit, centered=FALSE), survival, span=.1)
plot(survfit(simDataSurvival ~ Genetics.EP300, data=simDataFrame))

#' ### 1. RFX model
#+ simDataCoxRFXFit, cache=TRUE
simDataCoxRFXFit <- CoxRFX(simDataFrame[, names(coef(coxRFXFit))], simDataSurvival, groups = groups[whichRFXTD], sigma0 = 0.1)
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


#' ### Test interactions
X <- cbind(dataList$Genetics, dataList$Cytogenetics,dataList$Treatment,dataList$Clinical)
scope <- sub("(GeneGene|GeneTreat|CytoTreat|GeneCyto).","",colnames(dataFrame)[grep("(GeneGene|GeneTreat|CytoTreat|GeneCyto)", names(dataFrame))])
pInt <- TestInteractions(X, survival, scope)

simPInt <- TestInteractions(X, simSurvival, scope)

simDataPInt <- TestInteractions(simDataFrame[,groups %in% c("Clinical","Genetics","Treatment","Cytogenetics")], simDataSurvival, scope)
simDataPInt <- TestInteractions(simDataFrame[1:1000,groups %in% c("Clinical","Genetics","Treatment","Cytogenetics")], simDataSurvival[1:1000], scope)


set.seed(42)
v <- coxRFXFit$sigma2 * coxRFXFit$df[-8] / as.numeric(table(groups[whichRFXTD])-1)
#v["GeneGene"] <- v["GeneTreat"] <- v["CytoTreat"] <- v["GeneCyto"] <-  0.1
mainIdx <- groups %in% c("Clinical","Genetics","Treatment","Cytogenetics")
pairs <- GetPairs(names(simDataFrame)[mainIdx], scope)
#simulations <- list()
#for(size in c(100,1000, 10000)){
#	simulations[[as.character(size)]] <- list()
#	for(i in 1:3){
#		cat(".")
simCoef <- numeric(ncol(simDataFrame))
names(simCoef) <- colnames(simDataFrame)
w <- c(which(mainIdx), sample(which(!mainIdx), 500))
simCoef[w] <- rnorm(length(groups[w]), 0, sqrt(v[groups[w]]))	
s <- sample(10000, size)
simDataRisk <- (as.matrix(simDataFrame[s,]) %*% simCoef)
simDataRisk <- simDataRisk - mean(simDataRisk)
simDataSurvival <- SimSurvNonp(simDataRisk, basehaz(coxRFXFit, centered=FALSE), survival)
simDataPInt <- TestInteractions(simDataFrame[s,mainIdx], simDataSurvival, pairs)
#		simulations[[as.character(size)]][[i]] <- list(simCoef=simCoef, simDataRisk=simDataRisk, simDataSurvival=simDataSurvival, simDataPInt=simDataPInt, s=s, w=w)
#	}
#}

plot( simCoef[!mainIdx], simDataPInt$pWald, log='y', col=ifelse(simCoef[!mainIdx]==0,2,1))
idx <- simCoef[!mainIdx] == 0; qqplot(simDataPInt$pWald[idx], simDataPInt$pWald[!idx], log='xy', xlab="Coef == 0", ylab="Coef != 0")
abline(0,1)
title("QQ-plot")

which.min(simDataPInt$pWald)
simDataPInt[41,]
plot(survfit(simDataSurvival ~ Genetics.GATA2 + Genetics.EZH2, data=simDataFrame[s,] ), col=col1)
plot(survfit(simDataSurvival ~ Genetics.GATA2 + Genetics.TP53, data=simDataFrame[s,] ), col=col1)

plot( simCoef[!mainIdx], simDataPInt$pLR, log='y', col=ifelse(simCoef[!mainIdx]==0,2,1))
idx <- simCoef[!mainIdx] == 0; qqplot(simDataPInt$pLR[idx], simDataPInt$pLR[!idx], log='xy', xlab="Coef == 0", ylab="Coef != 0")
abline(0,1)
title("QQ-plot")

meanRisk <- sapply(pairs, function(p) mean(simDataRisk[simDataFrame[s,p[1]] * simDataFrame[s,p[2]]==1]) - simCoef[p[1]] - simCoef[p[2]])
plot(simCoef[!mainIdx], simDataPInt$coef, xlab="True coefficient", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(simCoef[!mainIdx][simCoef[!mainIdx]!=0], simDataPInt$coef[simCoef[!mainIdx]!=0], use="complete", method="spearman"),2)))
plot(meanRisk, simDataPInt$coef, xlab="Residual risk", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(meanRisk, simDataPInt$coef, use="complete", method="spearman"),2)))

load("temp.RData")
plot(simCoef[!mainIdx], simDataPIntAll$coef, xlab="True coefficient", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(simCoef[!mainIdx][simCoef[!mainIdx]!=0], simDataPIntAll$coef[simCoef[!mainIdx]!=0], use="complete", method="spearman"),2)))
plot(meanRisk, simDataPIntAll$coef, xlab="Residual risk", ylab="Est. coefficient")
legend("topleft", bty="n", paste("rho =", round(cor(meanRisk, simDataPIntAll$coef, use="complete", method="spearman"),2)))

plot(simCoef[!mainIdx],p.adjust(simDataPInt$pWald,"BH"), log="y", main="Marginal only", col=ifelse(simCoef[!mainIdx]==0,2,1))
abline(h=0.1)
plot(simCoef[!mainIdx],p.adjust(simDataPIntAll$pWald,"BH"), log="y", main="Incl. all mains", col=ifelse(simCoef[!mainIdx]==0,2,1))
abline(h=0.1)

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
simulate <- function(X, coxRFXFit, whichGene, treatEffect = -0.5, geneTreatEffect = .5, nSim = 1, whichMain = colnames(X)){
	treatment <- rbinom(nrow(X), 1, 0.5) ## assume randomize treatment
	Y <- cbind(X, Test = treatment)
	risk <- as.matrix(X) %*% coxRFXFit$coef + treatment * treatEffect +  X[,whichGene] * treatment * geneTreatEffect 
	simSurvival <-  SimSurvNonp(risk, basehaz(coxRFXFit, centered=FALSE), coxRFXFit$surv)
	TestInteractions(Y, simSurvival, list(c("Test", whichGene)), whichMain=whichMain)
}

#+ simGenes, cache=TRUE
set.seed(42)
effectSizes <- seq(-1,1,0.1)
simASXL1 <- lapply(c(1000,10000), function(nSim) sapply(effectSizes, function(e) sapply(1:10, function(i) simulate(simDataFrame[1:nSim,whichRFX], coxRFXFit, "Genetics.ASXL1", geneTreatEffect = e)[1,1])))
simNPM1 <- lapply(c(1000,10000), function(nSim) sapply(effectSizes, function(e) sapply(1:10, function(i) simulate(simDataFrame[1:nSim,whichRFX], coxRFXFit, "Genetics.NPM1", geneTreatEffect = e)[1,1])))

par(mfrow=c(1,2))
for(x in c("simASXL1","simNPM1")){
	tmp <- get(x)
	boxplot(tmp[[1]], log="y", names=effectSizes, border="grey", ylim=range(tmp[[2]]))
	boxplot(tmp[[2]], log="y", names=effectSizes, add=TRUE, col=NA)
	abline(h=0.05)
	abline(h=0.05 / length(pairs))
	title(xlab='Interaction effect size', ylab="P-value", main=x)
}


interactionEffects <- seq(-2,2,0.25)
nData <- c(100, 1000, 10000)
treatmentEffects <- seq(0,2,0.25)
nSim <- 10

#+ eval=FALSE
set.seed(42)
simResults <- lapply(nData, 
		function(n) mclapply(colnames(simDataFrame)[which(groups=="Genetics")], 
					function(g) lapply(treatmentEffects,
								function(t) sapply(interactionEffects, 
											function(e) sapply(1:nSim, 
														function(i) simulate(simDataFrame[1:n,whichRFX], coxRFXFit, g, treatEffect = t, geneTreatEffect = e)[1,1]))), mc.cores=14))

load("../simResults.RData")

sim <-  array(simResults$pWald[order(simResults$indeces)], dim = (c(length(nData), sum(groups=="Genetics"),length(treatmentEffects), length(interactionEffects), nSim)))
sim <- aperm(sim, 5:1)
#r <- array(unlist(simResults), dim = rev(c(length(nData), sum(groups=="Genetics"),length(treatmentEffects), length(effectSizes), nSim)))
x <- colMeans(simDataFrame[,groups=="Genetics"])
o <- order(x)
l <- (length(interactionEffects) +1)/2

plot(NA,NA, xlim=range(x), ylim=c(1e-40,1), log="y", xlab="Mutation frequency", ylab="P-value")
for(n in seq_along(nData))
	for(e in 1:l){
		lines(x=x[o],y=colMeans(sim[,2*e-1,,o,n], dims=2, na.rm=TRUE), col=brewer.pal(l,"Spectral")[e], lty=4-n, type="b", pch=n)
	}
legend("bottomleft", legend=interactionEffects[seq(1,17,2)], col = brewer.pal(l,"Spectral"), lty=1, bty='n')
legend("bottom", pch=1:3, legend=nData, bty="n")

loglowess <- function(x,y,...){
	l <- lowess(log(x),y,...)
	l$x <- exp(l$x)
	l
}


plot(NA,NA, xlim=c(1,2500), ylim=c(1e-21,1), log="xy", xlab="Mutation number x effect^2", ylab="P-value")
for(n in seq_along(nData))
	for(e in 1:l)
		points(x=x[o]*nData[n]*abs(interactionEffects[2*e-1])^2,y=(colMeans(sim[,2*e-1,,o,n], dims=2, na.rm=TRUE)), col=brewer.pal(9,"Spectral")[e], lty=4-n, pch=n)
xx <- x[o]*nData[n]*abs(interactionEffects[2*e-1])^2
lines(xx, 10^-(xx/100))


plot(NA,NA, xlim=range(x)+ 1e-3, ylim=c(0,1), xlab="Mutation frequency", ylab="Power", log="x")
for(n in seq_along(nData))
	for(e in 1:l){
		points(x[o],colMeans(sim[,2*e-1,,o,n] < 5e-2 / (length(x)), dims=2, na.rm=TRUE), col=brewer.pal(l,"Spectral")[e], lty=4-n)
		s <- smooth.spline(log(x[o]),colMeans(sim[,2*e-1,,o,n] < 5e-2 / (length(x)), dims=2, na.rm=TRUE), spar=.5)
		lines(exp(s$x),s$y, col=brewer.pal(11,"Spectral")[e], lty=4-n)
	}
legend("topleft", legend=interactionEffects[seq(1,length(interactionEffects),2)], col = brewer.pal(l,"Spectral"), lty=1, bty='n')

plot(NA,NA, xlim=c(1,2500), ylim=c(0,1), xlab="Mutation number x effect^2", ylab="Power", log="x")
for(n in seq_along(nData)[3])
	for(e in 1:l){
		points(x[o] *nData[n]* abs(interactionEffects[2*e-1])^2,colMeans(sim[,2*e-1,,o,n] < 5e-2 / (length(x)), dims=2, na.rm=TRUE), col=brewer.pal(l,"Spectral")[e], lty=4-n, pch=n)
		s <- smooth.spline(log(x[o]),colMeans(sim[,2*e-1,,o,n] < 5e-2 / (length(x)), dims=2, na.rm=TRUE), spar=.5)
		lines(exp(s$x)*nData[n]*abs(interactionEffects[2*e-1])^2,s$y, col=brewer.pal(11,"Spectral")[e], lty=4-n)
	}
legend("topleft", legend=interactionEffects[seq(1,length(interactionEffects),2)], col = brewer.pal(11,"Spectral"), lty=1, bty='n')


#' Effect size vs. mutation freq
lev <-  c(1,2,5) * rep(10^(-3:-1), each=3)
c <- cut(x,lev)
tmp <- lapply( seq_along(nData), function(n){
			sapply(levels(c), function(lev){
						sapply(1:l, function(e){
									mean(sim[,2*e-1,,which(c==lev),n, drop=FALSE] < 5e-2 / (length(x)), na.rm=TRUE)
								})
					})})

par(mfrow=c(2,2))
for(n in 1:3){
	plot(NA,NA, xlim=exp(c(-1,1)), ylim=range(log10(lev)), xlab="Hazard", ylab="Mutation frequency",  yaxt="n", log="x")
	.filled.contour(z=tmp[[n]], x = exp(interactionEffects[seq(1,length(interactionEffects),2)]), y=log10(lev[-1]),levels=seq(0.1,1.1,0.1), col=colorRampPalette(brewer.pal(9,"Reds"))(11))
	contour(z=tmp[[n]], x = exp(interactionEffects[seq(1,length(interactionEffects),2)]), y=log10(lev[-1]), levels=seq(0.1,1,0.1), col="black", add=TRUE, axes=FALSE)
	axis(2, at=log10(lev[-9]), labels=lev[-9])
	title(paste('n =' ,nData[n]))
	abline(v=c(1/seq(1,2.5,.5),seq(1,2.5,.5)), col="grey")
}

#' Interaction effect vs main effect
plot(colMeans(sim * X$nData * x[X$gene]), rep(interactionEffects, 9*56*3))

