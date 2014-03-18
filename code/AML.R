#' AML data analysis
#' =================

#+ echo=FALSE
options(width=120)
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' Libraries
library(RColorBrewer)
set1 <- brewer.pal(8, "Set1")

#' Clinical data
#' -------------
#' Loading
clinicalData <- read.table("../data/Clinical.txt", sep="\t", header=TRUE)
clinicalData <- clinicalData[order(clinicalData$PDID),]
clinicalData$ERDate <- as.Date(as.character(clinicalData$ERDate), "%d-%b-%y")
clinicalData$CR_date <- as.Date(as.character(clinicalData$CR_date), "%d-%b-%y")
clinicalData$TPL_date <- as.Date(as.character(clinicalData$TPL_date), "%d-%b-%y")
clinicalData$Date_LF <- as.Date(as.character(clinicalData$Date_LF), "%d-%b-%y")
clinicalData$Recurrence_date <- as.Date(as.character(clinicalData$Recurrence_date), "%d-%b-%y")
levels(clinicalData$Study) <- c( "tr98A" ,   "tr98B" ,   "tr0704")


#' Mutation data
#' -------------
mutationData = read.table("../data/Genetic.txt", sep="\t", header=TRUE)
mutationData$SAMPLE_NAME <- factor(as.character(mutationData$SAMPLE_NAME), levels = levels(clinicalData$PDID)) ## Refactor
oncogenics <- (table(mutationData[mutationData$Result %in% c("ONCOGENIC","POSSIBLE"),c("SAMPLE_NAME","GENE")]) > 0)+0

#' Compute gene-gene interactions
#+ interactions, cache=TRUE
interactions <- sapply(1:ncol(oncogenics), function(i) sapply(1:ncol(oncogenics), function(j) {f<- try(fisher.test(oncogenics[,i], oncogenics[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
odds <- sapply(1:ncol(oncogenics), function(i) sapply(1:ncol(oncogenics), function(j) {f<- try(fisher.test(oncogenics[,i] + .5, oncogenics[,j] +.5), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
diag(interactions) <- 0
diag(odds) <- 1
colnames(odds) <- rownames(odds) <- colnames(interactions) <- rownames(interactions) <- colnames(oncogenics)
odds[10^-abs(interactions) > 0.05] = 1
odds[odds<1e-3] = 1e-4
odds[odds>1e3] = 1e4
logodds=log10(odds)

#' Heatmap
#+ dev=c("png","pdf"), fig.width=8, fig.height=8, fig.args=list( pointsize = 8)
#pdf(paste(Sys.Date(),"-Interactions.pdf", sep=""), 4,4, pointsize = 8) ## HEATMAP OF ALL
par(bty="n", mgp = c(2,.5,0), mar=c(3,3,2,2)+.1, las=2, tcl=-.33)
ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
h = hclust(dist(interactions[ix,ix]))
o = h$order #c(h$order,(length(h$order) +1):ncol(interactions))
#image(x=1:ncol(interactions), y=1:nrow(interactions), interactions[o,o], col=brewer.pal(11,"BrBG"), breaks = c(-100, seq(-9,9,2),  100), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(interactions)+3), ylim=c(0, ncol(interactions)+3))
image(x=1:ncol(interactions), y=1:nrow(interactions), log10(odds[o,o]), col=brewer.pal(9,"BrBG"), breaks = c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(interactions)+3), ylim=c(0, ncol(interactions)+3))
mtext(side=1, at=1:ncol(interactions), colnames(interactions)[o], cex=.5, font=3)
mtext(side=2, at=1:ncol(interactions), colnames(interactions)[o], cex=.5, font=3)
abline(h = length(h$order)+.5, col="white", lwd=1)
abline(v = length(h$order)+.5, col="white", lwd=1)
abline(h=0:ncol(interactions)+.5, col="white", lwd=.5)
abline(v=0:ncol(interactions)+.5, col="white", lwd=.5)
w = arrayInd(which(p.adjust(10^-abs(interactions[o,o]), method="BH") < .1), rep(nrow(interactions),2))
points(w, pch=".", col="white")
w = arrayInd(which(p.adjust(10^-abs(interactions[o,o])) < .05), rep(nrow(interactions),2))
points(w, pch="*", col="white")
image(y = 1:8, x=rep(ncol(interactions),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
#axis(side = 4, at = 8:12 + .5, cex.axis=.5, tcl=-.15, label=10^-abs(seq(1,9,2)), las=1, lwd=.5)
#axis(side = 4, at = 1:5 + .5, cex.axis=.5, tcl=-.15, label=10^-abs(seq(-9,-1,2)), las=1, lwd=.5)
axis(side = 4, at = seq(1,7) + .5, cex.axis=.5, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
#lines(rep(ncol(interactions),2)+c(1,4), c(6,6)+.5, col="white")
#mtext(side=4, at=-1,  "P odds < 1", cex=.5, line=-.5)
#mtext(side=4, at=15,  "P odds > 1", cex=.5, line=-.5)
points(x=rep(ncol(interactions),2)+2.5, y=10:11, pch=c("*","."))
image(x=rep(ncol(interactions),2)+c(2,3), y=(12:13) -0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
mtext(side=4, at=12:10, c("P > 0.05", "BH < 0.1", "Bf. < 0.05"), cex=.5, line=0.2)
#dev.off()

#' Survival analyses
#' -----------------
library(survival)
all(rownames(oncogenics)==clinicalData$PDID)
survival <- Surv(clinicalData$efs, clinicalData$EFSSTAT)
#survival <- Surv(clinicalData$OS, clinicalData$Status)


#' ### Proportional hazards test
#+ fig.width=14, fig.height=5
coxModel <- coxph(survival ~ gender + AOD + Study, data=clinicalData) # Most basic
phTest <- cox.zph(coxModel)
phTest
par(mfrow=c(1,4), bty="n")
for(i in 1:4) plot(phTest[i])

source("~/Git/Projects/CoxHD/CoxHD/R/ecoxph.R")

#' #### A few functions
makeInteractions <- function(X,Y){
	do.call(cbind, lapply(X, `*`, Y))
}

makeInteger <- function(F){
	res <- as.data.frame(lapply(levels(F), `==`, F))
	colnames(res) <- levels(F)
	res + 0
}

trialArms <- data.frame(makeInteger(clinicalData$Study), 
		tr98B_ATRA = clinicalData$Study == "tr98B" & clinicalData$Randomisation_Arm == 2,
		tr0704_ATRA = clinicalData$Study == "tr0704" & clinicalData$Randomisation_Arm %in% c(3,4),
		tr0704_VPA = clinicalData$Study == "tr0704" & clinicalData$Randomisation_Arm %in% c(3,5),
		check.names=FALSE
) + 0

#' ## Prepating data
dataList <-list(Genetics = data.frame(oncogenics[,colSums(oncogenics)>0]),
		Cytogenetics = clinicalData[,43:78],
		Treatment = cbind(trialArms[,2:3], ATRA = trialArms$tr0704_ATRA | trialArms$tr98B_ATRA, VPA=trialArms$tr0704_VPA, TPL=clinicalData$Time_Diag_TPL < survival[,1] & !is.na(clinicalData$Time_Diag_TPL)),
		Clinical = cbind(clinicalData[, c("AOD","gender","Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH_","HB","platelet")], makeInteger(clinicalData$TypeAML)[,-1]))#,
		#MolRisk = makeInteger(clinicalData$M_Risk))
dataList$Genetics$CEBPA <- clinicalData$CEBPA > 0
dataList$Genetics$FLT3_ITD <- clinicalData$FLT3_ITD > 0
dataList$Genetics$FLT3_TKD <- clinicalData$FLT3_TKD > 0
dataList$Genetics$IDH2_p172 <- clinicalData$IDH2 == "172"
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Genetics = dataList$Genetics + 0
dataList$GeneGene <- makeInteractions(data.frame(oncogenics), data.frame(oncogenics))[,as.vector(upper.tri(matrix(0,ncol=ncol(oncogenics), nrow=ncol(oncogenics))))] ## TODO: Fixme
dataList$GeneGene <- dataList$GeneGene[,colSums(dataList$GeneGene)>0] 
dataList$GeneCyto <- makeInteractions(dataList$Genetics, dataList$Cytogenetics)
dataList$GeneCyto <- dataList$GeneCyto[,colSums(dataList$GeneCyto, na.rm=TRUE) > 0]
dataList$GeneTreatment <- makeInteractions(dataList$Genetics, dataList$Treatment)
dataList$GeneTreatment <- dataList$GeneTreatment[,colSums(dataList$GeneTreatment, na.rm=TRUE) > 0]

dataFrame <- do.call(cbind,dataList)
dataFrame <- StandardizeMagnitude(dataFrame)

groups <- factor(unlist(sapply(names(dataList), function(x) rep(x, ncol(dataList[[x]])))))

## Poor man's imputation
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
dataFrame <- as.data.frame(sapply(dataFrame, poorMansImpute))

#' #### Basic KM plot
#+ eval=FALSE, echo=FALSE
set1 <- brewer.pal(9, "Set1")
plot(survfit(survival ~ Genetics.TP53 + Treatment.TPL, data=dataFrame), col=rep(set1[1:2],each=2), lty=1:2, mark=16)


#' ### Fit static model
#+ coxRFXFit, warning=FALSE, cache=TRUE
## Pre-select based on univariate test
univP <- sapply(dataFrame, function(x){summary(coxph(survival~x))$logtest[3]})
whichRFX <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
## Fit Cox model
coxRFXFit <- CoxRFX(dataFrame[,whichRFX], survival, groups=groups[whichRFX], sigma0=0.1)

#' Coefficients
boxplot(coxRFXFit$coefficients ~ groups[whichRFX], las=2, border=set1, horizontal=TRUE)
abline(v=0)

#' ### Time-dependent survival effects (TPL)
#' Construct a time-dependent Surv() object by splitting TPL patients into pre and post
t <- clinicalData$Time_Diag_TPL
t[is.na(t)] <- Inf
e <- clinicalData$efs
tplIndex <-  t < e
survivalTD <-  Surv(time = rep(0, nrow(clinicalData)), time2=pmin(e, t), event=ifelse(tplIndex, 0, clinicalData$EFSSTAT) )
survivalTD <- rbind(survivalTD, 
		Surv(time=t[which(tplIndex)],
				time2=e[which(tplIndex)], 
				event=clinicalData$EFSSTAT[which(tplIndex)])
)
survivalTD = Surv(survivalTD[,1],survivalTD[,2],survivalTD[,3])
rm(e,t)
tplSplit <- c(1:nrow(clinicalData), which(tplIndex))

#' Data frame
dataFrameTD <- dataFrame[tplSplit,]
dataFrameTD[which(tplIndex), grep("TPL", colnames(dataFrameTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

patientTD <- c(clinicalData$PDID, clinicalData$PDID[which(tplIndex)])

#' Fit TD Cox RFX model
#+ coxRFXFitTD, cache=TRUE
whichRFX <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
coxRFXFitTD <- CoxRFX(dataFrameTD[,whichRFX], survivalTD, groups[whichRFX], sigma0=0.1)

plot(coef(coxRFXFit), coef(coxRFXFitTD), col=set1[groups[whichRFX]]) # Note the sign change for TPL..
abline(0,1)

boxplot(coef(coxRFXFitTD)~groups[whichRFX], border=set1, horizontal=TRUE, las=2)
abline(v=0)


#' ## Performance metrics
#+ fig.width=14, fig.height=5
for(g in levels(groups)){
	par(mar=c(7,4,2,0))
	x <- sort(coxRFXFitTD$coefficients[groups[whichRFX]==g])
	v <- diag(coxRFXFitTD$var)[groups[whichRFX]==g][ order(coxRFXFitTD$coefficients[groups[whichRFX]==g])]
	plot(x, las=2, pch=16, xaxt="n", main=g, cex = .5+pmin(3,-log10(pchisq((x-coxRFXFitTD$mu[g])^2/v, 1, lower.tail=FALSE))))
	segments(1:length(x),x - 2*sqrt(v) ,1:length(x), x+2*sqrt(v))
	axis(side=1, at = 1:length(x), labels = names(x), las=2, cex.axis=.6)
	abline(v=1:length(x), lty=3, col="grey")
	abline(h=coxRFXFitTD$mu[g])
	abline(h = coxRFXFitTD$mu[g] + c(-1,1)*sqrt(coxRFXFitTD$sigma2[g]), lty=2)
}

#' ### Partial risk contributions
partRisk <- PartialRisk(coxRFXFitTD)
varianceComponents <- rowSums(cov(partRisk, use="complete"))
varianceComponents
pie(abs(varianceComponents), col=brewer.pal(8,"Set1"))
title("Risk contributions")

#' #### Stars
#+ fig.width=12, fig.height=12
library(HilbertVis)
nStars <- 32
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
s <- sample(nrow(partRisk),nStars^2) #1:(nStars^2)
h <- hclust(dist(partRisk[s,]))
stars(exp(-(partRisk[s,][h$order,] - rep(colMeans(partRisk), each=nStars^2)))/2, scale=FALSE, locations=locations, key.loc=c(-1.5,-1.5), col.lines=rep(1,(nStars^2)), col.stars = brewer.pal(11,'RdBu')[cut(survivalTD[s,2][h$order], quantile(survivalTD[,2], seq(0,1,0.1), na.rm=TRUE))])
symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=3)

#' ### A few performance measures
#' Harrel's C
#library(Hmisc)
totalRisk <- rowSums(partRisk)
survConcordance( survivalTD~totalRisk)
predictiveRisk <- rowSums(partRisk[,-which(colnames(partRisk) %in% c("Treatment","GT"))])
survConcordance( survivalTD~predictiveRisk)

#' AUC
library(survAUC)
s <- survival
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRisk[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRisk[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' Risk quartiles
predictiveRiskGroups <-  cut(predictiveRisk[!is.na(predictiveRisk)], quantile(predictiveRisk, seq(0,1,0.25)), labels = c("very low","low","high","very high"))
survConcordance( survivalTD~as.numeric(predictiveRiskGroups))
plot(survfit(survivalTD[!is.na(predictiveRisk)] ~ predictiveRiskGroups), col=brewer.pal(4,"Spectral")[4:1], mark=16)
table(clinicalData$M_Risk, predictiveRiskGroups[1:nrow(clinicalData)])[c(2,3,4,1),]

#' Partial values of Harrel's C
c <- PartialC(coxRFXFitTD)
b <- barplot(c[1,], col=set1, las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' Cross-validation
#' ----------------
set.seed(42)
testIx <- sample(c(TRUE,FALSE), nrow(dataFrame), replace=TRUE, prob=c(0.66,0.34))
testIxTD <- testIx[tplSplit]

#' Pre-select based on univariate test
#+ warning=FALSE
p <- sapply(dataFrame, function(x){summary(coxph(survivalTD[testIxTD]~x[testIxTD]))$logtest[3]})
whichTrain <- which(p.adjust(p,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))

#' Fit Cox model
#+ coxRFXFitTDTrain, cache=TRUE
coxRFXFitTDTrain <- CoxRFX(dataFrameTD[testIxTD,whichTrain], survivalTD[testIxTD], groups=groups[whichTrain], sigma0 = 0.1)

#' Partial contributions
partialRiskTest <- PartialRisk(coxRFXFitTDTrain, newX=dataFrameTD[!testIxTD, whichTrain])
c <- PartialC(coxRFXFitTDTrain, newX = dataFrameTD[!testIxTD, whichTrain], newSurv = survivalTD[!testIxTD])
b <- barplot(c[1,], col=set1, las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' Overall
totalRiskTest <- rowSums(partialRiskTest)
survConcordance(survivalTD[!testIxTD]~totalRiskTest )

#' Compared to molecular risk
predictiveRiskTest <- rowSums(partialRiskTest[,-which(colnames(partRisk) %in% c("Treatment","GeneTreatment"))])
survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)

#barplot(c(CGP=survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)$concordance, MolecularRisk = survConcordance(survivalTD[!testIxTD] ~ c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[tplSplit][!testIxTD]])[[1]]))

s <- survivalTD[!testIxTD]
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTest[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskTest[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' Stability selection
#' --------------------
library(parallel)
library(glmnet)
set.seed(42)
source("/Users/mg14/Git/Projects/CoxHD/CoxHD/R/stacoxph.R")
source("/Users/mg14/Git/Projects/CoxHD/CoxHD/R/r_concave_tail.R")

#' Fit model
#+ stabCox, cache=TRUE
stabCox <- stacoxph(dataFrame[!is.na(survival),], survival[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

plot(stabCox)

selected <- which(stabCox$Pi > .8)
selected
coxFit <- coxph(survival[testIx] ~ ., data=dataFrame[testIx,][,names(selected)])
summary(coxFit)

totalRiskStabTest <- as.matrix(dataFrame[!testIx,][,names(selected)]) %*% coxFit$coefficients
survConcordance(survival[!testIx] ~ totalRiskStabTest)

plot(survfit(survival ~ FLT3_ITD + DNMT3A, data = dataList$Genetics), col=c("grey",brewer.pal(3,"Set1")))
legend("topright", bty="n",  lty=1, col=c("grey",brewer.pal(3,"Set1")), c("WT/WT","WT/DNMT3A-","FLT3-/WT","FLT3-/DNMT3A-"))

#' Probably effect of time-dep TPL:
plot(survfit(survival ~ `Genetics.TP53` + `Treatment.TPL`, data = dataFrame), col=c("grey",brewer.pal(3,"Set1")))
legend("topright", bty="n",  lty=1, col=c("grey",brewer.pal(3,"Set1")), c("WT/TPL-","WT/TPL+","TP53-/TPL-","TP53-/TPL+"))

#' Simulation study
#' ----------------
source("/Users/mg14/Git/Projects/CoxHD/CoxHD/R/functions.R")
simSurvNonp

#' ### Survival only
set.seed(42)
simSurvival <- simSurvNonp(rowSums(PartialRisk(coxRFXFit)), basehaz(coxRFXFit, centered=FALSE), survival, span=.1)
plot(survfit(simSurvival ~ 1))
survConcordance(simSurvival[testIx] ~ rowSums(PartialRisk(coxRFXFit))[testIx])
survConcordance(simSurvival[testIx][-na.action(coxFit)] ~ predict(coxFit))

#+ simCoxRFXFit, cache=TRUE
simCoxRFXFit <- CoxRFX(dataFrame[,whichRFX], simSurvival, groups=groups[whichRFX], sigma0 = 0.1)

boxplot(coef(simCoxRFXFit)~groups[whichRFX], border=set1, horizontal=TRUE, las=2)
abline(v=0)

#+ simStabCox, cache=TRUE
s <- simSurvival
s[s[,1]==0,1] <- .Machine$double.eps
simStabCox <- stacoxph(dataFrame[!is.na(survival),], s[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

#' Plots
plot(simStabCox)
plot(simStabCox$Pi, stabCox$Pi)

#' Selected variables
selected <- which(simStabCox$Pi > .8)
selected
simCoxFit <- coxph(simSurvival[testIx] ~ ., data=dataFrame[testIx,][,names(selected)])
summary(simCoxFit)

#' ### Survival and covariates
#' Simulate data 
#+ simData, cache=TRUE
set.seed(42)
simDataNonp
data <- data.frame(Clinical=dataList$Clinical, Genetics=dataList$Genetics, Cytogenetics=dataList$Cytogenetics, Treatment=dataList$Treatment)
simData <- simDataNonp(data, nData = nrow(data))
g <- sub("\\..+","", colnames(data))
simDataFrame <- data.frame(simData,
		GeneGene=makeInteractions(simData[,g=="Genetics"], simData[,g=="Genetics"])[,as.vector(upper.tri(matrix(0,ncol=sum(g=="Genetics"), nrow=sum(g=="Genetics"))))],
		GeneCyto=makeInteractions(simData[,g=="Genetics"], simData[,g=="Cytogenetics"]),
		GeneTreatment=makeInteractions(simData[,g=="Genetics"], simData[,g=="Treatment"])
)
colnames(simDataFrame) <- gsub("\\.(Genetics|Treatment|Cytogenetics)","", colnames(simDataFrame))
simDataFrame <- simDataFrame[,colnames(dataFrame)]
simDataFrame <- as.data.frame(sapply(simDataFrame, poorMansImpute))
simDataFrame <- StandardizeMagnitude(simDataFrame)

#simDataFrameTD <- simDataFrame[tplSplit,]
#simDataFrameTD[which(tplIndex), grep("TPL", colnames(simDataFrameTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

#' Coefficients
set.seed(42)
v <- coxRFXFit$sigma2 * coxRFXFit$df[-8] / (table(groups[whichRFX])-1)
simCoef <- numeric(ncol(simDataFrame))
names(simCoef) <- colnames(simDataFrame)
simCoef[whichRFX] <- rnorm(length(groups[whichRFX]), coxRFXFit$mu[groups[whichRFX]], sqrt(v[groups[whichRFX]]))

#' Survival
set.seed(42)
simDataRisk <- (as.matrix(simDataFrame) %*% simCoef)[,1] 
simDataSurvival <- simSurvNonp(simDataRisk, basehaz(coxRFXFit, centered=FALSE), survival, span=.1)
plot(survfit(simDataSurvival ~ Genetics.TP53, data=simDataFrame))

#' ### RFX model
#+ simDataCoxRFXFit, cache=TRUE
simDataCoxRFXFit <- CoxRFX(simDataFrame[, names(coef(coxRFXFit))], simDataSurvival, groups = groups[whichRFX], sigma0 = 0.1)
plot(coef(simDataCoxRFXFit), simCoef[whichRFX], col=set1[groups[whichRFX]])
plot(simDataCoxRFXFit$sigma2, v, col=set1, pch=19)


#' ### Stability selection
#+ simDataStabCox, cache=TRUE
set.seed(42)
s <- simDataSurvival
s[s[,1]==0,1] <- .Machine$double.eps
simDataStabCox <- stacoxph(simDataFrame[!is.na(survival),], s[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

#' Plots
plot(simDataStabCox)
plot(simDataStabCox$Pi, stabCox$Pi)
plot(simCoef, simDataStabCox$Pi)
plot(colMeans(simDataFrame), simDataStabCox$Pi, cex=sqrt(abs(simCoef))*2+.1, log="x")


#' Selected variables
selected <- which(simDataStabCox$Pi > .8)
selected
simDataCoxFit <- coxph(simDataSurvival[testIx] ~ ., data=simDataFrame[testIx,][,names(selected)])
summary(simDataCoxFit)


#' Germline polymorphisms
#' ----------------------
#+ eval=FALSE
load("/Volumes/mg14/subclones/snps.RData")
getSNPs <- function(samples){
	isCoding = function(subs){
		var = sub("(.*?)?(p\\.|$)","",info(subs)$VW, perl=TRUE)
		var[is.na(var)]=""
		x = sub("(^.)(.+)","\\1",var)
		y = sub("(.+)(.$)","\\2",var)
		return(x != y & x != "")
	}
	d <-  dir("/Volumes/nst_links/live/845/", pattern="WGA*")
	r <- t(sapply(samples, function(s){
						i<<-i+1
						cat(ifelse(i%%100 == 0, "\n","."))
				stem <<- grep(s, d, value=TRUE)[1]
				if(is.na(stem))
					return(rep(NA, ncol(oncogenics)))
				vcf <- readVcf(paste("/Volumes/nst_links/live/845/",stem,"/",stem,".cave.annot.vcf.gz", sep=""), "GRCh37")
				#genes <- sub("\\|.+","",info(vcf)$VD[info(vcf)$SNP & isCoding(vcf)])
				genes <- sub("\\|.+","",info(vcf)$VD[vcf %over% ensemblVariation & isCoding(vcf)])
				colnames(oncogenics) %in% genes
			}))
	colnames(r) <- colnames(oncogenics)
	r
}

#+ eval=FALSE
snp <- getSNPs(rownames(oncogenics))
save(snp, file="snp.RData")

#+ eval=FALSE
dataFrameTD <-  dir("/Volumes/nst_links/live/845/", pattern="WGA*")
allCaveOut <- sapply(rownames(oncogenics), function(s){
			i<<-i+1
			cat(ifelse(i%%100 == 0, "\n","."))
			stem <<- grep(s, dataFrameTD, value=TRUE)[1]
			if(is.na(stem))
				return(NA)
			readVcf(paste("/Volumes/nst_links/live/845/",stem,"/",stem,".cave.annot.vcf.gz", sep=""), "GRCh37")
		})
save(allCaveOut, file="allCaveOut.RData")
v <- snps[!is.na(snps$MF)]
s <- matrix(0, ncol = ncol(oncogenics), nrow=nrow(oncogenics), dimnames = dimnames(oncogenics))
i <- 1
for(vcf in allCaveOut){
			if(is.na(vcf))
				next
			genes <- sub("\\|.+","",info(vcf)$VD[vcf %over% v & isCoding(vcf)])
			s[i, colnames(s) %in% genes] <- 1
			i <- i+1
		}
		
#' Todo's
#' ------
#' * TODO: MI
#' * TODO: Latest data
#' * TODO: Simulation study