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
clinicalData <- read.table("../data/Ulm1.5_s_Clinical.txt", sep="\t", header=TRUE)
clinicalData <- clinicalData[order(clinicalData$PDID),]
clinicalData$ERDate <- as.Date(as.character(clinicalData$ERDate), "%d-%b-%y")
clinicalData$CR_date <- as.Date(as.character(clinicalData$CR_date), "%d-%b-%y")
clinicalData$TPL_date <- as.Date(as.character(clinicalData$TPL_date), "%d-%b-%y")
clinicalData$Date_LF <- as.Date(as.character(clinicalData$Date_LF), "%d-%b-%y")
clinicalData$Recurrence_date <- as.Date(as.character(clinicalData$Recurrence_date), "%d-%b-%y")
levels(clinicalData$Study) <- c( "tr98A" ,   "tr98B" ,   "tr0704")
clinicalData$VPA[is.na(clinicalData$VPA)] <- 0
clinicalData$ATRA_arm[is.na(clinicalData$ATRA_arm)] <- 0
colnames(clinicalData) <- gsub('\\.',"",colnames(clinicalData))

#' Mutation data
#' -------------
mutationData = read.table("../data/Ulm1.5_s_Genetic.txt", sep="\t", header=TRUE)
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

trialArms <- makeInteger(clinicalData$Study)

#' ## Prepating data
dataList <-list(Genetics = data.frame(oncogenics[,colSums(oncogenics)>0]),
		Cytogenetics = clinicalData[,46:68],
		Treatment = cbind(trialArms[,2:3], ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL=clinicalData$Time_Diag_TPL < survival[,1] & !is.na(clinicalData$Time_Diag_TPL)),
		Clinical = cbind(clinicalData[, c("AOD","gender","Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH_","HB","platelet")], makeInteger(clinicalData$TypeAML)[,-1]))#,
		#MolRisk = makeInteger(clinicalData$M_Risk))
dataList$Genetics$CEBPA <- clinicalData$CEBPA > 0
dataList$Genetics$FLT3 <- NULL
dataList$Genetics$FLT3_ITD <- clinicalData$ITD > 0
dataList$Genetics$FLT3_TKD <- clinicalData$TKD > 0
dataList$Genetics$IDH2_p172 <- clinicalData$IDH2 == "172"
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Genetics = dataList$Genetics + 0
dataList$GeneGene <- makeInteractions(data.frame(dataList$Genetics), data.frame(dataList$Genetics))[,as.vector(upper.tri(matrix(0,ncol=ncol(dataList$Genetics), nrow=ncol(dataList$Genetics))))]
dataList$GeneGene <- dataList$GeneGene[,colSums(dataList$GeneGene, na.rm=TRUE)>0] 
dataList$GeneCyto <- makeInteractions(dataList$Genetics, dataList$Cytogenetics)
dataList$GeneCyto <- dataList$GeneCyto[,colSums(dataList$GeneCyto, na.rm=TRUE) > 0]
dataList$GeneTreat <- makeInteractions(dataList$Genetics, dataList$Treatment)
dataList$GeneTreat <- dataList$GeneTreat[,colSums(dataList$GeneTreat, na.rm=TRUE) > 0]
dataList$CytoTreat <- makeInteractions(dataList$Cytogenetics, dataList$Treatment)
dataList$CytoTreat <- dataList$CytoTreat[,colSums(dataList$CytoTreat, na.rm=TRUE) > 0]

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
#whichRFX <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
whichRFX <- which(colSums(dataFrame)>5 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))

## Fit Cox model
coxRFXFit <- CoxRFX(dataFrame[,whichRFX], survival, groups=groups[whichRFX], sigma0=0.1)
#coxRFXFit <- CoxRFX(dataFrame[,whichRFX], survival, groups=groups[whichRFX], sigma0=0.01) ## reassuring that sigma doesn't shrink further


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
#whichRFXTD <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
whichRFXTD <- which(colSums(dataFrame)>5 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))

coxRFXFitTD <- CoxRFX(dataFrameTD[,whichRFXTD], survivalTD, groups[whichRFXTD], sigma0=0.1)

boxplot(coef(coxRFXFitTD)/log(2) ~ factor(coxRFXFitTD$groups, levels=levels(groups)[order(coxRFXFitTD$mu)]), col=set1, horizontal=TRUE, las=2)
abline(v=0)

plot(coef(coxRFXFit), coef(coxRFXFitTD), col=set1[groups[whichRFXTD]]) # Note the sign change for TPL..
abline(0,1)

#' ## Performance metrics
#+ fig.width=14, fig.height=5
for(g in levels(groups)){
	par(mar=c(7,4,2,0))
	x <- sort(coxRFXFitTD$coefficients[groups[whichRFXTD]==g])
	v <- diag(coxRFXFitTD$var)[groups[whichRFXTD]==g][ order(coxRFXFitTD$coefficients[groups[whichRFXTD]==g])]
	plot(x, las=2, pch=16, xaxt="n", main=g, cex = .5+pmin(3,-log10(pchisq((x-coxRFXFitTD$mu[g])^2/v, 1, lower.tail=FALSE))))
	segments(1:length(x),x - 2*sqrt(v) ,1:length(x), x+2*sqrt(v))
	axis(side=1, at = 1:length(x), labels = names(x), las=2, cex.axis=.6)
	abline(v=1:length(x), lty=3, col="grey")
	abline(h=coxRFXFitTD$mu[g])
	abline(h = coxRFXFitTD$mu[g] + c(-1,1)*sqrt(coxRFXFitTD$sigma2[g]), lty=2)
}

#' ### Partial risk contributions
partRiskTD <- PartialRisk(coxRFXFitTD)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponents <- diag(cov(partRiskTD, use="complete"))
varianceComponents
pie(abs(varianceComponents), col=brewer.pal(8,"Set1"))
title("Risk contributions")

#' #### Stars
#+ fig.width=12, fig.height=12
library(HilbertVis)
nStars <- 32
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
s <- sample(nrow(partRiskTD),nStars^2) #1:(nStars^2)
h <- hclust(dist(partRiskTD[s,]))
#stars(exp((partRisk[s,][h$order,] - rep(colMeans(partRisk), each=nStars^2)))/2, scale=FALSE, locations=locations, key.loc=c(0,-3), col.lines=rep(1,(nStars^2)), col.stars = brewer.pal(11,'RdBu')[cut(survivalTD[s,2][h$order], quantile(survivalTD[,2], seq(0,1,0.1), na.rm=TRUE))])
x <- partRiskTD - rep(colMeans(partRiskTD), each=nrow(partRiskTD))
x <- x/(2*sd(x)) + 1 
#x <- exp(x)
stars(x[s,][h$order,]/2, scale=FALSE, locations=locations, key.loc=c(0,-3), col.lines=rep(1,(nStars^2)), col.stars = (brewer.pal(11,'RdBu'))[cut(survivalTD[s,2][h$order], quantile(survivalTD[,2], seq(0,1,0.1), na.rm=TRUE))])
symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=1)

#' ### A few performance measures
#' Harrel's C
#library(Hmisc)
totalRiskTD <- rowSums(partRiskTD)
survConcordance( survivalTD~totalRiskTD)
predictiveRiskTD <- rowSums(partRiskTD[,-which(colnames(partRiskTD) %in% c("Treatment","GeneTreat","CytoTreat"))])
survConcordance( survivalTD~predictiveRiskTD)

#' AUC
library(survAUC)
s <- survival
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTD[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskTD[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' Risk quartiles
predictiveRiskGroups <-  cut(predictiveRiskTD[!is.na(predictiveRiskTD)], quantile(predictiveRiskTD, seq(0,1,0.25)), labels = c("very low","low","high","very high"))
survConcordance( survivalTD~as.numeric(predictiveRiskGroups))
plot(survfit(survivalTD[!is.na(predictiveRiskTD)] ~ predictiveRiskGroups), col=brewer.pal(4,"Spectral")[4:1], mark=16)
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
predictiveRiskTest <- rowSums(partialRiskTest[,-which(colnames(partRiskTD) %in% c("Treatment","GeneTreat","CytoTreat"))])
survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)

#barplot(c(CGP=survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)$concordance, MolecularRisk = survConcordance(survivalTD[!testIxTD] ~ c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[tplSplit][!testIxTD]])[[1]]))

s <- survival[!testIx]
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTest[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskTest[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' ### Recursive partitioning
r <- rpart(survival ~ ., data=dataFrame)
plot(r)
text(r)
survConcordance(na.omit(survival)~predict(r))


#' Stability selection
#' --------------------
library(parallel)
library(glmnet)
set.seed(42)
source("/Users/mg14/Git/Projects/CoxHD/CoxHD/R/stacoxph.R")
source("/Users/mg14/Git/Projects/CoxHD/CoxHD/R/r_concave_tail.R")

#' Fit model
#+ stabCox, cache=TRUE
stabCox <- StabCox(dataFrame[!is.na(survival),], survival[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

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
partRisk <- PartialRisk(coxRFXFit)
totalRisk <- rowSums(partRisk)
simSurvival <- simSurvNonp(totalRisk, basehaz(coxRFXFit, centered=FALSE), survival, span=.1)
plot(survfit(simSurvival ~ 1))
lines(survfit(survival ~ 1), col="red")
survConcordance(simSurvival[testIx] ~ totalRisk[testIx])
survConcordance(simSurvival[testIx][-na.action(coxFit)] ~ predict(coxFit))

simPredictedRisk <- PredictRiskMissing(simCoxRFXFit)
plot(totalRisk, simPredictedRisk[,1], xlab="True risk (simulated)",ylab="Estimated risk")
segments(totalRisk, simPredictedRisk[,1] + sqrt(simPredictedRisk[,2]), totalRisk, simPredictedRisk[,1] - sqrt(simPredictedRisk[,2]), col="#00000022")
abline(0,1, col="red")

par(mfrow=c(3,3))
p <- PartialRisk(simCoxRFXFit)
v <- PartialRiskVar(simCoxRFXFit)
for(i in 1:8){
	plot(partRisk[,i], p[,i], main=colnames(partRisk)[i], xlab="True risk component (simulated)",ylab="Estimated risk component")
	segments(partRisk[,i], p[,i] + sqrt(v[,i]), partRisk[,i], p[,i] - sqrt(v[,i]), col="#00000022")
	abline(0,1, col="red")
}

#+ simCoxRFXFit, cache=TRUE
simCoxRFXFit <- CoxRFX(dataFrame[,whichRFX], simSurvival, groups=groups[whichRFX], sigma0 = 0.1)

boxplot(coef(simCoxRFXFit)~groups[whichRFX], border=set1, horizontal=TRUE, las=2)
abline(v=0)

plot(coef(simCoxRFXFit), coef(coxRFXFit), col = set1[groups[whichRFX]])

#+ simStabCox, cache=TRUE
set.seed(42)
s <- simSurvival
s[s[,1]==0,1] <- .5 #Machine$double.eps
simStabCox <- StabCox(dataFrame[!is.na(survival),], s[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

#' Plots
plot(simStabCox)
plot(simStabCox$Pi, stabCox$Pi)

#' Selected variables
selected <- which(simStabCox$Pi > .8)
selected
simCoxFit <- coxph(simSurvival[testIx] ~ ., data=dataFrame[testIx,][,names(selected)])
summary(simCoxFit)

#' ### Stability selection with interactions
#' Only include interaction terms when main effects are present. 
#' On the other hand, the exclusive selection of a product term could mean that the most likely explanation is that the two main terms are zero and only the interaction is non-zero..
StabCoxInteractions

#+ stabCoxInt, cache=TRUE
set.seed(42)
simStabCoxInt <- StabCoxInteractions(dataFrame[!is.na(survival),groups %in% c("Genetics","Clinical","Treatment","Cytogenetics")], simSurvival[!is.na(survival)], bootstrap.samples=50)

s <- survival
s[s[,1]==0,1] <- .5 #Machine$double.eps
stabCoxInt <- StabCoxInteractions(dataFrame[!is.na(survival),groups %in% c("Genetics","Clinical","Treatment","Cytogenetics")], s[!is.na(survival)], bootstrap.samples=50)

s <- strsplit(gsub("(Genetics|Cytogenetics|Treatment|Clinical).","", names(stabCoxInt$Pi)), "\\.")
m <- sapply(s, function(x) grep( paste(paste(x, collapse="."),paste(rev(x), collapse="."), sep="|"), colnames(dataFrame))[1])

#' Plots
plot(stabCoxInt)
plot(stabCoxInt$Pi, stabCox$Pi[m])

plot(simStabCoxInt)
plot(stabCoxInt$Pi, simStabCoxInt$Pi[names(stabCoxInt$Pi)])



#' Selected variables
selected <- which(simStabCoxInt$Pi > .8)
selected
simCoxFit <- coxph(simSurvival[testIx] ~ ., data=simStabCoxInt[testIx,][,names(selected)])
summary(simCoxFit)

#' ### Survival and covariates
#' Simulate data 
#+ simData, cache=TRUE
set.seed(42)
simDataNonp
data <- data.frame(Clinical=dataList$Clinical, Genetics=dataList$Genetics, Cytogenetics=dataList$Cytogenetics, Treatment=dataList$Treatment)
#simData <- simDataNonp(data, nData = 10000, m=5)
#simData <- simDataNonp(data, nData = nrow(data))
#save(file="simData.RData", simData)
load(file="simData.RData")

g <- sub("\\..+","", colnames(data))
simDataFrame <- data.frame(simData,
		GeneGene=makeInteractions(simData[,g=="Genetics"], simData[,g=="Genetics"])[,as.vector(upper.tri(matrix(0,ncol=sum(g=="Genetics"), nrow=sum(g=="Genetics"))))],
		GeneCyto=makeInteractions(simData[,g=="Genetics"], simData[,g=="Cytogenetics"]),
		GeneTreat=makeInteractions(simData[,g=="Genetics"], simData[,g=="Treatment"]),
		CytoTreat=makeInteractions(simData[,g=="Cytogenetics"], simData[,g=="Treatment"])
)
colnames(simDataFrame) <- gsub("\\.(Genetics|Treatment|Cytogenetics)","", colnames(simDataFrame))
simDataFrame <- simDataFrame[,colnames(dataFrame)]
simDataFrame <- as.data.frame(sapply(simDataFrame, poorMansImpute))
simDataFrame <- StandardizeMagnitude(simDataFrame)

#simDataFrameTD <- simDataFrame[tplSplit,]
#simDataFrameTD[which(tplIndex), grep("TPL", colnames(simDataFrameTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

#' Coefficients
set.seed(42)
v <- coxRFXFit$sigma2 * coxRFXFit$df[-8] / as.numeric(table(groups[whichRFXTD])-1)
simCoef <- numeric(ncol(simDataFrame))
names(simCoef) <- colnames(simDataFrame)
simCoef[whichRFXTD] <- rnorm(length(groups[whichRFXTD]), coxRFXFit$mu[groups[whichRFXTD]], sqrt(v[groups[whichRFXTD]]))

#' Survival
set.seed(42)
simDataRisk <- (as.matrix(simDataFrame) %*% simCoef)[,1] 
simDataRisk <- simDataRisk - mean(simDataRisk)
simDataSurvival <- simSurvNonp(simDataRisk, basehaz(coxRFXFit, centered=FALSE), survival, span=.1)
plot(survfit(simDataSurvival ~ Genetics.EP300, data=simDataFrame))

#' ### RFX model
#+ simDataCoxRFXFit, cache=TRUE
simDataCoxRFXFit <- CoxRFX(simDataFrame[, names(coef(coxRFXFit))], simDataSurvival, groups = groups[whichRFXTD], sigma0 = 0.1)
plot(simCoef[whichRFXTD],coef(simDataCoxRFXFit), col=set1[groups[whichRFXTD]])
plot(sapply(split(simCoef[whichRFXTD], groups[whichRFXTD]), var),simDataCoxRFXFit$sigma2, col=set1, pch=19)
plot(sapply(split(simCoef[whichRFXTD], groups[whichRFXTD]), mean),simDataCoxRFXFit$mu, col=set1, pch=19)


#' ### Stability selection
#+ simDataStabCox, cache=TRUE
set.seed(42)
s <- simDataSurvival
s[s[,1]==0,1] <- .Machine$double.eps
simDataStabCox <- StabCox(simDataFrame[!is.na(survival),], s[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)
#save(simDataStabCox, file="simDataStabCox")

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


#' ### Test interactions
X <- cbind(dataList$Genetics, dataList$Cytogenetics,dataList$Treatment,dataList$Clinical)
scope <- sub("(GeneGene|GeneTreat|CytoTreat|GeneCyto).","",colnames(dataFrame)[grep("(GeneGene|GeneTreat|CytoTreat|GeneCyto)", names(dataFrame))])
pInt <- testInteractions(X, survival, scope)

simPInt <- testInteractions(simDataFrame[groups %in% c("Clinical","Genetics","Treatment","Cytogenetics")], simDataSurvival, scope)

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
		
#' TODO's
#' ------
#' * Uncertainty of RFX model
#'   - Covariance of estimates
#'   - Covariance of covariates
#' * MI
#' * Natural limit to C?
#' * Missing data
#' * Simulation study
#'   - How to quantify result