#' AML data analysis
#' =================

#+ Preliminaries, echo=FALSE
options(width=120)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('png','pdf'), fig.ext=c('png','pdf'), fig.width=4, fig.height=4, smallMar=TRUE)
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' #### Libraries
library(RColorBrewer)
set1 <- brewer.pal(8, "Set1")
col1 <- brewer.pal(9, "Set1")[c(1,4,3,5,7,8,2,9)]
source("../../CoxHD/CoxHD/R/ecoxph.R")
source("../../CoxHD/CoxHD/R/functions.R")
source("../../../../Projects/sandbox/mg14.R")


#' 1. Load data
#' -------------
#' ### Clinical
#' Loading
clinicalData <- read.table("../data/Ulm1.8_sm_Clinical.txt", sep="\t", header=TRUE, na.strings = "na", comment.char = "", quote="\"")
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
dim(clinicalData)

#' ### Mutation data
mutationData = read.table("../data/Ulm1.8_sm_Genetic.txt", sep="\t", header=TRUE, strip.white = TRUE)
mutationData$SAMPLE_NAME <- factor(as.character(mutationData$SAMPLE_NAME), levels = levels(clinicalData$PDID)) ## Refactor
mutationTable <- (table(mutationData[mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & mutationData$FINAL_CALL == "OK"  ,c("SAMPLE_NAME","GENE")]) > 0)+0
dim(mutationTable)

all(rownames(mutationTable)==clinicalData$PDID)

#' #### A few functions
makeInteractions
makeInteger

trialArms <- makeInteger(clinicalData$Study)

#' 2. Preparing data
#' ------------------
#' #### Survival object
library(survival)
survival <- Surv(clinicalData$efs, clinicalData$EFSSTAT) #EFS
tplIdx <- (clinicalData$Time_Diag_TPL < clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL)) | (is.na(clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL)))

#' #### All data as list
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
		Cytogenetics = clinicalData[,c(47:53,55:71)],
		Treatment = cbind(trialArms[,2:3], ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL=tplIdx),
		Clinical = cbind(clinicalData[, c("AOD","gender","Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH_","HB","platelet","Splenomegaly")], makeInteger(clinicalData$TypeAML)[,-1]))#,
#MolRisk = makeInteger(clinicalData$M_Risk))
#dataList$Genetics$CEBPA <-  clinicalData$CEBPA # encoded as 0,1,2
dataList$Genetics$CEBPA_mono <-  clinicalData$CEBPA == 1 # encoded as 0,1,2
dataList$Genetics$CEBPA_bi <-  clinicalData$CEBPA == 2 # encoded as 0,1,2
dataList$Genetics$CEBPA <- NULL
dataList$Genetics$FLT3 <- NULL
dataList$Genetics$FLT3_ITD <- clinicalData$ITD != "0"
dataList$Genetics$FLT3_TKD <- clinicalData$TKD != "0"
dataList$Genetics$IDH2_p172 <- clinicalData$IDH2 == "172"
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Genetics$MLL_PTD <- clinicalData$MLL_PTD
dataList$Cytogenetics$MLL_PTD <- NULL
dataList$Genetics = dataList$Genetics + 0
dataList$GeneGene <- makeInteractions(data.frame(dataList$Genetics), data.frame(dataList$Genetics))[,as.vector(upper.tri(matrix(0,ncol=ncol(dataList$Genetics), nrow=ncol(dataList$Genetics))))]
dataList$GeneGene <- dataList$GeneGene[,colSums(dataList$GeneGene, na.rm=TRUE)>0] 
dataList$GeneCyto <- makeInteractions(dataList$Genetics, dataList$Cytogenetics)
dataList$GeneCyto <- dataList$GeneCyto[,colSums(dataList$GeneCyto, na.rm=TRUE) > 0]
dataList$GeneTreat <- makeInteractions(dataList$Genetics, dataList$Treatment[c("TPL","ATRA","VPA")])
dataList$GeneTreat <- dataList$GeneTreat[,colSums(dataList$GeneTreat, na.rm=TRUE) > 0]
dataList$CytoTreat <- makeInteractions(dataList$Cytogenetics, dataList$Treatment[c("TPL","ATRA","VPA")])
dataList$CytoTreat <- dataList$CytoTreat[,colSums(dataList$CytoTreat, na.rm=TRUE) > 0]

#' #### Condensing to a data.frame
dataFrame <- do.call(cbind,dataList)
dataFrame <- StandardizeMagnitude(dataFrame)
names(dataFrame) <- unlist(sapply(dataList, names))
dim(dataFrame)

groups <- factor(unlist(sapply(names(dataList), function(x) rep(x, ncol(dataList[[x]])))))
table(groups)

#' Poor man's imputation
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
dataFrame <- as.data.frame(sapply(dataFrame, poorMansImpute))

#' Compute gene:gene interactions
#+ interactions, cache=TRUE, fig.width=6, fig.height=6
data <- cbind(dataList$Genetics, dataList$Cytogenetics)
data <- data[,colSums(data, na.rm=TRUE)>=16]
dim(data)
interactions <- sapply(1:ncol(data), function(i) sapply(1:ncol(data), function(j) {f<- try(fisher.test(data[,i], data[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
odds <- sapply(1:ncol(data), function(i) sapply(1:ncol(data), function(j) {f<- try(fisher.test(data[,i] + .5, data[,j] +.5), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
diag(interactions) <- 0
diag(odds) <- 1
colnames(odds) <- rownames(odds) <- colnames(interactions) <- rownames(interactions) <- colnames(data)
odds[10^-abs(interactions) > 0.05] = 1
odds[odds<1e-3] = 1e-4
odds[odds>1e3] = 1e4
logodds=log10(odds)

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


#' 3. Survival analyses
#' -----------------

#' ### 0. Proportional hazards test
#+ fig.width=14, fig.height=5
coxModel <- coxph(survival ~ gender + AOD + Study, data=clinicalData) # Most basic
phTest <- cox.zph(coxModel)
phTest
par(mfrow=c(1,4), bty="n")
for(i in 1:4) plot(phTest[i])


#' #### Basic KM plot
#+ eval=FALSE, echo=FALSE
plot(survfit(survival ~ TP53 + TPL, data=dataFrame), col=rep(set1[1:2],each=2), lty=1:2, mark=16)


#' ### 1. Random effects: static model
#+ coxRFXFit, warning=FALSE, cache=TRUE
## Pre-select based on univariate test
#univP <- sapply(dataFrame, function(x){summary(coxph(survival~x))$logtest[3]})
#whichRFX <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
mainGroups <- c("Genetics","Treatment","Clinical","Cytogenetics")
mainIdx <- groups %in% mainGroups
whichRFX <- which(colSums(dataFrame)>=8 | mainIdx) # ie, > 0.5%

## Fit Cox model
coxRFXFit <- CoxRFX(dataFrame[,whichRFX], survival, groups=groups[whichRFX], sigma0=0.1, nu=0)
#coxRFXFit <- CoxRFX(dataFrame[,whichRFX], survival, groups=groups[whichRFX], sigma0=0.01) ## reassuring that sigma doesn't shrink further


#' Coefficients
par(mar=c(5,7,1,1))
o <- order(coxRFXFit$mu)
boxplot(coef(coxRFXFit) ~ factor(coxRFXFit$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
abline(v=0)

#' ### 2. Random effects: Time-dependent model (TPL)
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

#' Construct data.frame
dataFrameTD <- dataFrame[tplSplit,]
dataFrameTD[which(tplIndex), grep("TPL", colnames(dataFrameTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

patientTD <- c(clinicalData$PDID, clinicalData$PDID[which(tplIndex)])

#' Fit TD Cox RFX model
#+ coxRFXFitTD, cache=TRUE
#whichRFXTD <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
whichRFXTD <- which(colSums(dataFrame)>=8 | mainIdx)

coxRFXFitTD <- CoxRFX(dataFrameTD[,whichRFXTD], survivalTD, groups[whichRFXTD], sigma0=0.1, nu=0)

par(mar=c(5,7,1,1))
o <- order(coxRFXFitTD$mu)
boxplot(coef(coxRFXFitTD)/log(2) ~ factor(coxRFXFitTD$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
abline(v=0)

plot(coef(coxRFXFit), coef(coxRFXFitTD), col=col1[groups[whichRFXTD]]) # Note the sign change for TPL..
abline(0,1)

#' #### Coefficient estimates (MAP)
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

#' #### Partial risk contributions
partRiskTD <- PartialRisk(coxRFXFitTD)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponents <- diag(cov(partRiskTD, use="complete"))
varianceComponents
pie(abs(varianceComponents), col=col1)
title("Risk contributions")

partRiskVar <- PartialRiskVar(coxRFXFitTD)
x=c(varianceComponents, Error=mean(rowSums(partRiskVar)))
pie(x, col=c(col1, "grey"), labels = paste(names(x), round(x/sum(x),2)))

m <- colMeans(partRiskVar)
barplot(varianceComponents +m, border=NA, col= colTrans(col1), las=2, xaxt="n", ylab="Risk variance component")
barplot(varianceComponents , border=NA, col= colTrans(col1,1), add=TRUE, xaxt="n", yaxt="n")
barplot(pmax(0,varianceComponents-m) , border=NA, col= colTrans(col1,0), add=TRUE, xaxt="n", yaxt="n") -> b
rotatedLabel(b, rep(0,8), names(varianceComponents), srt=45)

o <- order(varianceComponents)
stars(matrix(varianceComponents[o], nrow=1) +m, draw.segments=TRUE, col.segments=colTrans(col1)[o], scale=FALSE, col.lines=0, lty=0, labels="")
stars(matrix(varianceComponents[o], nrow=1) , draw.segments=TRUE, col.segments=colTrans(col1,1)[o], scale=FALSE, col.lines=0, lty=0, labels="", add=TRUE)
stars(matrix(pmax(0,varianceComponents[o] -m), nrow=1), draw.segments=TRUE, col.segments=col1[o], scale=FALSE, col.lines=0, lty=0, labels="", add=TRUE)



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

#' #### Harrel's C
#library(Hmisc)
totalRiskTD <- rowSums(partRiskTD)
survConcordance( survivalTD~totalRiskTD)
predictiveRiskTD <- rowSums(partRiskTD[,-which(colnames(partRiskTD) %in% c("Treatment","GeneTreat","CytoTreat"))])
survConcordance( survivalTD~predictiveRiskTD)

#' #### Genomic risk 
whichGenomic <- which(!colnames(partRiskTD) %in% c("Treatment","GeneTreat","CytoTreat","Clinical"))
genomicRiskTD <- rowSums(partRiskTD[,whichGenomic])
survConcordance( survivalTD~genomicRiskTD)

#' #### AUC
library(survAUC)
s <- survival
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], genomicRiskTD[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskTD[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)

#' #### Genomic risk groups
genomicRiskGroups <-  cut(genomicRiskTD[!is.na(genomicRiskTD)], quantile(genomicRiskTD, seq(0,1,0.25)), labels = c("very low","low","high","very high"))
survConcordance( survivalTD~as.numeric(genomicRiskGroups))
plot(survfit(survivalTD[!is.na(genomicRiskTD)] ~ genomicRiskGroups), col=brewer.pal(4,"Spectral")[4:1], mark=16)
table(clinicalData$M_Risk, genomicRiskGroups[1:nrow(clinicalData)])[c(2,3,4,1),]

#' Risk Plots
#+ genomicRisk, fig.width=6, fig.height=6
par(mfrow=c(4,1))
i <- 0
for(l in levels(clinicalData$M_Risk)[c(2,3,4,1)]){
	barplot(sapply(split(as.data.frame(partRiskTD[1:nrow(clinicalData),whichGenomic][clinicalData$M_Risk ==l,]), genomicRiskGroups[1:nrow(clinicalData)][clinicalData$M_Risk ==l] ), colMeans), beside=TRUE, legend=i==4, col=col1[whichGenomic], main=l, xlim=c(1,45))
	i <- i+1
}


#' Partial values of Harrel's C
c <- PartialC(coxRFXFitTD)
b <- barplot(c[1,], col=col1, las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' 4. Cross-validation
#' ----------------
#' Train only on 2/3 of the data
set.seed(42)
trainIdx <- sample(c(TRUE,FALSE), nrow(dataFrame), replace=TRUE, prob=c(0.66,0.34))
trainIdxTD <- trainIdx[tplSplit]

#' Pre-select based on univariate test
#+ warning=FALSE
#p <- sapply(dataFrame, function(x){summary(coxph(survivalTD[trainIdxTD]~x[trainIdxTD]))$logtest[3]})
whichTrain <- whichRFXTD# which(p.adjust(p,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))

#' ### 1. Static random effects model
#+ coxRFXFitTrain, cache=TRUE
coxRFXFitTrain <- CoxRFX(dataFrame[trainIdx,whichTrain], survival[trainIdx], groups=groups[whichTrain], sigma0 = 0.1, nu=0)

#' ### 2. Time-dependent random effects model
#+ coxRFXFitTDTrain, cache=TRUE
coxRFXFitTDTrain <- CoxRFX(dataFrameTD[trainIdxTD,whichTrain], survivalTD[trainIdxTD], groups=groups[whichTrain], sigma0 = 0.1, nu=0)

#' Partial contributions
partialRiskTest <- PartialRisk(coxRFXFitTDTrain, newX=dataFrameTD[!trainIdxTD, whichTrain])
c <- PartialC(coxRFXFitTDTrain, newX = dataFrameTD[!trainIdxTD, whichTrain], newSurv = survivalTD[!trainIdxTD])
b <- barplot(c[1,], col=col1, las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' Overall
totalRiskTest <- rowSums(partialRiskTest)
survConcordance(survivalTD[!trainIdxTD]~totalRiskTest )

#' Compared to molecular risk
predictiveRiskTest <- rowSums(partialRiskTest[,-which(colnames(partRiskTD) %in% c("Treatment","GeneTreat","CytoTreat"))])
survConcordance(survivalTD[!trainIdxTD] ~ predictiveRiskTest)

#barplot(c(CGP=survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)$concordance, MolecularRisk = survConcordance(survivalTD[!testIxTD] ~ c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[tplSplit][!testIxTD]])[[1]]))

s <- survival[!trainIdx]
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTest[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskTest[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' ### 3. Recursive partitioning & randomForestSRC
#+ tree, fig.width=5, fig.height=5
library(rpart)
library(randomForestSRC)
tree <- rpart(survival ~ ., data=dataFrame[mainIdx])
plot(tree)
text(tree)
survConcordance(na.omit(survival)~predict(tree))

#+ treeCV, fig.width=5, fig.height=5
treeCV <- rpart(survival[trainIdx] ~ ., data=dataFrame[trainIdx,mainIdx])
plot(treeCV)
text(treeCV)
survConcordance(survival[!trainIdx]~predict(treeCV, newdata=dataFrame[!trainIdx, mainIdx]))

#+ rForest, cache=TRUE
rForest <- rfsrc(Surv(time, status) ~.,data= cbind(time = survival[,1], status = survival[,2], dataFrame[,mainIdx]), ntree=100)
boxplot(rForest$importance ~ factor(as.character(groups[mainIdx])), border= col1, staplewex=0, pch=16, cex=0.75, ylab="RSF importance", lty=1, xaxt="n")
sapply(mainGroups, function(g) vimp(rForest, colnames(dataFrame)[which(groups==g)]))

survConcordance(na.omit(survival)~predict(rForest, importance="none")$predicted)

#+ rForestCV, cache=TRUE
rForestCV <- rfsrc(Surv(time, status) ~.,data= cbind(time = survival[,1], status = survival[,2], dataFrame[,mainIdx])[trainIdx,], ntree=100, importance="none")
p <- predict(rForestCV, newdata = dataFrame[!trainIdx,mainIdx], importance="none")
survConcordance(survival[!trainIdx]~ p$predicted)


#' ### 4. Stability selection: All terms
library(parallel)
library(glmnet)
source("../../CoxHD/CoxHD/R/stacoxph.R")
source("../../CoxHD/CoxHD/R/r_concave_tail.R")

#' Fit model
#+ stabCox, cache=TRUE
set.seed(42)
stabCox <- StabCox(dataFrame[!is.na(survival),], survival[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

plot(stabCox)

selected <- which(stabCox$Pi > .8)
selected
coxFit <- coxph(survival[trainIdx] ~ ., data=dataFrame[trainIdx,][,names(selected)])
summary(coxFit)

totalRiskStabTest <- as.matrix(dataFrame[!trainIdx,][,names(selected)]) %*% coxFit$coefficients
survConcordance(survival[!trainIdx] ~ totalRiskStabTest)


#' ### 5. Stability selection with interactions
#' Only include interaction terms when main effects are present. 
#' On the other hand, the exclusive selection of a product term could mean that the most likely explanation is that the two main terms are zero and only the interaction is non-zero..
StabCoxInteractions

#' #### Fit model to entire data set
#+ stabCoxInt, cache=TRUE
set.seed(42)
stabCoxInt <- StabCoxInteractions(dataFrame[!is.na(survival),groups %in% mainGroups], na.omit(survival), bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedInt <- names(which(stabCoxInt$Pi > 0.8))
selectedInt[21] <- "TPL:NPM1"

#' Concordance of the resulting static coxph model
coxFit <- coxph(survival[trainIdx] ~ ., data=dataFrame[trainIdx, selectedInt])
summary(coxFit)
survConcordance(survival[!trainIdx] ~ predict(coxFit, newdata = dataFrame[!trainIdx, selectedInt]))

#' Concordance of the corresponding time-dependent coxph model
coxFitTD <- coxph(survivalTD[trainIdxTD] ~ ., data=dataFrameTD[trainIdxTD, selectedInt])
summary(coxFitTD)
survConcordance(survivalTD[!trainIdxTD] ~ predict(coxFitTD, newdata = dataFrameTD[!trainIdxTD, selectedInt]))

#' KM plots of interaction terms
#+ stabCoxIntPlot
for(i in grep(":", selectedInt, value=TRUE)){
	j <- strsplit(i, ":")[[1]]
	plot(survfit(as.formula(paste("survivalTD ~ ", paste( j, collapse="+"))), data=dataFrameTD), col=set1)
	legend("topright", lty=1, col=set1, c("none", rev(j), "both"), bty="n")
}

#' Plots of the stability selection: No large difference between when interaction terms are enforce to only occur after main effects are included
plot(stabCoxInt)
plot(stabCoxInt$Pi, stabCox$Pi[match(names(stabCoxInt$Pi), names(stabCox$Pi))])

#' Partial risk for stabCoxInt
names(groups) <- colnames(dataFrame)
X <- as.matrix(dataFrameTD[selectedInt])
X <- X - rep(colMeans(X), each=nrow(X))
partRiskStabCoxInt <- sapply(levels(groups), function(l) {
			w <- groups[selectedInt] == l
			if(length(w)>0)
				X[,w, drop=FALSE] %*% coef(coxFitTD)[w]
			else
				rep(0, nrow(X))
		})
			
#' Genomic risk 
whichGenomicStabCoxInt <- which(!colnames(partRiskStabCoxInt) %in% c("Treatment","GeneTreat","CytoTreat","Clinical"))
genomicRiskStabCoxInt <- rowSums(partRiskStabCoxInt[,whichGenomicStabCoxInt])
survConcordance( survivalTD~genomicRiskStabCoxInt)

#' #### Genomic risk groups
genomicRiskGroupsStabCoxInt <-  cut(genomicRiskStabCoxInt[!is.na(genomicRiskStabCoxInt)], quantile(genomicRiskStabCoxInt, seq(0,1,0.25)), labels = c("very low","low","high","very high"))

survConcordance( survivalTD~as.numeric(genomicRiskGroupsStabCoxInt))
plot(survfit(survivalTD[!is.na(genomicRiskStabCoxInt)] ~ genomicRiskGroupsStabCoxInt), col=brewer.pal(4,"Spectral")[4:1], mark=16)
table(clinicalData$M_Risk, genomicRiskGroupsStabCoxInt[1:nrow(clinicalData)])[c(2,3,4,1),]

#' Risk Plots
#+ genomicRiskStabCoxInt, fig.width=6, fig.height=6
par(mfrow=c(4,1))
i <- 1
for(l in levels(clinicalData$M_Risk)[c(2,3,4,1)]){
	barplot(sapply(split(as.data.frame(partRiskStabCoxInt[1:nrow(clinicalData),whichGenomicStabCoxInt][clinicalData$M_Risk ==l,]), genomicRiskGroupsStabCoxInt[1:nrow(clinicalData)][clinicalData$M_Risk ==l] ), colMeans), beside=TRUE, legend=i==4, col=col1[whichGenomicStabCoxInt], main=l, xlim=c(1,30))
	i <- i+1
}

#' #### Fit model to the training subset only
#+ stabCoxIntCV, cache=TRUE
stabCoxIntCV <- StabCoxInteractions(dataFrame[!is.na(survival) & trainIdx,groups %in% mainGroups], na.omit(survival[trainIdx]), bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedIntCV <- names(which(stabCoxIntCV$Pi > 0.8))
survConcordance(survival[!trainIdx] ~ predict(coxph(survival[trainIdx] ~ ., data=dataFrame[trainIdx, selectedIntCV]), newdata = dataFrame[!trainIdx, selectedIntCV]))

tmp <- CoxRFX(dataFrame[trainIdx, selectedInt], survival[trainIdx], nu=0)
survConcordance(survival[!trainIdx] ~ as.matrix(dataFrame[!trainIdx, selectedInt]) %*% coef(tmp))


#' ### 6. Stepwise model selection
#' #### BIC
#+ coxBIC, cache=TRUE, warning=FALSE
c <- coxph(survivalTD[trainIdxTD] ~ 1, data=dataFrameTD[trainIdxTD,mainIdx])
coxBIC <- step(c, scope= as.formula(paste("survivalTD ~", paste(colnames(dataFrameTD)[mainIdx], collapse="+"))), k = log(sum(trainIdx)), trace=0)
summary(coxBIC)
survConcordance(survivalTD[!trainIdxTD] ~ predict(coxBIC, newdata=dataFrameTD[!trainIdxTD,mainIdx]))

#' #### AIC
#+ coxAIC, cache=TRUE, warning=FALSE
coxAIC <- step(c, scope= as.formula(paste("survivalTD ~", paste(colnames(dataFrameTD)[mainIdx], collapse="+"))), k = 2, trace=0)
summary(coxAIC)
survConcordance(survivalTD[!trainIdxTD] ~ predict(coxAIC, newdata=dataFrameTD[!trainIdxTD,mainIdx]))

#' ### 7. Summary of different models
#' #### Static models (may exaggerate TPL)
predictedRiskCV <- data.frame(
		stdRisk = c(4,1,2,3)[clinicalData$M_Risk[!trainIdx]],
		tree = predict(treeCV, newdata=dataFrame[!trainIdx, mainIdx]),
		rForest = predict(rForestCV, newdata = dataFrame[!trainIdx,mainIdx], importance="none")$predicted,
		coxRFX = as.matrix(dataFrame[!trainIdx,whichRFX]) %*% coef(coxRFXFitTrain),
		coxBIC = predict(coxBIC, newdata=dataFrame[!trainIdx,mainIdx]),
		coxAIC = predict(coxAIC, newdata=dataFrame[!trainIdx,mainIdx]),
		stabSel = predict(coxph(survival[trainIdx] ~ ., data=dataFrame[trainIdx, selectedIntCV]), newdata = dataFrame[!trainIdx, selectedIntCV])
)

#+ concordanceCV
concordanceCV <- sapply(predictedRiskCV, function(x) {c <- survConcordance(survival[!trainIdx] ~ x); c(c$concordance, c$std.err)})
concordanceCV
barplot(concordanceCV[1,], border=NA, col= set1, las=2, xaxt="n", ylab="Concordance") -> b
segments(b,concordanceCV[1,]-concordanceCV[2,],b,concordanceCV[1,]+concordanceCV[2,])
rotatedLabel(b, rep(0,length(b)), colnames(concordanceCV), srt=45)

#+ aucCV
aucCV <- sapply(predictedRiskCV, function(x) survivalROC(Stime=survival[!is.na(survival) & !trainIdx,1], status=survival[!is.na(survival) &  !trainIdx,2], marker = x[!is.na(survival[!trainIdx])], predict.time = 278, method="KM", cut.values=seq(-5,5,0.1))$AUC)
aucCV
barplot(aucCV, border=NA, col= set1, las=2, xaxt="n", ylab="AUC (median EFS)") -> b
rotatedLabel(b, rep(0,length(b)), names(aucCV), srt=45)


#' #### Time-dependent models
predictedRiskCVTD <- data.frame(
		stdRisk = c(4,1,2,3)[clinicalData$M_Risk[tplSplit][!trainIdxTD]],
		coxRFX = as.matrix(dataFrameTD[!trainIdxTD,whichRFX]) %*% coef(coxRFXFitTDTrain),
		coxBIC = predict(coxBIC, newdata=dataFrameTD[!trainIdxTD,mainIdx]),
		coxAIC = predict(coxAIC, newdata=dataFrameTD[!trainIdxTD,mainIdx]),
		stabSel = predict(coxph(survivalTD[trainIdxTD] ~ ., data=dataFrameTD[trainIdxTD, selectedIntCV]), newdata = dataFrameTD[!trainIdxTD, selectedIntCV])
)

#+ concordanceCVTD
concordanceCVTD <- sapply(predictedRiskCVTD, function(x) {c <- survConcordance(survivalTD[!trainIdxTD] ~ x); c(c$concordance, c$std.err)})
concordanceCVTD
barplot(concordanceCVTD[1,], border=NA, col= set1, las=2, xaxt="n", ylab="Concordance") -> b
segments(b,concordanceCVTD[1,]-concordanceCVTD[2,],b,concordanceCVTD[1,]+concordanceCVTD[2,])
rotatedLabel(b, rep(0,length(b)), colnames(concordanceCVTD), srt=45)

#' 5. Subsets and predictive power
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
						stbCx <- StabCox(dataFrame[!is.na(survival) & trn,mainIdx], survival[!is.na(survival) & trn], bootstrap.samples = 50, control="BH") ## Only main effects
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
#+ subsetOracle
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

#+ save
save(list = ls(), file=paste(Sys.Date(), "-AML.RData", sep=""))


#' 6. Germline polymorphisms
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
					return(rep(NA, ncol(mutationTable)))
				vcf <- readVcf(paste("/Volumes/nst_links/live/845/",stem,"/",stem,".cave.annot.vcf.gz", sep=""), "GRCh37")
				#genes <- sub("\\|.+","",info(vcf)$VD[info(vcf)$SNP & isCoding(vcf)])
				genes <- sub("\\|.+","",info(vcf)$VD[vcf %over% ensemblVariation & isCoding(vcf)])
				colnames(mutationTable) %in% genes
			}))
	colnames(r) <- colnames(mutationTable)
	r
}

#+ eval=FALSE
snp <- getSNPs(rownames(mutationTable))
save(snp, file="snp.RData")

#+ eval=FALSE
dataFrameTD <-  dir("/Volumes/nst_links/live/845/", pattern="WGA*")
allCaveOut <- sapply(rownames(mutationTable), function(s){
			i<<-i+1
			cat(ifelse(i%%100 == 0, "\n","."))
			stem <<- grep(s, dataFrameTD, value=TRUE)[1]
			if(is.na(stem))
				return(NA)
			readVcf(paste("/Volumes/nst_links/live/845/",stem,"/",stem,".cave.annot.vcf.gz", sep=""), "GRCh37")
		})
save(allCaveOut, file="allCaveOut.RData")
v <- snps[!is.na(snps$MF)]
s <- matrix(0, ncol = ncol(mutationTable), nrow=nrow(mutationTable), dimnames = dimnames(mutationTable))
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
