#' AML data analysis
#' =================

#+ echo=FALSE
options(width=120)
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' Libraries
library(RColorBrewer)

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
#+ cache=TRUE
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

source("~/Projects/CoxHD/CoxHD/R/ecoxph.R")

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

#' ### Anova with CoxRFX
#' Prepating data
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
dataList$GeneGene <- makeInteractions(data.frame(oncogenics), data.frame(oncogenics))[,as.vector(upper.tri(matrix(0,ncol=ncol(oncogenics), nrow=ncol(oncogenics))))] 
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


#' ### Fit model
#+ warning=FALSE
## Pre-select based on univariate test
univP <- sapply(dataFrame, function(x){summary(coxph(survival~x))$logtest[3]})
w <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))
## Fit Cox model
coxRFXFit <- CoxRFX(dataFrame[,w], survival, groups=groups[w], sigma0=0.1)

#' Coefficients
boxplot(coxRFXFit$coefficients ~ groups[w], las=2, border=set1)
abline(h=0)

#+ fig.width=14, fig.height=5
for(g in levels(groups)){
	par(mar=c(7,4,2,0))
	x <- sort(coxRFXFit$coefficients[groups[w]==g])
	v <- diag(coxRFXFit$var)[groups[w]==g][ order(coxRFXFit$coefficients[groups[w]==g])]
	plot(x, las=2, pch=16, xaxt="n", main=g, cex = .5+pmin(3,-log10(pchisq((x-coxRFXFit$mu[g])^2/v, 1, lower.tail=FALSE))))
	segments(1:length(x),x - 2*sqrt(v) ,1:length(x), x+2*sqrt(v))
	axis(side=1, at = 1:length(x), labels = names(x), las=2, cex.axis=.6)
	abline(v=1:length(x), lty=3, col="grey")
	abline(h=coxRFXFit$mu[g])
	abline(h = coxRFXFit$mu[g] + c(-1,1)*sqrt(coxRFXFit$sigma2[g]), lty=2)
}

#' ### Partial risk contributions
partRisk <- PartialRisk(coxRFXFit)
pie(rowSums(cov(partRisk, use="complete")), col=brewer.pal(8,"Set1"))
title("Risk contributions")

#' #### Stars
#+ fig.width=12, fig.height=12
library(HilbertVis)
nStars <- 32
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
h <- hclust(dist(partRisk[1:(nStars^2),]))
stars(exp(-(partRisk[1:(nStars^2),][h$order,] - rep(colMeans(partRisk), each=nStars^2)))/2, scale=FALSE, locations=locations, key.loc=c(-1.5,-1.5), col.lines=rep(1,(nStars^2)), col.stars = brewer.pal(11,'RdBu')[cut(survival[1:(nStars^2),1][h$order], quantile(survival[,1], seq(0,1,0.1), na.rm=TRUE))])
symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=3)

#' ### A few performance measures
#' Harrel's C
library(Hmisc)
totalRisk <- na.omit(rowSums(partRisk))
rcorr.cens(-totalRisk, survival)
predictiveRisk <- rowSums(partRisk[,-which(colnames(partRisk) %in% c("Treatment","GT"))])
rcorr.cens(-predictiveRisk, survival)

#' AUC
library(survAUC)
s <- survival
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRisk[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRisk[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' Risk quartiles
predictiveRiskGroups <-  cut(predictiveRisk[!is.na(predictiveRisk)], quantile(predictiveRisk, seq(0,1,0.25)), labels = c("very low","low","high","very high"))
rcorr.cens(-as.numeric(predictiveRiskGroups), survival)
plot(survfit(survival[!is.na(predictiveRisk)] ~ predictiveRiskGroups), col=brewer.pal(4,"Spectral"), mark=16)
table(clinicalData$M_Risk, predictiveRiskGroups)[c(2,3,4,1),]

#' Partical values of Harrel's C
barplot(PartialC(coxRFXFit), col=brewer.pal(8,"Set1"), las=2)
abline(h=.5)

#' ## Holding out a test set.
set.seed(42)
testIx <- sample(c(TRUE,FALSE), nrow(dataFrame), replace=TRUE, prob=c(0.66,0.34))

#' Pre-select based on univariate test
#+ warning=FALSE
p <- sapply(dataFrame, function(x){summary(coxph(survival[testIx]~x[testIx]))$logtest[3]})
w <- which(p.adjust(p,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))

#' Fit Cox model
coxRFXFitTrain <- CoxRFX(dataFrame[testIx,w], survival[testIx], groups=groups[w], sigma0 = 0.1)

#' Partial contributions
partialRiskTest <- PartialRisk(coxRFXFitTrain, newX=dataFrame[!testIx, w])
barplot(PartialC(coxRFXFitTrain, newX = dataFrame[!testIx, w], newSurv = survival[!testIx]), col=brewer.pal(8,"Set1"), las=2)
abline(h=.5)

#' Overall
totalRiskTest <- rowSums(partialRiskTest)
rcorr.cens(- totalRiskTest, survival[!testIx])

#' Compared to molecular risk
predictiveRiskTest <- rowSums(partialRiskTest[,-which(colnames(partRisk) %in% c("Treatment","GeneTreatment"))])
rcorr.cens(- predictiveRiskTest, survival[!testIx])

barplot(c(CGP=rcorr.cens(- predictiveRiskTest, survival[!testIx])[1], MolecularRisk = rcorr.cens(c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[!testIx]], survival[!testIx])[1]))

s <- survival[!testIx]
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTest[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskTest[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' Stabilitiy selection
#' --------------------
library(parallel)
library(glmnet)
set.seed(42)
source("/Users/mg14/Projects/CoxHD/CoxHD/R/stacoxph.R")
w <- TRUE
stabCox <- stacoxph(dataFrame[!is.na(survival),w], survival[!is.na(survival)], bootstrap.samples = 50, control="BH", level=.1)

plot(stabCox)

selected <- which(stabCox$Pi > .8)
selected
X <- as.matrix(dataFrame[testIx,names(selected)]) + 0
coxFit <- coxph(survival[testIx] ~ ., data=dataFrame[testIx,w][,names(selected)])
summary(coxFit)

totalRiskStabTest <- as.matrix(dataFrame[!testIx,w][,names(selected)]) %*% coxFit$coefficients
rcorr.cens(- totalRiskStabTest, survival[!testIx])

plot(survfit(survival ~ FLT3_ITD + DNMT3A, data = dataList$Genetics), col=c("grey",brewer.pal(3,"Set1")))
legend("topright", bty="n",  lty=1, col=c("grey",brewer.pal(3,"Set1")), c("WT/WT","WT/DNMT3A-","FLT3-/WT","FLT3-/DNMT3A-"))

#' Probably effect of time-dep TPL:
plot(survfit(survival ~ `Genetics.TP53` + `Treatment.TPL`, data = dataFrame), col=c("grey",brewer.pal(3,"Set1")))
legend("topright", bty="n",  lty=1, col=c("grey",brewer.pal(3,"Set1")), c("WT/TPL-","WT/TPL+","TP53-/TPL-","TP53-/TPL+"))

 
#' ### Time-dependent survival effects (TPL)
#' Construct a time-dependent Surv() object by splitting TPL patients into pre and post
#survWOTpl = Surv(time = rep(0, nrow(clinicalData)), time2=as.numeric(clinicalData$Date_LF -  clinicalData$ERDate), event=clinicalData$Status )
#survWOTpl[clinicalData$TPL == 1] = Surv(rep(0, sum(clinicalData$TPL)), as.numeric(clinicalData$TPL_date - clinicalData$ERDate)[clinicalData$TPL == 1], rep(0, sum(clinicalData$TPL)))
#t = rbind(survWOTpl, Surv(time = as.numeric(clinicalData$TPL_date - clinicalData$ERDate)[clinicalData$TPL == 1], time2 = as.numeric(clinicalData$Date_LF -  clinicalData$ERDate)[clinicalData$TPL == 1], clinicalData$Status[clinicalData$TPL == 1]))
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
#t = c(clinicalData$TPL, rep(1, sum(clinicalData$TPL)))
#t = c(rep(0,nrow(clinicalData)), rep(1, length(which(tplIndex))))
#plot(survfit(survivalTD ~ t), col=1:2)
#d = data.frame(risk=clinicalData$M_Risk, study=clinicalData$Study, I = survival[,1] > 257)
#d <- cbind(rbind(d, d[clinicalData$TPL==1,]), TPL=t==1)
#m <- factor(c(clinicalData$M_Risk, clinicalData$M_Risk[clinicalData$TPL==1]), labels=levels(clinicalData$M_Risk))
#coxph(survivalTD ~ TPL + I, data=d)
#
#s = Surv(time = rep(0, nrow(clinicalData)), time2=as.numeric(clinicalData$Date_LF -  clinicalData$ERDate), event=clinicalData$Status )
#plot(survfit(s ~ clinicalData$TPL + (survival[,1] > 257)), col=1:2, lty=rep(c(1:2), each=2))
#
#survAfterTpl <- Surv(time = as.numeric(clinicalData$Date_LF - clinicalData$TPL_date )[clinicalData$TPL == 1], clinicalData$Status[clinicalData$TPL == 1])
#
#s = Surv(c(survWOTpl[,2], survAfterTpl[,1]), c(survWOTpl[,3], survAfterTpl[,2]))
#
#trialArm = factor(paste(clinicalData$Study, clinicalData$Randomisation_Arm, sep="."))
#
#coxph(survival ~. +tr0704_ATRA:tr0704_VPA + clinicalData$AOD + clinicalData$gender+clinicalData$M_Risk, data=trialArms[,-1])		
#	
#simpleRisk <- data.frame()

#' TD Cox RFX model
dataFrameTD <- dataFrame[tplSplit,]
dataFrameTD$Treatment.TPL <- c(rep(0,nrow(clinicalData)), rep(1, length(which(tplIndex))))
w <- which(p.adjust(univP,"BH") < 0.1 | groups %in% c("Genetics","Treatment","Clinical","Cytogenetics"))

coxRFXFitTD <- CoxRFX(dataFrameTD[,w], survivalTD, groups[w], sigma0=0.1)

plot(coef(coxRFXFit), coef(coxRFXFitTD), col=set1[groups[w]]) # Note the sign change for TPL..
boxplot(coef(coxRFXFitTD)~groups[w], border=set1)

#' Germline polymorphisms
#' ----------------------
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