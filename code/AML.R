#' AML data analysis
#' =================

#+ Preliminaries, echo=FALSE
options(width=120)
pdf.options(pointsize=8)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' #### Libraries
library(RColorBrewer)
set1 <- brewer.pal(8, "Set1")
#col1 <- brewer.pal(9, "Set1")[c(1,6,4,3,5,7,8,2,9)]
library(CoxHD)
source("../../mg14/R/mg14.R")


#' 1. Load data
#' -------------
#' ### Clinical
#' Loading
clinicalData <- read.table("../data/Ulm1.17_MG_Clinical.txt", sep="\t", header=TRUE, na.strings = "na", comment.char = "", quote="\"")
clinicalData <- clinicalData[order(clinicalData$PDID),]
clinicalData$ERDate <- as.Date(as.character(clinicalData$ERDate), "%d-%b-%y")
clinicalData$CR_date <- as.Date(as.character(clinicalData$CR_date), "%d-%b-%y")
clinicalData$TPL_date <- as.Date(as.character(clinicalData$TPL_date), "%d-%b-%y")
clinicalData$Date_LF <- as.Date(as.character(clinicalData$Date_LF), "%d-%b-%y")
clinicalData$Recurrence_date <- as.Date(as.character(clinicalData$Recurrence_date), "%d-%b-%y")
levels(clinicalData$Study) <- c( "AMLSG0704" ,   "AMLHD98A" ,   "AMLHD98B")
clinicalData$VPA[is.na(clinicalData$VPA)] <- 0
clinicalData$ATRA_arm[is.na(clinicalData$ATRA_arm)] <- 0
colnames(clinicalData) <- gsub('\\.',"",colnames(clinicalData))
clinicalData <- clinicalData[!is.na(clinicalData$TypeAML),] ## remove unknown patients
clinicalData$PDID <- factor(as.character(clinicalData$PDID))
dim(clinicalData)

#' ### Mutation data
mutationData = read.table("../data/Ulm1.14_MG_Genetic.txt", sep="\t", header=TRUE, strip.white = TRUE)
mutationData$SAMPLE_NAME <- factor(as.character(mutationData$SAMPLE_NAME), levels = levels(clinicalData$PDID)) ## Refactor
mutationTable <- (table(mutationData[mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & mutationData$FINAL_CALL == "OK" ,c("SAMPLE_NAME","GENE")]) > 0)+0
dim(mutationTable)

all(rownames(mutationTable)==clinicalData$PDID)

#' #### A few functions
MakeInteractions
MakeInteger

trialArms <- MakeInteger(clinicalData$Study)

#' 2. Preparing data
#' ------------------
tplIdxEfs <- (clinicalData$Time_Diag_TPL < clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL)) | (is.na(clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL))) # Need to restrict to those with transplant prior to event..
tplIdxOs <- !is.na(clinicalData$TPL_type)

#' #### All data as list
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
		Cytogenetics = clinicalData[,50:74],
		Trial = trialArms[,1:2],
		Treatment = data.frame(ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL_efs=tplIdxEfs, TPL_os=tplIdxOs),
		Clinical = cbind(clinicalData[, c("AOD","gender","Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH","HB","platelet","Splenomegaly")], MakeInteger(clinicalData$TypeAML)[,-1]))#,
#MolRisk = makeInteger(clinicalData$M_Risk))
#dataList$Genetics$CEBPA <-  clinicalData$CEBPA # encoded as 0,1,2
dataList$Genetics$CEBPA_mono <-  clinicalData$CEBPA == 1 # encoded as 0,1,2
dataList$Genetics$CEBPA_bi <-  clinicalData$CEBPA == 2 # encoded as 0,1,2
dataList$Genetics$CEBPA <- NULL
dataList$Genetics$FLT3 <- NULL
dataList$Genetics$FLT3_ITD <- clinicalData$FLT3_ITD != "0"
dataList$Genetics$FLT3_TKD <- clinicalData$FLT3_TKD != "0"
dataList$Genetics$FLT3_other <- clinicalData$FLT3_other != "0"
dataList$Genetics$IDH2_p172 <- table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("172", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2_p140 <-  table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("140", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2 <- NULL
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Cytogenetics$MLL_PTD <- NULL
dataList$Genetics = dataList$Genetics + 0
dataList$GeneGene <- MakeInteractions(data.frame(dataList$Genetics), data.frame(dataList$Genetics))[,as.vector(upper.tri(matrix(0,ncol=ncol(dataList$Genetics), nrow=ncol(dataList$Genetics))))]
dataList$GeneGene <- dataList$GeneGene[,colSums(dataList$GeneGene, na.rm=TRUE)>0] 
dataList$CytoCyto <- MakeInteractions(dataList$Cytogenetics, dataList$Cytogenetics)[,sapply(1:ncol(dataList$Cytogenetics), `<`, 1:ncol(dataList$Cytogenetics))]
dataList$CytoCyto <- dataList$CytoCyto[, colSums(dataList$CytoCyto, na.rm=TRUE) > 0]
dataList$GeneCyto <- MakeInteractions(dataList$Genetics, dataList$Cytogenetics)
dataList$GeneCyto <- dataList$GeneCyto[,colSums(dataList$GeneCyto, na.rm=TRUE) > 0]
dataList$GeneTreat <- MakeInteractions(dataList$Genetics, dataList$Treatment)
dataList$GeneTreat <- dataList$GeneTreat[,colSums(dataList$GeneTreat, na.rm=TRUE) > 0]
dataList$CytoTreat <- MakeInteractions(dataList$Cytogenetics, dataList$Treatment)
dataList$CytoTreat <- dataList$CytoTreat[,colSums(dataList$CytoTreat, na.rm=TRUE) > 0]

#' #### Condensing to a data.frame
dataFrame <- do.call(cbind,dataList)
names(dataFrame) <- unlist(sapply(dataList, names))
dataFrame <- StandardizeMagnitude(dataFrame)
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
library(survival)

#' ### -1. Survival objects
#' #### Event-free survival
#' Static EFS
efs <- Surv(clinicalData$efs, clinicalData$EFSSTAT) #EFS

#' Construct a time-dependent Surv() object by splitting TPL patients into pre and post
t <- clinicalData$Time_Diag_TPL
t[is.na(t)] <- Inf
e <- clinicalData$efs
tplIndexEfs <-  t < e
efsTD <-  Surv(time = rep(0, nrow(clinicalData)), time2=pmin(e, t), event=ifelse(tplIndexEfs, 0, clinicalData$EFSSTAT) )
efsTD <- rbind(efsTD, 
		Surv(time=t[which(tplIndexEfs)],
				time2=e[which(tplIndexEfs)], 
				event=clinicalData$EFSSTAT[which(tplIndexEfs)])
)
efsTD = Surv(efsTD[,1],efsTD[,2],efsTD[,3])
rm(e,t)
tplSplitEfs <- c(1:nrow(clinicalData), which(tplIndexEfs))

#' Construct data.frame
dataFrameEfsTD <- dataFrame[tplSplitEfs,]
dataFrameEfsTD[which(tplIndexEfs), grep("TPL", colnames(dataFrameEfsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

patientTD <- c(clinicalData$PDID, clinicalData$PDID[which(tplIndexEfs)])


#' #### Overall survival
#' OS
os <- Surv(clinicalData$OS, clinicalData$Status) #OS
t <- clinicalData$Time_Diag_TPL
t[is.na(t)] <- Inf
o <- clinicalData$OS
tplIndexOs <-  t < o
osTD <-  Surv(time = rep(0, nrow(clinicalData)), time2=pmin(o, t), event=ifelse(tplIndexOs, 0, clinicalData$Status) )
osTD <- rbind(osTD, 
		Surv(time=t[which(tplIndexOs)],
				time2=o[which(tplIndexOs)], 
				event=clinicalData$Status[which(tplIndexOs)])
)
osTD = Surv(osTD[,1],osTD[,2],osTD[,3])
rm(o,t)
tplSplitOs <- c(1:nrow(clinicalData), which(tplIndexOs))

#' Construct data.frame
dataFrameOsTD <- dataFrame[tplSplitOs,]
dataFrameOsTD[which(tplIndexOs), grep("TPL", colnames(dataFrameOsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

#' Some indeces
mainGroups <- c("Genetics","Treatment","Clinical","Cytogenetics", "Trial")
mainIdx <- groups %in% mainGroups
efsIdx <- !grepl("TPL_os", colnames(dataFrame))
whichRFXEfs <- which((colSums(dataFrame)>=8 | mainIdx) & efsIdx) # ie, > 0.5%
mainIdxEfs <- mainIdx & efsIdx
osIdx <- !grepl("TPL", colnames(dataFrame)) ## Exclude TPL from OS analyses..
whichRFXOs <- which((colSums(dataFrame)>=8 | mainIdx) & osIdx) # ie, > 0.5%
mainIdxOs <- mainIdx & osIdx


#' ### 0. Proportional hazards test
#+ fig.width=14, fig.height=5
coxModel <- coxph(os ~ gender + AOD + Study, data=clinicalData) # Most basic
phTest <- cox.zph(coxModel)
phTest
par(mfrow=c(1,4), bty="n")
for(i in 1:4) plot(phTest[i])


#' #### Basic KM plots
#' TP53 and TPL
f <- formula("osTD ~ TP53 + TPL_os")
s <- survfit(f, data=dataFrameOsTD)
c <- coxph(f, data=dataFrameOsTD)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)


#' FLT ITD v TKD
f <- formula("os ~ FLT3_ITD + FLT3_TKD")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)


#' CEBPA mono v bi
f <- formula("os ~ CEBPA_mono + CEBPA_bi")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)

#' IDH2 p140 v p172
f <- formula("os ~ IDH2_p140 + IDH2_p172")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)

#' #### NONC and complex karyotyope
table(complex=clinicalData$complex, `#Recurr.Cyto`=rowSums(dataList$Cytogenetics))
summary(coxph(os ~ ., data=dataFrame[mainIdxOs][colSums(dataFrame[mainIdxOs]) > 10]))

NONC <- rowSums(cbind(dataList$Cytogenetics[names(dataList$Cytogenetics)!="complex"], dataList$Genetics), na.rm=TRUE)
summary(coxph(os ~ . + NONC, data=dataFrame[mainIdxOs][colSums(dataFrame[mainIdxOs]) > 10]))

#+ NONC, fig.width=3, fig.height=2
r=rainbow(9, s=1, v=0.66, start=0.2,end=1)
s <- survfit(os ~ pmin(NONC,8))
plot(s, col=r, xlab="Days", ylab="Survival")
legend('topright', bty='n', col=r,legend= 0:8, lty=1, title="# onc. mut.")

plot(0:8,summary(s)$table[,"median"], xlab="Number of oncogenic mutations", ylab="Median survival time", pch=19)
segments(0:8, summary(s)$table[,"0.95LCL"],0:8, sapply(summary(s)$table[,"0.95UCL"], function(x) ifelse(is.na(x), max(os[,1]), x)))

#' ### 1. Random effects: static model
#+ coxRFXFitEfs, warning=FALSE, cache=TRUE
## Fit Cox model
coxRFXFitEfs <- CoxRFX(dataFrame[,whichRFXEfs], efs, groups=groups[whichRFXEfs], sigma0=0.1, nu=0)
#+ coxRFXFitOs, warning=FALSE, cache=TRUE
coxRFXFitOs <- CoxRFX(dataFrame[,whichRFXOs], os, groups=groups[whichRFXOs], sigma0=0.1, nu=0)

#coxRFXFit <- CoxRFX(dataFrame[,whichRFX], survival, groups=groups[whichRFX], sigma0=0.01) ## reassuring that sigma doesn't shrink further

#' Coefficients
par(mar=c(5,7,1,1))
col1 <- c(brewer.pal(12, "Paired")[c(6,3,4,5,9,1,2,7,8)],"#999999")
names(col1) <- levels(groups)[order(toupper(levels(groups)))]
#o <- order(coxRFXFitEfs$mu)
#boxplot(coef(coxRFXFitEfs) ~ factor(coxRFXFitEfs$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
plot(coxRFXFitEfs, col=col1, order=order(coxRFXFitEfs$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

#' ### 2. Random effects: Time-dependent model (TPL)
#' Fit TD CoxRFX model
#+ coxRFXFitEfsTD, cache=TRUE
coxRFXFitEfsTD <- CoxRFX(dataFrameEfsTD[,whichRFXEfs], efsTD, groups[whichRFXEfs], sigma0=0.1, nu=0)
#+ coxRFXFitOsTD, cache=TRUE
#coxRFXFitOsTD <- CoxRFX(dataFrameOsTD[,whichRFXOs], osTD, groups[whichRFXOs], sigma0=0.1, nu=0)


par(mar=c(5,7,1,1))
#o <- order(coxRFXFitEfsTD$mu)
#boxplot(coef(coxRFXFitEfsTD)/log(2) ~ factor(coxRFXFitEfsTD$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
plot(coxRFXFitEfsTD, col=col1, order=order(coxRFXFitEfsTD$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

plot(coef(coxRFXFitEfs), coef(coxRFXFitEfsTD), col=col1[groups[whichRFXEfs]]) # Note the sign change for TPL..
abline(0,1)

par(mar=c(5,7,1,1))
#o <- order(coxRFXFitOsTD$mu)
#boxplot(coef(coxRFXFitOsTD)/log(2) ~ factor(coxRFXFitOsTD$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
plot(coxRFXFitOs, col=col1, order=order(coxRFXFitOs$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

#' #### Coefficient estimates (MAP)
#+ fig.width=14, fig.height=5
for(model in c("coxRFXFitEfsTD","coxRFXFitOs")){
	m <- get(model)
	for(g in levels(groups)){
		par(mar=c(7,4,2,0))
		x <- sort(m$coefficients[m$groups==g])
		v <- diag(m$var)[m$groups==g][ order(m$coefficients[m$groups==g])]
		plot(x, las=2, pch=16, xaxt="n", main=paste(model,g), cex = .5+pmin(3,-log10(pchisq((x-m$mu[g])^2/v, 1, lower.tail=FALSE))), xlab="", ylab="coefficient")
		segments(1:length(x),x - 2*sqrt(v) ,1:length(x), x+2*sqrt(v), col=col1[g])
		axis(side=1, at = 1:length(x), labels = names(x), las=2, cex.axis=.6)
		abline(v=1:length(x), lty=3, col="grey")
		abline(h=m$mu[g])
		abline(h = m$mu[g] + c(-1,1)*sqrt(m$sigma2[g]), lty=2)
	}
}

#' #### OS v EFS
#+ EFSvsOScoef, fig.width=5, fig.height=5
efs2Os <- match(names(coef(coxRFXFitOs)), sub("efs","os",names(coef(coxRFXFitEfsTD))))
plot(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs), pch=16, xlab="Coefficient EFS", ylab="Coefficient OS", col=col1[coxRFXFitEfsTD$groups])
segments((coef(coxRFXFitEfsTD) - sqrt(diag(coxRFXFitEfsTD$var)))[efs2Os], coef(coxRFXFitOs),(coef(coxRFXFitEfsTD)+ sqrt(diag(coxRFXFitEfsTD$var)))[efs2Os], coef(coxRFXFitOs),  col=paste(col1,"44",sep="")[coxRFXFitEfsTD$groups])
segments(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs) - sqrt(diag(coxRFXFitOs$var)),coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs)+sqrt(diag(coxRFXFitOs$var)),  col=paste(col1,"44",sep="")[coxRFXFitEfsTD$groups])
text(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs), ifelse(sqrt(coef(coxRFXFitEfsTD)[efs2Os]^2 + coef(coxRFXFitOs)^2) > 0.5, names(coef(coxRFXFitEfsTD))[efs2Os], ""), pos=3)
abline(0,1)
abline(h=0, lty=3)
abline(v=0, lty=3)
points(coxRFXFitEfsTD$mu, coxRFXFitOs$mu, bg=col1, pch=21, cex=2)
title(main=paste("rho =", round(cor(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs), use="c"), 2)))

#+ EFSvsOSsigma
plot(coxRFXFitEfsTD$sigma2, coxRFXFitOs$sigma2, col=col1, pch=19, xlab="Variance EFS", ylab="Variance OS")
text(coxRFXFitEfsTD$sigma2, coxRFXFitOs$sigma2, pos=3, levels(coxRFXFitEfsTD$groups))
abline(0,1)


#' #### Partial risk contributions
simpleGroupsEfs <- mergeLevels(coxRFXFitEfsTD$groups, mergeList=list(Genetics=c("Genetics","GeneGene"), Cytogenetics=c("Cytogenetics","CytoCyto","GeneCyto"), Treatment=c("Treatment","GeneTreat","CytoTreat")))
partRiskEfsTDsimple <- PartialRisk(coxRFXFitEfsTD, groups = simpleGroupsEfs)
partRiskEfsTD <- PartialRisk(coxRFXFitEfsTD)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponentsEfs <- diag(cov(partRiskEfsTDsimple, use="complete"))
varianceComponentsEfs
PlotVarianceComponents(coxRFXFitEfs, col=col1)
title("Risk contributions EFS")

PlotVarianceComponents(coxRFXFitEfsTD, col=col1)
title("Risk contributions EFS (time-dep TPL)")


partRiskVar <- PartialRiskVar(coxRFXFitEfsTD,  groups = simpleGroupsEfs)
x=c(varianceComponentsEfs, Error=mean(rowSums(partRiskVar)))
pie(x, col=c(col1[names(varianceComponentsEfs)], "grey"), labels = paste(names(x), round(x/sum(x),2)))

simpleGroupsOs <- mergeLevels(coxRFXFitOs$groups, mergeList=list(Genetics=c("Genetics","GeneGene"), Cytogenetics=c("Cytogenetics","CytoCyto","GeneCyto"), Treatment=c("Treatment","GeneTreat","CytoTreat")))
partRiskOsSimple <- PartialRisk(coxRFXFitOs, groups = simpleGroupsOs)
partRiskOs <- PartialRisk(coxRFXFitOs)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponentsOs <- diag(cov(partRiskOsSimple, use="complete"))
varianceComponentsOs
PlotVarianceComponents(coxRFXFitOs, col=col1)
title("Risk contributions OS")

PlotVarianceComponents(coxRFXFitOs, col=col1, groups=simpleGroupsOs)
title("Risk contributions OS")

#o <- order(varianceComponents)
#stars(matrix(varianceComponents[o], nrow=1) +m, draw.segments=TRUE, col.segments=colTrans(col1)[o], scale=FALSE, col.lines=0, lty=0, labels="")
#stars(matrix(varianceComponents[o], nrow=1) , draw.segments=TRUE, col.segments=colTrans(col1,1)[o], scale=FALSE, col.lines=0, lty=0, labels="", add=TRUE)
#stars(matrix(pmax(0,varianceComponents[o] -m), nrow=1), draw.segments=TRUE, col.segments=col1[o], scale=FALSE, col.lines=0, lty=0, labels="", add=TRUE)

#' #### Distribution of survival effects
#' Effect sizes
#+ effectSizesKM, fig.width=2, fig.height=2
H0 <- basehaz(coxRFXFitOs, centered = TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
i <- 1
x <- seq(0,2500, 10)
for(g in levels(groups)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(coxRFXFitOs$coef[groups[whichRFXOs] == g]), each=length(x))  )), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(sqrt(coxRFXFitOs$sigma2[g]) * c(-1,1), each=length(x))  + coxRFXFitOs$mu[g])), col=paste(col1[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( coxRFXFitOs$mu[g] )), col=col1[i], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) ), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}

#' Partial risk
#+ partialRiskKM, fig.width=2, fig.height=2
i <- 1
for(g in levels(simpleGroupsEfs)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	m <- mean(rowSums(partRiskOsSimple[,colnames(partRiskOsSimple)!=g]))
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(partRiskOsSimple[,g]), each=length(x)) +m)), col=paste(col1[g],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskOsSimple[,g], c(0.25,0.75)), each=length(x))+m)), col=paste(col1[g],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskOsSimple[,g], c(0.05,0.95)), each=length(x))+m)), col=paste(col1[g],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( median(partRiskOsSimple[,g])+m)), col=col1[g], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) *exp(+m)), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}
par(mfrow=c(1,1))

#' #### Confidence intervals
#' Estimate confidence intervals by parametric bootstrap. Note that the usual sample with replacement 
#' yields inconsistencies for the interaction terms due to the overdispersed correlations.
#+ parBoot, cache=TRUE
set.seed(42)
risk <- as.matrix(dataFrame[whichRFXOs]) %*% coxRFXFitOs$coefficients
risk <- risk - mean(risk)
parBoot <- mclapply(1:100, function(i) {
			s <- SimSurvNonp(risk, os)
			c <- CoxRFX(dataFrame[whichRFXOs], s, groups=groups[whichRFXOs], sigma0=0.1, nu=0)
			c$X <- NULL # set X to zero to save mem
			return(c)
		}, mc.cores=10)

#+ parBootPlot, fig.width=3, fig.height=2, cache=TRUE
boxplot(t(sapply(parBoot, `[[`, "sigma2")), border=col1, lty=1, pch=16, staplewex=0, ylab="sigma2", las=2)
abline(h=0, lty=3)
points(coxRFXFitOs$sigma2,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "mu")), border=col1, lty=1, pch=16, staplewex=0, ylab="mu", las=2)
abline(h=0, lty=3)
points(coxRFXFitOs$mu,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "df")), border=col1, lty=1, pch=16, staplewex=0, ylab="df", las=2)
abline(h=0, lty=3)
points(coxRFXFitOs$df,  pch=19)

v <- t(sapply(parBoot, VarianceComponents, newX=dataFrame[whichRFXOs]))
boxplot(v, border=col1, lty=1, pch=16, staplewex=0, ylab="variance comp.", las=2)
abline(h=0, lty=3)
points(VarianceComponents(coxRFXFitOs),  pch=19)

round(cov(v, 2))

rm(parBoot)

#' #### Large overview panel
#+ overviewRFX, fig.width=8, fig.height=8
layout(matrix(1:6, nrow=2), width=c(6,4,2),height=c(3,12))
par(bty="n", mar=c(1,4,3,2), mgp=c(2.5,1,0), tcl=-.25)
v <- VarianceComponents(coxRFXFitOs, type="rowSums")
w <- TRUE#coxRFXFitOsTD$groups %in% mainGroups
f <- factor(as.character(coxRFXFitOs$groups[w]), levels=levels(coxRFXFitOs$groups[w])[order(abs(v), decreasing=TRUE)])
o <- order(f,coxRFXFitOs$coefficients[w])
c <- coef(coxRFXFitOs)[w][o]
plot(c, type='p',  ylab="Coefficient", xaxt='n', col=col1[as.character(f)[o]], xaxs="i", lwd=2, pch=NA, ylim=c(-2,2))
abline(h=-1:1, col="grey")
#segments(seq_along(c), coxRFXFitOsTD$mu[as.character(f)[o]], seq_along(c), c, col=col1[as.character(f)[o]])
y <- seq(-2,2,l=100)
points(c, col=col1[as.character(f)[o]], pch=16, cex=.66)
par(xpd=NA)
for(l in levels(f)){
	x <- dnorm(y,coxRFXFitOs$mu[l], sqrt(coxRFXFitOs$sigma2[l]))
	x <- x/max(x)*15
	#polygon(x + (cumsum(table(f)) - table(f))[l], y, col=paste(col1[l],"44",sep=""), lty=0)
	lines(x + (cumsum(table(f)) - table(f))[l], y, col=col1[l], lwd=1)
	segments((cumsum(table(f)) - table(f))[l],coxRFXFitOs$mu[l],(cumsum(table(f)) - table(f))[l]+max(x),coxRFXFitOs$mu[l], col=col1[l])
}
mtext(side=3, levels(f), at= table(f)/2+ cumsum(c(0,table(f)[-nlevels(f)])), cex=.66)
par(bty="L", mar=c(4,4,1,2), mgp=c(2.5,1,0))
X <- as.matrix(dataFrame[whichRFXOs][w][o])
p <- order(X %*% c)
X <- X[p,]
x <- t(X)/apply(X,2,max)
#h <- hclust(dist(t(x)))
image(x=1:nrow(x),y=1:ncol(x),x*.9 + as.numeric(f[o])-1 + 1e-5, useRaster=TRUE, 
		col=sapply(col1[levels(f)], function(x) c(colorRampPalette(c("white",x))(11)[-1])), 
		breaks=seq(0,nlevels(f), 0.1),  ylab="Patients",xlab="Variable", xlim=c(0,nrow(x)), ylim=c(0,ncol(x)))

par(mar=c(4,.5,3,0.5), xpd=NA)
PlotVarianceComponents(coxRFXFitOs, col = col1, type='diag')
mtext(side=3, "Variance components",line=2, cex=.66)

xo <- 3
par(mar=c(4,0,1,0.5), xpd=NA, bty="n")
plot(NA,NA,ylim=c(1,nrow(X)), xlim=c(-xo/2, (nlevels(f)-.5)*xo), xaxt="n", yaxt="n", yaxs="i", xlab="Partial log hazard", ylab="n")
i <- 0
for(l in levels(f)){
	h <- X[,f[o]==l] %*% c[f[o]==l]
	h <- h - mean(h)
	#plot(h, seq_along(h), pch=NA, xlab="log hazard", yaxs="i", yaxt="n", ylab='', xlim=c(-2,2))
	segments(h + i*xo, seq_along(h), + i*xo, seq_along(h), col=col1[l])
	#mtext(side=3, l, cex=.66, line=.5)
	axis(side=1, at=-1:1 + xo*i, labels=-1:1)
	x <- seq(-1.5,1.5,l=100)
	lines(x +i*xo, dnorm(x, 0, sd(h))*100 /  dnorm(0, 0, sd(h)) + length(h)*1.01, col=col1[l])
	i <- i+1 
}

h <- X %*% c
h <- h - mean(h)

H0 <- basehaz(coxph(os ~ 1), centered=TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
#plot(survfit(osTD ~1))
#lines(0:5000, exp(-hazardDist(0:5000)),  col='red')
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
i <- 1
x <- seq(from=-4,to=4, l=512)
par(mar=c(1,.5,3,0.5), xpd=FALSE)
plot(x, invHazardDist(-log(l[2]) /exp(x) ), col='black', lty=c(2,1,2)[2], xlab="", ylab="Time", type='l',ylim=c(0,5000), xlim=range(h), xaxt="n")
abline(v=seq(-2,4,2), col="grey")
#for(i in seq_along(l)[-1])
#	lines(x, invHazardDist(-log(l[i]) /exp(x) ), col='black', lty=c(2,1,2)[i])
polygon(c(x, rev(x)), c(invHazardDist(-log(l[1]) /exp(x) ), rev(invHazardDist(-log(l[3]) /exp(x))) ), border=NA, col="#88888888")
#axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
#mtext(side=4, line=2.5, "Time")
mtext(side=3, at = log(-log(l)/hazardDist(par("usr")[4])), text=paste(100*l, "%", sep=""), cex=.66)
par(xpd=NA)
mtext(side=3, "Survival",line=2, cex=.66)
mtext(side=2, line=2.5, "Time",cex=.66)

par(mar=c(4,.5,1,0.5), xpd=FALSE)
plot(h, seq_along(h), pch=NA, xlab="Total log hazard", yaxs="i", yaxt="n", ylab='', xlim=range(h))
abline(v=seq(-2,4,2), col="grey")
segments(h, seq_along(h),0, seq_along(h))
#mtext(side=3, "Total log hazard", cex=.66, line=.5)
par(xpd=NA)
x <- seq(-3,5,l=100)
lines(x , dnorm(x, 0, sd(h))*100 /  dnorm(0, 0, sd(h)) + length(h)*1.01)


#' #### Stars
#+ stars, fig.width=12, fig.height=12
library(HilbertVis)
nStars <- 32
for(l in list(c("efsTD","partRiskEfsTD"),c("os","partRiskOs"),c("efsTD","partRiskEfsTDsimple"),c("os","partRiskOsSimple"))){
	t <- get(l[1])
	p <- get(l[2])
	locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
	s <- sample(nrow(p),nStars^2) #1:(nStars^2)
	h <- hclust(dist(p[s,]))
	x <- p - rep(colMeans(p), each=nrow(p))
	x <- x/(2*sd(x)) + 1
	if(l=="efsTD")
		c <- cut((t[s,2]-t[s,1])[h$order], quantile(t[,2]-t[,1], seq(0,1,0.1), na.rm=TRUE))
	else
		c <- cut(t[s,1][h$order], quantile(t[,1], seq(0,1,0.1), na.rm=TRUE))
	stars(x[s,][h$order,]/2, scale=FALSE, locations=locations, key.loc=c(0,-3), col.lines=rep(1,(nStars^2)), col.stars = (brewer.pal(11,'RdBu'))[c])
	symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=1)
}

#' #### Summary of stars
#' Overlaying all stars one appreciates the variance in each component.
#+ starsOverlay, fig.width=3, fig.height=3
x <- PartialRisk(coxRFXFitOs)
x <- x[1:nrow(dataFrame),]
x <- x - rep(colMeans(x), each=nrow(x))
t <- sd(x)
x <- x/(2*t) + 1 
q <- apply(x, 2, quantile, c(0.1,0.5,0.9))
rotation <- function(theta) matrix(c(cos(theta), sin(theta),  -sin(theta),cos(theta)), ncol=2)
j <- apply(x, 2, function(y) {v <- violinJitter(y)$y; v/max(v)})/5
c <- cos(seq(0,2*pi, l=ncol(x)+1))[-(ncol(x)+1)]
s <- sin(seq(0,2*pi, l=ncol(x)+1))[-(ncol(x)+1)]
par(bty="n", xpd=NA)
plot(NA,NA, xlim=c(-1,1)*max(x), ylim=c(-1,1)*max(x), xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:nrow(dataFrame)){
	polygon(x[i,]*c, x[i,]*s, col=NA, border="#DDDDDDAA")
}
symbols(rep(0,3),rep(0,3),circles=log(c(0.66,1,1.5))/(2*t)+1, inches=FALSE, add=TRUE, fg="white", lwd=c(1,2,1))
points(t(x) * c - t(j) *s, t(x) * s + t(j) * c, pch=16, cex=.25, col=col1[colnames(x)])
segments(q[1,]*c, q[1,]*s,q[3,]*c,q[3,]*s, lwd=2)
m <- apply(x,2,max)
for(i in seq_along(m))
	text(colnames(x)[i],x=(m[i] + strwidth(colnames(x)[i], units="user")/1.5)*c[i], y=(m[i] + strwidth(colnames(x)[i], units="user")/1.5)*s[i], srt = ((360/ncol(x) *(i-1)+90) %% 180) -90 )
text(log(c(0.66,1,1.5))/(2*t)+1, c(0,0,0)-.1,labels=c(0.66,1,1.5), cex=.33)

#' #### Harrel's C
#library(Hmisc)
#' EFS
totalRiskEfsTD <- rowSums(partRiskEfsTDsimple)
survConcordance( efsTD~totalRiskEfsTD)
predictiveRiskEfsTD <- rowSums(partRiskEfsTDsimple[,-which(colnames(partRiskEfsTDsimple) %in% c("Treatment","GeneTreat","CytoTreat","Trial"))])
survConcordance( efsTD~predictiveRiskEfsTD)
#' OS
totalRiskOs <- rowSums(partRiskOsSimple)
survConcordance( os~totalRiskOs)

#' #### Genomic risk 
whichGenomic <- setdiff(colnames(partRiskEfsTDsimple),c("Treatment","GeneTreat","CytoTreat","Clinical"))
#' EFS
genomicRiskEfsTD <- rowSums(partRiskEfsTDsimple[,whichGenomic])
survConcordance( efsTD~genomicRiskEfsTD)
#' OS
genomicRiskOs <- rowSums(partRiskOsSimple[,whichGenomic])
survConcordance( os~genomicRiskOs)


#' #### AUC
library(survAUC)
s <- efs
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], genomicRiskEfsTD[!is.na(s)], seq(0,1000,10)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskEfsTD[!is.na(s)], seq(0,1000,10)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)

#' #### Genomic risk groups
#' EFS
genomicRiskGroupsEfs <-  cut(genomicRiskEfsTD[!is.na(genomicRiskEfsTD)], quantile(genomicRiskEfsTD, seq(0,1,0.25)), labels = c("very low","low","high","very high"))
survConcordance( efsTD~as.numeric(genomicRiskGroupsEfs))
plot(survfit(efsTD[!is.na(genomicRiskEfsTD)] ~ genomicRiskGroupsEfs), col=brewer.pal(4,"Spectral")[4:1], mark=3, cex=.5)
table(clinicalData$M_Risk, genomicRiskGroupsEfs[1:nrow(clinicalData)])[c(2,3,4,1),]

#' OS
genomicRiskGroupsOs <-  cut(genomicRiskOs[!is.na(genomicRiskOs)], quantile(genomicRiskOs, seq(0,1,0.25)), labels = c("very low","low","high","very high"))
survConcordance( os~as.numeric(genomicRiskGroupsOs))
plot(survfit(os[!is.na(genomicRiskOs)] ~ genomicRiskGroupsOs), col=brewer.pal(4,"Spectral")[4:1], mark=3, cex=.5)
table(clinicalData$M_Risk, genomicRiskGroupsOs[1:nrow(clinicalData)])[c(2,3,4,1),]

#' Risk Plots
#+ genomicRisk, fig.width=6, fig.height=6
par(mfrow=c(4,1))
i <- 0
for(l in levels(clinicalData$M_Risk)[c(2,3,4,1)]){
	barplot(sapply(split(as.data.frame(partRiskOsSimple[1:nrow(clinicalData),whichGenomic][clinicalData$M_Risk ==l,]), genomicRiskGroupsOs[1:nrow(clinicalData)][clinicalData$M_Risk ==l] ), colMeans), beside=TRUE, legend=i==4, col=col1[whichGenomic], main=l, xlim=c(1,45))
	i <- i+1
}


#' Partial values of Harrel's C
c <- PartialC(coxRFXFitEfsTD, groups=simpleGroupsEfs)
b <- barplot(c[1,], col=col1[colnames(c)], las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)


#' 4. Cross-validation
#' ----------------
#' Train only on 2/3 of the data
set.seed(42)
trainIdx <- sample(c(TRUE,FALSE), nrow(dataFrame), replace=TRUE, prob=c(0.66,0.34))
trainIdxEfsTD <- trainIdx[tplSplitEfs]
trainIdxOsTD <- trainIdx[tplSplitOs]

#' ### 1. Static random effects model
#+ coxRFXFitTrain, cache=TRUE
coxRFXFitEfsTrain <- CoxRFX(dataFrame[trainIdx,whichRFXEfs], efs[trainIdx], groups=groups[whichRFXEfs], sigma0 = 0.1, nu=0)
coxRFXFitOsTrain <- CoxRFX(dataFrame[trainIdx,whichRFXOs], os[trainIdx], groups=groups[whichRFXOs], sigma0 = 0.1, nu=0)

#' ### 2. Time-dependent random effects model
#+ coxRFXFitEfsTDTrain, cache=TRUE
coxRFXFitEfsTDTrain <- CoxRFX(dataFrameEfsTD[trainIdxEfsTD,whichRFXEfs], efsTD[trainIdxEfsTD], groups=groups[whichRFXEfs], sigma0 = 0.1, nu=0)

#' Partial contributions
partialRiskOsTest <- PartialRisk(coxRFXFitOsTrain, newX=dataFrame[!trainIdx, whichRFXOs])
c <- PartialC(coxRFXFitOsTrain, newX = dataFrame[!trainIdx, whichRFXOs], newSurv = os[!trainIdx])
b <- barplot(c[1,], col=col1[colnames(c)], las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' Overall
totalRiskOsTest <- rowSums(partialRiskOsTest)
survConcordance(os[!trainIdx]~totalRiskOsTest )

#' Compared to molecular risk
predictiveRiskTest <- rowSums(partialRiskOsTest[,-which(colnames(partRiskOsSimple) %in% c("Treatment","GeneTreat","CytoTreat","Trial"))])
survConcordance(os[!trainIdx] ~ predictiveRiskTest)

#barplot(c(CGP=survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)$concordance, MolecularRisk = survConcordance(survivalTD[!testIxTD] ~ c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[tplSplit][!testIxTD]])[[1]]))

s <- efs[!trainIdx]
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTest[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskOsTest[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' ### 3. Recursive partitioning & randomForestSRC on OS
#+ tree, fig.width=5, fig.height=5
library(rpart)
library(randomForestSRC)
tree <- rpart(os ~ ., data=dataFrame[mainIdxOs & osIdx])
plot(tree)
text(tree)
survConcordance(na.omit(os)~predict(tree))

#+ treeCV, fig.width=5, fig.height=5
treeCV <- rpart(os[trainIdx] ~ ., data=dataFrame[trainIdx,mainIdxOs & osIdx])
plot(treeCV)
text(treeCV)
survConcordance(os[!trainIdx]~predict(treeCV, newdata=dataFrame[!trainIdx, mainIdxOs & osIdx]))

#+ rForest, cache=TRUE
rForest <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs & osIdx]), ntree=100)
boxplot(rForest$importance ~ factor(as.character(groups[mainIdxOs & osIdx])), border= col1, staplewex=0, pch=16, cex=0.75, ylab="RSF importance", lty=1, xaxt="n")
rForestVimp <- sapply(mainGroups, function(g) vimp(rForest, colnames(dataFrame)[which(groups==g)]))

survConcordance(na.omit(os)~predict(rForest, importance="none")$predicted)

#+ rForestCV, cache=TRUE
rForestCV <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs & osIdx])[trainIdx,], ntree=100, importance="none")
p <- predict(rForestCV, newdata = dataFrame[!trainIdx,mainIdxOs & osIdx], importance="none")
survConcordance(os[!trainIdx]~ p$predicted)


#' ### 4. Stability selection: All terms

#' Fit model
#+ coxCPSS, cache=TRUE
set.seed(42)
coxCPSSEfs <- CoxCPSS(dataFrame[!is.na(efs),efsIdx], efs[!is.na(efs)], bootstrap.samples = 50, control="BH", level=.1)

plot(coxCPSSEfs)

selectedEfs <- names(which(coxCPSSEfs$Pi > .8))
selectedEfs
coxFitEfs <- coxph(as.formula(paste("efs[trainIdx] ~ ", paste(selectedEfs, collapse="+"))), data=dataFrame[trainIdx, ])
summary(coxFitEfs)

totalRiskStabTest <- predict(coxFitEfs, newdata=dataFrame[!trainIdx,])
survConcordance(efs[!trainIdx] ~ totalRiskStabTest)


#' ### 5. Stability selection with interactions
#' Only include interaction terms when main effects are present. 
#' On the other hand, the exclusive selection of a product term could mean that the most likely explanation is that the two main terms are zero and only the interaction is non-zero..
CoxCPSSInteractions

#' #### Fit model to entire data set
#+ CoxCPSSInt, cache=TRUE
set.seed(42)
coxCPSSIntEfs <- CoxCPSSInteractions(dataFrame[!is.na(efs),groups %in% mainGroups & efsIdx], na.omit(efs), bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedIntEfs <- names(which(coxCPSSIntEfs$Pi > 0.8))
#selectedIntEfs[selectedIntEfs=="NPM1:TPL"] <- "TPL:NPM1"

#' Concordance of the resulting static coxph model
coxFitEfs <- coxph(as.formula(paste("efs[trainIdx] ~ ", paste(selectedIntEfs, collapse="+"))), data=dataFrame[trainIdx, ])
summary(coxFitEfs)
survConcordance(efs[!trainIdx] ~ predict(coxFitEfs, newdata = dataFrame[!trainIdx, ]))

#' Concordance of the corresponding time-dependent coxph model
coxFitEfsTD <- coxph(as.formula(paste("efsTD[trainIdxEfsTD] ~ ", paste(selectedIntEfs, collapse="+"))), data=dataFrameEfsTD[trainIdxEfsTD, ])
summary(coxFitEfsTD)
survConcordance(efsTD[!trainIdxEfsTD] ~ predict(coxFitEfsTD, newdata = dataFrameEfsTD[!trainIdxEfsTD, ]))

#' KM plots of interaction terms
#+ CoxCPSSIntPlot
for(i in grep(":", selectedIntEfs, value=TRUE)){
	j <- strsplit(i, ":")[[1]]
	plot(survfit(as.formula(paste("efsTD ~ ", paste( j, collapse="+"))), data=dataFrameEfsTD), col=set1)
	legend("topright", lty=1, col=set1, c("none", rev(j), "both"), bty="n")
}

#' Plots of the stability selection: No large difference between when interaction terms are enforce to only occur after main effects are included
plot(coxCPSSIntEfs)
plot(coxCPSSIntEfs$Pi, coxCPSSEfs$Pi[match(names(coxCPSSIntEfs$Pi), names(coxCPSSEfs$Pi))])

#' Partial risk for c
names(groups) <- colnames(dataFrame)
selectedIntEfs <- c(selectedIntEfs[selectedIntEfs %in% names(dataFrameEfsTD)],
		unlist(sapply(strsplit(selectedIntEfs[!selectedIntEfs %in% names(dataFrame)], ":"), function(x) paste(rev(x), collapse=":"))))
X <- as.matrix(dataFrameEfsTD[selectedIntEfs])
X <- X - rep(colMeans(X), each=nrow(X))
partRiskCoxCPSSInt <- sapply(levels(groups), function(l) {
			w <- groups[selectedIntEfs] == l
			if(length(w)>0)
				X[,w, drop=FALSE] %*% coef(coxFitEfsTD)[w]
			else
				rep(0, nrow(X))
		})
			
#' Genomic risk 
whichGenomicCoxCPSSInt <- which(!colnames(partRiskCoxCPSSInt) %in% c("Treatment","GeneTreat","CytoTreat","Clinical","Trial"))
genomicRiskCoxCPSSInt <- rowSums(partRiskCoxCPSSInt[,whichGenomicCoxCPSSInt])
survConcordance( efsTD~genomicRiskCoxCPSSInt)

#' #### Genomic risk groups
genomicRiskGroupsCoxCPSSInt <-  cut(genomicRiskCoxCPSSInt[!is.na(genomicRiskCoxCPSSInt)], quantile(genomicRiskCoxCPSSInt, seq(0,1,0.25)), labels = c("very low","low","high","very high"))

survConcordance( efsTD~as.numeric(genomicRiskGroupsCoxCPSSInt))
plot(survfit(efsTD[!is.na(genomicRiskCoxCPSSInt)] ~ genomicRiskGroupsCoxCPSSInt), col=brewer.pal(4,"Spectral")[4:1], mark=16)
table(clinicalData$M_Risk, genomicRiskGroupsCoxCPSSInt[1:nrow(clinicalData)])[c(2,3,4,1),]

#' Risk Plots
#+ genomicRiskCoxCPSSInt, fig.width=6, fig.height=6
par(mfrow=c(4,1))
i <- 1
for(l in levels(clinicalData$M_Risk)[c(2,3,4,1)]){
	barplot(sapply(split(as.data.frame(partRiskCoxCPSSInt[1:nrow(clinicalData),whichGenomicCoxCPSSInt][clinicalData$M_Risk ==l,]), genomicRiskGroupsCoxCPSSInt[1:nrow(clinicalData)][clinicalData$M_Risk ==l] ), colMeans), beside=TRUE, legend=i==4, col=col1[whichGenomicCoxCPSSInt], main=l, xlim=c(1,30))
	i <- i+1
}

#' #### Fit model to the training subset only
#+ CoxCPSSIntEfsCV, cache=TRUE
coxCPSSIntEfsCV <- CoxCPSSInteractions(dataFrame[!is.na(efs) & trainIdx,groups %in% mainGroups & efsIdx], na.omit(efs[trainIdx]), bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedIntEfsCV <- names(which(coxCPSSIntEfsCV$Pi > 0.8))
survConcordance(efs[!trainIdx] ~ predict(coxph(as.formula(paste("efs[trainIdx] ~", paste(selectedIntEfsCV, collapse="+"))), data=dataFrame[trainIdx, ]), newdata = dataFrame[!trainIdx, ]))

#tmp <- CoxRFX(dataFrame[trainIdx, selectedIntEfs], efs[trainIdx], nu=0)
#survConcordance(efs[!trainIdx] ~ as.matrix(dataFrame[!trainIdx, selectedIntEfs]) %*% coef(tmp))

#' #### CPSS on OS
#+ CoxCPSSIntOs, cache=TRUE
set.seed(42)
coxCPSSIntOs <- CoxCPSSInteractions(dataFrame[!is.na(os),groups %in% mainGroups & osIdx], na.omit(os), bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedIntOs <- names(which(coxCPSSIntOs$Pi > 0.8))

#' The corresponding coxph model (CV)
coxFitCPSSIntOsCV <- coxph(as.formula(paste("os[trainIdx] ~", paste(selectedIntOs, collapse="+"))), data=dataFrame[trainIdx, ])
summary(coxFitCPSSIntOsCV)
survConcordance(os[!trainIdx] ~ predict(coxFitCPSSIntOsCV, newdata = dataFrame[!trainIdx, ]))

#' The corresponding coxph model
coxFitCPSSIntOs <- coxph(as.formula(paste("os ~", paste(selectedIntOs, collapse="+"))), data=dataFrame)
riskCPSSIntOs <- predict(coxFitCPSSIntOs, newdata=dataFrame)
summary(coxFitCPSSIntOs)

#' With cross-val
#+ CoxCPSSIntOsCV, cache=TRUE
coxCPSSIntOsCV <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx,groups %in% mainGroups & osIdx], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedIntOsCV <- names(which(coxCPSSIntOsCV$Pi > 0.8))
coxFitCPSSIntOsCV <- coxph(as.formula(paste("os[trainIdx] ~", paste(selectedIntOsCV, collapse="+"))), data=dataFrame[trainIdx, ])
survConcordance(os[!trainIdx] ~ predict(coxFitCPSSIntOsCV, newdata = dataFrame[!trainIdx, ]))

#' Visualisation of risk
#+ riskCPSSIntOs, fig.width=3.5, fig.height=2.5
risk <- riskCPSSIntOs[dataFrame$TPL_os==0]
x <- seq(from=-4,to=4, l=512)
r <- sapply(levels(clinicalData$M_Risk)[c(2,3,4,1)], function(r){
			i <- clinicalData$M_Risk[tplSplitOs]==r
			d <- density(na.omit(risk[i]), from=-4,to=4)$y * mean(i, na.rm=TRUE)
		})
par(mar=c(4,4,3,4)+.1, bty="n")
plot(exp(x),rowSums(r), type='l', lty=0,xlab="Hazard", ylab="Prop. patients", log='x')
for(i in 1:4)
	polygon(exp(c(x, rev(x))), c(rowSums(r[,1:i, drop=FALSE]), rev(rowSums(cbind(0,r)[,1:i, drop=FALSE]))), col=set1[c(3,2,4,1)][i], border=NA)

H0 <- basehaz(coxFitCPSSIntOs, centered=TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
for(i in seq_along(l))
lines(exp(x), invHazardDist(-log(l[i]) /exp(x) )/10000, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
mtext(side=4, "Time", line=2.5)
mtext(side=3, at = -log(l)/hazardDist(par("usr")[4]*10000), text=paste(100*l, "% survive", sep=""))
legend("topright", levels(clinicalData$M_Risk)[c(2,3,4,1)], fill=set1[c(3,2,4,1)], bty="n", title="M risk")


#' ### 6. Stepwise model selection

#' #### BIC
#' Whole data set
#+ coxBIC, cache=TRUE, warning=FALSE
c <- coxph(os ~ 1, data=dataFrame[,mainIdxOs & osIdx])
scopeStep <- as.formula(paste("os ~", paste(colnames(dataFrame)[mainIdxOs& osIdx], collapse="+")))
coxBICOs <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
summary(coxBICOs)

#' Training subset
#+ coxBICTrain, cache=TRUE, warning=FALSE
c <- coxph(os[trainIdx] ~ 1, data=dataFrame[trainIdx,mainIdxOs& osIdx])
scopeStep <- as.formula(paste("os ~", paste(colnames(dataFrame)[mainIdxOs& osIdx], collapse="+")))
coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
summary(coxBICOsTrain)
survConcordance(os[!trainIdx] ~ predict(coxBICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]))

#' #### AIC
#' Whole data sets
#+ coxAIC, cache=TRUE, warning=FALSE
coxAICOs <- step(c, scope= scopeStep, k = 2, trace=0)
summary(coxAICOs)

#' Training subset
#+ coxAICTrain, cache=TRUE, warning=FALSE
coxAICOsTrain <- step(c, scope= scopeStep, k = 2, trace=0)
summary(coxAICOsTrain)
survConcordance(os[!trainIdx] ~ predict(coxAICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]))

#' ### 7. Summary of different models
#' #### Static models
predictedRiskCV <- data.frame(
		stdRisk = c(4,1,2,3)[clinicalData$M_Risk[!trainIdx]],
		tree = predict(treeCV, newdata=dataFrame[!trainIdx, mainIdxOs]),
		rForest = predict(rForestCV, newdata = dataFrame[!trainIdx,mainIdxOs], importance="none")$predicted,
		coxRFX = as.matrix(dataFrame[!trainIdx,whichRFXOs]) %*% coef(coxRFXFitOsTrain),
		coxBIC = predict(coxBICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]),
		coxAIC = predict(coxAICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]),
		coxCPSS = predict(coxFitCPSSIntOsCV, newdata = dataFrame[!trainIdx, ])
)

#+ concordanceCV
concordanceCV <- sapply(predictedRiskCV, function(x) {c <- survConcordance(os[!trainIdx] ~ x); c(c$concordance, c$std.err)})
concordanceCV
o <- order(concordanceCV[1,])
barplot(concordanceCV[1,o], border=NA, col= set1, las=2, xaxt="n", ylab="Concordance") -> b
segments(b,concordanceCV[1,o]-concordanceCV[2,o],b,concordanceCV[1,o]+concordanceCV[2,o])
rotatedLabel(b, rep(0,length(b)), colnames(concordanceCV)[o], srt=45)

#+ aucCV
library(survivalROC)
aucCV <- sapply(predictedRiskCV, function(x) survivalROC(Stime=efs[!is.na(efs) & !trainIdx,1], status=efs[!is.na(efs) &  !trainIdx,2], marker = x[!is.na(efs[!trainIdx])], predict.time = 278, method="KM", cut.values=seq(-5,5,0.1))$AUC)
aucCV
barplot(aucCV, border=NA, col= set1, las=2, xaxt="n", ylab="AUC (median EFS)") -> b
rotatedLabel(b, rep(0,length(b)), names(aucCV), srt=45)



#' 6. TCGA validation
#' ------------------
#' ### Load data
tcgaClinical <- read.table("../data/TCGA_clin.txt", sep="\t", header=TRUE)
tcgaGenetic <- read.table("../data/TCGA_gen.txt", sep="\t", header=TRUE)
tcgaGenetic$TCGA_ID <- factor(as.character(tcgaGenetic$TCGA_ID), levels = levels(tcgaClinical$TCGA_ID))
g <- as.character(tcgaGenetic$Hugo_Symbol)
g[tcgaGenetic$Hugo_Symbol=="FLT3" & tcgaGenetic$Variant_Type == 'INS'] <- "FLT3_ITD"
g[tcgaGenetic$Hugo_Symbol=="FLT3" & tcgaGenetic$Variant_Type == 'SNP'] <- "FLT3_TKD"
tcgaMutation <- (table(tcgaGenetic$TCGA_ID,g)) + 0
t <- data.frame(tcgaMutation[,]>0, CEBPA_mono = tcgaMutation[,"CEBPA"]==1,CEBPA_bi = tcgaMutation[,"CEBPA"]>1,tcgaClinical[,-c(1,2,4,5,6,12,24)], MakeInteger(tcgaClinical$TypeAML)) + 0
w <- grep("_10+$",colnames(dataFrame), value=TRUE)
f <- as.numeric(sub(".+_","",w))
n <- sub("_10+","",w)
f <- f[n %in% colnames(tcgaClinical)]
n <- n[n %in% colnames(tcgaClinical)]
t[n] <- t[n] / rep(f, each=nrow(t))
colnames(t)[match(n,colnames(t))] <- paste(n,f,sep="_")
rm(w,n,f,g)

tcgaData <- dataFrame[1:nrow(t),]
tcgaData[,] <- NA
w <- intersect(names(t), names(tcgaData))
tcgaData[w] <- t[w]
tcgaData$TPL_os <- tcgaData$TPL_efs <- NA
tcgaData[groups=="Genetics"][is.na(tcgaData[groups=="Genetics"])] <- 0
rm(t,w)
tcgaSurvival <- Surv(tcgaClinical$OS, tcgaClinical$Status)

#' ### Plot mutation frequencies
plot(colMeans(dataFrame[groups=="Genetics"])[colnames(tcgaMutation)], colMeans(tcgaMutation))
text(colMeans(dataFrame[groups=="Genetics"])[colnames(tcgaMutation)], colMeans(tcgaMutation), colnames(tcgaMutation))
cor(colMeans(dataFrame[groups=="Genetics"])[colnames(tcgaMutation)], colMeans(tcgaMutation), use='c')

#' ### NPM1 survival
plot(survfit(tcgaSurvival ~ NPM1, data=tcgaData), col=set1[1:2])
lines(survfit(os ~ NPM1, data=dataFrame), col=set1, lty=3,mark=NA)
legend("topright", col=c(set1[1:2],"black","black"), c("NPM1 wt", "NPM1 mut","TCGA","AML"), lty=c(1,1,1,3), bty='n')

#' ### Analyse risk
#' #### CoxRFX model and covariance based imputation
tcgaRiskRFXOs <- PredictRiskMissing(coxRFXFitOs, tcgaData[whichRFXOs])
survConcordance(tcgaSurvival ~ tcgaRiskRFXOs[,1])

#' #### CPSS model
tcgaDataImputed <- as.data.frame(ImputeXMissing(dataFrame[mainIdxOs], newX=tcgaData[mainIdxOs]))
tcgaRiskCPSSOs <- predict(coxFitCPSSIntOs, newdata=tcgaDataImputed)
survConcordance(tcgaSurvival ~ tcgaRiskCPSSOs)

#' Blind imputation (mean only)
f <- function(X) {X <- sapply(X, poorMansImpute);X[is.na(X)] <- 0; X}
survConcordance(tcgaSurvival ~ predict(coxFitCPSSIntOs, newdata=as.data.frame(f(tcgaData[mainIdxOs]))))

#' #### Cytogenetic risk
survConcordance(tcgaSurvival ~ c(3,1,2)[tcgaClinical$C_Risk])

#' #### PINA score
PINAOs <- function(X){
	coef <- c( NPM1=-1.2,
			FLT3_ITD=-.26,
			`NPM1:FLT3_ITD`=.89,
			CEBPA_bi=-1.3,
			wbc_log10=.57,
			age=0.044,
			ecog24=.4)
	x <- cbind(X[,colnames(X) %in% names(coef)], wbc_log10 = log10(100*1e3*pmax(X[,"wbc_100"], 0.001)), age = X[,"AOD_10"]*10, ecog24 = X[,"Performance_ECOG"]>=2)
	risk <- as.matrix(x[,names(coef)]) %*% coef
	group <- cut(risk, c(min(risk), 4,5.4, max(risk)), labels = c("low","int","high"))
	return(data.frame(risk, group))
}
pinaOs <- PINAOs(dataFrame)

#+ PINAos
nkIdx <- clinicalData$NK == 1
plot(survfit(os[nkIdx] ~ pinaOs[nkIdx,2]), col=rev(set1[1:3]))
survConcordance(os[nkIdx] ~ pinaOs[nkIdx,1])

#' Compared to CPSS (AML data)
survConcordance(os[nkIdx] ~ predict(coxFitCPSSIntOs, newdata=dataFrame)[nkIdx])
survConcordance(os[nkIdx&!trainIdx] ~ predict(coxFitCPSSIntOsCV, newdata=dataFrame)[nkIdx&!trainIdx])

#' And on TCGA data
tcgaPinaOs <- PINAOs(cbind(tcgaDataImputed, `NPM1:FLT3_ITD` = tcgaDataImputed[,"NPM1"]*tcgaDataImputed[,"FLT3_ITD"]))
tcgaNkIdx <- tcgaClinical$karyotype == "Normal"
survConcordance(tcgaSurvival[tcgaNkIdx] ~ tcgaPinaOs[tcgaNkIdx,1])
survConcordance(tcgaSurvival[tcgaNkIdx] ~ tcgaRiskCPSSOs[tcgaNkIdx])


#' #### Other models
tcgaRisk <- data.frame(
		stdRisk = c(3,1,2)[tcgaClinical$C_Risk],
		tree = predict(tree, newdata=tcgaDataImputed),
		rForest = predict(rForest, newdata = tcgaDataImputed, importance="none")$predicted,
		PINAos = tcgaPinaOs[,1],
		coxRFX = tcgaRiskRFXOs[,1],
		coxBIC = predict(coxBICOs, newdata=tcgaDataImputed),
		coxAIC = predict(coxAICOs, newdata=tcgaDataImputed),
		coxCPSS = tcgaRiskCPSSOs
)

#+ concordanceTCGA
tcgaConcordance <- sapply(tcgaRisk, function(x) {c <- survConcordance(tcgaSurvival ~ x); c(c$concordance, c$std.err)})
tcgaConcordance
o <- order(tcgaConcordance[1,])
barplot(tcgaConcordance[1,o], border=NA, col= set1, las=2, xaxt="n", ylab="Concordance") -> b
segments(b,tcgaConcordance[1,o]-tcgaConcordance[2,o],b,tcgaConcordance[1,o]+tcgaConcordance[2,o])
rotatedLabel(b, rep(0,length(b)), colnames(tcgaConcordance)[o], srt=45)

#+ kmTCGA
risk <- cut(tcgaRiskCPSSOs, quantile(tcgaRiskCPSSOs), labels=c("low","int1","int2","high"))
s <- survfit(tcgaSurvival ~ risk)
plot(s, col=set1[c(3,2,4,1)], mark=NA)
legend("topright", bty="n", rownames(summary(s)$table), col=set1[c(3,2,4,1)], lty=1)

#+ riskTCGA, fig.width=3.5, fig.height=2.5
risk <- tcgaRiskCPSSOs - mean(tcgaRiskCPSSOs)
x <- seq(from=-4,to=4, l=512)
r <- sapply(levels(clinicalData$M_Risk)[c(2,3,4,1)], function(r){
			i <- clinicalData$M_Risk[tplSplitOs]==r
			d <- density(na.omit(risk[i]), from=-4,to=4)$y * mean(i, na.rm=TRUE)
		})
par(mar=c(4,4,3,4)+.1, bty="n")
plot(exp(x),rowSums(r), type='l', lty=0,xlab="Hazard", ylab="Prop. patients", log='x')
for(i in 1:4)
	polygon(exp(c(x, rev(x))), c(rowSums(r[,1:i, drop=FALSE]), rev(rowSums(cbind(0,r)[,1:i, drop=FALSE]))), col=set1[c(3,2,4,1)][i], border=NA)

H0 <- basehaz(coxph(tcgaSurvival ~ risk), centered=TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
for(i in seq_along(l))
	lines(exp(x), invHazardDist(-log(l[i]) /exp(x) )/10000, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
mtext(side=4, "Time", line=2.5)
mtext(side=3, at = -log(l)/hazardDist(par("usr")[4]*10000), text=paste(100*l, "% survive", sep=""))
legend("topright", levels(clinicalData$M_Risk)[c(2,3,4,1)], fill=set1[c(3,2,4,1)], bty="n", title="M risk")


#' 7. Miscellania
#' --------------
#' ### 1. Test of hotspot mutations other than above
#+ coxRFXHotspots, cache=TRUE
h <- paste(mutationData$GENE, sub("([0-9])([A-Z])$","\\1",mutationData$AA_CHANGE), sep="_")
t <- table(h)
i <- mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & mutationData$FINAL_CALL == "OK" & h %in% names(t[t>=5])
hotspots <- data.frame(table(mutationData$SAMPLE_NAME[i], h[i])[,]) ## Table of hotspots
genes <- data.frame(mutationTable[,colSums(mutationTable)>0]) ## The original gene table
g <- c(rep("genes", ncol(genes)), rep("hotspots", ncol(hotspots))) ## groups
coxRFXHotspots <- CoxRFX(cbind(genes, hotspots), os, g, sigma0 = 0.1, nu=0, which.mu=1)
plot(coxRFXHotspots)

#' Basic CV
#+ coxRFXHotspotsCV
t <-  CoxRFX(cbind(genes, hotspots)[trainIdx,], os[trainIdx], g, sigma0 = 0.1, nu=0, which.mu=1)
survConcordance( os[!trainIdx] ~ (as.matrix(cbind(genes, hotspots)) %*% coef(t))[!trainIdx])
t0 <-  CoxRFX(genes[trainIdx,], os[trainIdx], sigma0 = 0.1, nu=0)
survConcordance( os[!trainIdx] ~ (as.matrix(genes) %*% coef(t0))[!trainIdx])

#' ### 2. Most frequent interaction terms

#' #### Second order
x <- cbind(mutationTable, dataList$Cytogenetics)
t <- unlist(sapply(2:ncol(x), function(j) sapply(1:(j-1), function(k) {s <- sum(x[,j]*x[,k], na.rm=TRUE); names(s) <- paste(colnames(x)[c(j,k)], collapse=":"); return(s)})))
head(sort(t, decreasing=TRUE), 50)

#' #### Third order
u <- unlist(sapply(3:ncol(x), function(i) sapply(2:(i-1), function(j) sapply(1:(j-1), function(k) {s <- sum(x[,i]*x[,j]*x[,k], na.rm=TRUE); names(s) <- paste(colnames(x)[c(i,j,k)], collapse=":"); return(s)}))))
head(sort(u, decreasing=TRUE), 50)

#' ### CPSS on OS for non WHO only
#+ CoxCPSSIntOsNonWHO, cache=TRUE
set.seed(42)
w <- !is.na(os) & clinicalData$WHOcat == "no" &! is.na( clinicalData$WHOcat)
coxCPSSIntOsNonWHO <- CoxCPSSInteractions(dataFrame[w,groups %in% mainGroups & osIdx], os[w], bootstrap.samples=50, scope = which(groups %in% c("Genetics","Cytogenetics","Treatment")))
selectedIntOsNonWHO <- names(which(coxCPSSIntOsNonWHO$Pi > 0.8))

#' The corresponding time-dep coxph model
coxFitIntOsNonWHO <- coxph(as.formula(paste("os ~", paste(selectedIntOsNonWHO, collapse="+"))), data=dataFrame)
summary(coxFitIntOsNonWHO)


#' ### 3. Very time-dependent model
#' Attempting to rescue TPL in OS.
#### More time-dependence
rdIdx <- clinicalData$ereignart %in% levels(clinicalData$ereignart)[3:4]
rdTime <- clinicalData$EFS
rdTime[!rdIdx] <- NA

d <- sapply(1:nrow(clinicalData), function(i){
			i <<- i
			t <- c(as.numeric(clinicalData[i,c("CR_date","TPL_date","Date_LF","Recurrence_date")]) - as.numeric(clinicalData$ERDate[i]), rdTime[i])
			o <- order(t, na.last=NA)
			r <- rank(t)
			r[is.na(t)] <- length(t)+1
			names(r) <- c("CR","TPL","LF", "REC","RD")
			tt <- c(0,t[o])
			if(length(o)==0)
				return(c(rep(NA,7),i))
			s <- cbind(time1 = tt[-length(tt)], time2=tt[-1], event=c(rep(0, length(o)-1), clinicalData$Status[i]), outer(0:(length(o)-1), r[-3], `>=`)+0, i=i)[diff(tt)!=0,]
			return(s)
		})
d <- as.data.frame(do.call("rbind",d))

#' A coxph model
summary(coxph(Surv(time1, time2, event) ~ ., data=d[,-8]))

#' With CPSS selected terms
summary(coxph(Surv(time1, time2, event) ~ ., data=cbind(d[,-8], dataFrame[d$i, grep("(TPL)|(:)", selectedIntOs, value=TRUE, invert = TRUE)])))

#' 8. Regression of blood counts
#' -----------------------------
#' ### Prepare data
#+ clinicalGlmnet, cache=TRUE
library(glmnet)
Y <- StandardizeMagnitude(dataList$Clinical)
X <- as.matrix(dataFrame[groups %in% c("Genetics","Cytogenetics")])
set.seed(42)
clinModels = lapply(Y, function(y){
			if (class(y) %in% c("numeric","integer")){
				if(all(y %in% c(0,1,NA)))
					cv.glmnet(X[!is.na(y),], na.omit(y), family = "binomial", alpha=1, standardize=FALSE, nfolds=5)
				else if(all(y %in% c(0:100,NA)))
					cv.glmnet(X[!is.na(y),], na.omit(y), family = "poisson", alpha=1, standardize=FALSE, nfolds=5)
				else
					cv.glmnet(X[!is.na(y),], na.omit(y), family = "gaussian", alpha=1, standardize=FALSE, nfolds=5)
			}
			else if (class(y)=="factor")
				cv.glmnet(X[!is.na(y),], na.omit(y), family="multinomial",  alpha=1, standardize=FALSE, nfolds=5)
		})

#' ### Plot
#+ clinicalGlmnetLolly, fig.width=7.5, fig.height=3.75
par(bty="n", mgp = c(2.5,.5,0), mar=c(3,4,2,4)+.1, las=2, tcl=-.25)
i = 1
n <- colnames(X)
annot <- 1 + grepl("^[A-Z]",n) 
names(annot) <- n
for(m in clinModels){
	plotcvnet(m, X, main=names(clinModels)[i],  col0="black", cex=1, simple.annot = annot, col=set1[c(3,2,4)], xlim=c(0.5,35.5))
	i = i+1
	legend("topright", col=c(set1[c(1,3)],"black")[c(1,3,2)], c(expression(paste("Explained variance ",R^2)), expression(paste("Lasso penalty ",lambda)), expression(paste("Model coefficient ", beta))), box.lty=0, bg="#FFFFFF33", pch=c(NA,NA,19), lty=c(1,1,NA), cex=.8, pt.cex = 1)
}

#' ### Heatmap of GLMs
j <- 0
z <- sapply(clinModels,function(x){ ## Creating matrix
			j <<- j+1
			w <- which.min(x$cvm)
			c <- x$glmnet.fit$beta[,w]
			yj <- sapply(c("Genetics","Cytogenetics"), function(i){
						w <- names(groups[groups == i])
						X[,w] %*% c[w] 
					})
			cj <- rowSums(cov(yj))
			y <- Y[,j] - x$glmnet.fit$a0[w]
			covj <- colMeans((y-mean(y))*(yj - rep(colMeans(yj), each=nrow(yj))))
			r2 <- cj
			R2 <- 1 - x$cvm[w]/x$cvm[1]
			c(c, NA,  r2/sum(r2)*R2, R2=R2)
		})
r <- sapply(clinModels, function(x) rank(apply(x$glmnet.fit$beta,1, function(y) which(y!=0)[1]), ties="min")) ## R
s <- sapply(clinModels, function(x) {
			w = which.min(x$cvm)
			w <- rev(which(x$cvm[1:w] > x$cvup[w]))[1] +1
			if(!is.na(w))
				sum(x$glmnet.fit$beta[,w]!=0)
			else
				0
		})
p <- sapply(clinModels, function(x) {w <- which.min(x$cvm); (1 - x$cvup[w]/x$cvm[1]) > 0 })
R2 <- z[nrow(z) - 2:0,]
R2[is.na(R2)] <- 0
z <- z[-(nrow(z) - 0:3),]
m <- nrow(z) 
z[1:m,] <- pmin(z[1:m,],0.1)
z[1:m,] <- pmax(z[1:m,],-.0999)
z[z==0] <- NA

#' Plot
#+ clinicalGlmnetHeatmap, fig.width=15, fig.height=3.75
layout(matrix(c(1,2),1,2), c(7.5,1.5), c(2,2), TRUE)
par(bty="n", mgp = c(3,.5,0), mar=c(4,10,2,0)+.1, las=1, tcl=-.25)
## Heatmap
w <- order(R2[nrow(R2),])
o <- order(-annot,rowSums(r, na.rm=TRUE))
image(y=1:ncol(z)-.5, x=1:nrow(z), z[o,w], breaks=c(-2,seq(-0.1,0.1,l=51)), col=c("grey",colorRampPalette(brewer.pal(9,"RdYlBu"))(50)), xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0,ncol(z)+1))
#abline(v=c(18.5,27.5), lwd=0.5)
rotatedLabel(y0=rep(0.5,nrow(z)), labels=gsub("_","/",rownames(z))[o], x0=1:nrow(z), font=c(rep(3,sum(groups=="Genetics")),rep(1,sum(groups=="Cytogenetics"))), col = set1[c(3,2,4)][annot], cex=0.9)
mtext(side=2, line=.2, text=colnames(z)[w], las=2, at=1:ncol(z)-.5)
text(y=rep(1:ncol(z)-.5, each=nrow(r)), x=rep(1:nrow(r), ncol(z)), r[o,w] * (0!=(z[o,w])), cex=0.66, font=ifelse(r[o,w] <= rep(s[ncol(r):1], each=nrow(r)), 2,1))
points(y=rep(1:ncol(z)-.5, each=nrow(r)), x=rep(1:nrow(r), ncol(z)), pch=ifelse(is.na(z[o,w]) | z[o,w]==0, ".",NA))
mtext(side=1, at=sum(groups=="Genetics")/2, "Genetics", col=set1[2], line=2.5 )
mtext(side=1, at=sum(groups=="Cytogenetics")/2 + sum(groups=="Genetics"), "Cytogenetics", col=set1[3], line=2.5 )
## Legends
mtext(side=3, "Model coefficients", at = 7, line=0.5)
clip(-10,50,0,15)
image(y=dim(z)[2] + c(0.5,1.5), x=1+ 1:7 , matrix(seq(-0.99,1,l=7), ncol=1), breaks=c(-2,seq(-1,1,l=51)), col=c("grey",colorRampPalette(brewer.pal(9,"RdYlBu"))(50)), xaxt="n", yaxt="n", xlab="", ylab="", add=TRUE)
text(y=dim(z)[2]+1, x=c(1,9), 0.1*c(-1,1))
points(y=11,x=5, pch=".")
rect(19.5,ncol(z) +.5,20.5, ncol(z) +1.5, lwd=0.5)
text(20,ncol(z) +1,1, cex=0.66)
text(21, ncol(z) +1, "LASSO rank", pos=4)
## Side panel
u <- par("usr")
par(bty="n", mgp = c(3,.5,0), mar=c(4,1,2,1)+.1, las=1, tcl=-.25)
plot(NA,NA, xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0,ncol(z)+1), xlim=c(0,.9), yaxs="i")
barplot(R2[1:2,w], border=NA, col=paste(set1[c(2,3,4)],"88",sep=""), horiz=TRUE, names.arg=rep(NA,ncol(R2)), width=0.95, space=0.0525, add=TRUE) -> b
points(R2[3,w]+0.1,b, pch=ifelse(p[w],"*",NA))
mtext(side=3, "Variance components", line=.5)
mtext(side=1, expression(paste("Explained variance ",R^2)), line=2.5)

#+ save
save(list = ls(), file=paste(Sys.Date(), "-AML.RData", sep=""))



		
#' 8. Germline polymorphisms
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
dataFrame <-  dir("/Volumes/nst_links/live/845/", pattern="WGA*")
allCaveOut <- sapply(rownames(mutationTable), function(s){
			i<<-i+1
			cat(ifelse(i%%100 == 0, "\n","."))
			stem <<- grep(s, dataFrame, value=TRUE)[1]
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
		
#' Session
#' -----
#+ sessionInfo, eval=TRUE
sessionInfo()
		
#' TODO's
#' ------
#' * Number of total interactions in each cat/exp v obs.
#' * Gene-wise contributions to reclassifaction
#' * Reclassification with OS
#' * Missing data
