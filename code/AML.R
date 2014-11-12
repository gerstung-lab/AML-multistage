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
levels(clinicalData$Study) <- c(`_07-04`="AMLSG0704" ,   `98A`="AMLHD98A" ,  `98B`="AMLHD98B")[levels(clinicalData$Study)]
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

#' 2. Preparing data
#' ------------------
tplIdxEfs <- (clinicalData$Time_Diag_TPL < clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL)) | (is.na(clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL))) # Need to restrict to those with transplant prior to event..
tplIdxOs <- !is.na(clinicalData$TPL_type)

#' #### All data as list
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
		Cytogenetics = clinicalData[,50:74],
		Nuisance = data.frame( MakeInteger(clinicalData$Study)[,1:2], Date=scale(as.numeric(clinicalData$ERDate), scale=FALSE)),
		Treatment = data.frame(ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL_efs=tplIdxEfs, TPL_os=tplIdxOs),
		Demographics = clinicalData[,c("AOD","gender")],
		Clinical = cbind(clinicalData[, c("Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH","HB","platelet","Splenomegaly")], MakeInteger(clinicalData$TypeAML)[,-1]))#,
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
dataList$GeneGene$`NPM1:FLT3_ITD:DNMT3A` <- (rowSums(dataList$Genetics[c('NPM1',"FLT3_ITD","DNMT3A")])==3)+0 ## Add NPM1:FLT3_ITD:DNMT3A product term as well
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

groups <- unlist(sapply(names(dataList), function(x) rep(x, ncol(dataList[[x]]))))
groups[grepl("^(t_)|(inv)", colnames(dataFrame)) &! grepl(":", colnames(dataFrame))] <- "Translocations"
groups <- factor(groups)
names(groups) <- colnames(dataFrame)
table(groups)

#' Poor man's imputation
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
dataFrame <- as.data.frame(sapply(dataFrame, poorMansImpute))

#' Compute gene:gene interactions
#+ interactions, cache=TRUE, fig.width=6, fig.height=6
genomicData <- cbind(dataList$Genetics, dataList$Cytogenetics)
genomicData <- (sapply(unique(sub("(_ITD)|(_TKD)|(_other)|(_mono)|(_bi)|(_p172)|(_p140)","",colnames(genomicData))), function(x) rowSums(genomicData[grep(x, colnames(genomicData))])) > 0)+0
genomicData <- genomicData[,colSums(genomicData, na.rm=TRUE)>=8]
genomicGroups <- factor(grepl("^[a-z]", colnames(genomicData)) + grepl("t_*[0-9M]", colnames(genomicData) ), labels=c("Genetics","CNA","BT"))
dim(genomicData)
logPInt <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(genomicData[,i], genomicData[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
odds <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(table(genomicData[,i], genomicData[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
pairs <- sapply(1:ncol(genomicData), function(i) colMeans(genomicData * genomicData[,i], na.rm=TRUE))
diag(logPInt) <- 0
diag(odds) <- 1
colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(genomicData)
odds[odds<1e-3] = 1e-4
odds[odds>1e3] = 1e4
odds[10^-abs(logPInt) > 0.05] = 1
logOdds=log10(odds)
diag(logPInt) <- NA

par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)
ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
d <- dist(t(genomicData[,ix]) + 10*as.numeric(genomicGroups))
h = hclust(d, method="average")
o = order(genomicGroups,-colSums(genomicData, na.rm=TRUE))#order(cmdscale(d, k=1))#h$order #c(h$order,(length(h$order) +1):ncol(interactions))
M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
M[lower.tri(M)] <- cut(logOdds[o,o][lower.tri(M)], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)
M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(genomicData), breaks=c(-1,0,5,10,20,50,100,200,600))) + 9 
image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(9,"BrBG"), c("white",brewer.pal(7,"Reds"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
l <- colnames(logPInt)[o]
mtext(side=1, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
mtext(side=2, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
abline(h=0:ncol(logPInt)+.5, col="white", lwd=.5)
abline(v=0:ncol(logPInt)+.5, col="white", lwd=.5)
P <- 10^-abs(logPInt[o,o])
P[upper.tri(P)] <- NA
w = arrayInd(which(p.adjust(P, method="BH") < .1), rep(nrow(logPInt),2))
points(w, pch=".", col="black")
w = arrayInd(which(p.adjust(P) < .05), rep(nrow(logPInt),2))
points(w, pch="*", col="black")
image(y = 1:9 +18, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=c("white",brewer.pal(7,"Reds")), add=TRUE)
axis(side = 4, at = seq(1,7) + 19, cex.axis=.66, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
mtext(side=4, at=28, "Frequency", las=2, line=-1,cex=.66)
image(y = 1:8 +5, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
axis(side = 4, at = seq(1,7) + 5.5, cex.axis=.66, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
mtext(side=4, at=14, "Odds ratio", las=2, line=-1,cex=.66)
mtext(side=4, at=4, "Significance", las=2, line=-1,cex=.66)
points(x=rep(ncol(logPInt),2)+2.5, y=1:2, pch=c("*","."))
image(x=rep(ncol(logPInt),2)+c(2,3), y=(2:3) +0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=.66, line=0.2)

#' Hotspots
#+ interactionsHotspot, cache=TRUE, fig.width=2.5, fig.height=6
par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)
h <- "(_ITD)|(_TKD)|(_other)|(_mono)|(_bi)|(_p172)|(_p140)"
k <- grep(h, colnames(dataList$Genetics), value=TRUE)
hotspots <- dataList$Genetics[k]
for(g in c("TET2","WT1")){
	t <- table(data.frame(mutationData$SAMPLE_NAME, mutationData$CONSEQUENCE=="misssense")[mutationData$Result=="ONCOGENIC" & mutationData$GENE==g,] )[,] > 0
	colnames(t) <- paste(g,c("_other","_misssense"), sep="")
	hotspots <- cbind(hotspots, t[,2:1]+0)
}
t <- table(data.frame(mutationData$SAMPLE_NAME, DNMT3A_hotspot=grepl("p.R882", mutationData$AA_CHANGE))[mutationData$Result=="ONCOGENIC" & mutationData$GENE=="DNMT3A",] )[,] > 0
colnames(t) <- c("DNMT3A_other","DNMT3A_p882")
hotspots <- cbind(hotspots, t[,2:1]+0)
t <- table(data.frame(mutationData$SAMPLE_NAME, grepl("p.Q61", mutationData$AA_CHANGE))[mutationData$Result=="ONCOGENIC" & mutationData$GENE=="NRAS",] )[,] > 0
colnames(t) <- c("NRAS_p61","NRAS_p12/13")
hotspots <- cbind(hotspots, t[,2:1]+0)
hotspots <- hotspots[order(colSums(genomicData, na.rm=TRUE)[sub("_.+","",colnames(hotspots))], decreasing=TRUE)]
godspots <- unique(sub("_.+","",colnames(hotspots)))

potspots <- sapply(godspots, function(i) {
								g <- grep(i, colnames(hotspots), value=TRUE)
								l <- length(g)
								h <- factor(rowSums(hotspots[g] * rep(2^(1:l)/2, each=nrow(hotspots))), levels=0:2^l)
								sapply(1:ncol(genomicData), function(j) {
											t <- table(h, genomicData[,j])
											f<- try(fisher.test(t[-1,]) , silent=TRUE) 
											if(class(f)=="try-error") 0 
											else f$p.val})
							})
oddspots <- sapply(hotspots, function(x) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(table(x, genomicData[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else if (f$p.val > 0.05) 1 else f$estimate} ))
rownames(oddspots) <- rownames(potspots) <- colnames(genomicData)
oddspots[oddspots<1e-3] = 1e-4
oddspots[oddspots>1e3] = 1e4
logOdds=log10(oddspots)
c <- cut(logOdds[o,], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)
M <-  matrix( as.numeric(c), ncol=ncol(oddspots), nrow=nrow(oddspots))
colnames(M) <- colnames(oddspots)
rownames(M) <- rownames(oddspots)[o]
for(s in godspots){
	M[s,grep(s, colnames(M))] <- NA
	potspots[s,s] <- NA
}
image(x=1:ncol(M), y=1:nrow(M), t(M), col=brewer.pal(9,"BrBG"), breaks=0:9, xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(M)+3), ylim=c(0, nrow(M)+3))
l <- rownames(oddspots)[o]
mtext(side=1, at=1:ncol(M), colnames(oddspots), cex=.66, font=3)
mtext(side=2, at=1:nrow(M), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
abline(h=0:ncol(M)+.5, col="white", lwd=.5)
abline(v=0:ncol(M)+.5, col="white", lwd=.5)
x0 <- cumsum(c(0,table(sub("_.+","",colnames(hotspots)))[godspots]))
x00 <-  x0[-length(x0)] + .5  +  table(sub("_.+","",colnames(hotspots)))[godspots]/2
abline(v=x0+.5, col="white", lwd=2)
i <- 1
for(q in list(p.adjust(potspots[o,], method="BH") < .1,p.adjust(potspots[o,]) < .05)){
	w = arrayInd(which(q), rep(nrow(potspots),2))
	x <- x00[w[,2]]
	rect(x0[w[,2]]+.5, w[,1]-.5,x0[w[,2]+1]+.5, w[,1]+.5, lwd=c(1,2)[i])
	points(x,w[,1], pch=c(".","*")[i], col="black")
	i <- 1+i
}
image(y = 1:8 +5, x=rep(ncol(M),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
axis(side = 4, at = seq(1,7) + 5.5, cex.axis=.66, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
mtext(side=4, at=14, "Odds ratio", las=2, line=-1,cex=.66)
mtext(side=4, at=4, "Significance", las=2, line=-1,cex=.66)
points(x=rep(ncol(M),2)+2.5, y=1:2, pch=c("*","."))
image(x=rep(ncol(M),2)+c(2,3), y=(2:3) +0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
rect(ncol(M)+2,.5,ncol(M)+3,1.5, lwd=2)
rect(ncol(M)+2,1.5,ncol(M)+3,2.5, lwd=1)
mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=.66, line=0.2)



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
t[is.na(t) | !clinicalData$TPL_Phase %in% "CR1" | !clinicalData$TPL_type %in% c("ALLO","FREMD") ] <- Inf ## Only allografts in CR1
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
mainGroups <- grep("[A-Z].+[A-Z]",levels(groups), invert=TRUE, value=TRUE)
mainGroups
mainIdx <- groups %in% mainGroups
efsIdx <- !grepl("TPL_os", colnames(dataFrame))
whichRFXEfs <- which((colSums(dataFrame)>=8 | mainIdx) & efsIdx) # ie, > 0.5%
mainIdxEfs <- mainIdx & efsIdx
osIdx <- !grepl("TPL", colnames(dataFrame)) ## Exclude TPL from OS analyses..
whichRFXOs <- which((colSums(dataFrame)>=8 | mainIdx) & osIdx) # ie, > 0.5%
mainIdxOs <- mainIdx & osIdx
osTDIdx <- !grepl("TPL_efs", colnames(dataFrame))
whichRFXOsTD <- which((colSums(dataFrame)>=8 | mainIdx) & osTDIdx) # ie, > 0.5%
mainIdxOsTD <- mainIdx & osTDIdx

#' ### 0. Proportional hazards test
#+ fig.width=14, fig.height=5
coxModel <- coxph(os ~ gender + AOD + Study, data=clinicalData) # Most basic
phTest <- cox.zph(coxModel)
phTest
par(mfrow=c(1,4), bty="n")
for(i in 1:4) plot(phTest[i])


#' #### Basic KM plots
#' TP53 and TPL
#+ kmTP53vTPL, fig.width=3, fig.height=2.5
f <- formula("osTD ~ TP53 + TPL_os")
s <- survfit(f, data=dataFrameOsTD)
c <- coxph(f, data=dataFrameOsTD)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)


#' FLT ITD v TKD
#+ kmITDvTKD, fig.width=3, fig.height=2.5
f <- formula("os ~ FLT3_ITD + FLT3_TKD")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)


#' CEBPA mono v bi
#+ kmMonoVBi, fig.width=3, fig.height=2.5
f <- formula("os ~ CEBPA_mono + CEBPA_bi")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)

#' IDH2 p140 v p172
#+ km140v172, fig.width=3, fig.height=2.5
f <- formula("os ~ IDH2_p140 + IDH2_p172")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)

#' NPM1, FLT3_ITD and DNMT3A
coxph(os ~ NPM1*FLT3_ITD*DNMT3A, data=dataFrame)
s <- survfit(os ~ NPM1+FLT3_ITD+DNMT3A, data=dataFrame)
plot(s, col=set1)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)


#' #### NONC and complex karyotyope
table(complex=clinicalData$complex, `#Recurr.Cyto`=rowSums(dataList$Cytogenetics))
summary(coxph(os ~ ., data=dataFrame[mainIdxOs][colSums(dataFrame[mainIdxOs]) > 10]))

NONC <- rowSums(cbind(dataList$Cytogenetics[names(dataList$Cytogenetics)!="complex"], dataList$Genetics), na.rm=TRUE)
summary(coxph(os ~ . + NONC, data=dataFrame[mainIdxOs][colSums(dataFrame[mainIdxOs]) > 10]))

#+ NONC, fig.width=2.5, fig.height=2.5
r=colorRampPalette(set1[c(3,2,4,1,5)])(9)
s <- survfit(os ~ pmin(NONC,8))
plot(s, col=r, xlab="Days", ylab="Survival", mark=NA)
legend('topright', bty='n', col=r,legend= 0:8, lty=1, title="# onc. mut.")

plot(0:8,summary(s)$table[,"median"], xlab="Number of oncogenic mutations", ylab="Median survival time", pch=19, col=r)
segments(0:8, summary(s)$table[,"0.95LCL"],0:8, sapply(summary(s)$table[,"0.95UCL"], function(x) ifelse(is.na(x), max(os[,1]), x)), col=r)

#' ### 1. Random effects: static model
#+ coxRFXFitEfs, warning=FALSE, cache=TRUE
## Fit Cox model
coxRFXFitEfs <- CoxRFX(dataFrame[,whichRFXEfs], efs, groups=groups[whichRFXEfs])
#+ coxRFXFitOs, warning=FALSE, cache=TRUE
coxRFXFitOs <- CoxRFX(dataFrame[,whichRFXOs], os, groups=groups[whichRFXOs])
#+ coxRFXFitOsMain, warning=FALSE, cache=TRUE
coxRFXFitOsMain <- CoxRFX(dataFrame[,mainIdxOs], os, groups=groups[mainIdxOs])


#' Coefficients
par(mar=c(5,7,1,1))
colGroups <- c(brewer.pal(12, "Paired")[c(6,3,4,5,12,9,1,2,7)],"#999999", brewer.pal(12, "Paired")[c(10,8)])
names(colGroups) <- levels(groups)[order(toupper(levels(groups)))]
#o <- order(coxRFXFitEfs$mu)
#boxplot(coef(coxRFXFitEfs) ~ factor(coxRFXFitEfs$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
plot(coxRFXFitEfs, col=colGroups, order=order(coxRFXFitEfs$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

#' ### 2. Random effects: Time-dependent model (TPL)
#' Fit TD CoxRFX model
#+ coxRFXFitEfsTD, cache=TRUE
coxRFXFitEfsTD <- CoxRFX(dataFrameEfsTD[,whichRFXEfs], efsTD, groups[whichRFXEfs])
#+ coxRFXFitOsTD, cache=TRUE
coxRFXFitOsTD <- CoxRFX(dataFrameOsTD[,whichRFXOsTD], osTD, groups[whichRFXOsTD])


par(mar=c(5,7,1,1))
#o <- order(coxRFXFitEfsTD$mu)
#boxplot(coef(coxRFXFitEfsTD)/log(2) ~ factor(coxRFXFitEfsTD$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
plot(coxRFXFitEfsTD, col=colGroups, order=order(coxRFXFitEfsTD$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

plot(coef(coxRFXFitEfs), coef(coxRFXFitEfsTD), col=colGroups[groups[whichRFXEfs]]) # Note the sign change for TPL..
abline(0,1)

par(mar=c(5,7,1,1))
#o <- order(coxRFXFitOsTD$mu)
#boxplot(coef(coxRFXFitOsTD)/log(2) ~ factor(coxRFXFitOsTD$groups, levels=levels(groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
plot(coxRFXFitOs, col=colGroups, order=order(coxRFXFitOs$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

##' #### Coefficient estimates (MAP)
##+ fig.width=14, fig.height=5
#for(model in c("coxRFXFitEfsTD","coxRFXFitOs")){
#	m <- get(model)
#	for(g in levels(groups)){
#		par(mar=c(7,4,2,0))
#		x <- sort(m$coefficients[m$groups==g])
#		v <- diag(m$var)[m$groups==g][ order(m$coefficients[m$groups==g])]
#		plot(x, las=2, pch=16, xaxt="n", main=paste(model,g), cex = .5+pmin(3,-log10(pchisq((x-m$mu[g])^2/v, 1, lower.tail=FALSE))), xlab="", ylab="coefficient")
#		segments(1:length(x),x - 2*sqrt(v) ,1:length(x), x+2*sqrt(v), col=colGroups[g])
#		axis(side=1, at = 1:length(x), labels = names(x), las=2, cex.axis=.6)
#		abline(v=1:length(x), lty=3, col="grey")
#		abline(h=m$mu[g])
#		abline(h = m$mu[g] + c(-1,1)*sqrt(m$sigma2[g]), lty=2)
#	}
#}

#' #### OS v EFS
#+ EFSvsOScoef, fig.width=4, fig.height=4
efs2Os <- match(names(coef(coxRFXFitOs)), sub("efs","os",names(coef(coxRFXFitEfsTD))))
plot(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs), pch=16, xlab="Coefficient EFS", ylab="Coefficient OS", col=colGroups[coxRFXFitEfsTD$groups])
segments((coef(coxRFXFitEfsTD) - sqrt(diag(coxRFXFitEfsTD$var)))[efs2Os], coef(coxRFXFitOs),(coef(coxRFXFitEfsTD)+ sqrt(diag(coxRFXFitEfsTD$var)))[efs2Os], coef(coxRFXFitOs),  col=paste(colGroups,"44",sep="")[coxRFXFitEfsTD$groups])
segments(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs) - sqrt(diag(coxRFXFitOs$var)),coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs)+sqrt(diag(coxRFXFitOs$var)),  col=paste(colGroups,"44",sep="")[coxRFXFitEfsTD$groups])
text(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs), ifelse(sqrt(coef(coxRFXFitEfsTD)[efs2Os]^2 + coef(coxRFXFitOs)^2) > 0.5, names(coef(coxRFXFitEfsTD))[efs2Os], ""), pos=3)
abline(0,1)
abline(h=0, lty=3)
abline(v=0, lty=3)
points(coxRFXFitEfsTD$mu, coxRFXFitOs$mu, bg=colGroups, pch=21, cex=2)
title(main=paste("rho =", round(cor(coef(coxRFXFitEfsTD)[efs2Os], coef(coxRFXFitOs), use="c"), 2)))

#+ EFSvsOSsigma
plot(coxRFXFitEfsTD$sigma2, coxRFXFitOs$sigma2, col=colGroups, pch=19, xlab="Variance EFS", ylab="Variance OS")
text(coxRFXFitEfsTD$sigma2, coxRFXFitOs$sigma2, pos=3, levels(coxRFXFitEfsTD$groups))
abline(0,1)


#' #### Partial risk contributions
simpleGroupsEfs <- mergeLevels(coxRFXFitEfsTD$groups, mergeList=list(Genetics=c("Genetics","GeneGene"), Cytogenetics=c("Cytogenetics","CytoCyto","GeneCyto"), Treatment=c("Treatment","GeneTreat","CytoTreat")))
partRiskEfsTDsimple <- PartialRisk(coxRFXFitEfsTD, groups = simpleGroupsEfs)
partRiskEfsTD <- PartialRisk(coxRFXFitEfsTD)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponentsEfs <- diag(cov(partRiskEfsTDsimple, use="complete"))
varianceComponentsEfs
PlotVarianceComponents(coxRFXFitEfs, col=colGroups)
title("Risk contributions EFS")

PlotVarianceComponents(coxRFXFitEfsTD, col=colGroups)
title("Risk contributions EFS (time-dep TPL)")


partRiskVar <- PartialRiskVar(coxRFXFitEfsTD,  groups = simpleGroupsEfs)
x=c(varianceComponentsEfs, Error=mean(rowSums(partRiskVar)))
pie(x, col=c(colGroups[names(varianceComponentsEfs)], "grey"), labels = paste(names(x), round(x/sum(x),2)))

simpleGroupsOs <- mergeLevels(coxRFXFitOs$groups, mergeList=list(Genetics=c("Genetics","GeneGene"), Cytogenetics=c("Cytogenetics","CytoCyto","GeneCyto"), Treatment=c("Treatment","GeneTreat","CytoTreat")))
partRiskOsMain <- PartialRisk(coxRFXFitOsMain)
partRiskOs <- PartialRisk(coxRFXFitOs)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponentsOs <- diag(cov(partRiskOsMain, use="complete"))
varianceComponentsOs
PlotVarianceComponents(coxRFXFitOs, col=colGroups)
title("Risk contributions OS")

PlotVarianceComponents(coxRFXFitOsMain, col=colGroups)
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
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(coxRFXFitOs$coef[groups[whichRFXOs] == g]), each=length(x))  )), col=paste(colGroups[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(sqrt(coxRFXFitOs$sigma2[g]) * c(-1,1), each=length(x))  + coxRFXFitOs$mu[g])), col=paste(colGroups[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( coxRFXFitOs$mu[g] )), col=colGroups[i], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) ), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}

#' Partial risk
#+ partialRiskKM, fig.width=2, fig.height=2
i <- 1
for(g in colnames(partRiskOsMain)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	m <- mean(rowSums(partRiskOsMain[,colnames(partRiskOsMain)!=g]))
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(partRiskOsMain[,g]), each=length(x)) +m)), col=paste(colGroups[g],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskOsMain[,g], c(0.25,0.75)), each=length(x))+m)), col=paste(colGroups[g],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskOsMain[,g], c(0.05,0.95)), each=length(x))+m)), col=paste(colGroups[g],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( median(partRiskOsMain[,g])+m)), col=colGroups[g], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) *exp(+m)), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}
par(mfrow=c(1,1))

#' #### Confidence intervals
#' Estimate confidence intervals by parametric bootstrap. Note that the usual sample with replacement 
#' yields inconsistencies for the interaction terms due to the overdispersed correlations.
#+ parBoot, cache=TRUE, eval=TRUE
set.seed(42)
risk <- as.matrix(dataFrame[mainIdxOs]) %*% coxRFXFitOsMain$coefficients
risk <- risk - mean(risk)
parBoot <- mclapply(1:100, function(i) {
			s <- SimSurvNonp(risk, os)
			c <- try(CoxRFX(dataFrame[mainIdxOs], s, groups=groups[mainIdxOs], sigma0=0.1, nu=0))
			if(class(c)=="try-error")
				return(s)
			c$X <- NULL # set X to zero to save mem
			return(c)
		}, mc.cores=10)

#' Plots of mean, sigma and df
#+ parBootPlots, fig.width=3, fig.height=2,  eval=TRUE
boxplot(t(sapply(parBoot, `[[`, "sigma2")), border=colGroups[names(parBoot[[1]]$sigma2)], lty=1, pch=16, staplewex=0, ylab="sigma2", las=2, log="y", ylim=c(1e-3,1))
abline(h=0, lty=3)
points(coxRFXFitOsMain$sigma2,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "mu")), border=colGroups[names(parBoot[[1]]$mu)], lty=1, pch=16, staplewex=0, ylab="mu", las=2)
abline(h=0, lty=3)
points(coxRFXFitOsMain$mu,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "df")), border=colGroups[names(parBoot[[1]]$mu)], lty=1, pch=16, staplewex=0, ylab="df", las=2)
abline(h=0, lty=3)
points(coxRFXFitOsMain$df,  pch=19)

#' Coefficients
#+ parBootSignif, fig.width=2.5, fig.height=2.5
v <- apply(sapply(parBoot, `[[`, "coefficients"), 1, var, na.rm=TRUE)
w <- diag(coxRFXFitOsMain$var) ## H^{-1}
w2 <- diag(coxRFXFitOsMain$var2) ## H^{-1} I H^{-1}
c <- coef(coxRFXFitOsMain)
plot(c^2/v, c^2/w, log="xy", xlab="Chi2 (bootstrap)", ylab="Chi2 (analyt.)", cex=.66)
par(xpd=NA)
points(c^2/v, c^2/w2, pch=16, cex=.7)
arrows(c^2/v, c^2/w, c^2/v,c^2/w2, length=0.05)
abline(0,1)
abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))

#' Table with significance
#+ parBootTable, results='asis'
library(xtable)
print(xtable(data.frame(group = groups[mainIdxOs],coef=round(c,4), sd = round(sqrt(w2),4), boot=sig2star(pchisq(c^2/v,1, lower.tail=FALSE)), var2=sig2star(pchisq(c^2/w2,1, lower.tail=FALSE)),var=sig2star(pchisq(c^2/w,1, lower.tail=FALSE)))),  type="html")

#' Volcano plot
#+ volcano, fig.width=2.5, fig.height=2.5
plot(c, c^2/w2, log='', col=colGroups[as.character(coxRFXFitOsMain$groups)], pch=19, ylab="Z-score",xlab="log hazard")
abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))
p <- pchisq(c^2/w2,1,lower.tail=FALSE) ## pvalues coxRFX
w <- which(p.adjust(p,"BH") < 0.1)
par(xpd=NA)
text(c[w], (c^2/w2)[w], names(c[w]), pos=3)

#' P-values and random model
#+ pValues, fig.width=2.5, fig.height=2.5, cache=TRUE
set.seed(42)
X <- apply(coxRFXFitOsMain$X, 2,sample) ## random covariates
coxRFXFitOsRain <- CoxRFX(X, os, groups=coxRFXFitOsMain$groups,  sigma0=0.1, nu=1) ## model
w2 <- diag(coxRFXFitOsRain$var2) 
c <- coef(coxRFXFitOsRain)
p2 <- pchisq(c^2/w2,1,lower.tail=FALSE)
plot(seq(0,1,l=length(p)+1)[-1],sort(p2), xlab="P-value (expected)", ylab="P-value (observed)", pch=16, col="grey")
abline(0,1)
points(seq(0,1,l=length(p)+1)[-1],sort(p), xlab="P-value (expected)", pch=16, col="black")
legend("topleft",bty="n", c("observed","randomised"), pch=16, col=c("black","grey"))



#+ parBootVarianceComp, fig.width=3, fig.height=2, cache=TRUE, eval=TRUE
v <- t(sapply(parBoot, function(x) {t <- try(VarianceComponents(x, newX=dataFrame[mainIdxOs])); if(class(t)=="try-error") rep(NA, nlevels(x$groups)+1) else t}))
boxplot(v, border=colGroups[colnames(v)], lty=1, pch=16, staplewex=0, ylab="variance comp.", las=2)
abline(h=0, lty=3)
points(VarianceComponents(coxRFXFitOsMain),  pch=19)

round(cov(v), 2)


rm(parBoot)

#' #### Large overview panel
#+ overviewRFX, fig.width=7, fig.height=5
for(model in c("coxRFXFitOs", "coxRFXFitOsMain")){
	layout(matrix(1:6, nrow=2), width=c(6,6,2),height=c(3,12))
	par(bty="n", mar=c(0,4,4,2), mgp=c(2,.5,0), tcl=-.25)
	v <- VarianceComponents(get(model), type="rowSums")
	w <- TRUE#get(model)TD$groups %in% mainGroups
	f <- factor(as.character(get(model)$groups[w]), levels=levels(get(model)$groups[w])[order(abs(v), decreasing=TRUE)])
	o <- order(f,get(model)$coefficients[w])
	c <- coef(get(model))[w][o]
	plot(c, type='p',  ylab="log hazard/variable", xaxt='n', col=colGroups[as.character(f)[o]], xaxs="i", lwd=2, pch=NA)
	y <- seq(-2,2,l=100)
	par(xpd=FALSE)
	abline(h=seq(-.5,.5,.25), col="grey", lty=c(1,3))
	colRamp <- sapply(colGroups[levels(f)], function(x) c(colorRampPalette(c("white",x))(11)[-1]))
#for(l in levels(f)){
#	x <- dnorm(y,get(model)$mu[l], sqrt(get(model)$sigma2[l]))
#	x <- x/max(x)*15
#	#polygon(x + (cumsum(table(f)) - table(f))[l], y, col=colRamp[1,l], lty=0)
#	#lines(x + (cumsum(table(f)) - table(f))[l], y, col=colRamp[1,l], lwd=1)
#	#segments((cumsum(table(f)) - table(f))[l],get(model)$mu[l],(cumsum(table(f)) - table(f))[l]+max(x),get(model)$mu[l], col=colRamp[8,l])
#	rect((cumsum(table(f)) - table(f))[l], get(model)$mu[l] - 2*sqrt(get(model)$sigma2[l]) , cumsum(table(f))[l] ,get(model)$mu[l]+ 2*sqrt(get(model)$sigma2[l]), col=colRamp[1,l], lty=0)
#	rect((cumsum(table(f)) - table(f))[l], get(model)$mu[l] - 1*sqrt(get(model)$sigma2[l]) , cumsum(table(f))[l] ,get(model)$mu[l]+ 1*sqrt(get(model)$sigma2[l]), col=colRamp[4,l], lty=0)	
#	segments((cumsum(table(f)) - table(f))[l], get(model)$mu[l]  , cumsum(table(f))[l] ,get(model)$mu[l], col=colRamp[8,l])	
#}
	par(xpd=NA)
	segments(seq_along(c), get(model)$mu[as.character(f)[o]], seq_along(c), c, col=colGroups[as.character(f)[o]])
#points(c, col=col1[as.character(f)[o]], pch=16, cex=.66)
#mtext(side=3, levels(f), at= table(f)/2+ cumsum(c(0,table(f)[-nlevels(f)])), cex=.66)
	rotatedLabel(table(f)/2+ cumsum(c(0,table(f)[-nlevels(f)])), y = rep(par("usr")[4]*.8, nlevels(f)), pos=3, labels=levels(f))
	par(bty="L", mar=c(4,4,1,2))
	X <- get(model)$X[,w][,o]
	p <- order(X %*% c)
	X <- X[p,]
	x <- t(X)/apply(X,2,max)
#h <- hclust(dist(t(x)))
	image(x=1:nrow(x),y=1:ncol(x),x*.9 + as.numeric(f[o])-1 + 1e-5, useRaster=TRUE, 
			col=colRamp, 
			breaks=seq(0,nlevels(f), 0.1),  ylab="Patients",xlab="Variable", xlim=c(0,nrow(x)), ylim=c(0,ncol(x)))
	
	par(mar=c(3,.5,1,0.5), xpd=NA)
	PlotVarianceComponents(get(model), col = colGroups)
	mtext(side=3, at = 0, "Variance components",line=0, cex=.66)
	
	xo <- 3
	par(mar=c(4,0,1,0.5), xpd=NA, bty="n")
	plot(NA,NA,ylim=c(1,nrow(X)), xlim=c(-xo/2, (nlevels(f)-.5)*xo), xaxt="n", yaxt="n", yaxs="i", xlab="Partial log hazard", ylab="")
	i <- 0
	for(l in levels(f)){
		h <- X[,f[o]==l] %*% c[f[o]==l]
		h <- h - mean(h)
		#plot(h, seq_along(h), pch=NA, xlab="log hazard", yaxs="i", yaxt="n", ylab='', xlim=c(-2,2))
		par(xpd=FALSE)
		abline(v=-1:1 + i*xo, col="grey")
		par(xpd=NA)
		segments(h + i*xo, seq_along(h), + i*xo, seq_along(h), col=colGroups[l])
		#mtext(side=3, l, cex=.66, line=.5)
		axis(side=1, at=-1:1 + xo*i, labels=-1:1)
		x <- seq(-1.5,1.5,l=101)
		lines(x +i*xo, dnorm(x, 0, sd(h))*100 /  dnorm(0, 0, sd(h)) + length(h)*1.01, col=colGroups[l])
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
	x <- seq(from=-3,to=3, l=512)
	par(mar=c(1,.5,3,1), xpd=FALSE)
	plot(x, invHazardDist(-log(l[2]) /exp(x) ), col='black', lty=c(2,1,2)[2], xlab="", ylab="Time", type='l',ylim=c(0,5000), xlim=c(min(h), 4), xaxt="n")
	abline(v=seq(-2,3,1), col="grey")
#for(i in seq_along(l)[-1])
#	lines(x, invHazardDist(-log(l[i]) /exp(x) ), col='black', lty=c(2,1,2)[i])
	polygon(c(x, rev(x)), c(invHazardDist(-log(l[1]) /exp(x) ), rev(invHazardDist(-log(l[3]) /exp(x))) ), border=NA, col="#88888888")
#axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
#mtext(side=4, line=2.5, "Time")
	mtext(side=3, at = log(-log(l)/hazardDist(par("usr")[4])), text=paste(100*l, "%", sep=""), cex=.66)
	par(xpd=NA)
	mtext(side=3, "Survival",line=2, cex=.66)
	mtext(side=2, line=2.5, "Time",cex=.66)
	q <- quantile(os[,1], seq(0,1,0.1), na.rm=TRUE)
	image(x=c(3.5,4), y = q, matrix(1:10, nrow=1), col=brewer.pal(10,'RdBu'), add=TRUE)
	
	
	par(mar=c(4,.5,1,1), xpd=FALSE)
	plot(h, seq_along(h), pch=NA, xlab="Total log hazard", yaxs="i", yaxt="n", ylab='', xlim=c(min(h), 4), xaxt='n')
	axis(side=1, at=-2:3)
	abline(v=seq(-2,3,1), col="grey")
	segments(h, seq_along(h),0, seq_along(h))
#mtext(side=3, "Total log hazard", cex=.66, line=.5)
	par(xpd=NA)	
	x <- seq(-3,5,l=100)
	lines(x , dnorm(x, 0, sd(h))*100 /  dnorm(0, 0, sd(h)) + length(h)*1.01)
	
	c <- cut(os[p,1], q)
#c[os[p,2]==0] <- NA
	image(x=c(3.5,4), y = c(0,seq_along(c)), matrix(as.numeric(c), nrow=1), col=brewer.pal(10,'RdBu'), add=TRUE, useRaster=TRUE)
	points(x=rep(4.1, sum(os[p,2]==0, na.rm=TRUE)), y=which(os[p,2]==0), pch=".")
}

#' #### Stars
#+ stars, fig.width=12, fig.height=12
set.seed(42)
library(HilbertVis)
nStars <- 32
for(l in list(c("efsTD","partRiskEfsTD"),c("os","partRiskOs"),c("efsTD","partRiskEfsTDsimple"),c("os","partRiskOsMain"))){
	t <- get(l[1])
	p <- get(l[2])
	p <- p[,colnames(p)!="Nuisance"]
	locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
	s <- sample(nrow(p),nStars^2) #1:(nStars^2)
	h <- hclust(dist(p[s,]))
	x <- p - rep(colMeans(p), each=nrow(p))
	x <- x/(2*sd(x)) + 1
	if(l[1]=="efsTD")
		c <- cut((t[s,2]-t[s,1])[h$order], quantile(t[,2]-t[,1], seq(0,1,0.1), na.rm=TRUE))
	else
		c <- cut(t[s,1][h$order], quantile(t[,1], seq(0,1,0.1), na.rm=TRUE))
	if(l[2]=="partRiskOsMain")
		x <- x[,c("Clinical","Demographics","Genetics","Cytogenetics","Translocations")]
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
points(t(x) * c - t(j) *s, t(x) * s + t(j) * c, pch=16, cex=.25, col=colGroups[colnames(x)])
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
totalRiskOs <- rowSums(partRiskOsMain)
survConcordance( os~totalRiskOs)


#' Partial values of Harrel's C
c <- PartialC(coxRFXFitEfsTD, groups=simpleGroupsEfs)
b <- barplot(c[1,], col=colGroups[colnames(c)], las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)


#' 4. Test/train splits
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
b <- barplot(c[1,], col=colGroups[colnames(c)], las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' Overall
totalRiskOsTest <- rowSums(partialRiskOsTest)
survConcordance(os[!trainIdx]~totalRiskOsTest )

#' Compared to molecular risk
predictiveRiskTest <- rowSums(partialRiskOsTest[,-which(colnames(partRiskOsMain) %in% c("Treatment","GeneTreat","CytoTreat","Trial"))])
survConcordance(os[!trainIdx] ~ predictiveRiskTest)

#barplot(c(CGP=survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)$concordance, MolecularRisk = survConcordance(survivalTD[!testIxTD] ~ c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[tplSplit][!testIxTD]])[[1]]))

s <- efs[!trainIdx]
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], predictiveRiskTest[!is.na(s)], seq(0,5000,100)))
plot(AUC.uno(s[!is.na(s)],s[!is.na(s)], totalRiskOsTest[!is.na(s)], seq(0,5000,100)), add=TRUE, col="black")
legend("bottomright",c("w/o treatment","w treatment"), col=2:1, lty=1)


#' ### 3. Recursive partitioning & randomForestSRC on OS
#+ tree, fig.width=3, fig.height=3
library(rpart)
library(randomForestSRC)
tree <- rpart(os ~ ., data=dataFrame[mainIdxOs & osIdx])
plot(tree)
text(tree)
survConcordance(na.omit(os)~predict(tree))

#+ treeCV, fig.width=3, fig.height=3
treeCV <- rpart(os[trainIdx] ~ ., data=dataFrame[trainIdx,mainIdxOs & osIdx])
plot(treeCV)
text(treeCV)
survConcordance(os[!trainIdx]~predict(treeCV, newdata=dataFrame[!trainIdx, mainIdxOs & osIdx]))

#+ rForest, cache=TRUE
rForest <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs & osIdx]), ntree=100)
boxplot(rForest$importance ~ factor(as.character(groups[mainIdxOs & osIdx])), border= colGroups, staplewex=0, pch=16, cex=0.75, ylab="RSF importance", lty=1, xaxt="n")
rForestVimp <- sapply(mainGroups, function(g) vimp(rForest, colnames(dataFrame)[which(groups==g)]))

survConcordance(na.omit(os)~predict(rForest, importance="none")$predicted)

#+ rForestCV, cache=TRUE
rForestCV <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs & osIdx])[trainIdx,], ntree=100, importance="none")
p <- predict(rForestCV, newdata = dataFrame[!trainIdx,mainIdxOs & osIdx], importance="none")
survConcordance(os[!trainIdx]~ p$predicted)


#' ### 4. Stability selection: All terms

#' Fit model
#+ coxCPSS, cache=TRUE, fig.width=3, fig.height=2.5

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
#+ CoxCPSSIntPlot, fig.width=3, fig.height=2.5
for(i in grep(":", selectedIntEfs, value=TRUE)){
	j <- strsplit(i, ":")[[1]]
	plot(survfit(as.formula(paste("efsTD ~ ", paste( j, collapse="+"))), data=dataFrameEfsTD), col=set1)
	legend("topright", lty=1, col=set1, c("none", rev(j), "both"), bty="n")
}

#' Plots of the stability selection: No large difference between when interaction terms are enforce to only occur after main effects are included
plot(coxCPSSIntEfs$Pi, coxCPSSEfs$Pi[match(names(coxCPSSIntEfs$Pi), names(coxCPSSEfs$Pi))])


			


#' #### Fit model to the training subset only
#+ CoxCPSSIntEfsCV, cache=TRUE
scope <- c("Genetics","Cytogenetics","Treatment","Translocations")
coxCPSSIntEfsCV <- CoxCPSSInteractions(dataFrame[!is.na(efs) & trainIdx,groups %in% mainGroups & efsIdx], na.omit(efs[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
selectedIntEfsCV <- names(which(coxCPSSIntEfsCV$Pi > 0.8))
survConcordance(efs[!trainIdx] ~ predict(coxph(as.formula(paste("efs[trainIdx] ~", paste(selectedIntEfsCV, collapse="+"))), data=dataFrame[trainIdx, ]), newdata = dataFrame[!trainIdx, ]))

#tmp <- CoxRFX(dataFrame[trainIdx, selectedIntEfs], efs[trainIdx], nu=0)
#survConcordance(efs[!trainIdx] ~ as.matrix(dataFrame[!trainIdx, selectedIntEfs]) %*% coef(tmp))

#' #### CPSS on OS
#+ CoxCPSSIntOs, cache=TRUE
set.seed(42)
coxCPSSIntOs <- CoxCPSSInteractions(dataFrame[!is.na(os),groups %in% mainGroups & osIdx], na.omit(os), bootstrap.samples=50, scope = which(groups %in% scope))
selectedIntOs <- names(which(coxCPSSIntOs$Pi > 0.8))
coxCPSSIntOs

#' The corresponding coxph model (CV)
coxFitCPSSIntOsCV <- coxph(as.formula(paste("os[trainIdx] ~", paste(selectedIntOs, collapse="+"))), data=dataFrame[trainIdx, ])
survConcordance(os[!trainIdx] ~ predict(coxFitCPSSIntOsCV, newdata = dataFrame[!trainIdx, ]))


#' With cross-val
#+ CoxCPSSIntOsCV, cache=TRUE
coxCPSSIntOsCV <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx,groups %in% mainGroups & osIdx], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
selectedIntOsCV <- names(which(coxCPSSIntOsCV$Pi > 0.8))
survConcordance(os[!trainIdx] ~ predict(coxCPSSIntOsCV, newdata = dataFrame[!trainIdx, ]))


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
		MRisk = c(4,1,3,2)[clinicalData$M_Risk[!trainIdx]],
		tree = predict(treeCV, newdata=dataFrame[!trainIdx, mainIdxOs]),
		rForest = predict(rForestCV, newdata = dataFrame[!trainIdx,mainIdxOs], importance="none")$predicted,
		coxRFX = as.matrix(dataFrame[!trainIdx,whichRFXOs]) %*% coef(coxRFXFitOsTrain),
		coxBIC = predict(coxBICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]),
		coxAIC = predict(coxAICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]),
		coxCPSS = predict(coxFitCPSSIntOsCV, newdata = dataFrame[!trainIdx, ])
)

#+ concordanceCV,fig.width=3, fig.height=2.5
concordanceCV <- sapply(predictedRiskCV, function(x) {c <- survConcordance(os[!trainIdx] ~ x); c(c$concordance, c$std.err)})
concordanceCV
o <- order(concordanceCV[1,])
barplot(concordanceCV[1,o], border=NA, col= set1[-6], las=2, xaxt="n", ylab="Concordance", ylim=c(0.5,0.75), xpd=FALSE) -> b
segments(b,concordanceCV[1,o]-concordanceCV[2,o],b,concordanceCV[1,o]+concordanceCV[2,o])
rotatedLabel(b, rep(0.49,length(b)), colnames(concordanceCV)[o], srt=45)

#+ aucCV, fig.width=3, fig.height=2.5
library(survivalROC)
library(survAUC)
#aucCV <- sapply(predictedRiskCV, function(x) survivalROC(Stime=os[!is.na(os) & !trainIdx,1], status=os[!is.na(os) &  !trainIdx,2], marker = scale(x[!is.na(os[!trainIdx])]), predict.time = 278, method="KM", cut.values=seq(-5,5,0.1))$AUC)
aucCV <- sapply(predictedRiskCV, function(x) AUC.uno(os[!is.na(os) & trainIdx], na.omit(os[!trainIdx][!is.na(x)]), scale(na.omit(x[!is.na(os[!trainIdx])])), c(90,365,1000))$auc)
o <- order(colMeans(aucCV))
barplot(aucCV[,o], border=1, col= rep(set1[-6],each=3), las=2, xaxt="n", ylab="AUC", beside=TRUE, density=c(NA, 48,24), ylim=c(0.5,0.75), xpd=FALSE) -> b
legend("topleft", bty="n", c("3mo","1yr","3yr"), fill='black', density=c(NA, 48,24))
rotatedLabel(b[seq(3, length(b), 3)], rep(0.49,length(predictedRiskCV)), names(predictedRiskCV)[o], srt=45)

#+ vennDiagram, fig.width=2.5, fig.height=2.5
library(VennDiagram)
grid.newpage()
pushViewport(viewport(w = .9, h = .9))
grid.draw(venn.diagram(list(CPSS=selectedIntOs, BIC=names(coef(coxBICOs)), AIC=names(coef(coxAICOs))), filename=NULL, lty=1, 
				col=set1[1:3], fill=set1[1:3], alpha=0.3, euler.d=TRUE, fontfamily="Helvetica", cat.fontfamily="Helvetica", euler.diagram=TRUE))

#' ### 8. Systematic cross-validation
#+ systematicCV, cache=TRUE
replicates <- 100 ## number of replicates
allConcordanceCV <- mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			c <- coxph(os[trainIdx] ~ 1, data=dataFrame[trainIdx,mainIdxOs])
			scopeStep <- as.formula(paste("os[trainIdx] ~", paste(colnames(dataFrame)[mainIdxOs], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxCPSSOsTrain <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx, mainIdxOs], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,mainIdxOs], os[trainIdx], groups=groups[mainIdxOs])
			rForestOsTrain <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs])[trainIdx,], ntree=100, importance="none")
			return(c(
					MRisk = survConcordance(os[!trainIdx]~c(4,1,3,2)[clinicalData$M_Risk[!trainIdx]])$concordance,
					BIC=survConcordance(os[!trainIdx]~predict(coxBICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]))$concordance,
					AIC=survConcordance(os[!trainIdx]~predict(coxAICOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]))$concordance,
					CPSS=survConcordance(os[!trainIdx]~predict(coxCPSSOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]))$concordance,
					#RFX=survConcordance(os[!tainIdx]~predict(coxRFXOsTrain, newdata=dataFrame[!trainIdx,mainIdxOs]))$concordance,
					RFX=survConcordance(os[!trainIdx]~as.matrix(dataFrame[!trainIdx,mainIdxOs]) %*% coef(coxRFXOsTrain))$concordance,
					rForest=survConcordance(os[!trainIdx]~predict(rForestOsTrain, newdata = dataFrame[!trainIdx,mainIdxOs & osIdx], importance="none")$predicted)$concordance
			))
		}, mc.cores=20)

allConcordanceCV <- do.call("rbind",allConcordanceCV)
colnames(allConcordanceCV) <- sub(".concordant","",colnames(allConcordanceCV))

#+ systematicCVboxplot, fig.width=2, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(t(apply(-allConcordanceCV,1,rank))),factor, levels=1:6)),table)
o <- order(colMeans(allConcordanceCV))
boxplot(allConcordanceCV[,o], notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:6, rep(par("usr")[3],6), colnames(allConcordanceCV)[o])

#+ systematicCVrank, fig.width=2, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
barplot(r[,o]/replicates, col=set1[c(3,2,4,1,5,7)], ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],6), colnames(allConcordanceCV)[o])
legend(par("usr")[1],1.5, fill=set1[c(3,2,4,1,5,7)], legend=1:6, bty="n", border=NA, horiz=TRUE, title="Rank")

#' Time-dependent models only
#+ systematicCVTD, cache=TRUE
replicates <- 100 ## number of replicates
allConcordanceCVTD <- do.call("rbind",mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrameOsTD)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			c <- coxph(osTD[trainIdx] ~ 1, data=dataFrameOsTD[trainIdx,mainIdxOsTD])
			scopeStep <- as.formula(paste("osTD[trainIdx] ~", paste(colnames(dataFrameOsTD)[mainIdxOsTD], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxRFXOsTrain <- CoxRFX(dataFrameOsTD[trainIdx,mainIdxOsTD], osTD[trainIdx], groups=groups[mainIdxOsTD])
			return(c(
							BIC=survConcordance(osTD[!trainIdx]~predict(coxBICOsTrain, newdata=dataFrameOsTD[!trainIdx,mainIdxOsTD]))$concordance,
							AIC=survConcordance(osTD[!trainIdx]~predict(coxAICOsTrain, newdata=dataFrameOsTD[!trainIdx,mainIdxOsTD]))$concordance,
							RFX=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,mainIdxOsTD]) %*% coef(coxRFXOsTrain))$concordance
					))
		}, mc.cores=20))

#' RFX models with different interaction terms
#+ systematicCVRFX, cache=TRUE
replicates <- 100 ## number of replicates
whichRFXOsTDGG <- which((colSums(dataFrame)>=8 | mainIdxOsTD) & osTDIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%
allConcordanceCVRFX <- do.call("rbind",mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrameOsTD)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			coxRFXOsMain <- CoxRFX(dataFrameOsTD[trainIdx,mainIdxOsTD], osTD[trainIdx], groups=groups[mainIdxOsTD])
			coxRFXOsGG <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTDGG], osTD[trainIdx], groups=groups[whichRFXOsTDGG])
			coxRFXOsAll <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTD], osTD[trainIdx], groups=groups[whichRFXOsTD])
			return(c(
							Main=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,mainIdxOsTD]) %*% coef(coxRFXOsMain))$concordance,
							GeneGene=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTDGG]) %*% coef(coxRFXOsGG))$concordance,
							AllInt=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTD]) %*% coef(coxRFXOsAll))$concordance	
							))
		}, mc.cores=20))

#' ### 6. Models on genomics only
#+ genomicModels, cache=TRUE
genomicsIndex <- groups %in% c("Genetics","Cytogenetics", "Translocations")
coxRFXOsGenomics <- CoxRFX(dataFrame[genomicsIndex], os, groups=groups[genomicsIndex])
coxCPSSOsGenomics <- CoxCPSSInteractions(dataFrame[genomicsIndex], os, bootstrap.samples=500)
coxCPSSOsGenomics

#' #### 5 fold cross validation
#+ genomicModelsCV, cache=TRUE
set.seed(42)
cvIndex <- sample(1:5,nrow(dataFrame), replace=TRUE)
cvRFXOsGenomics <- sapply(1:50, function(x){
			i <- 1
			cvIndex <- sample(1:5,nrow(dataFrame), replace=TRUE)
			fit <- CoxRFX(dataFrame[cvIndex!=i,genomicsIndex], os[cvIndex!=i], groups=groups[genomicsIndex])
			survConcordance(os[cvIndex==i] ~ as.matrix(dataFrame[cvIndex==i,genomicsIndex]) %*% coef(fit))$concordance
		}) 
cvCPSSOsGenomics <- sapply(1:10, function(i){
			i <- 1
			cvIndex <- sample(1:5,nrow(dataFrame), replace=TRUE)
			fit <- CoxCPSS(dataFrame[cvIndex!=i,genomicsIndex], os[cvIndex!=i], control="BH")
			survConcordance(os[cvIndex==i] ~ predict(fit, newdata=dataFrame[cvIndex==i,genomicsIndex]))$concordance
		}) 

#' #### Genomic risk groups
#+ groupsCPSS, fig.width=3, fig.height=2.5
genomicRiskCoxCPSSInt <- predict(coxCPSSOsGenomics)
genomicRiskGroupsCoxCPSSInt <-  cut(genomicRiskCoxCPSSInt[!is.na(genomicRiskCoxCPSSInt)], quantile(genomicRiskCoxCPSSInt, seq(0,1,0.25)), labels = c("very low","low","high","very high"), include.lowest=TRUE)
survConcordance( os~as.numeric(genomicRiskGroupsCoxCPSSInt))
plot(survfit(os[!is.na(genomicRiskCoxCPSSInt)] ~ genomicRiskGroupsCoxCPSSInt), col=set1[c(3,2,4,1)])
table(clinicalData$M_Risk, genomicRiskGroupsCoxCPSSInt[1:nrow(clinicalData)])[c(2,4,3,1),]

#' Partial risk for c
names(groups) <- colnames(dataFrame)
X <- model.matrix(coxCPSSOsGenomics$coxph$formula, data=dataFrame)[,-1]
#X <- X - rep(colMeans(X), each=nrow(X))
partRiskCoxCPSSInt <- sapply(unique(groups[colnames(X)]), function(l) {
			w <- groups[colnames(X)] == l
			if(length(w)>0)
				X[,w, drop=FALSE] %*% coef(coxCPSSOsGenomics$coxph)[w]
			else
				rep(0, nrow(X))
		})
colnames(partRiskCoxCPSSInt) <- as.character(unique(groups[colnames(X)]))

#' Risk Plots
#+ genomicRiskCoxCPSSInt, fig.width=6, fig.height=6
par(mfrow=c(4,1))
i <- 1
for(l in levels(clinicalData$M_Risk)[c(2,3,4,1)]){
	barplot(sapply(split(as.data.frame(partRiskCoxCPSSInt[1:nrow(clinicalData),][clinicalData$M_Risk ==l,]), genomicRiskGroupsCoxCPSSInt[1:nrow(clinicalData)][clinicalData$M_Risk ==l] ), colMeans), beside=TRUE, legend=i==4, col=colGroups[1:3], main=l, xlim=c(1,30))
	i <- i+1
}


#' Visualisation of risk
#+ riskCPSSIntOs, fig.width=3, fig.height=2.5
p <- predict(coxCPSSOsGenomics, newdata=dataFrame)
x <- seq(from=-4,to=4, l=512)
r <- sapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(r){
			i <- clinicalData$M_Risk==r
			d <- density(na.omit(p[i]), from=-4,to=4)$y * mean(i, na.rm=TRUE)
		})
par(mar=c(4,4,3,4)+.1, bty="n")
plot(exp(x),rowSums(r), type='l', lty=0,xlab="Hazard", ylab="Prop. patients", log='x')
for(i in 1:4)
	polygon(exp(c(x, rev(x))), c(rowSums(r[,1:i, drop=FALSE]), rev(rowSums(cbind(0,r)[,1:i, drop=FALSE]))), col=set1[c(3,2,4,1)][i], border=NA)

H0 <- basehaz(coxCPSSOsGenomics$coxph, centered=TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
s <- par("usr")[4]/5000
for(i in seq_along(l))
	lines(exp(x), invHazardDist(-log(l[i]) /exp(x) )*s, col='black', lty=c(2,1,2)[i])
p <- pretty(c(0,5000))
axis(side=4, at=p*s, labels=p)
mtext(side=4, "Time", line=2.5)
mtext(side=3, at = -log(l)/hazardDist(par("usr")[4]/s), text=paste(100*l, "% survive", sep=""))
legend("topright", levels(clinicalData$M_Risk)[c(2,4,3,1)], fill=set1[c(3,2,4,1)], bty="n", title="M risk")

boxplot(predict(coxCPSSOsGenomics, newdata=dataFrame) ~ clinicalData$M_Risk, horizontal=TRUE)

#' #### BIC
#+ coxBICOsGenomcs, cache=TRUE, warning=FALSE
c <- coxph(os ~ 1, data=dataFrame[,genomicsIndex])
scopeStep <- as.formula(paste("os ~", paste(colnames(dataFrame)[genomicsIndex], collapse="+")))
coxBICOsGenomics <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
summary(coxBICOsGenomics)


#' #### AIC
#+ coxAICOsGenomcs, cache=TRUE, warning=FALSE
coxAICOsGenomics <- step(c, scope= scopeStep, k = 2, trace=0)
summary(coxAICOsGenomics)

#+ vennDiagramGenomics, fig.width=2.5, fig.height=2.5
library(VennDiagram)
grid.newpage()
pushViewport(viewport(w = .9, h = .9))
grid.draw(venn.diagram(list(CPSS=names(coef(coxCPSSOsGenomics$coxph)), BIC=names(coef(coxBICOsGenomics)), AIC=names(coef(coxAICOsGenomics))), filename=NULL, lty=1, 
				col=set1[1:3], fill=set1[1:3], alpha=0.3, euler.d=TRUE, fontfamily="Helvetica", cat.fontfamily="Helvetica", euler.diagram=TRUE))



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
t <- data.frame(tcgaMutation[,]>0, CEBPA_mono = tcgaMutation[,"CEBPA"]==1,CEBPA_bi = tcgaMutation[,"CEBPA"]>1,tcgaClinical[,-c(1,2,4,5,6,13,25)], MakeInteger(tcgaClinical$TypeAML)) + 0
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
#+ tcgaFreq, fig.width=2.5, fig.height=2.5
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
tcgaRiskCPSSOs <- predict(coxCPSSIntOs, newdata=tcgaDataImputed)
survConcordance(tcgaSurvival ~ tcgaRiskCPSSOs)

#' Blind imputation (mean only)
f <- function(X) {X <- sapply(X, poorMansImpute);X[is.na(X)] <- 0; X}
survConcordance(tcgaSurvival ~ predict(coxCPSSIntOs, newdata=as.data.frame(f(tcgaData[mainIdxOs]))))

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

#+ PINAos, fig.width=3, fig.height=2.5
nkIdx <- clinicalData$NK == 1
plot(survfit(os[nkIdx] ~ pinaOs[nkIdx,2]), col=rev(set1[1:3]))
survConcordance(os[nkIdx] ~ pinaOs[nkIdx,1])

#' Compared to CPSS (AML data)
survConcordance(os[nkIdx] ~ predict(coxCPSSIntOs, newdata=dataFrame)[nkIdx])
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

#+ concordanceTCGA, fig.width=3, fig.height=2.5
tcgaConcordance <- sapply(tcgaRisk, function(x) {c <- survConcordance(tcgaSurvival ~ x); c(c$concordance, c$std.err)})
tcgaConcordance
o <- order(tcgaConcordance[1,])
barplot(tcgaConcordance[1,o], border=NA, col= set1[-6], las=2, xaxt="n", ylab="Concordance", ylim=c(0.5,0.75), xpd=FALSE) -> b
segments(b,tcgaConcordance[1,o]-tcgaConcordance[2,o],b,tcgaConcordance[1,o]+tcgaConcordance[2,o])
rotatedLabel(b, rep(0.49,length(b)), colnames(tcgaConcordance)[o], srt=45)

#+ aucTCGA, fig.width=3, fig.height=2.5
library(survivalROC)
#aucCV <- sapply(predictedRiskCV, function(x) survivalROC(Stime=os[!is.na(os) & !trainIdx,1], status=os[!is.na(os) &  !trainIdx,2], marker = scale(x[!is.na(os[!trainIdx])]), predict.time = 278, method="KM", cut.values=seq(-5,5,0.1))$AUC)
tcgaAUC <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], c(90,365,1000))$auc)
o <- order(colMeans(tcgaAUC))
barplot(tcgaAUC[,o], border=1, col= rep(set1[-6],each=3), las=2, xaxt="n", ylab="AUC", beside=TRUE, density=c(NA, 48,24), ylim=c(0.5,0.85), xpd=FALSE) -> b
legend("topleft", bty="n", c("3mo","1yr","3yr"), fill='black', density=c(NA, 48,24))
rotatedLabel(b[seq(3, length(b), 3)], rep(0.49,length(tcgaRisk)), names(tcgaRisk)[o], srt=45)


#+ kmTCGA, fig.width=3, fig.height=2.5
risk <- cut(tcgaRiskCPSSOs, quantile(tcgaRiskCPSSOs), labels=c("1st Q","2nd Q","3rd Q","4th Q"))
s <- survfit(tcgaSurvival ~ risk)
plot(s, col=set1[c(3,2,4,1)], mark=NA)
legend("topright", bty="n", rownames(summary(s)$table), col=set1[c(3,2,4,1)], lty=1)

#+ riskTCGA, fig.width=3, fig.height=2.5
risk <- tcgaRiskCPSSOs - mean(tcgaRiskCPSSOs)
x <- seq(from=-4,to=4, l=512)
r <- sapply(levels(tcgaClinical$C_Risk)[c(2,3,1)], function(r){
			i <- tcgaClinical$C_Risk==r
			d <- density(na.omit(risk[i]), from=-4,to=4)$y * mean(i, na.rm=TRUE)
		})
par(mar=c(4,4,3,4)+.1, bty="n")
plot(exp(x),rowSums(r), type='l', lty=0,xlab="Hazard", ylab="Prop. patients", log='x', ylim=c(0,.55))
for(i in 1:3)
	polygon(exp(c(x, rev(x))), c(rowSums(r[,1:i, drop=FALSE]), rev(rowSums(cbind(0,r)[,1:i, drop=FALSE]))), col=set1[c(3,2,1)][i], border=NA)

H0 <- basehaz(coxph(tcgaSurvival ~ risk), centered=TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
for(i in seq_along(l))
	lines(exp(x), pmax(0,invHazardDist(-log(l[i]) /exp(x) ))/10000, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
mtext(side=4, "Time", line=2.5)
mtext(side=3, at = -log(l)/hazardDist(par("usr")[4]*10000), text=paste(100*l, "% survive", sep=""))
legend("topright", levels(tcgaClinical$C_Risk)[c(2,3,1)], fill=set1[c(3,2,1)], bty="n", title="M risk")

#' #### Genomic models
tcgaRiskGenomics <- data.frame(
		stdRisk = c(3,1,2)[tcgaClinical$C_Risk],
		coxBIC = predict(coxBICOsGenomics, newdata=tcgaDataImputed),
		coxAIC = predict(coxAICOsGenomics, newdata=tcgaDataImputed),
		coxCPSS = predict(coxCPSSOsGenomics, newdata=tcgaDataImputed)
)

#+ concordanceTCGAGenomics, fig.width=2.5, fig.height=2.5
tcgaConcordanceGenomics <- sapply(tcgaRiskGenomics, function(x) {c <- survConcordance(tcgaSurvival ~ x); c(c$concordance, c$std.err)})
tcgaConcordanceGenomics
o <- order(tcgaConcordanceGenomics[1,])
barplot(tcgaConcordanceGenomics[1,o], border=NA, col= set1[-6], las=2, xaxt="n", ylab="Concordance", ylim=c(0.5,0.75), xpd=FALSE) -> b
segments(b,tcgaConcordanceGenomics[1,o]-tcgaConcordanceGenomics[2,o],b,tcgaConcordanceGenomics[1,o]+tcgaConcordanceGenomics[2,o])
rotatedLabel(b, rep(0.49,length(b)), colnames(tcgaConcordanceGenomics)[o], srt=45)

#+ kmTCGAGenomics, fig.width=3, fig.height=2.5
risk <- cut(tcgaRiskGenomics$coxCPSS, quantile(tcgaRiskGenomics$coxCPSS), labels=c("1st Q","2nd Q","3rd Q","4th Q"))
s <- survfit(tcgaSurvival ~ risk)
plot(s, col=set1[c(3,2,4,1)], mark=NA)
legend("topright", bty="n", rownames(summary(s)$table), col=set1[c(3,2,4,1)], lty=1)


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
coxCPSSIntOsNonWHO <- CoxCPSSInteractions(dataFrame[w,groups %in% mainGroups & osIdx], os[w], bootstrap.samples=50, scope = which(groups %in% scope))
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

tdData <- sapply(1:nrow(clinicalData), function(i){
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
tdData <- as.data.frame(do.call("rbind",tdData))
tdData$Age <- dataFrame$AOD_10[tdData$i] + tdData$time1/3650

#' A coxph model
summary(coxph(Surv(time1, time2, event) ~ . + TPL:CR + REC:TPL , data=tdData[,-8]))

#' With CPSS selected terms
summary(coxph(Surv(time1, time2, event) ~ ., data=cbind(tdData[,-8], dataFrame[tdData$i, grep("(TPL)|(:)", selectedIntOs, value=TRUE, invert = TRUE)])))

#' ### 4. Splice factor AML
whichSplice <- c("SFRS2","SF3B1","U2AF1","ZRSR2","U2AF2", "SF1","SF3A1")

#' Frequency
#+ spliceBar, fig.width=3, fig.height=2.5
colSums(dataFrame[whichSplice])
barplot(colSums(dataFrame[whichSplice]), col=set1, ylab="# Cases")
spliceSum <- rowSums(genes[,whichSplice])
table(spliceSum)
spliceAny <- spliceSum > 0
spliceFactor <- factor(apply(dataFrame[whichSplice], 1, function(x) {if(sum(x)==0) "none" else if(sum(x)==2) "multiple" else whichSplice[x==1]}))
dataFrame[spliceSum==2, whichSplice]

#' Overlaps
table(spliceAny, clinicalData$WHOcat)
f <- factor(paste(spliceAny,clinicalData$WHOcat)[!is.na(clinicalData$WHOcat)], levels=c("TRUE no","FALSE no" , "FALSE yes" ,"TRUE yes"  ))
pie(table(f), col=set1[c(3,1,1,3)], density=c(NA,NA,36,36), angle=45, init.angle=90, labels=table(f))
legend("topleft", bty="n", fill=c(set1[c(3,1)],"black"), density=c(NA,NA,36), c("SF mutations", "no SF mutations", "rec. WHO class"))

table(spliceAny, dataFrame$sAML)
fisher.test(table(spliceAny, dataFrame$sAML))
cor(dataFrame[selectedIntOs], spliceAny)

#' Age
#+ spliceAge, fig.width=1.5, fig.height=2.5
boxplot(dataFrame$AOD_10*10 ~ spliceSum, ylab ="Age", xlab="SF mutations", staplewex=0, lty=1)
summary(glm(dataFrame$AOD_10*10 ~ spliceAny + clinicalData$WHOcat))

#+ spliceBlast, fig.width=1.5, fig.height=2.5
boxplot(dataFrame$BM_Blasts_100 * 100 ~ spliceAny, ylab ="% Bone marrow blasts", xlab="SF mutations", staplewex=0, lty=1)
summary(glm(dataFrame$BM_Blasts_100 * 100 ~ spliceAny + clinicalData$WHOcat))


#' Survival
#+ spliceKM, fig.width=3, fig.height=2.5
coxph(os ~ .,data=dataFrame[whichSplice[-6]])
coxph(os ~ ., data=cbind(dataFrame[selectedIntOs], SF = spliceAny+0))
plot(survfit(os ~ spliceAny), mark=NA, col=c("grey","black"))
i <- 1
for(w in whichSplice){
	lines(survfit(os ~ 1, subset=dataFrame[w]==1), mark=NA, col=set1[i], conf.int=FALSE)
	i <- i+1
}
legend("topright", bty="n", col=c("grey", set1[seq_along(whichSplice)], "black"), legend = c("none", whichSplice,"any"), lty=1)

#' CPSS
#+ coxCPSSOsSplice, cache=TRUE
set.seed(42)
coxCPSSOsSplice <- CoxCPSSInteractions(cbind(dataFrame[mainIdxOs], SF=spliceAny+0), os, scope = which(groups %in% scope))
coxCPSSOsSplice

#' TCGA
#+ tcgaSF, fig.width=3, fig.height=2.5
tcgaAnySplice <- rowSums(tcgaDataImputed[whichSplice]>0)
plot(survfit(tcgaSurvival ~ tcgaAnySplice), mark=NA, col=c("grey","black"))
i <- 1
for(w in whichSplice){
	if(sum(tcgaDataImputed[w]==1)>0)
		lines(survfit(tcgaSurvival ~ 1, subset=tcgaDataImputed[w]==1), mark=NA, col=set1[i], conf.int=FALSE)
	i <- i+1
}
legend("topright", bty="n", col=c("grey", set1[seq_along(whichSplice)], "black"), legend = c("none", whichSplice,"any"), lty=1)

coxph(tcgaSurvival ~ tcgaAnySplice)

survConcordance(tcgaSurvival ~ predict(coxCPSSOsSplice, newdata=cbind(tcgaDataImputed, SF = tcgaAnySplice)))
tcgaConcordance

#' ### 5. Predicted variance
#' #### Cumulative risk with interaction terms
#+ cumRiskGeneGene, fig.width=3, fig.height=2.5
X <- dataFrame[groups=="Genetics"]
genRisk <- X * rep(coef(coxRFXFitOs)[coxRFXFitOs$groups=="Genetics"], each=nrow(dataFrame))
geneGeneRisk <- dataFrame[whichRFXOs][coxRFXFitOs$groups=="GeneGene"] * rep(coef(coxRFXFitOs)[coxRFXFitOs$groups=="GeneGene"], each=nrow(dataFrame)) 
o <- order(colSums(X), decreasing=TRUE)
genRisk <- genRisk[,o]
cumRisk <- apply(genRisk, 1, cumsum)
cumVar <- apply(cumRisk,1,var)
s <- strsplit(colnames(geneGeneRisk), ":")
n <- colnames(genRisk)
r <- 1:ncol(X)
names(r)  <- names(X)[o]
i <- sapply(s, function(i){
			max(r[i])
		})
cumRiskInt <- t(cumRisk) + sapply(1:ncol(genRisk), function(x){
			rowSums(geneGeneRisk[,which(i<=x), drop=FALSE])
		})
cumVarInt <- apply(cumRiskInt,2,var)
plot(cumsum(colMeans(X)[o]), cumVar, xlim=c(0,4), ylim=c(0,.3), xlab="Mean number of driver mutations/patient", ylab="Predicted variance")
a1 <- coxRFXFitOs$sigma2["Genetics"] + coxRFXFitOs$mu["Genetics"]^2
a2 <- var(coef(coxRFXFitOs)[coxRFXFitOs$groups=="Genetics"])+ coxRFXFitOs$mu["Genetics"]^2
abline(0, a1, lty=2)
abline(0, a2, lty=3)
points(3.7, 3.7*a2, col='red')
points(3.7, 3.7*a1, col='red', pch=19)
points(cumsum(colMeans(X)[o]), cumVarInt, pch=19)

#' #### Cumulative risk with main terms
#+ cumRiskGenes, fig.width=3, fig.height=2.5
X <- dataFrame[groups=="Genetics"]
genRisk <- X * rep(coef(coxRFXFitOsMain)[coxRFXFitOsMain$groups=="Genetics"], each=nrow(dataFrame))
o <- order(colSums(X), decreasing=TRUE)
genRisk <- genRisk[,o]
cumRisk <- apply(genRisk, 1, cumsum)
cumVar <- apply(cumRisk,1,var)
plot(c(0,cumsum(colMeans(X)[o])), c(0,cumVar), xlim=c(0,4), ylim=c(0,.3), xlab="Mean number of driver mutations/patient", ylab="Predicted variance", pch=19, type='b')
a1 <- coxRFXFitOsMain$sigma2["Genetics"] + coxRFXFitOsMain$mu["Genetics"]^2
a2 <- mean(coef(coxRFXFitOsMain)[coxRFXFitOsMain$groups=="Genetics"]^2)
abline(0, a1, lty=2)
abline(0, a2, lty=3)
points(3.7, 3.7*a2, col='red')
points(3.7, 3.7*a1, col='red', pch=19)
axis(side=1, at=c(cumsum(colMeans(X)[o])[c(1:10,15,20,30,40,50)], 3.7), labels=c(1:10,15,20,30,40,50, "Exome (TCGA)"), tcl=.5, las=2, hadj=0, mgp = c(-1.5,-0.5,0))
mtext(side=1, line=-3, "Genes sequenced", las=1)

#' ### 6. NPM1 FLT3_ITD interaction
#' There is a reported interaction term for NPM1:FLT3_ITD, which can be explained by covarying factors:
#+ intNPM1vFLT3, fig.width=3, fig.height=2.5
m <- coxph(os ~  NPM1 + FLT3_ITD + NPM1:FLT3_ITD, data=dataFrame); m
c <- coef(m)
m <- coxph(os ~  NPM1 + FLT3_ITD, data=dataFrame); m
d <- coef(m)
f <- factor(paste(c("-","NPM1")[dataFrame$NPM1+1],c("-","FLT3_ITD")[dataFrame$FLT3_ITD+1], sep=" / "))
s <- survfit(formula(m), data=dataFrame)
plot(s, col=set1)
legend("topright", col=set1, lty=1, rownames(summary(s)$table), bty='n')

#+ intNPM1vFLT3box, fig.width=2, fig.height=2.5
par(mar=c(6,4,1,1))
boxplot(predict(coxRFXFitOsMain) ~ f , staplewex=0, lty=1, border=brewer.pal(4,"Pastel1"), ylab="log hazard", xaxt="n", ylim=c(-2.5,3.5), pch=".")
rotatedLabel(1:4, rep(min(predict(coxRFXFitOsMain)),4), levels(f))
points(c(0, c[2], c[1], sum(c[1:3])), col=set1, pch=19)
points(c(0, d[2], d[1], sum(d)), col=set1)
arrows(4,sum(d),4,sum(c), length=.2, col=set1[4])
legend("topright", pch=c(46,1,19), col=c("grey","red","red"), c("all variables","NPM1 + FLT3 only", "NPM1 + FLT3 + interaction"), bty="n" )

#' ### 7. BT model
copyNumbers = cbind(dataList$Cytogenetics[grep(c("minus|plus|mono"), colnames(dataList$Cytogenetics))], clinicalData$gender)
copyNumbers$minus7 <- (copyNumbers$minus7 | copyNumbers$minus7q) +0
copyNumbers$minus7q <- NULL
for(i in 1:ncol(copyNumbers)){
	if(grepl("plus", colnames(copyNumbers)[i]))
		copyNumbers[,i] = copyNumbers[,i] * 3
}
copyNumbers[copyNumbers==0 | is.na(copyNumbers)] = 2
colnames(copyNumbers) = c(5,7,8,9,12,13,17,18,20,21,22,"Y",11,4,"X")
rownames(copyNumbers) <- clinicalData$PDID
copyNumbers$Y <- copyNumbers$Y - c(2,1)[copyNumbers$X]

cn = sapply(1:nrow(mutationData), function(i) {c=copyNumbers[mutationData$SAMPLE_NAME[i],match(mutationData$CHR[i], colnames(copyNumbers))]; if(length(c)==0) 2 else c})
mutationData$VAF <- as.numeric(as.character(mutationData$VAF))
mutationData$TUM_DEPTH <- as.numeric(as.character(mutationData$TUM_DEPTH))

vafCn <- mutationData$VAF/100*cn ## Approx mutant cell fraction
vafCn[which(vafCn > 1.25)] <- mutationData$VAF[which(vafCn > 1.25)] ## Wrong CN adjust
vafCn[vafCn > 1] <- 1 ## Random
vaf <- mutationData$VAF/100 
precedence <- matrix(0, ncol=ncol(mutationTable), nrow=ncol(mutationTable), dimnames=list(colnames(mutationTable), colnames(mutationTable)))
plist <- list()
ix=mutationData$GENE %in% colnames(precedence) & mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & !mutationData$CONSEQUENCE %in% c("ITD","PTD","deletion","insertion", "deletion_frameshift","insertion_frameshift","na")
for(s in clinicalData$PDID){
	l <- list()
	for(i in which(mutationData$SAMPLE_NAME==s & ix))
		for(j in which(mutationData$SAMPLE_NAME==s & ix)){
			if(!is.na(cn[i]) & !is.na(cn[j]) & i!=j){
				m <- round(matrix(c(
										vafCn[i]*mutationData$TUM_DEPTH[i],
										mutationData$TUM_DEPTH[i]-vafCn[i]*mutationData$TUM_DEPTH[i], 
										vafCn[j]*mutationData$TUM_DEPTH[j],
										mutationData$TUM_DEPTH[j]-vafCn[j]*mutationData$TUM_DEPTH[j]),
								ncol=2))
				f <- try(fisher.test(m, alternative="greater")$p.value< 0.05 , silent=TRUE) ## Fisher test
				if(class(f)!="try-error")
					if(f & vafCn[i] >= 1 - vafCn[j]){ ## Pidgeonhole
						precedence[as.character(mutationData$GENE[i]),as.character(mutationData$GENE[j])] <- precedence[as.character(mutationData$GENE[i]),as.character(mutationData$GENE[j])] + 1
						l <- c(l, list(c(as.character(mutationData$GENE[i]),as.character(mutationData$GENE[j]))))		
					}
				
			}
		}
	plist[[s]] <- l
}

makeDesign <- function(I) {
	w <- which(lower.tri(I), arr.ind=TRUE)
	x <- matrix(0, nrow(w), nrow(I))
	for(i in 1:nrow(w)){
		x[i,w[i,1]] <- 1
		x[i,w[i,2]] <- -1
	}
	return(x)
}

btModel <- function(I){
	y <- cbind(I[lower.tri(I)], t(I)[lower.tri(I)])
	x <- makeDesign(I = I)
	glm.fit(x=x[,-1],y=y, family=binomial())
}

nCasesGene <- table(factor(unlist(sapply(plist, function(x) unique(unlist(x)))), levels=colnames(precedence)))
w <- which(nCasesGene > 5)

fit <- btModel(precedence[w,w]+.01)
c <- c(0,coef(fit))
names(c) <- colnames(precedence)[w]
o <- rank(c)
v <- pmin(2,sqrt(c(0,diag(chol2inv(fit$qr$qr)))))

#+ bradleyTerryAML, fig.width=4, fig.height=2.5
l <- names(c)
m <- paste("n=",nCasesGene[w], sep="")
plot(-c, o, xlab="Relative time", yaxt="n", pch=19, col="grey", ylab="", xlim=range(-c+3*c(-v,v)))
segments(-c-v, o,-c+v,o, col="grey")
text(-c-v ,o,l, font=3, pos=2)
text(-c+v ,o,m, font=1, pos=4)

#sapply(split(mutationData$CHR, mutationData$GENE) ,unique)

#+ subcloneAML, fig.width=2.5, fig.height=2.5
t <- table(sapply(plist, length)>0)
pie(t, labels=paste(t, c("clonal/NA","polyclonal")))

#+ vafAML_NPM1_DNMT3A_FLT3, fig.width=2, fig.height=2.5
vaf <- t(sapply(split(mutationData[,c("VAF","GENE")], mutationData$SAMPLE_NAME), function(x) x$VAF[match(levels(mutationData$GENE), x$GENE)]))
colnames(vaf) <- levels(mutationData$GENE)
boxplot(vaf[,c("NPM1","DNMT3A","FLT3")], staplewex=0, pch=16, lty=1, ylab="VAF")
pastel1 <- brewer.pal(9,"Pastel1")
segments(2,vaf[, "DNMT3A"],1,vaf[, "NPM1"], col=pastel1[1])
segments(2,vaf[, "DNMT3A"],3,vaf[, "FLT3"], col=pastel1[2])
segments(1,vaf[, "NPM1"],3,vaf[, "FLT3"], col=pastel1[3])
m <- colMeans(vaf[!is.na(rowSums(vaf[,c("NPM1","DNMT3A","FLT3")])),c("NPM1","DNMT3A","FLT3")])
segments(2,m["DNMT3A"],1,m["NPM1"], col=set1[1],lwd=2)
segments(2,m["DNMT3A"],3,m["FLT3"], col=set1[2], lwd=2)
segments(1,m["NPM1"],3,m["FLT3"], col=set1[3], lwd=2)
precedence[c("NPM1","DNMT3A","FLT3"),c("NPM1","DNMT3A","FLT3")]

f <- sapply(plist, function(x) if("NPM1" %in% sapply(x,  `[`,1) ) "clonal" else  if("NPM1" %in% sapply(x, `[`,2))  "subclonal" else "absent")
f[dataFrame$NPM1==1 & f=="absent"] <- "n/a"

#+ vafAML_KRAS_NRAS, fig.width=1.5, fig.height=2.5
boxplot(vaf[,c("KRAS","NRAS")], staplewex=0, pch=16, lty=1, ylab="VAF")
pastel1 <- brewer.pal(9,"Pastel1")
segments(1,vaf[, "KRAS"],2,vaf[, "NRAS"], col=pastel1[1])
m <- colMeans(vaf[!is.na(rowSums(vaf[,c("KRAS","NRAS")])),c("KRAS","NRAS")])
segments(1,m["KRAS"],2,m["NRAS"], col=set1[1], lwd=2)
precedence[c("KRAS","NRAS"),c("KRAS","NRAS")]

#' 8. Regression of blood counts
#' -----------------------------
#' ### Prepare data
#+ clinicalGlmnet, cache=TRUE
library(glmnet)
Y <- StandardizeMagnitude(cbind(dataList$Clinical, dataList$Demographics))
X <- as.matrix(dataFrame[groups %in% c("Genetics","Cytogenetics","Translocations")])
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
	plotcvnet(m, X, main=names(clinModels)[i],  col0="black", cex=1, simple.annot = annot, col=pastel1[c(3,2,4)], xlim=c(0.5,35.5))
	i = i+1
	legend("topright", col=c(pastel1[c(1,3)],"black")[c(1,3,2)], c(expression(paste("Explained variance ",R^2)), expression(paste("Lasso penalty ",lambda)), expression(paste("Model coefficient ", beta))), box.lty=0, bg="#FFFFFF33", pch=c(NA,NA,19), lty=c(1,1,NA), cex=.8, pt.cex = 1)
}

#' ### Heatmap of GLMs
j <- 0
z <- sapply(clinModels,function(x){ ## Creating matrix
			j <<- j+1
			w <- which.min(x$cvm)
			c <- x$glmnet.fit$beta[,w]
			yj <- sapply(c("Genetics","Cytogenetics","Translocations"), function(i){
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
rotatedLabel(y0=rep(0.5,nrow(z)), labels=gsub("_","/",rownames(z))[o], x0=1:nrow(z), font=c(rep(3,sum(groups=="Genetics")),rep(1,sum(groups%in%c("Cytogenetics","Translocations")))), col = pastel1[c(3,2,4)][annot], cex=0.9)
mtext(side=2, line=.2, text=colnames(z)[w], las=2, at=1:ncol(z)-.5)
text(y=rep(1:ncol(z)-.5, each=nrow(r)), x=rep(1:nrow(r), ncol(z)), r[o,w] * (0!=(z[o,w])), cex=0.66, font=ifelse(r[o,w] <= rep(s[ncol(r):1], each=nrow(r)), 2,1))
points(y=rep(1:ncol(z)-.5, each=nrow(r)), x=rep(1:nrow(r), ncol(z)), pch=ifelse(is.na(z[o,w]) | z[o,w]==0, ".",NA))
mtext(side=1, at=sum(groups=="Genetics")/2, "Genetics", col=pastel1[2], line=2.5 )
mtext(side=1, at=sum(groups=="Cytogenetics")/2 + sum(groups=="Genetics"), "Cytogenetics", col=pastel1[3], line=2.5 )
mtext(side=1, at=sum(groups=="Tranlocations")/2 + sum(groups%in%c("Genetics","Translocations")), "Translocations", col=pastel1[3], line=2.5 )
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
barplot(R2[1:2,w], border=NA, col=paste(pastel1[c(2,3,4)],"88",sep=""), horiz=TRUE, names.arg=rep(NA,ncol(R2)), width=0.95, space=0.0525, add=TRUE) -> b
points(R2[3,w]+0.1,b, pch=ifelse(p[w],"*",NA))
mtext(side=3, "Variance components", line=.5)
mtext(side=1, expression(paste("Explained variance ",R^2)), line=2.5)

#+ save
#save(list = ls(), file=paste(Sys.Date(), "-AML.RData", sep=""))


#' 8. Data for web tool
#' --------------------
#' ### Prelim
#' Times for allografts pre and post relapse, after 1CR only
alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD") # only allografts
alloTime1CR <- clinicalData$Time_1CR_TPL + .5 # +.5 to make > 0
alloTime1CR[!alloIdx] <- NA
alloTimeRel <- clinicalData$TPL_date - clinicalData$Recurrence_date + .5
alloTimeRel[!alloIdx] <- NA
alloTimeRel[alloTimeRel <0] <- NA

#' Create data frames for each phase
w <- mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Nuisance"
t <- clinicalData$Recurrence_date
t[is.na(t)] <- as.Date(1e6, origin="2000-01-01")
cirData <- MakeTimeDependent(dataFrame[w], timeEvent=alloTime1CR, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=!is.na(clinicalData$Recurrence_date)+0)
nrmData <- MakeTimeDependent(dataFrame[w], timeEvent=alloTime1CR, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=is.na(clinicalData$Recurrence_date) & clinicalData$Status)
#prsData <- makeTimeDependent(dataFrame[w], timeTpl=alloTimeRel, timeSurv=as.numeric(clinicalData$Date_LF - clinicalData$Recurrence_date)+1, event=clinicalData$Status)
prsData <- MakeTimeDependent(dataFrame[w], timeEvent=alloTime1CR, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), timeStart = as.numeric(clinicalData$Recurrence_date- clinicalData$CR_date), status=clinicalData$Status)
cirData$transplant1CR <- cirData$event
cirData$event <- NULL
cirData$transplantRel <- 0
nrmData$transplant1CR <- nrmData$event
nrmData$event <- NULL
nrmData$transplantRel <- 0
prsData$transplant1CR <- nrmData$transplant1CR[1:nrow(dataFrame)][prsData$index]
prsData$transplantRel <- prsData$event
prsData$event <- NULL

#' Fit models
crGroups <- c(as.character(groups[w]), "Treatment","Treatment")
names(crGroups) <- c(names(dataFrame)[w],"transplant1CR","transplantRel")
coxRFXNrmTD <- CoxRFX(nrmData[names(crGroups)], Surv(nrmData$time1, nrmData$time2, nrmData$status), groups=crGroups)
coxRFXNrmTD$coefficients["transplantRel"] <- 0
prsData$time1[!is.na(prsData$time1)] <- 0
coxRFXPrsTD <-  CoxRFX(prsData[names(crGroups)], Surv(prsData$time1, prsData$time2, prsData$status), groups=crGroups, sigma0 = 0.1, nu=0.1)
coxRFXCirTD <-  CoxRFX(cirData[names(crGroups)], Surv(cirData$time1, cirData$time2, cirData$status), groups=crGroups)
coxRFXCirTD$coefficients["transplantRel"] <- 0

#' Plots
#+ cirSplits, fig.width=2.5, fig.height=2.5
par(mar=c(3,3,2,1), bty="n", mgp=c(2,.5,0))
riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)
r <- coxRFXCirTD$X %*% coef(coxRFXCirTD) - cirData$transplant1CR * coef(coxRFXCirTD)["transplant1CR"]
f <- function(x) 1-x
plot(survfit(Surv(time1, time2, event) ~ clinicalData$M_Risk[cirData$index], data=cirData), col=riskCol, ylab="CIR", xlab="Time after CR", main="Molecular risk groups, all cases", fun=f , ylim=c(0,1))
legend("bottomright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)
Q <- numeric(nrow(cirData))
for(l in levels(clinicalData$M_Risk)){
	w <- which(clinicalData$M_Risk[cirData$index]==l)
	q <- cut(r[w], quantile(r[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
	Q[w] <- q
	plot(survfit(Surv(time1, time2, status) ~ q + transplant1CR, data=cirData[w,]), col=rep(sapply(2:0,function(x) colTrans(riskCol[l],x)), each=2), lty=c(1,3), ylab="CIR", main=paste(l, "terciles"), mark=NA, xlab="Time after CR", fun=f, ylim=c(0,1))
	legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[l])
}

#' Prevalence of risk factors
p <- lapply(levels(clinicalData$M_Risk), function(l) {
			w <- which(clinicalData$M_Risk==l)
			q <- cut(r[w], quantile(r[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
			sapply(split(cirData[w, names(crGroups)], q), colMeans)
		})
names(p) <- levels(clinicalData$M_Risk)

#+relapseFactors, fig.width=5,fig.heigh=4
par(mfrow=c(4,1), xpd=NA)
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	t <- t((p[[l]]) * coef(coxRFXCirTD))[,crGroups != "Treatment"]
	z <- (coef(coxRFXCirTD)/sqrt(diag(coxRFXCirTD$var2)))[crGroups != "Treatment"]
	o <- order(z)
	w <- c(1:15,ncol(t)-14:0)
	b <- barplot(t[,o][,w], las=2, col=sapply(2:0,function(x) colTrans(riskCol[l],x)), beside=TRUE, ylim=c(-.5,.5), names.arg=rep("", length(w)))
	rotatedLabel(b[2,],pmin(0,apply(t,2,min)[o][w]), colnames(t)[o][w])
	s <- matrix(rep(sqrt(diag(coxRFXCirTD$var2)[crGroups != "Treatment"]), each=3) * t/rep(coef(coxRFXCirTD)[crGroups != "Treatment"], each=3), nrow=3)[,o][,w]
	segments(b[1,], (colMeans(coxRFXCirTD$X)*coef(coxRFXCirTD))[crGroups != "Treatment"][o][w] ,b[3,], (colMeans(coxRFXCirTD$X)*coef(coxRFXCirTD))[crGroups != "Treatment"][o][w])
	segments(b,t[,o][,w]-s, b,t[,o][,w]+s)
}

p <- as.data.frame(PartialRisk(coxRFXCirTD)[1:nrow(clinicalData),])
s <- do.call("rbind",lapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(l) {
			w <- which(clinicalData$M_Risk==l)
			q <- cut(r[w], quantile(r[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
			t(sapply(split(p[w, ], q), colMeans) +.5 - colMeans(p))
		}))

#+relapseStars, fig.width=3,fig.heigh=3
c <- sapply(2:0, function(t) sapply(riskCol[c(2,4,3,1)], function(c) colTrans(c,t)))
g <- expand.grid(1:3,1:4-1)*3
stars(2*s[,c("Clinical","Demographics","Genetics","Cytogenetics","Translocations","Treatment")], scale=FALSE, col.stars=t(c), key.loc = c(13,0), locations=g, labels=NA)
symbols(g[,1], g[,2], circles=rep(1,12), inches=FALSE, add=TRUE)
text(1, 0:3*3, names(riskCol[c(2,4,3,1)]), pos=2)
text(1:3*3, 11, c("Low","Intermediate","High"), pos=3)

#' Find prototypes
prototypes <- sapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(l) sapply(1:3, function(i){
						d <- dist(t(t(coxRFXCirTD$X[which(clinicalData$M_Risk[cirData$index]==l & Q==i &! is.na(clinicalData$CR_date[cirData$index])), ]) ))
						as.numeric(names(which.min(rowMeans(as.matrix(d), na.rm=TRUE))))
					}))

c <- sapply(2:0, function(t) sapply(riskCol[c(2,4,3,1)], function(c) colTrans(c,t)))
g <- expand.grid(1:3,1:4-1)*3
stars(2*t(t(p[prototypes,])- colMeans(p))[,c("Clinical","Demographics","Genetics","Cytogenetics","Translocations","Treatment")] +1, scale=FALSE, col.stars=t(c), key.loc = c(13,0), locations=g, labels=NA)
symbols(g[,1], g[,2], circles=rep(1,12), inches=FALSE, add=TRUE)
text(1, 0:3*3, names(riskCol[c(2,4,3,1)]), pos=2)
text(1:3*3, 11, c("Low","Intermediate","High"), pos=3)


#' Who achieves CR?
c <- as.numeric(clinicalData$CR_date - clinicalData$ERDate)
e <- c < os[,1]
e[is.na(e)] <- 0
c[is.na(c) &! is.na(os[,1])] <- os[is.na(c) &! is.na(os[,1]),1]
cr <- Surv(time=pmin(c, os[,1]), event = e)
coxRFXCr <- CoxRFX(dataFrame[mainIdxOs &! grepl("TPL", names(dataFrame))], cr, groups=groups[mainIdxOs &! grepl("TPL", names(dataFrame))])

#' Four plots comparing different intervals
par(mfrow=c(2,2), xpd=FALSE)
PlotVarianceComponents(coxRFXCr, col=colGroups)
title(main="CR")
PlotVarianceComponents(coxRFXCirTD, col=colGroups)
title(main="CIR")
PlotVarianceComponents(coxRFXNrmTD, col=colGroups)
title(main="NRM")
PlotVarianceComponents(coxRFXPrsTD, col=colGroups)
title(main="PRS")



#' 9. Germline polymorphisms
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
#' * Gene-wise contributions to reclassifaction
