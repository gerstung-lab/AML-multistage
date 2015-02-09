#' AML data analysis
#' =================
#' Code run on 
system("hostname -f", intern=TRUE)
#' at
Sys.time()
#' using
#+ run, eval=FALSE, echo=TRUE
library(knitr)
spin("../../code/AML.R")


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
library(CoxHD)
library(mg14)
set1 <- brewer.pal(8, "Set1")

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

#' ### 1. Survival objects
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

tplIndexEfs <- (clinicalData$Time_Diag_TPL < clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL)) | (is.na(clinicalData$efs & !is.na(clinicalData$Time_Diag_TPL))) # Need to restrict to those with transplant prior to event..

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
osYr <- os
osYr[,1] <- osYr[,1]/365
osYrTD <- osTD
osYrTD[,1] <- osYrTD[,1]/365

#' ### 2. Covariates
#' #### All data as list
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
		Cytogenetics = clinicalData[,50:74],
		Nuisance = data.frame( MakeInteger(clinicalData$Study)[,1:2], Date=scale(as.numeric(clinicalData$ERDate), scale=FALSE), MissingCyto=is.na(clinicalData$t_15_17)+0),
		Treatment = data.frame(ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL_efs=tplIndexEfs, TPL_os=tplIndexOs),
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
groups[grepl("^(t_)|(inv)", colnames(dataFrame)) &! grepl(":", colnames(dataFrame))] <- "BT"
groups[groups=="Cytogenetics"] <- "CNA"
groups <- factor(groups)
names(groups) <- colnames(dataFrame)
table(groups)

#' Poor man's imputation by column means
#+ dataFrame, cache=TRUE
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
dataFrame <- as.data.frame(sapply(dataFrame, poorMansImpute))

#' ### Gene:Gene interactions
#+ interactions, cache=TRUE, fig.width=6, fig.height=6
genomicData <- cbind(dataList$Genetics, dataList$Cytogenetics)
genomicData <- (sapply(unique(sub("(_ITD)|(_TKD)|(_other)|(_mono)|(_bi)|(_p172)|(_p140)","",colnames(genomicData))), function(x) rowSums(genomicData[grep(paste(x,"($|_.+)",sep=""), colnames(genomicData))])) > 0)+0
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
colnames(t) <- c("NRAS_p12/13","NRAS_p61")
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
f <- function(t) t[1]/t[2] * t[4]/t[3]
oddspots <- sapply(hotspots, function(x) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(table(x, genomicData[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else if (f$p.val > 0.05) 1 else f$estimate} ))
#oddspots <- sapply(hotspots, function(x) sapply(1:ncol(genomicData), function(j) f(table(x, genomicData[,j])+.5)))#if (f$p.val > 0.05) 1 else f$estimate} ))
rownames(oddspots) <- rownames(potspots) <- colnames(genomicData)
oddspots[oddspots<1e-3] = 1e-4
oddspots[oddspots>1e3] = 1e4
logOdds=log10(oddspots)
#logOdds[] <- 0
#w <- which(potspots[o,]<0.05, arr.ind=TRUE)
#for(i in 1:nrow(w))
#	logOdds[rownames(w)[i], grep(paste(colnames(potspots)[w[i,2]], "($|_.+)",sep=""),colnames(logOdds))] <- log10(oddspots)[rownames(w)[i], grep(paste(colnames(potspots)[w[i,2]], "($|_.+)",sep=""),colnames(logOdds))]
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
for(q in list(potspots[o,]<0.05,p.adjust( potspots[o,], method="BH") < .1,p.adjust(potspots[o,]) < .05)){
	w = arrayInd(which(q), rep(nrow(potspots),2))
	x <- x00[w[,2]]
	rect(x0[w[,2]]+.5, w[,1]-.5,x0[w[,2]+1]+.5, w[,1]+.5, lwd=c(NA,1,2)[i])
	points(x,w[,1], pch=c("",".","*")[i], col="black")
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

#' Construct data.frame for EFS
dataFrameEfsTD <- dataFrame[tplSplitEfs,]
dataFrameEfsTD[which(tplIndexEfs), grep("TPL", colnames(dataFrameEfsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

patientTD <- c(clinicalData$PDID, clinicalData$PDID[which(tplIndexEfs)])

#' Construct data.frame for OS
#+ dataFrameOsTD, cache=TRUE
dataFrameOsTD <- dataFrame[tplSplitOs,]
dataFrameOsTD[which(tplIndexOs), grep("TPL", colnames(dataFrameOsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

#' Some indeces
#+ indeces, cache=TRUE
mainGroups <- grep("[A-Z][a-z]+[A-Z]",levels(groups), invert=TRUE, value=TRUE)
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
whichRFXOsGG <- which((colSums(dataFrame)>=8 | mainIdxOs) & osIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%


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
f <- formula("osYrTD ~ TP53 + TPL_os")
s <- survfit(f, data=dataFrameOsTD)
c <- coxph(f, data=dataFrameOsTD)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)


#' FLT ITD v TKD
#+ kmITDvTKD, fig.width=3, fig.height=2.5
f <- formula("osYr ~ FLT3_ITD + FLT3_TKD")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)


#' CEBPA mono v bi
#+ kmMonoVBi, fig.width=3, fig.height=2.5
f <- formula("osYr ~ CEBPA_mono + CEBPA_bi")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)

#' IDH2 p140 v p172
#+ km140v172, fig.width=3, fig.height=2.5
f <- formula("osYr ~ IDH2_p140 + IDH2_p172")
s <- survfit(f, data=dataList$Genetics)
c <- coxph(f, data=dataList$Genetics)
summary(c)
p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)

#' NPM1, FLT3_ITD and DNMT3A
coxph(osYr ~ NPM1*FLT3_ITD*DNMT3A, data=dataFrame)
s <- survfit(os ~ NPM1+FLT3_ITD+DNMT3A, data=dataFrame)
plot(s, col=set1)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)


#' #### NONC and complex karyotyope
table(complex=clinicalData$complex, `#Recurr.Cyto`=rowSums(dataList$Cytogenetics))
summary(coxph(osYr ~ ., data=dataFrame[mainIdxOs][colSums(dataFrame[mainIdxOs]) > 10]))

NONC <- rowSums(cbind(dataList$Cytogenetics[names(dataList$Cytogenetics)!="complex"], dataList$Genetics), na.rm=TRUE)
summary(coxph(osYr ~ . + NONC, data=dataFrame[mainIdxOs][colSums(dataFrame[mainIdxOs]) > 10]))

#+ NONC, fig.width=2.5, fig.height=2
r=colorRampPalette(set1[c(3,2,4,1,5)])(9)
s <- survfit(osYr ~ pmin(NONC,8))
plot(s, col=r, xlab="Years", ylab="Survival", mark=NA)
legend('topright', bty='n', col=r,legend= 0:8, lty=1, title="# onc. mut.", ncol=3)

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
colGroups <- c(brewer.pal(12, "Paired")[c(10)],brewer.pal(12, "Paired")[c(6,4,3,5,12,9,1,2,7)],"#999999", brewer.pal(12, "Paired")[c(8)])
names(colGroups) <- levels(groups)[order(toupper(levels(groups)))]
plot(coxRFXFitEfs, col=colGroups, order=order(coxRFXFitEfs$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

#' ### 2. Random effects: Time-dependent model (TPL)
#' Fit TD CoxRFX model
#+ coxRFXFitEfsTD, cache=TRUE
coxRFXFitEfsTD <- CoxRFX(dataFrameEfsTD[,whichRFXEfs], efsTD, groups[whichRFXEfs])
#+ coxRFXFitOsTD, cache=TRUE
coxRFXFitOsTD <- CoxRFX(dataFrameOsTD[,whichRFXOsTD], osTD, groups[whichRFXOsTD])

#+ coxRFXFitOsTDMain, cache=TRUE
coxRFXFitOsTDMain <- CoxRFX(dataFrameOsTD[,mainIdxOsTD], osTD, groups[mainIdxOsTD])

#+ coxRFXFitOsTDGG, cache=TRUE
whichRFXOsTDGG <- which((colSums(dataFrame)>=8 | mainIdxOsTD) & osTDIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%
coxRFXFitOsTDGG <- CoxRFX(dataFrameOsTD[,whichRFXOsTDGG], osTD, groups[whichRFXOsTDGG])
coxRFXFitOsTDGGc <- CoxRFX(dataFrameOsTD[,whichRFXOsTDGG], osTD, groups[whichRFXOsTDGG], which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 


par(mar=c(5,7,1,1))
plot(coxRFXFitEfsTD, col=colGroups, order=order(coxRFXFitEfsTD$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)

par(mar=c(5,7,1,1))
plot(coxRFXFitOs, col=colGroups, order=order(coxRFXFitOs$mu), xlim=c(-2,2), xaxt="n", xlab="Hazard ratio")
axis(side=1, at= log(c(0.1,0.2,0.5,1,2,5)), labels = c(0.1,0.2,0.5,1,2,5))
abline(v=0)


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
partRiskOsMain <- PartialRisk(coxRFXFitOsMain)
partRiskOs <- PartialRisk(coxRFXFitOs)
#varianceComponents <- rowSums(cov(partRiskTD, use="complete"))
varianceComponentsOs <- diag(cov(partRiskOsMain, use="complete"))
varianceComponentsOs
PlotVarianceComponents(coxRFXFitOs, col=colGroups)
title("Risk contributions OS")

PlotVarianceComponents(coxRFXFitOsMain, col=colGroups)
title("Risk contributions OS")

PlotVarianceComponents(coxRFXFitOsTDGGc, col=colGroups)
title("Risk contributions OS (time-dep)")

#' #### Distribution of survival effects
#' Effect sizes
#+ effectSizesKM, fig.width=1.5, fig.height=1.5
H0 <- basehaz(coxRFXFitOsTDGGc, centered = TRUE)
hazardDist <- splinefun(H0$time/365, H0$hazard, method="monoH.FC") ## in years
i <- 1
x <- seq(0,7, 0.01)
for(g in levels(groups)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Years")
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(coxRFXFitOs$coef[groups[whichRFXOs] == g]), each=length(x))  )), col=paste(colGroups[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(sqrt(coxRFXFitOs$sigma2[g]) * c(-1,1), each=length(x))  + coxRFXFitOs$mu[g])), col=paste(colGroups[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( coxRFXFitOs$mu[g] )), col=colGroups[i], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) ), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}

#' Partial risk
#+ partialRiskKM, fig.width=1.5, fig.height=1.5
i <- 1
partRiskOsTDGGc <- PartialRisk(coxRFXFitOsTDGGc)
for(g in colnames(partRiskOsTDGGc)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Years")
	m <- 0 #mean(rowSums(partRiskOsTDGGc[,colnames(partRiskOsTDGGc)!=g]))
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(partRiskOsTDGGc[,g]), each=length(x)) +m)), col=paste(colGroups[g],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskOsTDGGc[,g], c(0.25,0.75)), each=length(x))+m)), col=paste(colGroups[g],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskOsTDGGc[,g], c(0.05,0.95)), each=length(x))+m)), col=paste(colGroups[g],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( median(partRiskOsTDGGc[,g])+m)), col=colGroups[g], type="l", lwd=2)	
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
risk <- as.matrix(dataFrame[whichRFXOsTDGG]) %*% coxRFXFitOsTDGGc$coefficients
risk <- risk - mean(risk)
parBoot <- mclapply(1:100, function(i) {
			s <- SimSurvNonp(risk, os)
			c <- try(CoxRFX(dataFrame[whichRFXOsTDGG], s, groups=groups[whichRFXOsTDGG], sigma0=0.1, nu=0))
			if(class(c)=="try-error")
				return(s)
			c$Z <- NULL # set X to zero to save mem
			return(c)
		}, mc.cores=10)

#' Plots of mean, sigma and df
#+ parBootPlots, fig.width=3, fig.height=2,  eval=TRUE
boxplot(t(sapply(parBoot, `[[`, "sigma2")), border=colGroups[names(parBoot[[1]]$sigma2)], lty=1, pch=16, staplewex=0, ylab="sigma2", las=2, log="y", ylim=c(1e-3,1))
abline(h=0, lty=3)
points(coxRFXFitOsTDGGc$sigma2,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "mu")), border=colGroups[names(parBoot[[1]]$mu)], lty=1, pch=16, staplewex=0, ylab="mu", las=2)
abline(h=0, lty=3)
points(coxRFXFitOsTDGGc$mu,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "df")), border=colGroups[names(parBoot[[1]]$mu)], lty=1, pch=16, staplewex=0, ylab="df", las=2)
abline(h=0, lty=3)
points(coxRFXFitOsTDGGc$df,  pch=19)

#' Coefficients
#+ parBootSignif, fig.width=2.5, fig.height=2.5
v <- apply(sapply(parBoot, `[[`, "coefficients"), 1, var, na.rm=TRUE)
w <- diag(coxRFXFitOsTDGGc$var) ## H^{-1}
w2 <- diag(coxRFXFitOsTDGGc$var2) ## H^{-1} I H^{-1}
c <- coef(coxRFXFitOsTDGGc)
plot(c^2/v, c^2/w, log="xy", xlab="Chi2 (bootstrap)", ylab="Chi2 (analyt.)", cex=.66)
par(xpd=NA)
points(c^2/v, c^2/w2, pch=16, cex=.7)
arrows(c^2/v, c^2/w, c^2/v,c^2/w2, length=0.05)
abline(0,1)
abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))

#' Table with significance
#+ parBootTable, results='asis'
library(xtable)
print(xtable(data.frame(group = groups[whichRFXOsTDGG],coef=round(c,4), sd = round(sqrt(w2),4), boot=sig2star(pchisq(c^2/v,1, lower.tail=FALSE)), var2=sig2star(pchisq(c^2/w2,1, lower.tail=FALSE)),var=sig2star(pchisq(c^2/w,1, lower.tail=FALSE)))),  type="html")

#' Volcano plot
#+ volcanoGGc, fig.width=3.3, fig.height=3
par(mar=c(3,3,1,3)+.1,  bty="n", mgp=c(2,.5,0))
i <- coxRFXFitOsTDGGc$groups %in% c("Genetics", "CNA","BT","GeneGene","Treatment")#apply(coxRFXFitOsTDGGc$Z,2,min) == 0 & apply(coxRFXFitOsTDGGc$Z,2,max) == 1
plot(c, c^2/w2, log='', col=paste(colGroups[as.character(coxRFXFitOsTDGGc$groups)],"BB", sep=""), pch=ifelse(i,16,16), ylab="Chi2",xlab="log hazard", cex=ifelse(i, sqrt(colMeans(coxRFXFitOsTDGGc$Z[!rev(duplicated(rev(tplSplitOs))),])*50),1), xlim=range(c*1.2))
#abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))
p <- pchisq(c^2/w2,1,lower.tail=FALSE) ## pvalues coxRFX
w <- which(p.adjust(p,"BY") < 0.1)
points(c[w], (c^2/w2)[w],  pch=1, cex=ifelse(i[w], sqrt(colMeans(coxRFXFitOsTDGGc$Z[!rev(duplicated(rev(tplSplitOs))),w])*50),1))
w <- which(p.adjust(p,"bonf") < 0.05)
par(xpd=NA)
text(c[w], (c^2/w2)[w], names(c[w]), pos=3)
pp <- pretty(pchisq(par("usr")[3:4],1,lower.tail=FALSE, log=TRUE)/log(10))
q <- qchisq(pp*log(10),1,lower.tail=FALSE, log.p=TRUE)
axis(side=4, at = q[q<par("usr")[4]], labels=pp[q<par("usr")[4]])
mtext(side=4, "log10 P-value", line=1.5)
u <- par("usr")
f <- c(0.01,0.05,0.1,0.2,0.5)
s <- sqrt(f*50)
legend("topright",legend=f, pch=16, pt.cex=s, bty='n', col="grey")

#' P-values and random model
#+ pValuesMain, fig.width=2.5, fig.height=2.5, cache=TRUE
set.seed(42)
Z <- apply(coxRFXFitOsTDGGc$Z, 2,sample)[1:nrow(dataFrame),] ## random covariates
coxRFXFitOsRain <- CoxRFX(Z, os, groups=coxRFXFitOsTDGGc$groups, nu=1) ## model
w2 <- diag(coxRFXFitOsRain$var2) 
c <- coef(coxRFXFitOsRain)
p2 <- pchisq(c^2/w2,1,lower.tail=FALSE)
plot(seq(0,1,l=length(p2)+1)[-1],sort(p2), xlab="P-value (expected)", ylab="P-value (observed)", pch=16, col="grey")
abline(0,1)
points(seq(0,1,l=length(p)+1)[-1],sort(p), pch=16)
legend("topleft",bty="n", c("observed","randomised"), pch=16, col=c("black","grey"))



#+ parBootVarianceComp, fig.width=3, fig.height=2, cache=TRUE, eval=TRUE
v <- t(sapply(parBoot, function(x) {t <- try(VarianceComponents(x, newZ=dataFrame[whichRFXOsTDGG])); if(class(t)=="try-error") rep(NA, nlevels(x$groups)+1) else t}))
boxplot(v, border=colGroups[colnames(v)], lty=1, pch=16, staplewex=0, ylab="variance comp.", las=2)
abline(h=0, lty=3)
points(VarianceComponents(coxRFXFitOsTDGGc),  pch=19)

round(cov(v), 2)


rm(parBoot)

#' Volcano plot
#+ volcanoAll, fig.width=2.5, fig.height=2.5
c <- coef(coxRFXFitOsTD)
v <- diag(coxRFXFitOsTD$var2)
plot(c, -pchisq(c^2/v, 1,lower.tail = FALSE, log=TRUE), log='', col=colGroups[as.character(coxRFXFitOsTD$groups)], pch=19, ylab="Chi2",xlab="log hazard")
abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))
p <- pchisq(c^2/v,1,lower.tail=FALSE) ## pvalues coxRFX
w <- which(p.adjust(p,"BH") < 0.05)
par(xpd=NA)
text(c[w], (c^2/v)[w], names(c[w]), pos=3)

#' #### Large overview panel
#+ overviewRFX, fig.width=6, fig.height=5
for(model in c("coxRFXFitOs", "coxRFXFitOsMain","coxRFXFitOsTDGGc")){
	layout(matrix(1:6, nrow=2), width=c(4,5,2),height=c(6,12))
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
	par(xpd=NA)
	segments(seq_along(c), get(model)$mu[as.character(f)[o]], seq_along(c), c, col=colGroups[as.character(f)[o]])
	rotatedLabel(table(f)/2+ cumsum(c(0,table(f)[-nlevels(f)])), y = rep(par("usr")[4]*.8, nlevels(f)), pos=3, labels=levels(f))
	par(bty="L", mar=c(4,4,1,2))
	r <- if(grepl("TD", model)) which(rev(!duplicated(rev(tplSplitOs)))) else 1:nrow(dataFrame) ## select final observation in TD case
	X <- get(model)$Z[r,w][,o]
	p <- order(X %*% c)
	X <- X[p,]
	x <- t(X)/apply(X,2,max)

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
	
	H0 <- basehaz(coxph(osYr ~ 1), centered=TRUE)
	hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
	invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
	l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
	i <- 1
	x <- seq(from=-3,to=3, l=512)
	par(mar=c(1,.5,5,3), xpd=FALSE)
	plot(x, invHazardDist(-log(l[2]) /exp(x) ), col='black', lty=c(2,1,2)[2], xlab="", ylab="Time", type='l',ylim=c(0,10), xlim=c(min(h), 4), xaxt="n", yaxt="n")
	abline(v=seq(-2,3,1), col="grey")
	polygon(c(x, rev(x)), c(invHazardDist(-log(l[1]) /exp(x) ), rev(invHazardDist(-log(l[3]) /exp(x))) ), border=NA, col="#88888888")
	rotatedLabel(x0 = log(-log(l)/hazardDist(par("usr")[4])), y0 = rep(par("usr")[4], length(labels)), pos=3, labels=paste(100*l, "%", sep=""))
	par(xpd=NA)
	mtext(side=3, "Survival",line=2, cex=.66)
	mtext(side=2, line=2.5, "Time",cex=.66)
	s <- get(model)$surv[r]
	s[,(ncol(s):1)[-1]] <- s[,(ncol(s):1)[-1]]/365
	q <- quantile(s[,ncol(s)-1], seq(0,1,0.1), na.rm=TRUE)
	image(x=c(3.5,4), y = q, matrix(1:10, nrow=1), col=brewer.pal(10,'RdBu'), add=TRUE)
	axis(side=4, at=axTicks(side=4))
	mtext(side=4, "Years", line=2.5, cex=.66)
	
	
	par(mar=c(4,.5,1,3), xpd=FALSE)
	plot(h, seq_along(h), pch=NA, xlab="Total log hazard", yaxs="i", yaxt="n", ylab='', xlim=c(min(h), 4), xaxt='n')
	axis(side=1, at=-2:3)
	abline(v=seq(-2,3,1), col="grey")
	segments(h, seq_along(h),0, seq_along(h))
	par(xpd=NA)	
	x <- seq(-3,5,l=100)
	lines(x , dnorm(x, 0, sd(h))*100 /  dnorm(0, 0, sd(h)) + length(h)*1.01)
	
	c <- cut( s[p,ncol(s)-1], q)
	#c[s[p,"status"]==0] <- NA 
	image(x=c(3.5,4), y = c(0,seq_along(c)), matrix(as.numeric(c), nrow=1), col=brewer.pal(10,'RdBu'), add=TRUE, useRaster=TRUE)
	points(x=rep(4.1, sum(s[p,"status"]==0, na.rm=TRUE)), y=which(s[p,"status"]==0), pch=".")
}

#' #### Stars
#+ stars, fig.width=12, fig.height=12
set.seed(42)
library(HilbertVis)
nStars <- 32
for(l in c("coxRFXFitOs","coxRFXFitOsMain","coxRFXFitOsTDGGc")){
	t <- os#get(l)$surv
	p <- PartialRisk(get(l),  newZ=dataFrame[, colnames(get(l)$Z)])
	p <- p[,colnames(p)!="Nuisance"]
	locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
	s <- sample(nrow(p),nStars^2) #1:(nStars^2)
	h <- hclust(dist(p[s,]))
	x <- p - rep(colMeans(p), each=nrow(p))
	x <- x/(2*sd(x)) + 1
		c <- cut(t[s,1][h$order], quantile(t[,1], seq(0,1,0.1), na.rm=TRUE))
	if(l=="coxRFXFitOsTDGGc")
		x <- x[,c("Clinical","Demographics","Genetics","GeneGene","CNA","BT","Treatment")]
	stars(x[s,][h$order,]/2, scale=FALSE, locations=locations, key.loc=c(0,-3), col.lines=rep(1,(nStars^2)), col.stars = (brewer.pal(11,'RdBu'))[c])
	symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=1)
	title(main=l)
}

#' #### Selected stars
#+ patientStars, fig.width=4, fig.height=1.5
patients <- c(
		which(dataFrame$`TP53`==1 & dataFrame$complex==1 & os[,1] < 300 & os[,2]==1)[1],
		which(dataFrame$`NPM1:FLT3_ITD:DNMT3A`==1 & os[,1] < 300 & os[,2]==1)[1],
		which(dataFrame$SFRS2==1 & clinicalData$WHOcat=='no' & os[,2]==1)[1],
		which(dataFrame$NPM1==1 & dataFrame$FLT3_ITD==0 & dataFrame$DNMT3A==0 & os[,1] > 2000)[1],
		which(dataFrame$t_15_17==1 & os[,1] > 2000)[1]
)
genotype <- apply(dataFrame[groups %in% c("BT","CNA","Genetics")]==1, 1,function(x) paste(names(which(x)), collapse=";"))

t <- os
p <- PartialRisk(coxRFXFitOsTDGGc, newZ=dataFrame[, whichRFXOsTDGG])
p <- p[,colnames(p)!="Nuisance"]
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
h <- hclust(dist(p[s,]))
x <- p - rep(colMeans(p), each=nrow(p))
x <- x/(2*sd(x)) + 1
c <- cut(t[patients,1], quantile(t[,1], seq(0,1,0.1), na.rm=TRUE))
x <- x[patients,c("Clinical","Demographics","Genetics","GeneGene","CNA","BT","Treatment")]
locations <- expand.grid(seq_along(patients)* 1.5, 1)
stars(x/2, scale=FALSE, locations=locations, key.loc=NA, col.lines=rep(1,(nStars^2)), col.stars = (brewer.pal(11,'RdBu'))[c])
symbols(locations[,1], locations[,2], circles=rep(.5,length(patients)), inches=FALSE, fg="grey", add=TRUE, lty=1)
text(locations[,1], locations[,2]-1,labels=clinicalData$PDID[patients], pos=1)
l <- apply(dataFrame[patients,c("gender","AOD_10","TPL_os","wbc_100")], 1,paste, collapse=";")
par(xpd=NA)
text(locations[,1], locations[,2]+1,labels=paste(gsub(";","\n",genotype[patients]),l, paste(round(os[patients,1],2), osYr[patients,2]), sep="\n"), pos=3)


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

#' OS
survConcordance(osTD~coxRFXFitOsTDGGc$linear.predictors)

#' Partial values of Harrel's C
c <- PartialC(coxRFXFitOsTDGGc)
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
coxRFXFitOsTrain <- CoxRFX(dataFrame[trainIdx,whichRFXOs], os[trainIdx], groups=groups[whichRFXOs])

#' ### 2. Time-dependent random effects model
#+ coxRFXFitOfsTDTrain, cache=TRUE
coxRFXFitOsTDTrain <- CoxRFX(dataFrameOsTD[trainIdxOsTD,whichRFXOsTD], osTD[trainIdxOsTD], groups=groups[whichRFXOsTD])

#' Partial contributions
partialRiskOsTest <- PartialRisk(coxRFXFitOsTrain, newZ=dataFrame[!trainIdx, whichRFXOs])
c <- PartialC(coxRFXFitOsTrain, newZ = dataFrame[!trainIdx, whichRFXOs], newSurv = os[!trainIdx])
b <- barplot(c[1,], col=colGroups[colnames(c)], las=2)
segments(b, c[1,]-c[2,], b, c[1,]+c[2,])
abline(h=.5)

#' Overall
totalRiskOsTest <- rowSums(partialRiskOsTest)
survConcordance(os[!trainIdx]~totalRiskOsTest )

#' Compared to molecular risk
predictiveRiskTest <- rowSums(partialRiskOsTest[,-which(colnames(partRiskOsMain) %in% c("Treatment","GeneTreat","CytoTreat","Nuisance"))])
survConcordance(os[!trainIdx] ~ predictiveRiskTest)

#barplot(c(CGP=survConcordance(survivalTD[!testIxTD] ~ predictiveRiskTest)$concordance, MolecularRisk = survConcordance(survivalTD[!testIxTD] ~ c(Favourable=1, Adverse=4, `inter-1`=2, `inter-2`=3)[clinicalData$M_Risk[tplSplit][!testIxTD]])[[1]]))

library(survivalROC)
library(survAUC)
s <- os[!trainIdx]
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



#' ### 5. Stability selection with interactions
#' Only include interaction terms when main effects are present. 
#' On the other hand, the exclusive selection of a product term could mean that the most likely explanation is that the two main terms are zero and only the interaction is non-zero..
CoxCPSSInteractions


#' #### CPSS on OS
#+ CoxCPSSIntOs, cache=TRUE
set.seed(42)
scope <- c("Genetics","CNA","Treatment","BT")
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
#' #### Static models
#+allModelsCV, cache=TRUE, cache.lazy=FALSE
replicates <- 100 ## number of replicates
allModelsCV <- mclapply(1:replicates, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			c <- coxph(os[trainIdx] ~ 1, data=dataFrame[trainIdx,mainIdxOs])
			scopeStep <- as.formula(paste("os[trainIdx] ~", paste(colnames(dataFrame)[mainIdxOs], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxCPSSOsTrain <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx, mainIdxOs], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,mainIdxOs], os[trainIdx], groups=groups[mainIdxOs])
			coxRFXOsTrain$Z <- NULL
			coxRFXOsGGc <- CoxRFX(dataFrame[trainIdx,whichRFXOsGG], os[trainIdx], groups=groups[whichRFXOsGG], which.mu=mainGroups)
			coxRFXOsGGc$Z <- NULL
			rForestOsTrain <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs])[trainIdx,], ntree=100, importance="none")
			return(list(
							BIC=coxBICOsTrain,
							AIC=coxAICOsTrain,
							CPSS=coxCPSSOsTrain,
							RFX=coxRFXOsTrain,
							RFXgg=coxRFXOsGGc,
							rForest=rForestOsTrain
					))
		}, mc.cores=10)

#' #### Predictions

#' All predictions
#+ allModelsCvPredictions, cache=TRUE
predictAllModels <- function(x, newdata){
	if("rfsrc" %in% class(x)){
		predict(x, newdata, importance="none")$predicted
	}else{
		predict(x, newdata)
	}
}

allModelsCvPredictions <- mclapply(seq_along(allModelsCV), function(foo){
			set.seed(foo)
			x <- allModelsCV[[foo]]
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[!trainIdx]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrame[!trainIdx,])
							}))
		}, mc.cores=10)

#' Harrel's C
#+ allModelsCv-C
foo <- 1
allModelsCvC <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x, 2 , function(p){						
						survConcordance(osYr[!trainIdx,] ~ p)$concordance
					})
		})
apply(allModelsCvC,1,quantile)
boxplot(t(allModelsCvC), ylim=c(0,0.5), notch=TRUE, ylab="Concordance", border=colModels, las=2, lty=1, pch=16, staplewex=0)

#+ allModelsCvBoxplot, fig.width=2, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(t(apply(-allModelsCvC,2,rank))),factor, levels=1:7)),table)
o <- order(apply(allModelsCvC,1,median))
boxplot(t(allModelsCvC[o,]), notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:7, rep(par("usr")[3],7), rownames(allModelsCvC)[o])

#+ allModelsCvRank, fig.width=2, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
barplot(r[,o]/replicates, col=c(set1[c(3,2,4,1,5,7)],"grey"), ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],6), colnames(allModelsCvC)[o])
legend(par("usr")[1],1.5, fill=c(set1[c(3,2,4,1,5,7)],"grey"), legend=1:6, bty="n", border=NA, horiz=TRUE, title="Rank")


#' Brier scores
#+ allModelsCv-Brier
foo <- 1
allModelsCvBrier<- sapply(allModelsCV, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			sapply(x, function(y){
						p <- predictAllModels(y, newdata=dataFrame)
						a <- predErr(Surv.rsp = osYr[trainIdx,], Surv.rsp.new = osYr[!trainIdx,], lp=p[trainIdx], lpnew = p[!trainIdx], times= c(90,365,1000)/365, type="brier")$error
					})
		})
apply(allModelsCvBrier,1,quantile)
rownames(allModelsCvBrier) <- paste(rep(names(allModelsCV[[1]]), each=3), rep(c(90,365,1000), length(allModelsCV[[1]])))
colModels <- c("#888888", set1[c(2,1,4,3,5,7)])
boxplot(t(allModelsCvBrier)[,rep(0:5*3, 3) + rep(1:3, each=6)],notch=TRUE, ylab="Brier score", border=rep(colModels[-1],3), las=2, lty=1, pch=16, staplewex=0)


#' GHCI
#+ allModelsCv-GHCI
allModelsCvGHCI<- sapply(allModelsCvPredictions, function(x){
			apply(x[,2:6], 2 , function(p){
						p <- GHCI(lpnew = na.omit(p))
					})
		})
apply(allModelsCvGHCI,1,quantile)
boxplot(t(allModelsCvGHCI),notch=TRUE, ylab="GH C", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)


#' OXS R2 estimates
#+ allModelsCv-OXS
foo <- 1
allModelsCvOXS <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x[,2:6], 2 , function(p){						
						a <- OXS(osYr[!trainIdx,], p, rep(0,length(p)))
					})
		})
apply(allModelsCvOXS,1,quantile)
boxplot(t(allModelsCvOXS), ylim=c(0,0.5), notch=TRUE, ylab="OXS R2", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)

#' Nagelk R2 estimates
#+ allModelsCv-Nagelk
foo <- 1
allModelsCvNagelk <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x[,2:6], 2 , function(p){						
						a <- Nagelk(osYr[!trainIdx,], p, rep(0,length(p)))
					})
		})
apply(allModelsCvNagelk,1,quantile)
boxplot(t(allModelsCvNagelk), ylim=c(0,0.4), notch=TRUE, ylab="Nagelk's R2", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)

#' UnoC
#+ allModelsCv-UnoC
foo <- 1
allModelsCvUnoC<- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x,2, function(p){
						a <- UnoC(Surv.rsp = osYr[trainIdx,], Surv.rsp.new = osYr[!trainIdx,][!is.na(p)],  lpnew = na.omit(p), time=5)
					})
		})
apply(allModelsCvUnoC,1,quantile)
boxplot(t(allModelsCvUnoC), notch=TRUE,  ylab="Uno's C", border=colModels, lty=1, pch=16, staplewex=0)

#' AUC UNO
#+ allModelsCv-AUCuno
t <- seq(0.1,5,0.1) #times
allModelsCvAuc <- sapply(seq_along(allModelsCvPredictions), function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			apply(allModelsCvPredictions[[foo]], 2, function(p){
						AUC.uno(osYr[trainIdx,], osYr[!trainIdx, ][!is.na(p)], scale(na.omit(p)), t)$auc
					})
		})
allModelsCvAuc <- array(allModelsCvAuc, dim=c(length(t),ncol(allModelsCvPredictions[[1]]),length(allModelsCvPredictions)))
plot(NA,NA, xlab="Years",ylab="AUC", xlim=range(t), ylim=c(0.5,0.8))
for(i in 1:dim(allModelsCvAuc)[2]){
	lines(t,rowMeans(allModelsCvAuc, dims=2)[,i], type='l', new=i==1, col=colModels[i])
}
legend("bottomright", colnames(allModelsCvPredictions[[1]]), bty="n", lty=1, col=colModels)


#' Wisdom of the Krauts?
foo <- 1
allModelsCvCKraut <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			r <- rowMeans(apply(x, 2 , rank))
			survConcordance(osYr[!trainIdx,] ~ r)$concordance
		})
quantile(allModelsCvCKraut)
boxplot(cbind(t(allModelsCvC),allModelsCvCKraut), notch=TRUE, ylab="Concordance", border=c(colModels,1), las=2, lty=1, pch=16, staplewex=0)

ranks <- apply(apply(-cbind(t(allModelsCvC),kraut=allModelsCvCKraut),1,rank, ties.method="random"),1,function(x) table(factor(x, levels=1:8)))
ranks <- ranks[,order(1:8 %*% ranks)]

#' Clean up.. 
rm(allModelsCV)

#' #### Time-dep models

#+ allModelsCVTD, cache=TRUE
replicates <- 100 ## number of replicates
allModelsCVTD <- mclapply(1:replicates, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			c <- coxph(osTD[trainIdx] ~ 1, data=dataFrameOsTD[trainIdx,mainIdxOsTD])
			scopeStep <- as.formula(paste("osTD[trainIdx] ~", paste(colnames(dataFrameOsTD)[mainIdxOsTD], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxRFXOsTrain <- CoxRFX(dataFrameOsTD[trainIdx,mainIdxOsTD], osTD[trainIdx], groups=groups[mainIdxOsTD])
			coxRFXOsTrain$Z <- NULL
			coxRFXOsGGc <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTDGG], osTD[trainIdx], groups=groups[whichRFXOsTDGG], which.mu=mainGroups)
			coxRFXOsGGc$Z <- NULL
			
			return(list(
							BIC=coxBICOsTrain,
							AIC=coxAICOsTrain,
							RFX=coxRFXOsTrain,
							RFXgg=coxRFXOsGGc							))
		}, mc.cores=10)


#' All predictions
#+ allModelsCvTdPredictions, cache=TRUE
allModelsCvTdPredictions <- mclapply(seq_along(allModelsCVTD), function(foo){
			set.seed(foo)
			x <- allModelsCVTD[[foo]]
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )[tplSplitOs]!=1
			cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[tplSplitOs][!trainIdx]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrameOsTD[!trainIdx,])
							}))
		}, mc.cores=10)

#' Harrel's C
#+ allModelsCvTd-C
foo <- 1
allModelsCvTdC <- sapply(allModelsCvTdPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x, 2 , function(p){						
						survConcordance(osYrTD[!trainIdx,] ~ p)$concordance
					})
		})
apply(allModelsCvTdC,1,quantile)
boxplot(t(allModelsCvTdC), ylim=c(0,0.5), notch=TRUE, ylab="Concordance", border=colModels, las=2, lty=1, pch=16, staplewex=0)


#+ allModelsCvTdCBoxplot, fig.width=1.5, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(t(apply(-allModelsCvTdC,2,rank))),factor, levels=1:6)),table)
o <- order(apply(allModelsCvTdC,1,median))
boxplot(t(allModelsCvTdC[o,]), notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:nrow(allModelsCvTdC), rep(par("usr")[3],ncol(allModelsCvTdC)), rownames(allModelsCvTdC)[o])

#+ allModelsCvTdCRank, fig.width=1.5, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
barplot(r[,o]/replicates, col=set1[c(3,2,4,1,5,7)][1:ncol(allModelsCvTdC)], ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],ncol(allModelsCvTdC)), colnames(allModelsCvTdC)[o])
legend(par("usr")[1],1.5, fill=set1[c(3,2,4,1,5,7)][1:nrow(allModelsCvTdC)], legend=1:nrow(allModelsCvTdC), bty="n", border=NA, horiz=TRUE, title="Rank")

#' #### Different RFX models

#' RFX models with different interaction terms
#+ allModelsCvRfx, cache=TRUE
replicates <- 100 ## number of replicates
allModelsCvRfxC <- do.call("rbind",mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrameOsTD)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			coxRFXOsMain <- CoxRFX(dataFrameOsTD[trainIdx,mainIdxOsTD], osTD[trainIdx], groups=groups[mainIdxOsTD])
			coxRFXOsGG <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTDGG], osTD[trainIdx], groups=groups[whichRFXOsTDGG])
			coxRFXOsGGc <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTDGG], osTD[trainIdx], groups=groups[whichRFXOsTDGG], which.mu=mainGroups)
			coxRFXOsAll <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTD], osTD[trainIdx], groups=groups[whichRFXOsTD])
			coxRFXOsAllc <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTD], osTD[trainIdx], groups=groups[whichRFXOsTD], which.mu=mainGroups)
			return(c(
							Main=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,mainIdxOsTD]) %*% coef(coxRFXOsMain))$concordance,
							GeneGene=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTDGG]) %*% coef(coxRFXOsGG))$concordance,
							GeneGeneCentred=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTDGG]) %*% coef(coxRFXOsGGc))$concordance,
							AllInt=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTD]) %*% coef(coxRFXOsAll))$concordance,	
							AllIntCentred=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTD]) %*% coef(coxRFXOsAllc))$concordance	
			))
		}, mc.cores=10))
colnames(allModelsCvRfxC) <- sub(".concordant","",colnames(allModelsCvRfxC))

#+ allModelsCvRfxBoxplot, fig.width=2, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(round(t(apply(-allModelsCvRfxC,1,rank)))),factor, levels=1:6)),table)
o <- order(colMeans(allModelsCvRfxC))
boxplot(allModelsCvRfxC[,o], notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:ncol(allModelsCvRfxC), rep(par("usr")[3],ncol(allModelsCvRfxC)), colnames(allModelsCvRfxC)[o])

#+ allModelsCvRfxRank, fig.width=2, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
barplot(r[,o]/replicates, col=set1[c(3,2,4,1,5,7)][1:ncol(allModelsCvRfxC)], ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],ncol(allModelsCvRfxC)), colnames(allModelsCvRfxC)[o])
legend(par("usr")[1],1.5, fill=set1[c(3,2,4,1,5,7)][1:ncol(allModelsCvRfxC)], legend=1:ncol(allModelsCvRfxC), bty="n", border=NA, horiz=TRUE, title="Rank")



#' #### Inter-study CV
#+ allModelsTrial, cache=TRUE
allModelsTrial <- mclapply(levels(clinicalData$Study), function(foo){
			#set.seed(foo)
			trainIdx <- clinicalData$Study != foo 
			c <- coxph(os[trainIdx] ~ 1, data=dataFrame[trainIdx,mainIdxOs])
			scopeStep <- as.formula(paste("os[trainIdx] ~", paste(colnames(dataFrame)[mainIdxOs], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxCPSSOsTrain <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx, mainIdxOs], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
			w <- colnames(dataFrame[mainIdxOs])
			w <- setdiff(w, names(which(colSums(dataFrame[trainIdx,w])==0)))
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,w], os[trainIdx], groups=groups[w], nu = if(foo=="AMLSG0704") 1 else 0) # add prior for 0704 (just one group member)
			coxRFXOsTrain$Z <- NULL
			w <- whichRFXOsGG
			w <- setdiff(w, which(colSums(dataFrame[trainIdx,w])==0))
			coxRFXOsGGc <- CoxRFX(dataFrame[trainIdx,w], os[trainIdx], groups=groups[w], which.mu=mainGroups, nu = if(foo=="AMLSG0704") 1 else 0)
			coxRFXOsGGc$Z <- NULL
			rForestOsTrain <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs])[trainIdx,], ntree=100, importance="none")
			return(list(
							BIC=coxBICOsTrain,
							AIC=coxAICOsTrain,
							CPSS=coxCPSSOsTrain,
							RFX=coxRFXOsTrain,
							RFXgg=coxRFXOsGGc,
							rForest=rForestOsTrain
					))
		}, mc.cores=3)
names(allModelsTrial) <- levels(clinicalData$Study)

allModelsTrialPredictions <- mclapply(names(allModelsTrial), function(foo){
			x <- allModelsTrial[[foo]]
			trainIdx <- clinicalData$Study != foo
			cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[!trainIdx]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrame[!trainIdx,])
							}))
		}, mc.cores=10)
names(allModelsTrialPredictions) <- names(allModelsTrial)

allModelsTrialC <- sapply(names(allModelsTrial), function(foo){
			trainIdx <- clinicalData$Study != foo
			apply(allModelsTrialPredictions[[foo]], 2 , function(p){						
						survConcordance(osYr[!trainIdx,] ~ p)$concordance
					})
		})

allModelsTrialC


		
#' ### 6. Models on genomics only
#+ genomicModels, cache=TRUE
genomicsIndex <- groups %in% c("Genetics","CNA", "BT")
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
tcgaSurvival <- Surv(tcgaClinical$OS/365, tcgaClinical$Status)

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
tcgaRiskRFXOs <- PredictRiskMissing(coxRFXFitOsTDGGc, tcgaData[whichRFXOsTDGG])
survConcordance(tcgaSurvival ~ tcgaRiskRFXOs[,1])

#' #### CPSS model
tcgaDataImputed <- as.data.frame(ImputeMissing(dataFrame[mainIdxOs], newX=tcgaData[mainIdxOs]))
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
tcgaAUC <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], c(90,365,1000)/365)$auc)
tcgaAUCi <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], sort(na.omit(tcgaSurvival[,1])))$iauc)
o <- order(colMeans(tcgaAUC))
barplot(tcgaAUC[,o], border=1, col= rep(c("grey",set1[-6]),each=3), las=2, xaxt="n", ylab="AUC", beside=TRUE, density=c(NA, 48,24), ylim=c(0.5,0.85), xpd=FALSE) -> b
legend("topleft", bty="n", c("3mo","1yr","3yr"), fill='black', density=c(NA, 48,24))
rotatedLabel(b[seq(3, length(b), 3)], rep(0.49,length(tcgaRisk)), names(tcgaRisk)[o], srt=45)


#+ kmTCGA, fig.width=3, fig.height=2.5
risk <- cut(tcgaRiskCPSSOs, quantile(tcgaRiskCPSSOs), labels=c("1st Q","2nd Q","3rd Q","4th Q"))
s <- survfit(tcgaSurvival ~ risk)
plot(s, col=set1[c(3,2,4,1)], mark=NA, xlab="Years", ylab="Survival")
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
	lines(exp(x), pmax(0,invHazardDist(-log(l[i]) /exp(x) ))/10000*365, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000*365)
mtext(side=4, "Time", line=2.5)
mtext(side=3, at = -log(l)/hazardDist(par("usr")[4]*10000*365), text=paste(100*l, "% survive", sep=""))
legend("topright", levels(tcgaClinical$C_Risk)[c(2,3,1)], fill=set1[c(3,2,1)], bty="n", title="M risk")

#' #### CV revisited
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
o <- order(colMeans(allModelsCvTdC))
boxplot(allModelsCvTdC[,o], notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n", border="black")
segments(1:6+.05,tcgaConcordance[1,c(1,3,7,6,8,5)]-tcgaConcordance[2,c(1,3,7,6,8,5)],1:6+.05,tcgaConcordance[1,c(1,3,7,6,8,5)], col='red')
segments(1:6+.05,tcgaConcordance[1,c(1,3,7,6,8,5)],1:6+.05,tcgaConcordance[1,c(1,3,7,6,8,5)]+tcgaConcordance[2,c(1,3,7,6,8,5)], col='red')
points(1:6+.05,tcgaConcordance[1,c(1,3,7,6,8,5)], col='red', pch=16, cex=2)
rotatedLabel(1:6, rep(par("usr")[3],6), colnames(allModelsCvTdC)[o])
legend("bottomright", c("CV x100", "TCGA +/- sd"), lty=c(1,1), bty="n", col=c(1,2), pch=c(22,16))


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
coxRFXHotspots <- CoxRFX(cbind(genes, hotspots), os, g,  which.mu=1)
plot(coxRFXHotspots)

#' Basic CV
#+ coxRFXHotspotsCV
t <-  CoxRFX(cbind(genes, hotspots)[trainIdx,], os[trainIdx], g, which.mu=1)
survConcordance( os[!trainIdx] ~ (as.matrix(cbind(genes, hotspots)) %*% coef(t))[!trainIdx])
t0 <-  CoxRFX(genes[trainIdx,], os[trainIdx],)
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

#' ### Interpolations
#' #### 1. Subsampling of patients
#+ subsetPatients, cache=TRUE
library(survivalROC)
set.seed(42)
subsets <- seq(100,1500,100)
subsetPatients <- lapply(subsets, function(s){
			mclapply(1:ceiling(50000/(1540-s)), function(foo){
						set.seed(s*foo)
						trn <- 1:nrow(dataFrame) %in% sample(nrow(dataFrame), s)
						tst <-  !trn 
						fit <- CoxRFX(dataFrameOsTD[tplSplitOs[trn], whichRFXOsTDGG], osTD[tplSplitOs[trn]], groups[whichRFXOsTDGG], which.mu=mainGroups, nu = 0.1)
						C <- survConcordance(osTD[tplSplitOs[tst]]~predict(fit, newdata=dataFrameOsTD[tplSplitOs[tst], whichRFXOsTDGG]))
						ROC <- survivalROC(Stime=os[!is.na(os) & tst,1], status=os[!is.na(os) & tst,2], marker = predict(fit, newdata=dataFrame[tst, whichRFXOsTDGG]), predict.time = 850, method="KM", cut.values=seq(-5,5,0.1))
						list(C, ROC, trn, tst, coef(fit))}, mc.cores=20)
		})

#+ subsetPatientsPlot, fig.width=2, fig.height=2
#pdf("subsetConcordance.pdf", 2.5,2.5, pointsize=8)
col1 <- colorRampPalette(set1[c(3,2,4,1,5)])(length(subsets))
plot(NA,NA, xlim=c(0,1),ylim=c(0,1), xlab="FPR",ylab="TPR")
abline(0,1, lty=3)
for(i in seq_along(subsets)){
	x <- sapply(subsetPatients[[i]], function(x) x[[2]]$FP)
	y <- sapply(subsetPatients[[i]], function(x) x[[2]]$TP)
	lines(rowMeans(x),rowMeans(y), col=col1[i], type="l")
}
#legend("bottomright", legend=rev(subsets), lty=1, col=col1[5:1], bty="n")

rangeplot <- function(x, y, col = 1, pch = 19, lty = 1, ylim=range(y),...){
	plot(x, colMeans(y), col = col, pch=pch, ylim = ylim, ..., xaxt="n")
	points(jitter(rep(x, each=nrow(y))),y,pch=1, col=rep(col, each=nrow(y)), cex=.2) 
	lines(x, colMeans(y), lwd=2)
	lines(x, colMeans(y) + 2*apply(y, 2, sd)/sqrt(nrow(y)))
	lines(x, colMeans(y) - 2*apply(y, 2, sd)/sqrt(nrow(y)))
	axis(at = x, labels=x, side=1)
	#segments(x,apply(y,2,min),x,apply(y,2,max), col=col, lty = lty)
}

rangeplot2 <- function(x, y, col = 1, pch = 19, lty = 1, ylim=range(unlist(y)),...){
	plot(x, sapply(y, mean), col = col, pch=pch, ylim = ylim, ..., xaxt="n")
	points(jitter(unlist(sapply(seq_along(y), function(i) rep(x[i], length(y[[i]]))))),unlist(y),pch=1, col=unlist(sapply(seq_along(y), function(i) rep(col[i], length(y[[i]])))), cex=.2) 
	lines(x, sapply(y, mean), lwd=2)
	lines(x, sapply(y, mean) + 2*sapply(y,  sd)/sqrt(sapply(y,length)))
	lines(x, sapply(y, mean) - 2*sapply(y,  sd)/sqrt(sapply(y,length)))
	axis(at = x, labels=x, side=1)
	#segments(x,apply(y,2,min),x,apply(y,2,max), col=col, lty = lty)
}

rangeplot3 <- function(x, y, col = 1, pch = 19, lty = 1, ylim=range(unlist(y)),...){
	plot(x, sapply(y, mean), col = col, pch=pch, ylim = ylim, ...)
	#points(jitter(unlist(sapply(seq_along(y), function(i) rep(x[i], length(y[[i]]))))),unlist(y),pch=1, col=unlist(sapply(seq_along(y), function(i) rep(col[i], length(y[[i]])))), cex=.2) 
	#lines(x, sapply(y, mean), lwd=2)
	s <- sapply(y,  sd)/sqrt(sapply(y,length)) 
	m <- sapply(y, mean)
	segments(x, m+s*2, x, m-s*2, col=col)
	#axis(at = x, labels=x, side=1)
	#segments(x,apply(y,2,min),x,apply(y,2,max), col=col, lty = lty)
}

rangeplot2(x=subsets, y = sapply(subsetPatients, function(x) sapply(x, function(y) y[[2]]$AUC)) , col=col1, xlab="Cohort", ylab="AUC", lty=1, ylim=c(0.7,0.85))
rangeplot2(x=subsets, y = sapply(subsetPatients, function(x) sapply(x, function(y) y[[1]]$concordance)) , col=col1, xlab="Cohort", ylab="Concordance", lty=1, ylim=c(0.65,.75), log='')

rangeplot3(x=subsets, y = sapply(subsetPatients, function(x) sapply(x, function(y) y[[1]]$concordance)) , col=col1, xlab="Cohort", ylab="Concordance", lty=1, ylim=c(0.67,.73), log='')


#lines(x=subsets, y = concordanceFromVariance(sapply(subsetPatients, function(x) {
#					mean(sapply(x, function(y) {
#										h <-  var(as.matrix(dataFrameOsTD[tplSplitOs[y[[3]]],whichRFXOsTDGG]) %*% y[[5]])
#									}))
#				})) , col=1, xlab="Cohort", ylab="Concordance", ylim=c(0.65,.75))
#

#' #### 2. Subsampling of genes
#+ subsetGenes, cache=TRUE
set.seed(42)
subsets <- seq(5,55,5)
genes <- names(whichRFXOsTDGG[groups=="Genetics"])
subsetGenes <- lapply(subsets, function(s){
			mclapply(1:100, function(foo){
						g <- sample(genes, s)
						ix <- !grepl(paste(g,collapse="|"), names(whichRFXOsTDGG))
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
						testIdx <-  !trainIdx 
						fit <- CoxRFX(dataFrameOsTD[tplSplitOs[trainIdx], whichRFXOsTDGG[ix]], osTD[tplSplitOs[trainIdx]], groups[whichRFXOsTDGG[ix]], which.mu=mainGroups, nu = 0.1)
						C <- survConcordance(osTD[tplSplitOs[testIdx]]~predict(fit, newdata=dataFrameOsTD[tplSplitOs[testIdx], whichRFXOsTDGG[ix]]))
						ROC <- survivalROC(Stime=os[!is.na(os) & testIdx,1], status=os[!is.na(os) & testIdx,2], marker = predict(fit, newdata=dataFrame[testIdx, whichRFXOsTDGG[ix]]), predict.time = 850, method="KM", cut.values=seq(-5,5,0.1))
						fit <- CoxRFX(dataFrameOsTD[, whichRFXOsTDGG[ix]], osTD, groups[whichRFXOsTDGG[ix]], which.mu=mainGroups, nu = 0.1)
						S <- cov(PartialRisk(fit))
						list(C, ROC, S, trainIdx, testIdx, ix, mean(rowMeans(dataFrame[setdiff(genes,g)])))
					}, mc.cores=20)
		})

#+ subsetGenesPlot, fig.width=2, fig.height=2
plot(sapply(subsetGenes, function(x) sapply(x, function(y) y[[7]]*sum(y[[6]][1:58]))), sapply(subsetGenes, function(x) sapply(x, function(y) y[[1]]$concordance)), xlab="Mean no. of drivers", ylab="Concordance")
plot(sapply(subsetGenes, function(x) sapply(x, function(y) y[[7]]*sum(y[[6]][1:58]))), sapply(subsetGenes, function(x) sapply(x, function(y) {t <- try(sum(y[[3]][c("Genetics"),c("Genetics")])); ifelse(class(t)=="try-error",NA,t)})), xlab="Mean no. of drivers", ylab="Variance")
plot(sapply(subsetGenes, function(x) sapply(x, function(y) y[[7]]*sum(y[[6]][1:58]))), sapply(subsetGenes, function(x) sapply(x, function(y) {t <- try(sum(y[[3]][c("GeneGene"),c("GeneGene")])); ifelse(class(t)=="try-error",NA,t)})), xlab="Mean no. of drivers", ylab="Variance")
plot(sapply(subsetGenes, function(x) sapply(x, function(y) y[[7]]*sum(y[[6]][1:58]))), sapply(subsetGenes, function(x) sapply(x, function(y) {t <- try(sum(y[[3]][c("Genetics","GeneGene"),c("Genetics","GeneGene")])); ifelse(class(t)=="try-error",NA,t)})), xlab="Mean no. of drivers", ylab="Variance")

#+ subsetGenesPlotTCGA, fig.width=2.5, fig.height=2.5
plot(sapply(subsetGenes, function(x) sapply(x, function(y) y[[7]]*sum(y[[6]][1:58]))), sapply(subsetGenes, function(x) sapply(x, function(y) {t <- try(sum(y[[3]][c("Genetics","GeneGene"),c("Genetics","GeneGene")])); ifelse(class(t)=="try-error",NA,t)})), xlab="Mean no. of drivers", ylab=expression(paste(Var,"[",h[g],"]")), xlim=c(0,3.8), ylim=c(0,.35), pch=16, col=c("#00000044"))
x <- c(0,3.7)
s <- coxRFXFitOsTDGGc$sigma2["Genetics"]
segments(c(2.3, 3.7), rep(par("usr")[3],2), c(2.3, 3.7), c(2.3, 3.7) * s, col="grey")
segments( rep(par("usr")[1],2),  c(2.3, 3.7) * s, c(2.3, 3.7), c(2.3, 3.7) * s, col="grey")
lines(x, x*s, col="red")
par(xpd=NA)
axis(at=c(2.3, 3.7), labels=c("111 genes", "TCGA (exome)"), tcl=0.5, side=1, mgp=c(-2.5,-2,0))

#' 8. Extrapolations
#' -----------------
#' ### 1. Generate new data
#' Simulate data using multiple imputation.
#+ simData, cache=TRUE
set.seed(42)
SimDataNonp
simData <- SimDataNonp(dataFrame[mainIdxOsTD], nData = 10000, m=10)
names(simData) <- names(dataFrame[mainIdxOsTD])

#' Merge into data.frame
#+ simDataFrame, cache=TRUE
set.seed(42)
g <- groups[mainIdxOsTD]
for(w in which(colSums(simData,na.rm=TRUE) == 0))
	simData[[w]] <- rbinom(nrow(simData),1,mean(dataFrame[mainIdxOsTD][,w]))
all(colSums(simData,na.rm=TRUE) != 0)
simDataFrame <- cbind(simData,
		MakeInteractions(simData[,g=="Genetics"], simData[,g=="Genetics"])[,as.vector(upper.tri(matrix(0,ncol=sum(g=="Genetics"), nrow=sum(g=="Genetics"))))])
for(n in unique(which(is.na(simDataFrame), arr.ind = TRUE)[,2]))
	simDataFrame[[n]] <- poorMansImpute(simDataFrame[[n]])
simDataFrame <- StandardizeMagnitude(simDataFrame)
simDataFrame <- simDataFrame[,colnames(simDataFrame)  %in% names(whichRFXOsTDGG) | colSums(simDataFrame)>=8]
simDataFrame$`NPM1:FLT3_ITD:DNMT3A` <- simDataFrame$NPM1 * simDataFrame$FLT3_ITD * simDataFrame$DNMT3A
dim(simDataFrame)

#' #### Basie simulations
simGroups <- factor(c(as.character(g), rep("GeneGene", ncol(simDataFrame)-length(g))))
names(simGroups) <- colnames(simDataFrame)
simCoef <- CoxHD:::SimCoef(coxRFXFitOsTDGGc, groups = simGroups)

simRisk <- as.matrix(simDataFrame[names(whichRFXOsTDGG)]) %*% simCoef[names(whichRFXOsTDGG)]
simSurv <- SimSurvNonp(simRisk, os)

survConcordance(simSurv ~ simRisk)

#' #### Save output
save(coxRFXFitOsTDGGc, whichRFXOsTDGG, simDataFrame, simGroups, os, mainGroups, file="sim2Data.RData")

#' ### 2. Simulation code
#' The following code is run on the farm
#+ farmulations, cache=FALSE
read_chunk('Farmulations2.R')

#' ### 3. Analysis
#' Read files
files <- dir("simRFX", pattern="Farmulations\\[1-1000\\]*", full.names = TRUE)
tmp <- new.env()
load(files[1], envir = tmp)

#' #### P-values
#' Plot the P-values as a function of Np\beta^2.
#+ pVarSchoenfeld, fig.width=2, fig.height=2, cache=TRUE
w <- groups[whichRFXOsTDGG] %in% c("Genetics","BT","CNA", "GeneGene") ## Which groups
psi <- mean(os[,2]) ## Fraction of uncensored observations
plot(colSums(simDataFrame[names(whichRFXOsTDGG[w])]) * tmp$simCoef[whichRFXOsTDGG[w]]^2 , CoxHD:::WaldTest( tmp$fit10000)$p[w], log="yx", pch=NA, xlab=expression(psi *N *p *beta^2), ylab="P-value", ylim=c(1e-50,1))
for(f in files[1:50]){
	load(f, envir = tmp)
	points(psi*colSums(simDataFrame[names(whichRFXOsTDGG[w])]) * tmp$simCoef[names(whichRFXOsTDGG[w])]^2 , CoxHD:::WaldTest( tmp$fit10000)$p[w],  col=colGroups[as.character(groups)[whichRFXOsTDGG[w]]], pch=1, cex=.5)
	points(psi*colSums(simDataFrame[tmp$w1000, names(whichRFXOsTDGG[w])]) * tmp$simCoef[names(whichRFXOsTDGG[w])]^2 , CoxHD:::WaldTest( tmp$fit1000)$p[w], col=colGroups[as.character(groups)[whichRFXOsTDGG[w]]], pch=2, cex=.5)
	if(tmp$fit100$iter[1] < 50) ## Exclude simulations without convergence
		points(psi*colSums(simDataFrame[tmp$w100, names(whichRFXOsTDGG[w])]) * tmp$simCoef[names(whichRFXOsTDGG[w])]^2 ,CoxHD:::WaldTest( tmp$fit100)$p[w],  col=colGroups[as.character(groups)[whichRFXOsTDGG[w]]], pch=3, cex=.5)
}
legend("bottomleft", lty=c(0,1),pch=c(1,NA), c("Simulations","Schoenfeld"), bty="n")
x <- 10^seq(-4,4,0.1)
lines(x, pnorm(sqrt(x), lower.tail = FALSE))


#' #### Power
#' The theoretical power (according to Schoenfeld/Schmoor is)
power <- function(beta, N, p, psi=0.5, alpha=0.05){
	pnorm(sqrt(N*psi*beta^2*p*(1-p))-qnorm(1-alpha/2))
}

#' Plot for observed cases and overlay a few usual suspects
#+ power1540, fig.width=3, fig.height=3
x <- seq(-2,2,0.01)
y <- 10^seq(-4,0,0.01)
colLevels <- colorRampPalette(brewer.pal(9, "Reds")[-(1:2)])(11)
g <- c("BT","CNA","Genetics","GeneGene")
xObs <- matrix(exp(rep(coxRFXFitOsTDGGc$mu[g], each=2) + c(-1,1) * rep(sqrt(coxRFXFitOsTDGGc$sigma2[g]),each=2)), nrow=2) ## Mean log haz +/- sd
yObsQ <- sapply(split(colMeans(dataFrameOsTD[whichRFXOsTDGG]), groups[whichRFXOsTDGG]),quantile, c(0.05,0.5,0.95))[,g] ## 5,50,95% frequency quantiles

contour(outer(x,y,function(x,y) power(x,1540,y)), x=exp(x),y=y, log='xy', xlab="Hazard ratio", ylab="Mutation frequency", main="N=1540", col=colLevels)
rect(xObs[1,],yObsQ[1,],xObs[2,],yObsQ[3,], border = colGroups[c("BT","CNA","Genetics","GeneGene")])
#segments(exp(coxRFXFitOsTDGGc$mu[g]),yObsQ[1,],exp(coxRFXFitOsTDGGc$mu[g]),yObsQ[3,], col = colGroups[g])
#segments(xObs[1,],yObsQ[2,],xObs[2,],yObsQ[2,], col = colGroups[g])

effects <- c("NPM1","TP53","inv3_t3_3","t_15_17","inv16_t16_16","CEBPA_bi","FLT3_ITD","complex","NPM1:FLT3_ITD:DNMT3A") ## A few interesting variables
points(exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), col=colGroups[as.character(groups[effects])], pch=19)
text(labels=effects,exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), pos=ifelse(sign(coef(coxRFXFitOsTDGGc)[effects])==1,4,2))
legend("bottom", lty=c(1,NA,NA,NA,NA,NA),pch=c(NA,19,22,22,22,22), c("Power","Selected variables", paste("Dist.", g)), col=c(colLevels[10], "black", colGroups[g]), bty="n", ncol=2)

#' Compared to other cohort sizes
#+ power100-10000, fig.width=1.5, fig.height=1.5
for(N in c(100,1000,10000)){
	contour(outer(x,y,function(x,y) power(x,N,y)), x=exp(x),y=y, log='xy', xlab="Hazard ratio", ylab="Mutation frequency", main=paste("N=",N,sep=""), col=colLevels, drawlabels=FALSE)
	rect(xObs[1,],yObsQ[1,],xObs[2,],yObsQ[3,], border = colGroups[g])
}

#' #### Concordance
#+ concordance100-1000, fig.width=2, fig.height=2, cache=TRUE
C <- sapply(files[1:500], function(f){
			load(f)
			c(`100`=survConcordance(SimSurvNonp(simRisk[w100], os)~fit100$linear.predictors)$concordance,
					`1000`=survConcordance(SimSurvNonp(simRisk[w1000], os)~fit1000$linear.predictors)$concordance,
					`10000`=survConcordance(SimSurvNonp(simRisk, os)~fit10000$linear.predictors)$concordance,
					Truth=survConcordance(SimSurvNonp(simRisk, os)~simRisk)$concordance)
		})
boxplot(t(C), staplewex=0, pch=16, lty=1, ylab="", ylab="Concordance", xaxt="n")
rotatedLabel(labels=(sub(".concordant","", rownames(C))))
abline(h=CoxHD:::concordanceFromVariance(var(tmp$fit10000$linear.predictors)))

#' #### Cohort size
#+ cohort, fig.width=2.5, fig.height=2.5
par(mar=c(3,3,1,1), bty='n', mgp=c(2,0.5,0))
cohort <- function(beta, p, psi=0.5, alpha=0.05, power=0.5){
	(qnorm(1-alpha/2) + qnorm(1-power) )^2 / (beta^2 * psi * p * (1-p))
}
x <- seq(-2,2, 0.01)
y <- 10^seq(-3,0, 0.01)
contour(outer(x,y,function(x,y) cohort(x,y, alpha=0.05/100)), x=exp(x),y=y, log='xy', xlab="Hazard ratio", ylab="Mutation frequency",  col=colLevels, levels=c(10,20,50,100,200,500,1000,2000,5000,10000,20000))
rect(xObs[1,],yObsQ[1,],xObs[2,],yObsQ[3,], border = colGroups[c("BT","CNA","Genetics","GeneGene")])
effects <- c("NPM1","TP53","inv3_t3_3","t_15_17","inv16_t16_16","CEBPA_bi","FLT3_ITD","complex","NPM1:FLT3_ITD:DNMT3A") ## A few interesting variables
points(exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), col=colGroups[as.character(groups[effects])], pch=19)
text(labels=effects,exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), pos=ifelse(sign(coef(coxRFXFitOsTDGGc)[effects])==1,4,2))
#legend("bottom", lty=c(1,NA,NA,NA,NA,NA),pch=c(NA,19,22,22,22,22), c("Power","Selected variables", paste("Dist.", g)), col=c(colLevels[10], "black", colGroups[g]), bty="n", ncol=2)

#' #### Number of cases
#+ cases, fig.width=2.2, fig.height=2
par(mar=c(3,5,1,1), bty='n', mgp=c(2.5,0.5,0))
cases <- function(beta, alpha=0.05, power=0.5, p = 1e-2, psi=0.5){
	(qnorm(1-alpha/2) + qnorm(1-power) )^2 / (beta^2 * (1-p) * psi) 
}
x <- seq(-1,1,0.01)

x0 <- log(c(0.01,0.02,0.05,0.1,0.2,0.5,1)+1)
plot(exp(x), cases(x, alpha=5e-2), log='yx', type='l', xlab="Hazard ratio", ylab="Minimal number of cases", las=1)
#lines(exp(x), cases(x, alpha=1e-2),  type='l', lty=2)
lines(exp(x), cases(x, alpha=1e-3),  type='l', lty=3)
segments(exp(x0), par("usr")[3],exp(x0),cases(x0, alpha=5e-2), col='grey')
segments(exp(x[1]), cases(x0, alpha=5e-2),exp(x0),cases(x0, alpha=5e-2), col='grey')
axis(side=2, at=cases(x0, alpha=5e-2), labels=exp(x0), tcl=.5, line=0, las=2, mgp=c(-2.5,-.5,0), hadj=0)
axis(side=1, at=c(seq(0.1,3,0.1)), labels=rep("",30), tcl=-.2, line=0, las=2)
axis(side=2, at=rep(c(1:10), 4) * 10^rep(1:4, each=10), labels=rep("",40), tcl=-.2, line=0, las=2)
legend("topright", legend=c("P < 0.05 *","P < 0.001 ***"), lty=c(1,3), bty="n")


#' 9. Regression of blood counts
#' -----------------------------
#' ### Prepare data
#+ clinicalGlmnet, cache=TRUE
library(glmnet)
Y <- StandardizeMagnitude(cbind(dataList$Clinical, dataList$Demographics))
X <- as.matrix(dataFrame[groups %in% c("Genetics","CNA","BT")])
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
annot <- 3 - 2 * grepl("^[A-Z]",n) 
annot[grep("(^t_)|(_t)",n)] <- 2
annot <- factor(annot, labels = c("Genetics","BT","CNA"))
names(annot) <- n
for(m in clinModels){
	plotcvnet(m, X, main=names(clinModels)[i],  col0="black", cex=1, simple.annot = annot, col=colGroups[levels(annot)], xlim=c(0.5,35.5))
	i = i+1
	legend("topright", col=c(pastel1[c(1,3)],"black")[c(1,3,2)], c(expression(paste("Explained variance ",R^2)), expression(paste("Lasso penalty ",lambda)), expression(paste("Model coefficient ", beta))), box.lty=0, bg="#FFFFFF33", pch=c(NA,NA,19), lty=c(1,1,NA), cex=.8, pt.cex = 1)
}

#' ### Heatmap of GLMs
j <- 0
z <- sapply(clinModels,function(x){ ## Creating matrix
			j <<- j+1
			w <- which.min(x$cvm)
			c <- x$glmnet.fit$beta[,w]
			yj <- sapply(c("Genetics","CNA","BT"), function(i){
						w <- names(groups[groups == i])
						X[,w] %*% c[w] 
					})
			cj <- rowSums(cov(yj))
			y <- Y[,j] - x$glmnet.fit$a0[w]
			covj <- colMeans((y-mean(y))*(yj - rep(colMeans(yj), each=nrow(yj))))
			r2 <- cj
			R2 <- 1 - x$cvm[w]/x$cvm[1]
			c(c, NA,  r2, R2=R2)
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
R2 <- z[nrow(z) - 3:0,]
R2[is.na(R2)] <- 0
z <- z[-(nrow(z) - 0:4),]
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
o <- order(annot,rowSums(r, na.rm=TRUE))
image(y=1:ncol(z)-.5, x=1:nrow(z), z[o,w], breaks=c(-2,seq(-0.1,0.1,l=51)), col=c("grey",colorRampPalette(brewer.pal(9,"RdYlBu"))(50)), xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0,ncol(z)+1))
#abline(v=c(18.5,27.5), lwd=0.5)
rotatedLabel(y0=rep(0.5,nrow(z)), labels=gsub("_","/",rownames(z))[o], x0=1:nrow(z), font=c(rep(3,sum(groups=="Genetics")),rep(1,sum(groups%in%c("CNA","BT")))), col = colGroups[as.character(annot[o])], cex=0.9)
mtext(side=2, line=.2, text=colnames(z)[w], las=2, at=1:ncol(z)-.5)
text(y=rep(1:ncol(z)-.5, each=nrow(r)), x=rep(1:nrow(r), ncol(z)), r[o,w] * (0!=(z[o,w])), cex=0.66, font=ifelse(r[o,w] <= rep(s[ncol(r):1], each=nrow(r)), 2,1))
points(y=rep(1:ncol(z)-.5, each=nrow(r)), x=rep(1:nrow(r), ncol(z)), pch=ifelse(is.na(z[o,w]) | z[o,w]==0, ".",NA))
mtext(side=1, at=sum(groups=="Genetics")/2, "Genetics", col=colGroups["Genetics"], line=2.5 )
mtext(side=1, at=sum(groups=="CNA")/2 + sum(groups=="Genetics"), "CNA", col=colGroups["CNA"], line=2.5 )
mtext(side=1, at=sum(groups=="Tranlocations")/2 + sum(groups%in%c("Genetics","BT")), "BT", col=colGroups["BT"], line=2.5 )
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
barplot(R2[1:3,w] * rep(R2[4,w] / colSums(R2[1:3,w]), each=3), border=NA, col=colGroups[levels(annot)], horiz=TRUE, names.arg=rep(NA,ncol(R2)), width=0.95, space=0.0525, add=TRUE) -> b
points(R2[4,w]+0.1,b, pch=ifelse(p[w],"*",NA))
mtext(side=3, "Variance components", line=.5)
mtext(side=1, expression(paste("Explained variance ",R^2)), line=2.5)

#+ save
#save(list = ls(), file=paste(Sys.Date(), "-AML.RData", sep=""))


#' 10. Data for web tool
#' --------------------
#' ### Prelim
#' Times for allografts pre and post relapse, after 1CR only
alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD") # only allografts
alloTimeCR1 <- clinicalData$Time_1CR_TPL + .5 # +.5 to make > 0
alloTimeCR1[!alloIdx | (clinicalData$TPL_date < clinicalData$Recurrence_date & !clinicalData$TPL_Phase %in% c("CR1","RD"))] <- NA

#' Create data frames for each phase
whichRFXCirTD <- whichRFXOsTDGG[grep("TPL",names(whichRFXOsTDGG), invert=TRUE)] #mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Nuisance"
t <- clinicalData$Recurrence_date
t[is.na(t)] <- as.Date(1e6, origin="2000-01-01")
cirData <- MakeTimeDependent(dataFrame[whichRFXCirTD], timeEvent=alloTimeCR1, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=!is.na(clinicalData$Recurrence_date)+0)
cirData$transplantCR1 <- cirData$event
cirData$event <- NULL
cirData$transplantRel <- 0
nrmData <- MakeTimeDependent(dataFrame[whichRFXCirTD], timeEvent=alloTimeCR1, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=is.na(clinicalData$Recurrence_date) & clinicalData$Status)
nrmData$transplantCR1 <- nrmData$event
nrmData$event <- NULL
nrmData$transplantRel <- 0
i <- !is.na(clinicalData$Recurrence_date)
#prsData <- makeTimeDependent(dataFrame[w], timeTpl=alloTimeRel, timeSurv=as.numeric(clinicalData$Date_LF - clinicalData$Recurrence_date)+1, event=clinicalData$Status)
prsData <- MakeTimeDependent(dataFrame[i,whichRFXCirTD], timeEvent=alloTimeCR1[i], timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date)[i], timeStart = as.numeric(clinicalData$Recurrence_date- clinicalData$CR_date)[i], status=clinicalData$Status[i])
prsData$transplantCR1 <- nrmData$transplantCR1[1:nrow(dataFrame)][prsData$index]
prsData$transplantRel <- prsData$event
prsData$event <- NULL

#' Fit models
crGroups <- c(as.character(groups[whichRFXCirTD]), "Treatment","Treatment")
names(crGroups) <- c(names(dataFrame)[whichRFXCirTD],"transplantCR1","transplantRel")
coxRFXNrmTD <- CoxRFX(nrmData[names(crGroups)], Surv(nrmData$time1, nrmData$time2, nrmData$status), groups=crGroups, which.mu = mainGroups)
coxRFXNrmTD$coefficients["transplantRel"] <- 0
#prsData$time1[!is.na(prsData$time1)] <- 0
prsData$time1[prsData$time1 >= prsData$time2] <- NA
coxRFXPrsTD <-  CoxRFX(prsData[names(crGroups)], Surv(prsData$time1, prsData$time2, prsData$status), groups=crGroups, nu=1, which.mu = mainGroups)
coxRFXCirTD <-  CoxRFX(cirData[names(crGroups)], Surv(cirData$time1, cirData$time2, cirData$status), groups=crGroups, which.mu = mainGroups)
coxRFXCirTD$coefficients["transplantRel"] <- 0

osData <- MakeTimeDependent(dataFrame[whichRFXCirTD], timeEvent=alloTimeCR1, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
osData$transplantCR1 <- osData$event
osData$transplantRel <- osData$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date | clinicalData$TPL_Phase != "CR1")  
osData$transplantCR1[osData$index %in% w] <- 0
#osData$transplantRel[!osData$index %in% w] <- 0
osData$transplantRel <- 0

coxRFXOsCR <- CoxRFX(osData[names(crGroups)], Surv(osData$time1, osData$time2, osData$status), groups=crGroups, which.mu = mainGroups)

#save(coxRFXCirTD, coxRFXNrmTD, coxRFXPrsTD, coxRFXOsCR, nrmData, cirData, prsData, osData, osCR, crGroups, file="../../../Projects/sandbox/relapse/predict.RData")

#' Prediction of OS and Cross-validation
#+concordanceCIRcv, cache=TRUE
PredictOS <- function(coxRFXNrmTD, coxRFXCirTD, coxRFXPrsTD, data, x =365){
	getS <- function(coxRFX, data) {		
		coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
		data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data))])
		r <- PredictRiskMissing(coxRFX, data, var="var2")
		H0 <- basehaz(coxRFX, centered = FALSE)
		hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
		x <- 0:5000
		S <- exp(-hazardDist(x))
		return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
	}
	hazCir <- getS(coxRFX = coxRFXCirTD, data = data)
	hazNrm <- getS(coxRFX = coxRFXNrmTD, data = data)
	hazPrs <- getS(coxRFX = coxRFXPrsTD, data = data)
	
	hazCir$Sadj <- cumsum(c(1,diff(hazCir$S)) * hazNrm$S)
	hazNrm$Sadj <- cumsum(c(1,diff(hazNrm$S)) * hazCir$S)
	
	w <- hazCir$x == x
	nrs <- hazNrm$Sadj[w]^exp(hazNrm$r[,1])
	nrsUp <- hazNrm$Sadj[w]^exp(hazNrm$r[,1] + 2*sqrt(hazNrm$r[,2]))
	nrsLo <- hazNrm$Sadj[w]^exp(hazNrm$r[,1] - 2*sqrt(hazNrm$r[,2]))
	
	cir <- hazCir$Sadj[w]^exp(hazCir$r[,1])
	cirUp <- hazCir$Sadj[w]^exp(hazCir$r[,1] + 2*sqrt(hazCir$r[,2]))
	cirLo <- hazCir$Sadj[w]^exp(hazCir$r[,1] - 2*sqrt(hazCir$r[,2]))
	
	rs <- 1 - (1-cir) * (1-hazPrs$S[w]^exp(hazPrs$r[,1]))
	rsUp <- 1 - (1-cirUp) * (1-hazPrs$S[w]^exp(hazPrs$r[,1] + 2* sqrt(hazPrs$r[,2])))
	rsLo <- 1 - (1-cirLo) * (1-hazPrs$S[w]^exp(hazPrs$r[,1] - 2* sqrt(hazPrs$r[,2])))
	
	return(data.frame(os=rs*nrs, osUp = rsUp*nrsUp, osLo = rsLo*nrsLo, rs=rs, cir=cir, nrs=nrs))
}

d <- MakeTimeDependent(dataFrame[whichRFXCirTD], timeEvent=alloTimeCR1, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
d$transplantCR1 <- d$event
d$transplantRel <- d$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date)  
d$transplantCR1[d$index %in% w] <- 0
d$transplantRel[!d$index %in% w] <- 0

replicates <- 100 ## number of replicates
concordanceCIRcv <- lapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			mclapply(1:replicates, function(foo){
						set.seed(foo)
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 4/5
						dNrm <- nrmData[nrmData$index %in% which(trainIdx),names(g)]
						sNrm <- Surv(nrmData$time1, nrmData$time2, nrmData$status)[nrmData$index %in% which(trainIdx)]
						coxRFXNrmTD <- CoxRFX(dNrm, sNrm, groups=g, nu=1, which.mu = mainGroups)
						coxRFXNrmTD$coefficients["transplantRel"] <- 0
						dPrs <- prsData[prsData$index %in% which(trainIdx), names(g)]
						sPrs <- Surv(prsData$time1, prsData$time2, prsData$status)[prsData$index %in% which(trainIdx)]
						coxRFXPrsTD <-  CoxRFX(dPrs, sPrs, groups=g, nu=1, which.mu = mainGroups)
						dCir <- cirData[cirData$index %in% which(trainIdx), names(g)]
						sCir <- Surv(cirData$time1, cirData$time2, cirData$status)[cirData$index %in% which(trainIdx)]
						coxRFXCirTD <-  CoxRFX(dCir, sCir, groups=g, which.mu = mainGroups)
						coxRFXCirTD$coefficients["transplantRel"] <- 0
						dOs <- osData[osData$index %in% which(trainIdx), names(g)]
						sOs <- Surv(osData$time1, osData$time2, osData$status)[osData$index %in% which(trainIdx)]
						coxRFXOsCR <- CoxRFX(dOs, sOs, groups=g, which.mu = mainGroups)
						
						allRisk365 <- PredictOS(coxRFXNrmTD = coxRFXNrmTD, coxRFXPrsTD = coxRFXPrsTD, coxRFXCirTD = coxRFXCirTD, d, 365)
						allRisk1000 <- PredictOS(coxRFXNrmTD = coxRFXNrmTD, coxRFXPrsTD = coxRFXPrsTD, coxRFXCirTD = coxRFXCirTD, d, 1000)
						
						p365 <- -allRisk365[,1]
						p1000 <-  -allRisk1000[,1]
						pCIR <- as.matrix(cirData[names(g)]) %*% coef(coxRFXCirTD)
						pPRS <- as.matrix(prsData[names(g)]) %*% coef(coxRFXPrsTD)
						pNRM <- as.matrix(nrmData[names(g)]) %*% coef(coxRFXNrmTD)
						pOS <- as.matrix(osData[names(g)]) %*% coef(coxRFXOsCR)
						
						C <- c(
								CIRrfx = survConcordance(Surv(time1, time2, status)~ pCIR, data=cirData, subset = cirData$index %in% which(!trainIdx) )$concordance,
								PRSrfx = survConcordance(Surv(time1, time2, status) ~ pPRS, data=prsData, subset=prsData$index %in% which(!trainIdx) )$concordance,
								NRMrfx = survConcordance(Surv(time1, time2, status)~  pNRM, data=nrmData, subset=nrmData$index %in% which(!trainIdx) )$concordance,
								OSrfx = survConcordance(Surv(time1, time2, status) ~ pOS, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
								OS365 = survConcordance(Surv(time1, time2, status) ~ p365, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
								OS1000 = survConcordance(Surv(time1,time2, status) ~ p1000, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance
						)
						
						coef <- cbind(CIRrfx=coef(coxRFXCirTD), PRSrfx=coef(coxRFXPrsTD), NRMrfx=coef(coxRFXNrmTD),  OSrfx=coef(coxRFXOsCR))
						
						return(list(C=C, coef=coef, allRisk365=allRisk365, allRisk1000=allRisk1000))
					}, mc.cores=10)
		})

riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)

#' Sources of mortality
#+ mortality, fig.width=3, fig.height=3
i <- 1
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	plot(NA,NA, ylim=c(0,1),  xlab="Years", ylab="Mortality", xlim=c(0,10), yaxs='i', xaxs='i')
	abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	#abline(v=seq(1,9), col='lightgrey')
	lines(survfit(Surv(time1/365, time2/365, status) ~ clinicalData$M_Risk[osData$index], data=osData, subset=clinicalData$M_Risk[osData$index]==l), col=riskCol[l],fun=function(x) 1-x ,mark=NA, lty=1, conf.int=FALSE)
	rsKM <- survfit(Surv(time1/365, time2/365, status) ~ 1, data=osData, subset=  clinicalData$M_Risk[osData$index]==l)
	nrsKM <- survfit(Surv(time1/365, time2/365, status) ~ 1, data=nrmData, subset=  clinicalData$M_Risk[nrmData$index]==l)
	
	rsCR <- cumsum(c(1,diff(rsKM$surv)) * splinefun(nrsKM$time, nrsKM$surv, method="monoH.FC")(rsKM$time))
	nrsCR <- cumsum(c(1,diff(nrsKM$surv)) * splinefun(rsKM$time, rsKM$surv, method="monoH.FC")(nrsKM$time))
	
	
	lines(rsKM$time, 1-rsCR, col=riskCol[l], lty=2, type='s')
	lines(nrsKM$time, 1-nrsCR, col=riskCol[l], lty=3, type='s')
	if(i ==1)
	legend(ifelse(i<=3,"topleft","bottomright"), c("total","relapse","non-rel"), lty=c(1,2,3), col="black", box.lty = 0, bg="white")
	i <- i+1
	mtext(l, side=3, font=2)
}

f <- function(x) 1-x
plot(survfit(Surv(time1/365, time2/365, status) ~ clinicalData$M_Risk[cirData$index], data=cirData), col=riskCol, ylab="CIR", xlab="Time after CR", main="Molecular risk groups, all cases", fun=f , ylim=c(0,1))
legend("bottomright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)

#' CIR Plots
#+ cirSplits, fig.width=3, fig.height=3
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
r <- coxRFXCirTD$Z %*% coef(coxRFXCirTD) - cirData$transplantCR1 * coef(coxRFXCirTD)["transplantCR1"]
Q <- numeric(nrow(cirData))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	w <- which(clinicalData$M_Risk[cirData$index]==l)
	q <- cut(r[w], quantile(r[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
	Q[w] <- q
	plot(NA,NA,  ylab="CIR", main=paste(l, "terciles"),  xlab="Years after CR", ylim=c(0,1), xlim=c(0,10), xaxs="i", yaxs="i")
	abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	fit <- survfit(Surv(time1/365, time2/365, status) ~ q + transplant1CR, data=cirData[w,])
	## adjust for competing risk (NRM)
	i <- c(0,diff(fit$surv))
	s <- split(fit$surv, cumsum(i>0)) # split into strata
	t <- split(fit$time, cumsum(i>0))
	fit$surv <- unlist(sapply(seq_along(s), function(i) cumsum(c(1,diff(s[[i]])) * splinefun(nrsKM$time, nrsKM$surv, method="monoH.FC")(t[[i]])))) #adjust increments by nrs KM est
	lines(fit, col=rep(sapply(2:0,function(x) colTrans(riskCol[l],x)), each=2), lty=c(1,3), mark=NA, xlab="Time after CR", fun=f)
	legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[l])
}

#' OS Plots
#+ osSplits, fig.width=3, fig.height=3
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
r <- coxRFXOsCR$Z %*% coef(coxRFXOsCR) - osData$transplantCR1 * coef(coxRFXOsCR)["transplantCR1"]
Q <- numeric(nrow(osData))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	w <- which(clinicalData$M_Risk[osData$index]==l)
	q <- cut(r[w], quantile(r[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
	Q[w] <- q
	plot(NA,NA,  ylab="OS", main=paste(l, "terciles"),  xlab="Years after CR", ylim=c(0,1), xlim=c(0,10), xaxs="i", yaxs="i")
	abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	fit <- survfit(Surv(time1/365, time2/365, status) ~ q + transplantCR1, data=osData[w,])
	lines(fit, col=rep(sapply(2:0,function(x) colTrans(riskCol[l],x)), each=2), lty=c(1,3), mark=NA, xlab="Time after CR")
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
	segments(b[1,], (colMeans(coxRFXCirTD$Z)*coef(coxRFXCirTD))[crGroups != "Treatment"][o][w] ,b[3,], (colMeans(coxRFXCirTD$Z)*coef(coxRFXCirTD))[crGroups != "Treatment"][o][w])
	segments(b,t[,o][,w]-s, b,t[,o][,w]+s)
}

#p <- as.data.frame(PartialRisk(coxRFXCirTD)[1:nrow(clinicalData),])
p <- as.data.frame(PartialRisk(coxRFXOsCR)[1:nrow(clinicalData),])
s <- do.call("rbind",lapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(l) {
			w <- which(clinicalData$M_Risk==l)
			q <- cut(r[w], quantile(r[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
			t(sapply(split(p[w, ], q), colMeans) +.5 - colMeans(p))
		}))

#+relapseStars, fig.width=3,fig.heigh=3
c <- sapply(2:0, function(t) sapply(riskCol[c(2,4,3,1)], function(c) colTrans(c,t)))
g <- expand.grid(1:3,1:4-1)*3
stars(2*s[,c("Clinical","Demographics","Genetics","GeneGene","CNA","BT","Treatment")], scale=FALSE, col.stars=t(c), key.loc = c(13,0), locations=g, labels=NA)
symbols(g[,1], g[,2], circles=rep(1,12), inches=FALSE, add=TRUE)
text(1, 0:3*3, names(riskCol[c(2,4,3,1)]), pos=2)
text(1:3*3, 11, c("Best","Intermediate","Worst"), pos=3)

#' Find prototypes
prototypes <- sapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(l) sapply(1:3, function(i){
						d <- dist(t(t(coxRFXCirTD$Z[which(clinicalData$M_Risk[cirData$index]==l & Q==i &! is.na(clinicalData$CR_date[cirData$index])), ]) ))
						#d <- dist(t(t(coxRFXOsCR$Z[which(clinicalData$M_Risk[cirData$index]==l & Q==i &! is.na(clinicalData$CR_date[cirData$index])), ]) ))
						as.numeric(names(which.min(rowMeans(as.matrix(d), na.rm=TRUE))))
					}))

c <- sapply(2:0, function(t) sapply(riskCol[c(2,4,3,1)], function(c) colTrans(c,t)))
g <- expand.grid(1:3,1:4-1)*3
stars(2*t(t(p[prototypes,])- colMeans(p))[,c("Clinical","Demographics","Genetics","CNA","BT","Treatment")] +1, scale=FALSE, col.stars=t(c), key.loc = c(13,0), locations=g, labels=NA)
symbols(g[,1], g[,2], circles=rep(1,12), inches=FALSE, add=TRUE)
text(1, 0:3*3, names(riskCol[c(2,4,3,1)]), pos=2)
text(1:3*3, 11, c("Low","Intermediate","High"), pos=3)

#+ starsCR, fig.width=4, fig.height=4
s <- p - rep(colMeans(p), each=nrow(p))
w <- sapply(split(1:1540, paste(clinicalData$M_Risk, Q[1:1540])), `[`, 1:12)
w <- w[,!grepl("NA", colnames(w))][,c(4:6,10:12,7:9,1:3)]
l <- stars(s[w,c("Clinical","Demographics","Genetics","GeneGene","CNA","BT","Treatment")] + .5, scale=FALSE, col.stars = mapply(function(i,j) {t <- try(c[i,j]); if(class(t)=="try-error") NA else t}, as.character(clinicalData$M_Risk[w]),Q[w]), labels="")
symbols(l[,1],l[,2], circles=rep(0.5, nrow(l)), inches=FALSE,add=TRUE)

#' Who achieves CR?
c <- as.numeric(clinicalData$CR_date - clinicalData$ERDate)
e <- c < os[,1]
e[is.na(e)] <- 0
c[is.na(c) &! is.na(os[,1])] <- os[is.na(c) &! is.na(os[,1]),1]
cr <- Surv(time=pmin(c, os[,1]), event = e)
coxRFXCr <- CoxRFX(dataFrame[whichRFXOsTDGG], cr, groups=groups[whichRFXOsTDGG], which.mu=mainGroups)

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



#' 10. Germline polymorphisms
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
