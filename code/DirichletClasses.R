#' AML classification using Dirichlet Processes
#' =======================
#' <script type="text/javascript">
#'     $(document).ready(function() {
#'         $('#posteriorMedian').DataTable();
#'     } );
#' </script>
#' 
#' Code run on 
options(markdown.HTML.header = "tmp.html")
system("hostname -f", intern=TRUE)
Sys.time()
getwd()
#' using
#+ run, eval=FALSE, echo=TRUE
library(knitr)
spin("DirichletClasses.R")


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

#' #### Libraries and data
library(CoxHD) # library(devtools); install_github("mg14/CoxHD/CoxHD")
library(mg14) # library(devtools); install_github("mg14/mg14")
library(hdp)
set1 <- brewer.pal(8, "Set1")
load("AML-2015-09-03.RData")


#' #### Preparation of data
#set.seed(42)
#binomialImpute <- function(x) {x[is.na(x)] <- rbinom(sum(is.na(x)), 1, mean(x, na.rm=TRUE)); return(x)}
d <- as.matrix(dataFrame[groups %in% c("Genetics","Fusions","CNA")])
d[!d %in% c(0,1)] <- NA
#genotypesImputed <- Matrix(apply(d,2, binomialImpute))
genotypesImputed <- Matrix(round(ImputeMissing(d)))

#' #### Setting up hdp
#+ DPsetup
n <- ncol(genotypesImputed)
shape <- 1
invscale <- 1
hdp <- hdp_init(ppindex=0, #index of the parent DP for initial DP
		cpindex=1, #index of alphaa and alphab for initial DP
		hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
		alphaa=shape,
		alphab=invscale)

hdp <- hdp_adddp(hdp,
		numdp=nrow(genotypesImputed), # one DP for every sample in that cancer type
		pp=1, # parent DP for group i is the i-th+1 overall DP because of the grandparent at position 1
		cp=1) # index of alphaa and alphab for each DP

# Assign the data from each patient to a child DP
hdp <- hdp_setdata(hdp = hdp, dpindex=1:nrow(genotypesImputed)+1, data=genotypesImputed)

# Activate the DPs with specified number of classes (signatures)
hdp <- dp_activate(hdp, 1:(nrow(genotypesImputed)+1), 5)

#' DP parameters
#+ DPparam
burnin <- 5000
postsamples <- 10000
spacebw <- 20
cpsamples <- 10

#' Run DP
#+ DPrun, cache=TRUE, cache.lazy=FALSE
set.seed(42)
output <- hdp_posterior(hdp, #activated hdp structure
		numburnin=burnin,
		numsample=postsamples,
		numspace=spacebw,
		doconparam=cpsamples)


#' DP checks
#+ DPchecks, cache=TRUE
plot(output$lik, type='l'); abline(v=burnin, lty=3)
plot(output$numclass, type='l')


#' AML classes
#+ posteriorMerged, cache=TRUE
posteriorMerged <- hdp_extract_signatures(output, prop.explained=0.99, cos.merge=0.95)

#+ classes, results='asis', cache=TRUE
#posteriorMeans <- Reduce("+",posteriorMerged$sigs_qq)/length(posteriorMerged$sigs_qq)
posteriorSamples <- array(unlist(posteriorMerged$sigs_qq), dim=c(dim(posteriorMerged$sigs_qq[[1]]), length(posteriorMerged$sigs_qq)))
rownames(posteriorSamples) <- colnames(genotypesImputed)
colnames(posteriorSamples) <- 1:ncol(posteriorSamples) -1
posteriorMeans <- rowMeans(posteriorSamples, dim=2)
posteriorQuantiles <- apply(posteriorSamples, 1:2, quantile, c(0.025,.5,0.975))
posteriorMode <- apply(posteriorSamples, 1:2, function(x) {t <- table(x); as.numeric(names(t)[which.max(t)])})
kable(posteriorQuantiles[2,,], "html", table.attr = 'id="posteriorMedian"') # Posterior median


#' Most prevalent lesions
genes <- apply(posteriorMeans, 2, function(x) paste(ifelse(x>10,rownames(posteriorMeans),"")[order(x, decreasing = TRUE)[1:5]], collapse=";"))
genes <- gsub(";+$","",genes)
genes

#' Assignment from posterior samples
#+ DPassignment, cache=TRUE, fig.width=8
library(RColorBrewer)
col <- c(brewer.pal(9,"Set1")[c(9,1:8)], brewer.pal(8,"Dark2"))
posteriorProbability <- apply(sapply(posteriorMerged$sigs_nd_by_dp, colMeans)[,-1],2,function(x) (x+.Machine$double.eps)/sum(x+.Machine$double.eps))
o <- order(apply(posteriorProbability,2,which.max))
barplot(posteriorProbability[,o], col=col, border=NA, ylab="Probability", xlab="Patient")
data.frame(Prob=rowMeans(posteriorProbability), genes)

#' ### Classes
#+ DPclass
dpClass <- factor(apply(posteriorProbability, 2, which.max)-1)
table(dpClass)
plot(seq(0,1,l=ncol(posteriorProbability)),sort(apply(posteriorProbability,2,max)), type='l', ylim=c(0,1) , xlab="Fraction of patients", ylab="Assignment probability")
boxplot(apply(posteriorProbability,2,max) ~ dpClass, col=col, ylab="Probability", xlab="Class")

#+ DPbar, fig.width=8
par(mar=c(6,3,1,1)+.1, cex=.8)
o <- order(colSums(genotypesImputed), decreasing=TRUE)
driverPrevalence <- t(sapply(split(as.data.frame(as.matrix(genotypesImputed)), dpClass), colSums)[o,])
b <- barplot(driverPrevalence, col=col, las=2, legend=TRUE, border=NA, args.legend=list(border=NA), names.arg=rep("", ncol(genotypesImputed)))
abline(h=seq(100,500,100), col="white")
rotatedLabel(b, labels=colnames(genotypesImputed)[o])

#' Driver signatures
#+ DPbars, fig.width=8, fig.height=2
par(mar=c(6,3,1,1)+.1, cex=.8)
t <- table(dpClass)
i <- 0; for(c in levels(dpClass)){i <- 1+i 
	b <- barplot(posteriorQuantiles[2,o,c]/t[i], col=col[i], las=2, legend=FALSE, border=NA,  names.arg=rep("", ncol(genotypesImputed)), main=paste("Class", i-1,t[i], "Patients","f =",round(t[i]/1540,2)), ylim=c(0,1))
	segments(b, posteriorQuantiles[1,o,c]/t[i], b, posteriorQuantiles[2,o,c]/t[i], col="white")
	segments(b, posteriorQuantiles[2,o,c]/t[i], b, posteriorQuantiles[3,o,c]/t[i], col=col[i])
	rotatedLabel(b, labels=colnames(genotypesImputed)[o])
}
#' ### Clinical associations
#+ DPclinical
boxplot(clinicalData$wbc ~ factor(dpClass), log='y', xlab="Class",ylab="wbc", col=col)
boxplot(clinicalData$BM_Blasts ~ factor(dpClass), xlab="Class",ylab="Blast %", col=col)
boxplot(clinicalData$AOD ~ factor(dpClass), xlab="Class",ylab="Age", col=col)
boxplot(rowSums(genotypesImputed) ~ factor(dpClass), xlab="Class",ylab="# Mutations", col=col)

#' ### Survival
#' Simple coxph
#+ DPsurvival
plot(survfit(os ~ dpClass), col=col)
legend("topright", legend = levels(dpClass), col=col, lty=1)
kable(summary(survfit(os ~ dpClass))$table)
summary(coxph(os ~ dpClass))

#' Risk variance
#+ RFX, cache=TRUE
library(CoxHD)
dataFrameOsTD <- dataFrame[tplSplitOs,]
dataFrameOsTD[which(tplIndexOs), grep("TPL", colnames(dataFrameOsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 
mainGroups <- grep("[A-Z][a-z]+[A-Z]",levels(groups), invert=TRUE, value=TRUE)
mainIdx <- groups %in% mainGroups
osTDIdx <- !grepl("TPL_efs", colnames(dataFrame))
mainIdxOsTD <- mainIdx & osTDIdx
whichRFXOsTDGG <- which((colSums(dataFrame)>=8 | mainIdxOsTD) & osTDIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%

coxRFXFitOsTDGGc <- CoxRFX(dataFrameOsTD[,whichRFXOsTDGG], osTD, groups[whichRFXOsTDGG], which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 
coxRFXFitOsTDGGc

d <- cbind(dataFrameOsTD[,whichRFXOsTDGG],DP=t(posteriorProbability)[tplSplitOs,-1])
coxRFXFitOsTDGGcDP <- CoxRFX(d, osTD, c(as.character(groups[whichRFXOsTDGG]),rep("DP", nlevels(dpClass)-1)), which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 
coxRFXFitOsTDGGcDP

PlotVarianceComponents(coxRFXFitOsTDGGcDP, col=col)
round(cov(PartialRisk(coxRFXFitOsTDGGcDP)),2)

#' ### Phylogeny
library(ape)
plot(nj(dist(t(posteriorMeans/(rep(rowSums(posteriorProbability), each=nrow(posteriorMeans)))>.1))))

#' ### Gene:Gene interactions
#' #### Population based
#+ computeInteractions, cache=TRUE
geneToClass <- factor(apply(posteriorMeans, 1,which.max) -1)
getOdds <- function(x) {
	f <- sapply(1:ncol(x), 
			function(i) sapply(1:ncol(x), 
						function(j) {
							if(j <= i) return(c(NA,NA))
							f<- try(fisher.test(x[,i], x[,j]), silent=TRUE)
							if(class(f)=="try-error") c(0,NA) 
							else if(f$estimate>1) c(-log10(f$p.val),f$estimate)
							else c(log10(f$p.val), f$estimate)}
				), 
			simplify="array")
	for(i in 1:2)
		f[i,,][upper.tri(f[i,,])] <- t(f[i,,])[upper.tri(f[i,,])]
	return(f)
}
f <- getOdds(genotypesImputed)
logPInt <- f[1,,]
odds <- f[2,,]
pairs <- sapply(1:ncol(genotypesImputed), function(i) colMeans(genotypesImputed * genotypesImputed[,i], na.rm=TRUE))
diag(logPInt) <- 0
diag(odds) <- 1
colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(genotypesImputed)
odds[odds<1e-3] = 1e-4
odds[odds>1e3] = 1e4

#+ plotInteractions, fig.width=6, fig.height=6
odds[10^-abs(logPInt) > 0.1] = 1
logOdds=log10(odds)
diag(logPInt) <- NA
par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)
ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
o = order(geneToClass, -colSums(genotypesImputed))
M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
M[lower.tri(M)] <- cut(logOdds[o,o][lower.tri(M)], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)
M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(genotypesImputed), breaks=c(-1,0,5,10,20,50,100,200,600))) + 9 
image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(9,"BrBG"), c("white",brewer.pal(7,"Reds"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
l <- colnames(logPInt)[o]
mtext(side=1, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=col[geneToClass[o]])
mtext(side=2, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=col[geneToClass[o]])
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
t <- c(0,table(geneToClass))
s <- cumsum(t)+.5
abline(h=s[-length(s)], col=col, lty=3)
abline(v=s[-length(s)], col=col, lty=3)
rect(s[-1],s[-1], s[-length(s)], s[-length(s)], border=col)


#' Expected heatmap
#+ heatmapExpected, fig.width=6, fig.height=6, cache=TRUE
set.seed(42)
t <- table(dpClass)
pp <- t(t(posteriorMeans)/as.numeric(t))
expectedOdds <- sapply(colnames(genotypesImputed), function(j){
			sapply(colnames(genotypesImputed), function(i){
						if(i==j) return(0)
						P <- Reduce("+",lapply(seq_along(t), function(k) {t[k] * (pp[i,k] * c(1,-1) + c(0,1)) %o% c(pp[j,k]* c(1,-1) + c(0,1))}))/sum(t)
						#res <- round(log10(M[1,1]*M[2,2]/M[1,2]/M[2,1]))
						#if( sum(M[,1]) * sum(M[1,])/sum(M) < 5 & res < 0) res <- 0
						#return(res)
						M <- matrix(rmultinom(1,sum(t), P), ncol=2)
						f <- fisher.test(round(M))
						res <- pmin(pmax(round(log10(f$estimate)),-4),4)
						if(f$p.value > 0.05)
							res <- 0
						return(res)
					})
		})

par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)
o = order(geneToClass, -colSums(genotypesImputed))
image(x=1:ncol(expectedOdds), y=1:nrow(expectedOdds), expectedOdds[o,o], col=brewer.pal(9,"BrBG"), breaks=c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(expectedOdds)+3), ylim=c(0, ncol(expectedOdds)+3))
l <- colnames(expectedOdds)[o]
mtext(side=1, at=1:ncol(expectedOdds), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=col[geneToClass[o]])
mtext(side=2, at=1:ncol(expectedOdds), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=col[geneToClass[o]])
abline(h=0:ncol(expectedOdds)+.5, col="white", lwd=.5)
abline(v=0:ncol(expectedOdds)+.5, col="white", lwd=.5)
t <- c(0,table(geneToClass))
s <- cumsum(t)+.5
rect(s[-1],s[-1], s[-length(s)], s[-length(s)], border=col)

#' #### Per class
#+ plotInteractionsClass, cache=TRUE, fig.width=6, fig.height=6
for(cls in levels(dpClass)){
	w <- dpClass == cls
	f <- getOdds(genotypesImputed[w,])
	logPInt <- f[1,,]
	odds <- f[2,,]
	pairs <- sapply(1:ncol(genotypesImputed), function(i) colMeans(genotypesImputed[w,] * genotypesImputed[w,i], na.rm=TRUE))
	diag(logPInt) <- 0
	diag(odds) <- 1
	colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(genotypesImputed)
	odds[odds<1e-3] = 1e-4
	odds[odds>1e3] = 1e4
	
	odds[10^-abs(logPInt) > 0.1] = 1
	logOdds=log10(odds)
	diag(logPInt) <- NA
	par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)
	ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
	o = order(geneToClass, -colSums(genotypesImputed))
	M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
	M[lower.tri(M)] <- cut(logOdds[o,o][lower.tri(M)], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)
	M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(genotypesImputed[w,]), breaks=c(-1,0,5,10,20,50,100,200,600))) + 9 
	image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(9,"BrBG"), c("white",brewer.pal(7,"Reds"))), breaks=0:17, xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
	l <- colnames(logPInt)[o]
	mtext(side=1, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=col[geneToClass[o]])
	mtext(side=2, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=col[geneToClass[o]])
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
	t <- c(0,table(geneToClass))
	s <- cumsum(t)+.5
	rect(s[-1],s[-1], s[-length(s)], s[-length(s)], border=col)
	title(main=paste0("Class ",cls,": ", genes[cls]))
}

#' Alternatively: Naive Bayes assignment
sigBayes <- function(genotype, sigs){
	dSig <- function(sig,genotype){
		dmultinom(genotype, prob=sig/sum(sig))
	}
	lik <- apply(sigs,2, dSig, genotype)
	lik/sum(lik)
}

naiveBayes <- t(apply(genotypesImputed,1,sigBayes,posteriorMeans))

#' Save
names(dpClass) <- clinicalData$PDID
save(dpClass, posteriorMeans, posteriorQuantiles, posteriorProbability, file='dpClass.RData')

#' ## Curated classification
c <- read.table("../data/reduced_classes.txt", header=TRUE, sep="\t")
curatedClass <- c$ReductionClass
names(curatedClass) <- c$Sample
rm(c)

library(clue)
t <- table(dpClass, curatedClass)
s <- solve_LSAP(t, maximum=TRUE)
t[,c(s, setdiff(1:ncol(t),s))]

#' ### Associations

X <- as.matrix(MakeInteger(curatedClass))
Z <- as.matrix(genotypesImputed)
Y <- as.matrix(dataFrame[groups %in% c("Clinical","Demographics")])

cv.glm <- function(x,y, fold=5, family="gaussian"){
	cvIdx <- sample(1:nrow(x)%% fold +1 )
	p <- numeric(length(y))
	for(i in 1:fold){
		p[cvIdx==i] <- predict(glm(y ~ ., data=x, subset=cvIdx!=i, family=family), newdata=x[cvIdx==i,])
	}
	if(family=="gaussian"){
		m <- mean((p-y)^2)
		s <- sd(sapply(1:fold, function(i) mean((p-y)[cvIdx==i]^2)))/sqrt(fold)
		return(c(avg=m, sd=s))
	}
	else if(family=="binomial"){
		m <- performance(prediction(p, y), "auc")@y.values[[1]]
		s <- sd(sapply(1:fold, function(i)  performance(prediction(p[cvIdx==i], y[cvIdx==i]), "auc")@y.values[[1]]))/sqrt(fold)
		return(c(avg=m, sd=s))
	}
}

#' #### White counts
#+ wbc_box
y <- log(clinicalData$wbc)
boxplot(y ~ factor(curatedClass), ylab="log wbc", las=2)
#+ wbc_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="White counts")


#' #### Bone marrow blasts
#+ BM_box
y <- car::logit(clinicalData$BM_Blasts/100 )
boxplot(y ~ factor(curatedClass), ylab="logit BM blasts", las=2)
#+ BM_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="BM blasts")

#' #### PB blasts
#+ PB_box
y <- car::logit(clinicalData$PB_Blasts/100 )
boxplot(y ~ factor(curatedClass), ylab="logit PB blasts", las=2)
#+ PB_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="PB blasts")

#' #### Age
#+ Age_box
y <- clinicalData$AOD
boxplot(y ~ factor(curatedClass), ylab="Age", las=2)
#+ Age_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="Age")

#' #### LDH
#+ LDH_box
y <- log(clinicalData$LDH)
boxplot(y ~ factor(curatedClass), ylab="log LDH", las=2)
#+ LDH_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="LDH")

#' #### Platelets
#+ platelets_box
y <- log(clinicalData$platelet)
boxplot(y ~ factor(curatedClass), ylab="log platelets", las=2)
#+ platelets_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="Platelets")

#' #### HB
#+ HB_box
y <- log(clinicalData$HB)
boxplot(y ~ factor(curatedClass), ylab="log HB", las=2)
#+ HB_bar, fig.width=1.5
set.seed(42)
mse0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), na.omit(y))
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), type="mse", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
mse1 <- c(min(g$cvm), g$cvsd[which.min(g$cvm)])
v <- var(y, use='c')
a <- (v-c(subtypes=mse0[1], `subtypes+genomics`=min(mse1[1])))/v*100
barplot(rbind(a,100-a), ylab="Explained variance (%)", names=rep("",2), ylim=c(0,100)) -> b
segments(b,(v-c(mse0[1]-mse0[2],mse1[1]-mse1[2]))*100/v,b,(v-c(mse0[1]+mse0[2], mse1[1]+mse1[2]))*100/v)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="HB")

oddsPlot <- function(t){
	plot(x=rep(1:2, ncol(t)), t[], ylab="Number of cases", xlab= names(dimnames(t))[1], xaxt="n", log='y')
	mtext(side=1,at=c(1,2), text=paste0(rownames(t), ", n=", rowSums(t)), las=1, pch=16)
	segments(1,t[1,],2, t[2,])
	p <- sapply(1:ncol(t), function(i) {f <- fisher.test(cbind(rowSums(t[,-i]), t[,i])); c(p.value=f$p.value, OR=f$estimate[1], f$conf.int)})
	mtext(side=4, at=t[2,], text=paste0(colnames(t), ", OR=", format(p[2,],digits=1), "(", apply(round(p[3:4,],2),2,paste, collapse="-"),")", sig2star(p[1,])))
}

#' #### Splenomegaly
#+ Splenomegaly_bar, fig.width=1.5
library(ROCR)
set.seed(42)
y <- clinicalData$Splenomegaly
table(Splenomegaly=y,factor(curatedClass))
auc0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), y[!is.na(y)], family="binomial")
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), family='binomial',type="auc", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
auc1 <-  c(max(g$cvm), g$cvsd[which.max(g$cvm)])
a <- c(subtypes=auc0[1], `subtypes+genomics`=auc1[1])*100
barplot(rbind(a,150-a)-50, ylab="AUC (%)", main="Splenomegaly",  offset=50, names=rep("",2),  ylim=c(50,100))
segments(b,c(auc0[1]-auc0[2],auc1[1]-auc1[2])*100,b,c(auc0[1]+auc0[2], auc1[1]+auc1[2])*100)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="Splenomegaly")


#' #### Gender
#+ Gender_bar, fig.width=1.5
set.seed(42)
y <- clinicalData$gender -1
table(Gender=factor(y, labels=c('male','female')),factor(curatedClass))
auc0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), y[!is.na(y)], family="binomial")
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), family='binomial',type="auc", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
auc1 <-  c(max(g$cvm), g$cvsd[which.max(g$cvm)])
a <- c(subtypes=auc0[1], `subtypes+genomics`=auc1[1])*100
barplot(rbind(a,150-a)-50, ylab="AUC (%)",   offset=50, names=rep("",2),  ylim=c(50,100))
segments(b,c(auc0[1]-auc0[2],auc1[1]-auc1[2])*100,b,c(auc0[1]+auc0[2], auc1[1]+auc1[2])*100)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="Gender")

#' #### CR
#+ CR_bar, fig.width=1.5
set.seed(42)
y <- !is.na(clinicalData$CR_date)
y[is.na(clinicalData$CR_date) & clinicalData$OS==0] <- NA
table(CR=y,factor(curatedClass))
auc0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), y[!is.na(y)], family="binomial")
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), family='binomial',type="auc", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
auc1 <-  c(max(g$cvm), g$cvsd[which.max(g$cvm)])
a <- c(subtypes=auc0[1], `subtypes+genomics`=auc1[1])*100
barplot(rbind(a,150-a)-50, ylab="AUC (%)",  offset=50, names=rep("",2),  ylim=c(50,100))
segments(b,c(auc0[1]-auc0[2],auc1[1]-auc1[2])*100,b,c(auc0[1]+auc0[2], auc1[1]+auc1[2])*100)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="Complete remission")

#' #### OS
#+ OS_bar, fig.width=1.5
set.seed(42)
y <- os[1:1540,2]
y[os[1:1540,1] < 3 * 365 & os[1:1540,2]==0] <- NA
table(OS=factor(y,labels=c("alive","dead")),factor(curatedClass))
auc0 <- cv.glm(as.data.frame(X[!is.na(y),-1]), y[!is.na(y)], family="binomial")
g <- cv.glmnet(cbind(X,Z)[!is.na(y),], na.omit(y), family='binomial',type="auc", alpha=1, penalty.factor=c(rep(0,ncol(X)), rep(1,ncol(Z))))
auc1 <-  c(max(g$cvm), g$cvsd[which.max(g$cvm)])
a <- c(subtypes=auc0[1], `subtypes+genomics`=auc1[1])*100
barplot(rbind(a,150-a)-50, ylab="AUC (%)",  offset=50, names=rep("",2),  ylim=c(50,100))
segments(b,c(auc0[1]-auc0[2],auc1[1]-auc1[2])*100,b,c(auc0[1]+auc0[2], auc1[1]+auc1[2])*100)
rotatedLabel(b, labels=c("subtypes","subtypes+genomics"))
title(main="Overall survival at 3yr")



#' ## Session
devtools::session_info()
