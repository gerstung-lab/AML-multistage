# TODO: Add comment
# 
# Author: mg14
###############################################################################


#' Dirichlet Classification
#' =======================

#' Code run on 
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
load("AML.RData")


#' #### Preparation of data
#set.seed(42)
#binomialImpute <- function(x) {x[is.na(x)] <- rbinom(sum(is.na(x)), 1, mean(x, na.rm=TRUE)); return(x)}
d <- as.matrix(dataFrame[groups %in% c("Genetics","BT","CNA")])
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
#+ DPchecks
plot_hdp_lik(output$lik, burnin)
plot_hdp_numclass_raw(output$numclass)
plot_hdp_data_assigned(output$classqq, output$numclass)


#' AML classes
#+ classes, results='asis', cache=TRUE
posteriorMerged <- hdp_extract_signatures(output, prop.explained=0.99, cos.merge=0.95)
posteriorMeans <- Reduce("+",posteriorMerged$sigs_qq)/length(posteriorMerged$sigs_qq)
rownames(posteriorMeans) <- colnames(genotypesImputed)
kable(round(posteriorMeans,1))

#' Most prevalent lesions
genes <- apply(posteriorMeans, 2, function(x) paste(ifelse(x>10,rownames(posteriorMeans),"")[order(x, decreasing = TRUE)[1:5]], collapse=";"))
genes

#' Assignment from posterior samples
#+ DPassignment, fig.width=8
library(RColorBrewer)
col <- c(brewer.pal(9,"Set1")[c(9,1:8)], brewer.pal(8,"Dark2"))
posteriorProbability <- apply(sapply(posteriorMerged$sigs_nd_by_dp, colMeans)[,-1],2,function(x) (x+.Machine$double.eps)/sum(x+.Machine$double.eps))
o <- order(apply(posteriorProbability,2,which.max))
barplot(posteriorProbability[,o], col=col, border=NA, ylab="Probability", xlab="Patient")
data.frame(Prob=rowMeans(posteriorProbability), genes)

#' Classes
#+ DPclass
dpClass <- factor(apply(posteriorProbability, 2, which.max)-1)
table(dpClass)
plot(seq(0,1,l=ncol(posteriorProbability)),sort(apply(posteriorProbability,2,max)), type='l', ylim=c(0,1) , xlab="Fraction of patients", ylab="Assignment probability")
boxplot(apply(posteriorProbability,2,max) ~ dpClass, col=col, ylab="Probability", xlab="Class")

#+ DPbar, fig.width=8
o <- order(colSums(genotypesImputed), decreasing=TRUE)
barplot(t(sapply(split(as.data.frame(as.matrix(genotypesImputed)), dpClass), colSums)[o,]), col=col, las=2, legend=TRUE, border=NA, args.legend=list(border=NA))
abline(h=seq(100,500,100), col="white")

#' Clinical associations
#+ DPclinical
boxplot(clinicalData$wbc ~ factor(dpClass), log='y', xlab="Class",ylab="wbc", col=col)
boxplot(clinicalData$BM_Blasts ~ factor(dpClass), xlab="Class",ylab="Blast %", col=col)
boxplot(clinicalData$AOD ~ factor(dpClass), xlab="Class",ylab="Age", col=col)
boxplot(rowSums(genotypesImputed) ~ factor(dpClass), xlab="Class",ylab="# Mutations", col=col)

#' Survival
#+ DPsurvival
plot(survfit(os ~ dpClass), col=col)
legend("topright", legend = levels(dpClass), col=col, lty=1)
kable(summary(survfit(os ~ dpClass))$table)
summary(coxph(os ~ dpClass))

#' Phylogeny
library(ape)
plot(nj(dist(t(posteriorMeans/(rep(rowSums(posteriorProbability), each=nrow(posteriorMeans)))>.1))))

#' Alternatively: Naive Bayes assignment
sigBayes <- function(genotype, sigs){
	dSig <- function(sig,genotype){
		dmultinom(genotype, prob=sig/sum(sig))
	}
	lik <- apply(sigs,2, dSig, genotype)
	lik/sum(lik)
}

naiveBayes <- t(apply(genotypesImputed,1,sigBayes,posteriorMeans))

#' ## Session
sessionInfo()
