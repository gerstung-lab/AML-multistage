# TODO: Add comment
# 
# Author: mg14
###############################################################################

mutationData = read.table("data/AML2_prelim_subs_for_peter2.txt", sep="\t", header=FALSE)
oncogenics <- (table(mutationData[,c(3,19)]) > 0)+0
rownames(oncogenics) <- sub("WGA_","",rownames(oncogenics))

interactions <- sapply(1:ncol(oncogenics), function(i) sapply(1:ncol(oncogenics), function(j) {f<- try(fisher.test(oncogenics[,i], oncogenics[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
odds <- sapply(1:ncol(oncogenics), function(i) sapply(1:ncol(oncogenics), function(j) {f<- try(fisher.test(oncogenics[,i] + .5, oncogenics[,j] +.5), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
diag(interactions) <- 0
diag(odds) <- 1
colnames(odds) <- rownames(odds) <- colnames(interactions) <- rownames(interactions) <- colnames(oncogenics)
odds[10^-abs(interactions) > 0.05] = 1
odds[odds<1e-3] = 1e-4
odds[odds>1e3] = 1e4
logodds=log10(odds)

pdf(paste(Sys.Date(),"-Interactions.pdf", sep=""), 4,4, pointsize = 8) ## HEATMAP OF ALL
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
dev.off()

clinicalData <- read.table("data/AML.clin.prelim4.txt", sep="\t", header=TRUE, na.strings = "na")
clinicalData$BM.blasts <- as.numeric(as.character(clinicalData$BM.blasts))

library(survival)
surv <- Surv(clinicalData$efs, clinicalData$efsstat)

m <- match(rownames(oncogenics), clinicalData$PDID)
source("~/Projects/hdcox/hdcox/R/ecoxph.R")


makeInteractions <- function(X,Y){
	do.call(cbind, lapply(X, `*`, Y))
}

makeInteger <- function(F){
	res <- do.call(cbind,lapply(levels(F), `==`, F))
	colnames(res) <- levels(F)
	res
}

library(RColorBrewer)
set1 <- brewer.pal(9, "Set1")
plot(survfit(surv[m] ~ TP53 + TPL, data=X), col=rep(set1[1:2],each=2), lty=1:2, mark=16)


X <- list(G = data.frame(oncogenics), 
		T = cbind(makeInteger(clinicalData$Study[m]), clinicalData$TPL[m]),
		C = clinicalData[m, c("Age", "gender","BM.blasts","wbc","LDH_")])
X$GG <- makeInteractions(data.frame(oncogenics), data.frame(oncogenics))
X$GG <- X$GG[,-seq(1,ncol(X$GG), ncol(oncogenics)+1)] 
X$GG <- X$GG[,colSums(X$GG)>0] 
X$GT <- makeInteractions(X$G, X$T)
X$GT <- X$GT[,colSums(X$GT, na.rm=TRUE) > 0]


c <- coxRFX(do.call(cbind,X), surv[m], groups=unlist(sapply(names(X), function(x) rep(x, ncol(X[[x]])))))

library(Hmisc)
rcorr.cens(-rowSums(partialRisk(c)), surv[m])
risk <- rowSums(partialRisk(c))
AUC.uno(surv[m][!is.na(risk)], surv[m][!is.na(risk)], risk[!is.na(risk)], 5000)

plot(survfit(surv[m][!is.na(risk)] ~ cut(risk[!is.na(risk)], 5)), col=rep(set1[1:2],each=2), lty=1:2, mark=16)

barplot(partialC(c))

set.seed(42)
testIx <- sample(c(TRUE,FALSE), length(m), replace=TRUE, prob=c(0.66,0.34))
d <- coxRFX(do.call(cbind,X)[testIx,], surv[m][testIx], groups=unlist(sapply(names(X), function(x) rep(x, ncol(X[[x]])))))

p <- as.matrix(do.call(cbind,X)[!testIx,]) %*% d$coefficients
rcorr.cens(- na.omit(p), na.omit(surv[m][!testIx]))
AUC.uno(na.omit(surv[m]), na.omit(surv[m][!testIx]), na.omit(p), 1000)




