# TODO: Add comment
# 
# Author: mg14
###############################################################################

data <- cbind(dataList$Genetics, dataList$Cytogenetics)
t <- table(rowSums(data, na.rm=TRUE))
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

colRamp <- colorRampPalette(set1[c(3,2,4,1,5)])

t <- table(rowSums(genomicData, na.rm=TRUE))
pdf("nOnc.pdf",4,3,pointsize=8)
par(mar=c(3,3,.2,.2)+.1, mgp=c(2,0.5,0), bty="n")
barplot(t/1540,col=colRamp(length(t)), xlab="Oncogenic mutations", ylab="Relative frequency", border=NA)
abline(h=seq(0.05,0.25,0.05), col="white")
dev.off()

pdf("nOncSurv.pdf",3,2.5, pointsize=8)
par(mar=c(3,3,0,0)+.1, mgp=c(2,0.5,0), bty="n")
plot(survfit(rfs ~ rowSums(data, na.rm=TRUE)), col=colRamp(11), xlab="Days", ylab="RFS", mark=NA)
legend("topright", legend=0:10, lty=1, bty="n",  col=colRamp(11), cex=0.66, title="# Onc")
dev.off()


pdf("allKM.pdf",3,2.5, pointsize=8)
par(mar=c(3,3,0,0)+.1, mgp=c(2,0.5,0), bty="n")
for(j in colnames(data)){
	plot(survfit(survival ~ data[[j]]), col=c("grey", set1[1]), main=j, xlab="Days", ylab="EFS")
}
dev.off()

pdf("pie.pdf", 3,3, pointsize=8)
partRiskVar <- PartialRiskVar(coxRFXFitTD)
x=c(varianceComponents, Residual=mean(rowSums(partRiskVar)))
pie(x, col=c(col1, "grey"), labels = paste(names(x), round(x/sum(x),2)))
dev.off()


par(mar=c(5,7,1,1))
o <- order(coxRFXFitTDMain$mu)
boxplot(coef(coxRFXFitTDMain)/log(2) ~ factor(coxRFXFitTDMain$groups, levels=levels(coxRFXFitTDMain$groups)[o]), border=col1[o], horizontal=TRUE, las=1, lty=1, pch=16, cex=.66, staplewex=0)
abline(v=0)

pdf("coefHaz.pdf",3.2,2.5, pointsize=8)
par(mar=c(3,6,1,3.5),mgp=c(2,0.5,0), bty="n")
p <- PartialRisk(coxRFXFit)
r <- rowSums(p)
H0 <- basehaz(coxRFXFit, centered = FALSE)
hazardDistNrm <- splinefun(H0$time, H0$hazard, method="monoH.FC")
x <- 10^seq(1,log10(2500), 0.1)
plot(x, exp(-hazardDistNrm(x)), log="x", ylim=c(0,1), type="l")
xx <- seq(250, 2000, 250)
o <- order(coxRFXFit$mu)
boxplot(exp(-exp(coef(coxRFXFit))) / exp(hazardDistNrm(xx[order(o)[coxRFXFit$groups]]))~ factor(coxRFXFit$groups, levels=levels(groups)[o]), at=xx,border=col1[o], horizontal=FALSE, las=1, lty=1, pch=NA, cex=.66, staplewex=0, add=TRUE, xaxt="n", yaxt="n", boxwex=200, col="NA")

boxplot(-coef(coxRFXFit) + ~ factor(coxRFXFit$groups, levels=levels(groups)[o]), border=col1[o], horizontal=FALSE, las=1, lty=1, pch=NA, cex=.66, staplewex=0, ylab="", ylim=c(-1,1), yaxt="n", labels="", xaxt="n")
abline(h=0)
x <- seq(-1,1,0.01)
p <- pretty(0:2000)
axis(side=1, at = p/250+.5,labels=p ,las=1)
dev.off()

par(mfrow=c(3,3))
i <- 1
for(g in levels(groups)){
	plot( x, exp(-hazardDistNrm(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	polygon( c(x, rev(x)), exp(-hazardDistNrm(c(x,rev(x)))*exp( rep(range(coxRFXFit$coef[groups[whichRFX] == g]), each=length(x)) + mean(r) )), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDistNrm(c(x,rev(x)))*exp( rep(sqrt(coxRFXFit$sigma2[g]) * c(-1,1), each=length(x)) + mean(r) + coxRFXFit$mu[g])), col=paste(col1[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDistNrm(x)*exp( coxRFXFit$mu[g] + mean(r))), col=col1[i], type="l", lwd=2)	
	lines( x, exp(-hazardDistNrm(x) * exp(mean(r))), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}

par(mfrow=c(3,3))
i <- 1
for(g in levels(groups)){
	plot( x, exp(-hazardDistNrm(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	m <- mean(rowSums(p[,colnames(p)!=g]))
	polygon( c(x, rev(x)), exp(-hazardDistNrm(c(x,rev(x)))*exp( rep(range(partRiskTD[,g]), each=length(x)) +m)), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDistNrm(c(x,rev(x)))*exp( rep(quantile(partRiskTD[,g], c(0.25,0.75)), each=length(x))+m)), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDistNrm(c(x,rev(x)))*exp( rep(quantile(partRiskTD[,g], c(0.05,0.95)), each=length(x))+m)), col=paste(col1[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDistNrm(x)*exp( median(partRiskTD[,g])+m)), col=col1[i], type="l", lwd=2)	
	lines( x, exp(-hazardDistNrm(x) *exp(+m)), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}


pdf("barplot.pdf", 7,2.5, pointsize=8)
par(mar=c(5,4,1,1))
w <- c("Genetics","CNA","BT")
s <- sort(colSums(dataFrame[groups %in% w]), d=TRUE)
barplot(s/nrow(dataFrame), col=colGroups[as.character(groups[names(s)])], las=2, border=NA, ylab="Relative frequency", names=NA)
rotatedLabel(.Last.value, rep(0, length(s)), names(s), cex=.5)
legend("topright", fill=colGroups[w], w, bty="n", border=NA)
abline(h=seq(0.05,0.25,0.05), col="lightgrey")
dev.off()

load("../cv.RData")
cv <- as.data.frame(tmp)
colnames(cv) <- grep(":", colnames(coxRFXFitTD$X), value=TRUE, invert=TRUE)
tmp <- as.matrix(cv)
o <- order(colSums(cv[-9,]))
barplot(as.matrix(cv[,o]), col=c(col1,"grey"), las=2, border=NA, space=0)
1/(exp(0) +1)

plot(rForest$importance, sum(varianceComponents) - colSums(tmp[-9,]), log="xy")
text(rForest$importance, sum(varianceComponents) - colSums(tmp[-9,]), colnames(tmp))


barplot(rbind(colSums(pmin(c,0)),pmax(c,0)), col=c("white", col1[colnames(c)]), border=NA, width=0.5, space=c(0.4,rep(0.7/0.5,8)), density=36, angle=-45) -> b
barplot(pmin(c,0), col=col1[colnames(c)], add=TRUE, density=36, border=NA, width=.5, space=.7/.5 , axes=FALSE, names.arg=rep("", n))
barplot(colSums(c), col=col1[colnames(c)], add=TRUE, width=.5, space=.7/.5 , axes=FALSE,names.arg=rep("", n))

layout(matrix(1:2, ncol=2), c(2,1), c(1))
image(1:9,1:9,c, col=colorRampPalette(rev(brewer.pal(11,"PRGn")))(100), breaks=seq(-.35,.35,l=101), xaxt='n', yaxt='n')
segments(1:9 - .4,rep(1:9, each=9) - .4*c, 1:9+.4, rep(1:9, each=9)+.4*c)
#text(1:9, rep(1:9, each=9), round(signif(c,2),3))
u <- par("usr")
barplot(colSums(c), horiz=TRUE, col=col1[colnames(c)], width=0.8, space=c(.8,rep(.2/.8,8)), ylim=u[1:2], yaxs="i")

MakeLego <- function(c, bricks=10){
	n <- ncol(c)
	m <- nrow(c)
	im <- matrix(0, ncol=bricks*n, nrow=bricks*m)
	for(i in 1:m)
		for(j in 1:n){
			l <- round(abs(100*c[i,j]))
			if(l==0)
				next
			s <- rep(sign(c[i,j]), l)
			ii <- 0:(length(s)-1) %/% bricks + 1
			jj <- 0:(length(s)-1) %% bricks + 1
#		if(sign(c[i,j])==1)
#			jj <- jj + m/2
#		else
#			jj <- m/2 - jj +1
			for(k in 1: length(s)){
				im[(i-1)*bricks + ii[k],(j-1)*bricks+jj[k]] <- s[k]
			}
		}
	return(im)
}
im <- MakeLego(c)
is <- MakeLego(matrix(rowSums(pmax(c,0)), ncol=1)) + 2* MakeLego(matrix(rowSums(pmin(c,0)), ncol=1))
it <- MakeLego(matrix(sum(c), nrow=1))


n <- ncol(c)
m <- nrow(c)

PlotLego <- function(im, m = nrow(im)/10, n=ncol(im)/10){
	x <- seq(0,n, l=ncol(im)+1)
	y <- seq(0,m, l=nrow(im)+1)
	plot(NA,NA, xlim=c(0,n), ylim=c(0,m), xaxt='n', yaxt='n',xlab='',ylab='', xaxs='i', yaxs='i')
	u <- par("usr")
	rect(u[1],u[3],u[2],u[4], col=grey(.95), border=NA)
	abline(v=1:(n-1), col='white')
#abline(h=0:(n-1)+.5, lty=2)
	abline(h=1:(m-1), col='white')
	for(j in seq_along(x)[-1])
		for(i in seq_along(y)[-1]){
			if(im[i-1,j-1]==0)
				next
			rect(x[j-1],y[i-1], x[j],y[i], col=ifelse(im[i-1,j-1]==-1, brewer.pal(3,"RdGy")[1], "black"), lwd=.5, border='white')
		}
}

layout(matrix(c(1,2,4,3), ncol=2, byrow=TRUE), c(9,1), c(9,1))
par(bty='n', mar=c(1,1,1,1))
#image(seq(0,n, l=ncol(im)+1), seq(0,n, l=nrow(im)+1),im, col=c("white",grey(.9),'darkgrey'), xaxt='n', yaxt='n', xlab='', ylab='n')
PlotLego(im)
par(mar=c(1,0,1,0.5))
#image(seq(0,1, l=nrow(is)+1), seq(0,n, l=ncol(is)+1),is, col=brewer.pal(3,"RdGy"), xaxt='n', yaxt='n', xlab='', ylab='n')
#abline(h=1:(n-1))
PlotLego(is)
par(mar=c(1,0,0,0.5))
image(seq(0,1, l=nrow(it)+1), seq(0,1, l=ncol(it)+1),it, col=brewer.pal(3,"RdBu"), xaxt='n', yaxt='n', xlab='', ylab='n', breaks=c(-2:1))

#barplot(colSums(c), horiz=TRUE, col=col1[colnames(c)], width=0.8, space=c(.2,rep(.2/.8,8)), ylim=u[1:2], yaxs="i", names.arg=rep('',n))


m <- split(1:nrow(dataFrameOsTD), tplSplitOs)
set.seed(42)
b <- mclapply(1:10, function(i) {
			s <- sample(1:nrow(dataFrame), replace=TRUE)
			u <- unlist(m[s]); 
			w <- which((colSums(dataFrame[s,])>=8 | mainIdx) & osIdx) # ie, > 0.5%
			CoxRFX(dataFrameOsTD[u,w], osTD[u], groups=groups[w], sigma0=0.1, nu=0)
		}, mc.cores=10)

set.seed(42)
s <- sample(1:nrow(dataFrame), replace=TRUE)
u <- unlist(m[s]); 
w <- whichRFXOs

c <- CoxRFX(dataFrameOsTD[u,w], osTD[u], groups=groups[w], sigma0=0.1, nu=0)

r <- coxRFXFitOsTD$X %*% coxRFXFitOsTD$coefficients
simSurvival <- SimSurvNonp(risk=r, H0=basehaz(coxph(osTD ~ r), centered=FALSE), surv=os)

simSurvivalTD <-  Surv(c(rep(0,nrow(dataFrame)),simSurvival[tplSplitOs[(nrow(dataFrame)+1):nrow(dataFrameOsTD)],1]),
		c(simSurvival[1:nrow(dataFrame),1], simSurvival[(nrow(dataFrame)+1):nrow(dataFrameOsTD),1]),
		c(simSurvival[1:nrow(dataFrame),2], 1-simSurvival[(nrow(dataFrame)+1):nrow(dataFrameOsTD),2]+2))



u <- unique(unlist(m[s]))
cc <- CoxRFX(dataFrameOsTD[u,whichRFXOs], osTD[u], groups=groups[whichRFXOs], sigma0=0.1, nu=0)


intIx <- clinicalData$Study == "tr98A" & !clinicalData$ereignart %in% c("Death_in_Induktion","RD_after_Ind1","RD_or_PR_after_Ind2") & (rowSums(clinicalData[,c("inv3_t3_3","minus5_5q","minus7","minus7q", "mono12_12p_abn12p" , "mono17_17p_abn17p","complex","t_8_21" )])==0 | clinicalData$inv16_t16_16)
highIx <- clinicalData$Study == "tr98A" & ! clinicalData$ereignart == "Death_in_Induktion" & (clinicalData$ereignart %in% c("RD_after_Ind1","RD_or_PR_after_Ind2") | (rowSums(clinicalData[,c("inv3_t3_3","minus5_5q","minus7","minus7q", "mono12_12p_abn12p" , "mono17_17p_abn17p","complex" )])>0 & !clinicalData$inv16_t16_16 &! clinicalData$t_8_21))
lowIx <- clinicalData$Study == "tr98A" & !clinicalData$ereignart %in% c("Death_in_Induktion","RD_after_Ind1","RD_or_PR_after_Ind2") & clinicalData$t_8_21

r <- factor(I(totalRiskOsTD - coxRFXFitOsTD$coefficients["TPL_os"] * dataFrameOsTD[,"TPL_os"]>1.6), labels=c("low","high"))[intIx[tplSplitOs]]
s <- survfit(osTD[intIx[tplSplitOs]] ~ TPL_os + r, data=dataFrameOsTD[intIx[tplSplitOs],])
plot(s, col=set1)
coxRFXOsTDTrainNon98Int <- CoxRFX(dataFrameOsTD[!intIx[tplSplitOs],whichRFXOs], osTD[!intIx[tplSplitOs]], groups[whichRFXOs], sigma0=0.1, nu=0)
riskCoxRFXOsTDTrainNon98Int <- as.matrix(dataFrameOsTD[tplSplitOs,whichRFXOs]) %*% coef(coxRFXOsTDTrainNon98Int)

r <- riskCoxRFXOsTDTrainNon98Int - dataFrameOsTD$TPL_os * coxRFXOsTDTrainNon98Int$coefficients["TPL_os"]

s <- survfit(osTD[intIx[tplSplitOs]] ~ TPL + highrisk, data=data.frame(TPL =  dataFrameOsTD$TPL_os, highrisk=r>median(r, na.rm=T))[intIx[tplSplitOs],])
plot(s, col=set1[1:4])
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
coxph(osTD[intIx[tplSplitOs]] ~ TPL_allo + highrisk , data=data.frame(TPL_allo = clinicalData$TPL_type[tplSplitOs] %in% c("ALLO"), highrisk=r>median(r))[intIx[tplSplitOs],])



#### More time-dependence
nonLethalEventIdx <- clinicalData$ereignart %in% levels(clinicalData$ereignart)[3:5]
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

summary(coxph(Surv(time1, time2, event) ~ ., data=cbind(d[,-8], dataFrame[d$i, grep("(TPL)|(:)", selectedIntOs, value=TRUE, invert = TRUE)])))


### Risk viz
risk <- riskCPSSIntOsTD[dataFrameOsTD$TPL_os==0]
x <- seq(from=-4,to=4, l=512)
r <- sapply(levels(clinicalData$M_Risk)[c(2,3,4,1)], function(r){
			i <- clinicalData$M_Risk[tplSplitOs]==r
			d <- density(na.omit(risk[i]), from=-4,to=4)$y * mean(i, na.rm=TRUE)
		})
plot(x,rowSums(r), type='l', lty=1,xlab="log hazard", ylab="Prop. patients", col=set1[1], lwd=2)
polygon(x,rowSums(r), col=colTrans(set1[1]), border=NA)


H0 <- basehaz(coxFitCPSSIntOsTD, centered=TRUE)
hazardDistNrm <- splinefun(H0$time, H0$hazard, method="monoH.FC")
#plot(survfit(osTD ~1))
#lines(0:5000, exp(-hazardDist(0:5000)),  col='red')
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
for(i in seq_along(l))
lines(x, invHazardDist(-log(l[i]) /exp(x) )/10000, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
mtext(side=4, line=2.5, "Time")
mtext(side=3, at = log(-log(l)/hazardDistNrm(par("usr")[4]*10000)), text=paste(100*l, "% survive", sep=""))


concordanceFromVariance <- function(x) {
	.cfv <- function(x){
		stopifnot(x >= 0)
		if(x == 0)
			return(.5)			
		f <- function(x, sigma2) dnorm(x, sd=sqrt(2*sigma2)) / (1 +  exp(- x))
		2*integrate(f, 0, 10*sqrt(x), sigma2=x)$value
	}
	sapply(x, .cfv)
}

x <- 10^seq(-4,3,0.1)

c <- sapply(x, concordanceFromVariance)
plot(x,c, type='l', xlab="Variance", ylab="Concordance", log='x')

C <- CoxRFX(dataFrame[groups=="Genetics"], os, sigma0=0.1)

o <- order(colMeans(dataFrame[groups=="Genetics"]), decreasing=TRUE)

R <- t(t(as.matrix(dataFrame[groups=="Genetics"][o])) * C$coefficients[o])
CR <- apply(apply(R,1,cumsum),1,var)
M <- cumsum(colMeans(dataFrame[groups=="Genetics"][o]))
plot(M,CR, xlab="Mutations/patient", ylab="Explained var.", xlim=c(0,3.7), ylim=c(0,3.7*C$sigma2))
points(3.7, 3.7*C$sigma2, col="red", cex=2)
abline(0,C$sigma2)

cToV <- splinefun(c(0.5,c), c(0,x), method="monoH.FC")

x <- PartialRisk(coxRFXFitOsTD)
x <- x[1:nrow(dataFrame),]
x <- x - rep(colMeans(x), each=nrow(x))
t <- sd(x)
x <- x/(2*t) + 1 


d <- sapply(1:ncol(x), function(i) density(x[,i], bw=diff(range(x[,i]/50)))[1:2])
q <- apply(x, 2, quantile, c(0.1,0.5,0.9))
sig <- apply(x,2, function(x) mean(x)+ (-1:1) * sd(x))
rotation <- function(theta) matrix(c(cos(theta), sin(theta),  -sin(theta),cos(theta)), ncol=2)
j <- apply(x, 2, function(y) {v <- violinJitter(y)$y; v/max(v)})/5

c <- cos(seq(0,2*pi, l=10))[-10]
s <- sin(seq(0,2*pi, l=10))[-10]
plot(t(x) * c - t(j) *s, t(x) * s + t(j) * c, pch=16, cex=.5, col=col1)


pdf("tmp.pdf", fillOddEven=TRUE)
par(bty="n", xpd=NA)
plot(NA,NA, xlim=c(-1,1)*max(x), ylim=c(-1,1)*max(x), xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:nrow(dataFrame)){
	polygon(x[i,]*c, x[i,]*s, col=NA, border="#DDDDDDAA")
}
symbols(rep(0,3),rep(0,3),circles=log(c(0.66,1,1.5))/(2*t)+1, inches=FALSE, add=TRUE, fg="white", lwd=c(1,2,1))
#polygon(t(q[c(3,1),c(1:9,1)])*c[c(1:9,1)], t(q[c(3,1), c(1:9,1)])*s[c(1:9,1)], col='#88888888')
points(t(x) * c - t(j) *s, t(x) * s + t(j) * c, pch=16, cex=.25, col=col1[colnames(x)])
#for(i in nrow(q):1){
#	polygon(q[i,]*c, q[i,]*s, col=NA, lty=1, lwd=c(1,2,1)[i])
#}
segments(q[1,]*c, q[1,]*s,q[3,]*c,q[3,]*s, lwd=2)
#points(q[2,]*c, q[2,]*s, pch=16, cex=.25, col=col1[colnames(x)])
m <- apply(x,2,max)
for(i in seq_along(m))
text(colnames(x)[i],x=(m[i] + strwidth(colnames(x)[i], units="user")/1.5)*c[i], y=(m[i] + strwidth(colnames(x)[i], units="user")/1.5)*s[i], srt = ((360/9 *(i-1)+90) %% 180) -90 )
text(log(c(0.66,1,1.5))/(2*t)+1, c(0,0,0)-.1,labels=c(0.66,1,1.5), cex=.33)
dev.off()
#stars(q[,,drop=FALSE], scale=FALSE, locations=matrix(0,ncol=2,nrow=3), add=TRUE)

for(i in 1:ncol(d))
polygon(cbind(d[,i]$x, d[,i]$y/max(d[,i]$y*4)) %*% t(rotation((i-1)*2*pi/9)), col=col1[i])


pdf("Schematic.pdf",8,8, pointsize=8)
layout(matrix(1:12, nrow=2), width=c(8,1,1,1,1,2),height=c(1,4))
par(bty="n", mar=c(1,4,3,2), mgp=c(2.5,1,0))
w <- TRUE#coxRFXFitOsTD$groups %in% mainGroups
f <- factor(as.character(coxRFXFitOsTD$groups[w]))
o <- order(f,coxRFXFitOsTD$coefficients[w])
c <- coef(coxRFXFitOsTD)[w][o]
plot(c, type='p',  ylab="Coefficient", xaxt='n', col=col1[as.character(f)[o]], xaxs="i", lwd=2, pch=NA)
segments(seq_along(c), coxRFXFitOsTD$mu[as.character(f)[o]], seq_along(c), c, col=col1[as.character(f)[o]])
mtext(side=3, levels(f), at= table(f)/2+ cumsum(c(0,table(f)[-4])))
par(bty="n", mar=c(4,4,1,2), mgp=c(2.5,1,0))
x <- as.matrix(dataFrame[w][o])
p <- order(x %*% c)
x <- t(x)/apply(x,2,max)
h <- hclust(dist(t(x)))
image(x[,p]*.9 + as.numeric(f[o])-1, useRaster=TRUE, 
		col=c(brewer.pal(9,"Reds"),  brewer.pal(9,"Greens")[c(1,1:9)], brewer.pal(9,"Blues")[c(1,1:9)], brewer.pal(9,"Oranges")[c(1,1:9)],"white"), 
		breaks=seq(0,4, 0.1), xaxt="n",yaxt="n", ylab="Patients",xlab="Variable")
x <- as.matrix(dataFrame[w][p,o])
for(l in levels(f)){
	plot.new()
	par(mar=c(4,0,1,0.5), xpd=NA)
	h <- x[,f[o]==l] %*% c[f[o]==l]
	h <- h - mean(h)
	plot(h, seq_along(h), pch=NA, xlab="log hazard", yaxs="i", yaxt="n", ylab='', xlim=c(-2,2))
	segments(h, seq_along(h),0, seq_along(h), col=col1[l])
	mtext(side=3, l, cex=.66, line=.5)
}

h <- x %*% c
h <- h - mean(h)

H0 <- basehaz(coxph(os ~ 1), centered=TRUE)
hazardDistNrm <- splinefun(H0$time, H0$hazard, method="monoH.FC")
#plot(survfit(osTD ~1))
#lines(0:5000, exp(-hazardDist(0:5000)),  col='red')
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
i <- 1
x <- seq(from=-4,to=4, l=512)
par(mar=c(1,.5,3,0.5), xpd=FALSE)
plot(x, invHazardDist(-log(l[i]) /exp(x) ), col='black', lty=c(2,1,2)[i], xlab="", ylab="Time", type='l',ylim=c(0,5000), xlim=range(h), xaxt="n")
abline(v=seq(-2,4,2), col="grey")
for(i in seq_along(l)[-1])
	lines(x, invHazardDist(-log(l[i]) /exp(x) ), col='black', lty=c(2,1,2)[i])
#axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
#mtext(side=4, line=2.5, "Time")
mtext(side=3, at = log(-log(l)/hazardDistNrm(par("usr")[4])), text=paste(100*l, "%", sep=""), cex=.66)
par(xpd=NA)
mtext(side=3, "Survival",line=2, cex=.66)
mtext(side=2, line=2.5, "Time",cex=.66)

par(mar=c(4,.5,1,0.5), xpd=FALSE)
plot(h, seq_along(h), pch=NA, xlab="log hazard", yaxs="i", yaxt="n", ylab='', xlim=range(h))
abline(v=seq(-2,4,2), col="grey")
segments(h, seq_along(h),0, seq_along(h))
mtext(side=3, "Total log hazard", cex=.66, line=.5)
dev.off()



set.seed(42)
risk <- as.matrix(dataFrame[whichRFXOs]) %*% coxRFXFitOs$coefficients
risk <- risk - mean(risk)
parBoot <- mclapply(1:100, function(i) {
			s <- SimSurvNonp(risk, os)
			c <- try(CoxRFX(dataFrame[whichRFXOs], s, groups=groups[whichRFXOs], sigma0=0.1, nu=0))
			if(is(c, "try-error"))
				return(s)
			c$X <- NULL # set X to zero to save mem
			return(c)
		}, mc.cores=10)


singularities <- sapply(strsplit("67 68 72 76 78 81 86 88 89 308 310 312 313 314 315 320 322 324 325 326 327 331 332 333 334 340 341 344 346 347 350 355 357 365"," "), as.numeric)[,1]

singularities <- sapply(strsplit("228 229 230 231 232 233 234 235 236 392 393"," "), as.numeric)[,1]

#### Cumulative risk genetic
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
plot(cumsum(colMeans(X)[o]), cumVar, xlim=c(0,4), ylim=c(0,.3), xlab="Mean number of genetic drivers", ylab="Predicted variance")
a1 <- coxRFXFitOs$sigma2["Genetics"] + coxRFXFitOs$mu["Genetics"]^2
a2 <- var(coef(coxRFXFitOs)[coxRFXFitOs$groups=="Genetics"])+ coxRFXFitOs$mu["Genetics"]^2
abline(0, a1, lty=2)
abline(0, a2, lty=3)
points(3.7, 3.7*a2, col='red')
points(3.7, 3.7*a1, col='red', pch=19)
points(cumsum(colMeans(X)[o]), cumVarInt, pch=19)


f <- function(fit){M <- matrix(unlist(sapply(split(fit$coefficients, fit$groups), function(x) c(rep(0, ncol(fit$X)),x)))[-(1:ncol(fit$X))], ncol=nlevels(fit$groups))
	as.matrix(fit$X)[,order(fit$groups)] %*% M
} 


power <- function(x,y){x*y / (x/(1-x)*y/(1-y)+x/(1-x)+y/(1-y)+1)}

t <- clinicalData$Recurrence_date - clinicalData$CR_date
d <- clinicalData$Date_LF - clinicalData$CR_date
t[is.na(t)] <- d[is.na(t)]
cir <- Surv(time=as.numeric(t), event=!is.na(clinicalData$Recurrence_date)+0)
lrm <- Surv(time=as.numeric(d), event=!is.na(clinicalData$Recurrence_date) & clinicalData$Status)
nrm <- Surv(time=as.numeric(pmin(d,t)), event=is.na(clinicalData$Recurrence_date) & clinicalData$Status)


plot(survfit(cir ~ clinicalData$Family_donnor,subset=clinicalData$M_Risk %in% c("Inter-1","Inter-2")),  col=1:2, log='', mark.time=FALSE, fun=function(y) 1-y, ylim=c(0,1), xlab="Time", ylab="CIR", conf.int=TRUE )
plot(survfit(lrm ~ clinicalData$Family_donnor,subset=clinicalData$M_Risk %in% c("Inter-1","Inter-2")),  col=1:2, log='', mark.time=FALSE, fun=function(y) 1-y, ylim=c(0,1), xlab="Time", ylab="LRM", conf.int=TRUE )
plot(survfit(nrm ~ clinicalData$Family_donnor,subset=clinicalData$M_Risk %in% c("Inter-1","Inter-2")),  col=1:2, log='', mark.time=FALSE, fun=function(y) 1-y, ylim=c(0,1), xlab="Time", ylab="NRM" )

plot(survfit(cir ~ c,subset=clinicalData$M_Risk %in% c("Inter-1","Inter-2")),  log='', mark.time=FALSE, fun=function(y) 1-y, ylim=c(0,1), xlab="Time", ylab="CIR", conf.int=FALSE)
lines(survfit(nrm ~ clinicalData$Family_donnor),  col=2, mark.time=FALSE, fun=function(y) 1-y)

plot(survfit(nrm ~ c,subset=clinicalData$M_Risk %in% c("Inter-1","Inter-2")),  log='', mark.time=FALSE, fun=function(y) 1-y, ylim=c(0,1), xlab="Time", ylab="NRM", conf.int=FALSE, col=set1)

c <- cut(coxRFXFitOs$linear.predictors, seq(-3,4))

plot(survfit(cir ~ 1, data=dataList$Cytogenetics), mark.time=FALSE, ylim=c(0,1), xlab="Time", ylab="Fraction", conf.int=TRUE )
lines(survfit(nrm ~ 1), col='red', mark=NA)

d <- cbind(dataFrame[mainIdx &! grepl("TPL",names(dataFrame))], TPL=clinicalData$Family_donnor)
fitLrm <- CoxCPSSInteractions(d[!is.na(lrm) &!is.na(d$TPL),], lrm[!is.na(lrm)&!is.na(d$TPL)])
fitNrm <- CoxCPSSInteractions(d[!is.na(nrm)&!is.na(d$TPL),], nrm[!is.na(nrm)&!is.na(d$TPL)])

coxLrm <- coxph(lrm ~ ., data=dataFrame[names(which(fitLrm$Pi > .8))])
coxNrm <- coxph(nrm ~ ., data=dataFrame[c(names(which(fitNrm$Pi > .8)),"Family_donor")])

H0Nrm <- basehaz(coxph(nrm~1), centered=TRUE)
l <- range(coxNrm$linear.predictor)
plot(H0Nrm$time, exp(-H0Nrm$hazard), type='l', ylim=c(0,1))
for(i in seq_along(l))
	lines(H0Nrm$time, exp(-H0Nrm$hazard * exp(l[i])), col='black', lty=2)

H0Lrm <- basehaz(coxph(lrm~1), centered=TRUE)
l <- range(coxNrm$linear.predictor)
lines(H0Lrm$time, exp(-H0Lrm$hazard), type='l', ylim=c(0,1), col='red')
for(i in seq_along(l))
	lines(H0Lrm$time, exp(-H0Lrm$hazard * exp(l[i])), lty=2, col='red')

par(mfrow=c(5,5))
for(i in order(coxLrm$linear.predictor, decreasing = FALSE)[1:25]){
	plot(H0Nrm$time, exp(-H0Nrm$hazard * exp(coxNrm$linear.predictor[i] - dataFrame$Family_donor[i]*coxNrm$coefficients["Family_donor"])), type='l', ylim=c(0,1))
	lines(H0Nrm$time, exp(-H0Nrm$hazard * exp(coxNrm$linear.predictor[i] - (dataFrame$Family_donor[i]-1)*(coxNrm$coefficients["Family_donor"]))), lty=2)
	lines(H0Lrm$time, exp(-H0Lrm$hazard * exp(coxLrm$linear.predictor[i] - dataFrame$Family_donor[i]*coxLrm$coefficients["Family_donor"])), type='l', col='red')
	lines(H0Lrm$time, exp(-H0Lrm$hazard * exp(coxLrm$linear.predictor[i] - (dataFrame$Family_donor[i]-1)*(coxLrm$coefficients["Family_donor"]))), lty=2, col='red')
}

par(mfrow=c(1,1))
plot(coxLrm$linear.predictor, coxNrm$linear.predictor)
cor(coxLrm$linear.predictor, coxNrm$linear.predictor)


## TPL as TD
coxNrmTD <- coxph(splitSurv(nrm, clinicalData$Time_1CR_TPL)~ ., data= data.frame(TPL=c(rep(0,nrow(nrm)), rep(1, length(which(clinicalData$Time_1CR_TPL < nrm[,1]))))))
coxLrmTD <- coxph(splitSurv(lrm, clinicalData$Time_1CR_TPL)~ ., data= data.frame(TPL=c(rep(0,nrow(lrm)), rep(1, length(which(clinicalData$Time_1CR_TPL < lrm[,1]))))))

alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD")
alloTime1CR <- clinicalData$Time_1CR_TPL
alloTime1CR[!alloIdx] <- NA
coxLrmTD <- coxph(splitSurv(lrm, alloTime1CR)~ ., data= data.frame(TPL=c(rep(0,nrow(lrm)), rep(1, length(which(alloTime1CR < lrm[,1]))))))

lrmTD <- splitSurv(lrm, alloTime1CR)
i <- mainIdx & !grepl("TPL", names(dataFrame))& groups!="Trial"
coxRFXLrmTD <- CoxRFX(cbind(dataFrame[splitIndex(lrm, alloTime1CR) ,i],TPL=c(rep(0,nrow(lrm)), rep(1, length(which(alloTime1CR < lrm[,1]))))), lrmTD, groups=c(as.character(groups[i]), "Treatment"))

nrmTD <- splitSurv(nrm, alloTime1CR)
i <- mainIdx & !grepl("TPL", names(dataFrame))#& groups!="Trial"
nrmData <- cbind(dataFrame[splitIndex(nrm, alloTime1CR) ,i],TPL=c(rep(0,nrow(nrm)), rep(1, length(which(alloTime1CR < nrm[,1])))))
coxRFXNrmTD <- CoxRFX(nrmData, nrmTD, groups=c(as.character(groups[i]), "Treatment"))

i <- mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Trial"
w <- which(clinicalData$M_Risk %in% c("Inter-1","Inter-2") & !is.na(clinicalData$Family_donnor) & clinicalData$Study == "AMLHD98A")
coxRFXLrm <- CoxRFX(cbind(dataFrame[w ,i],TPL=clinicalData$Family_donnor[w]), lrm[w], groups=c(as.character(groups[i]), "Treatment"), max.iter=100)
coxRFXNrm <- CoxRFX(cbind(dataFrame[w ,i],TPL=clinicalData$Family_donnor[w]), nrm[w], groups=c(as.character(groups[i]), "Treatment"), max.iter=100)


survProb <- function(model, logHR, time){
	H0 <- basehaz(model, centered=FALSE)
	hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
	exp(-hazardDist(time) * exp(logHR))
}

intSurvProbDiff <- function(model, logHR0, logHR1, time){
	H0 <- basehaz(model, centered=FALSE)
	hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
	lambda0 <- hazardDist(0:time)
	sapply(seq_along(logHR0), function(i)
				sum(exp(-lambda0 * exp(logHR0[i])) - exp(-lambda0 * exp(logHR1)[i]))
	)
}

x <- survProb(coxRFXLrmTD, coxRFXLrmTD$X[1:nrow(lrm),] %*% coef(coxRFXLrmTD) + coef(coxRFXLrmTD)["TPL"], 1000 ) - survProb(coxRFXLrmTD, coxRFXLrmTD$X[1:nrow(lrm),] %*% coef(coxRFXLrmTD), 1000 )
y <- survProb(coxRFXNrmTD, coxRFXNrmTD$X[1:nrow(nrm),] %*% coef(coxRFXNrmTD) + coef(coxRFXNrmTD)["TPL"], 1000 ) - survProb(coxRFXNrmTD, coxRFXNrmTD$X[1:nrow(nrm),] %*% coef(coxRFXNrmTD), 1000 )

x <- survProb(coxRFXLrm, coxRFXLrmTD$X[1:nrow(lrm),] %*% coef(coxRFXLrm) + coef(coxRFXLrm)["TPL"], 1000 ) - survProb(coxRFXLrm, coxRFXLrmTD$X[1:nrow(lrm),] %*% coef(coxRFXLrm), 1000 )
y <- survProb(coxRFXNrm, coxRFXNrmTD$X[1:nrow(nrm),] %*% coef(coxRFXNrm) + coef(coxRFXNrm)["TPL"], 1000 ) - survProb(coxRFXNrm, coxRFXNrmTD$X[1:nrow(nrm),] %*% coef(coxRFXNrm), 1000 )


plot(x,y,xlab="dLRM", ylab="dNRM")
abline(0,-1)
table(x> -y)

# Train
x <- intSurvProbDiff(coxRFXLrm, coxRFXLrmTD$X[1:nrow(lrm),] %*% coef(coxRFXLrm) + coef(coxRFXLrm)["TPL"], coxRFXLrmTD$X[1:nrow(lrm),] %*% coef(coxRFXLrm), 5000)
y <- intSurvProbDiff(coxRFXNrm, coxRFXNrmTD$X[1:nrow(nrm),] %*% coef(coxRFXNrm) + coef(coxRFXNrmTD)["TPL"], coxRFXNrmTD$X[1:nrow(nrm),] %*% coef(coxRFXNrm), 5000)

# Predict
w <- which(clinicalData$M_Risk %in% c("Inter-1","Inter-2") & !is.na(clinicalData$Family_donnor) & clinicalData$Study == "AMLSG0704")
d <- as.matrix(cbind(dataFrame[w,i], TPL=0))
x <- intSurvProbDiff(coxRFXLrm, d %*% coef(coxRFXLrm) + coef(coxRFXLrm)["TPL"], d %*% coef(coxRFXLrm), 5000)
y <- intSurvProbDiff(coxRFXNrm, d %*% coef(coxRFXNrm) + coef(coxRFXNrmTD)["TPL"], d %*% coef(coxRFXNrm), 5000)
survfit(lrm[w]~ (x> -y) + clinicalData$TPL_o[w])
plot(survfit(lrm[w]~ (x> -y) + clinicalData$TPL_o[w]), col=1:4)


H0LrmTD <- basehaz(coxLrmTD, centered=FALSE)

plot(survfit(splitSurv(nrm, clinicalData$Time_1CR_TPL)~ c(rep(0,nrow(nrm)), rep(1, length(which(clinicalData$Time_1CR_TPL < nrm[,1]))))), col=1:2, xlab="Days after CR", ylab="NRM")
legend("bottomright", bty="n", lty=1, col=1:2, c("no TPL","TPL"))

plot(survfit(splitSurv(lrm, clinicalData$Time_1CR_TPL)~ c(rep(0,nrow(lrm)), rep(1, length(which(clinicalData$Time_1CR_TPL < lrm[,1]))))), col=1:2, xlab="Days after CR", ylab="Relapse-free fraction")
legend("topright", bty="n", lty=1, col=1:2, c("no TPL","TPL"))

plot(survfit(splitSurv(lrm, alloTime1CR)~ c(rep(0,nrow(lrm)), rep(1, length(which(alloTime1CR < lrm[,1])))) + clinicalData$M_Risk[splitIndex(lrm, alloTime1CR)]), col=set1[1:4], xlab="Days after CR", ylab="Relapse-free fraction", mark=NA, lty=rep(1:2, each=4))
legend("topright", bty="n", lty=1:2,  c("no TPL","TPL"))

par(mfrow=c(3,3))
for(h in  seq(-1,1,0.25)){
plot(H0LrmTD$time,  exp(-H0LrmTD$hazard * exp(h )), type='l', col=1, ylim=c(0,1))
lines(H0LrmTD$time,  exp(-H0LrmTD$hazard * exp(h + coxRFXLrm$coef["TPL"] )), type='l', col=2)
title(main=h)
}

par(mfrow=c(2,2))

lrmTD <- splitSurv(lrm, alloTime1CR)
plot(survfit(lrmTD  ~ c(rep(0,nrow(lrm)), rep(1, sum(alloTime1CR < lrm[,1], na.rm=TRUE)))), col=1:2, ylab="Leukaemia-related survival")
legend("topright", bty="n", lty=1, col=1:2, legend=paste(c("no TPL (", "TPL ("), sum(!is.na(lrm)) * c(1,0) + c(-1,1) * sum(alloTime1CR < lrm[,1], na.rm=TRUE), c(")",")"), sep=""))

nrmTD <- splitSurv(nrm, alloTime1CR)
plot(survfit(nrmTD  ~ c(rep(0,nrow(nrm)), rep(1, sum(alloTime1CR < nrm[,1], na.rm=TRUE)))), col=1:2, ylab="Non-relapse mortality", fun = function(x) 1-x)
legend("topright", bty="n", lty=1, col=1:2, legend=paste(c("no TPL (", "TPL ("), sum(!is.na(nrm)) * c(1,0)+ c(-1,1) * sum(alloTime1CR < nrm[,1], na.rm=TRUE), c(")",")"), sep=""))

prs <- Surv(time=as.numeric(clinicalData$Date_LF - clinicalData$Recurrence_date), event=clinicalData$Status)
#plot(survfit(prm ~ 1))
t <- clinicalData$TPL_date - clinicalData$Recurrence_date
t[t < 0] <- NA
prd <- data.frame(postRelapseTPL=c(rep(0, nrow(prs)), rep(1,sum(t < prs[,1], na.rm=TRUE))))
i <- (clinicalData$TPL_date < clinicalData$Recurrence_date) 
i[is.na(clinicalData$TPL_date) &! is.na(clinicalData$Recurrence_date)] <- 0
s <- survfit(splitSurv(prs, as.numeric(t))  ~ prd$postRelapseTPL + i[splitIndex(prs, as.numeric(t))])
plot(s, col=1:3, ylab="Post-relapse survival")
coxph(splitSurv(prs, as.numeric(t)) ~ ., data=prd)
legend("topright", bty="n", lty=1, col=1:2, legend=paste(c("no TPL (", "TPL ("), sum(!is.na(prs)) * c(1,0) + c(-1,1) * sum(alloTime1CR < prs[,1], na.rm=TRUE), c(")",")"), sep=""))


pcm <- Surv(time=as.numeric(clinicalData$Date_LF - clinicalData$CR_date), event=clinicalData$Status)
#plot(survfit(pcm ~ 1))
t <- clinicalData$TPL_date - clinicalData$CR_date
t[t < 0] <- NA
pcd <- data.frame(postRemissionTPL=c(rep(0, nrow(pcm)), rep(1,sum(t < pcm[,1], na.rm=TRUE))))
plot(survfit(splitSurv(pcm, as.numeric(t))  ~ pcd$postRemissionTPL), col=1:2, ylab="Overall survival")
coxph(splitSurv(pcm, as.numeric(t)) ~ ., data=pcd)
legend("topright", bty="n", lty=1, col=1:2, legend=paste(c("no TPL (", "TPL ("), sum(!is.na(pcm)) * c(1,0) + c(-1,1) * sum(alloTime1CR < pcm[,1], na.rm=TRUE), c(")",")"), sep=""))



preRelTpl <- clinicalData$TPL_date < clinicalData$Recurrence_date
plot(survfit(Surv(time1, time2, event) ~ preRelTpl[tdData$i], data=tdData), col=1:2)


# Achieval of CR v death

c <- as.numeric(clinicalData$CR_date - clinicalData$ERDate)
e <- c < os[,1]
e[is.na(e)] <- 0
c[is.na(c) &! is.na(os[,1])] <- os[is.na(c) &! is.na(os[,1]),1]
cr <- Surv(time=pmin(c, os[,1]), event = e)
PlotVarianceComponents(CoxRFX(dataFrame[mainIdxOs &! grepl("TPL", names(dataFrame))], cr, groups=groups[mainIdxOs &! grepl("TPL", names(dataFrame))]))

### SF AML

plot(survfit(lrm ~ U2AF1, data=dataFrame), col=1:2)

survdiff(lrm ~ U2AF1, data=dataFrame)
survdiff(os ~ U2AF1, data=dataFrame)

whichSplice <- c("U2AF1","SFRS2","SF3B1","ZRSR2","U2AF2","SF3A1", "SF1")
colSums(dataFrame[whichSplice])
spliceSum <- rowSums(dataList$Genetics[,whichSplice])
table(spliceSum)
anySplice <- spliceSum > 0
cor(dataFrame[selectedIntOs], anySplice)
boxplot(dataFrame$AOD_10 ~ spliceSum, ylab ="Age/10", xlab="# SF mutations")
coxph(os ~ .,data=dataFrame[whichSplice[-6]])
coxph(os ~ anySplice)
plot(survfit(os ~ anySplice))

p <- colMeans(genes[,whichSplice])
p0 <- prod(1-p)
p1 <- sum(p * p0 %o% (1/(1-p)))
p2 <- (p * p0 %o% (1/(1-p)))

fit <- CoxCPSSInteractions(cbind(dataFrame[mainIdx &! grepl("TPL",names(dataFrame))], splice=anySplice>0), os)

boxplot(dataFrame$BM_Blasts_100 ~ spliceSum, ylab ="Age/10", xlab="# SF mutations")


spliceCases <- unique(as.character(mutationData$SAMPLE_NAME[mutationData$GENE %in% whichSplice])) 

pdf(width=8,height=8)
t <- get(l[1])
p <- get(l[2])
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
s <- sample(nrow(p),nStars^2) #1:(nStars^2)
h <- hclust(dist(p[s,]))
x <- p - rep(colMeans(p), each=nrow(p))
x <- x/(2*sd(x)) + 1
stars(x[s,][h$order,]/2, scale=FALSE, locations=locations, key.loc=c(0,-3), col.lines=rep(1,(nStars^2)), col.stars = set1[clinicalData$M_Risk[s][h$order]])
symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=1)
dev.off()

library(VennDiagram)
grid.newpage()
pushViewport(viewport(w = .9, h = .9))
grid.draw(venn.diagram(list(CPSS=selectedIntOs, BIC=names(coef(coxBICOs)), AIC=names(coef(coxAICOs))), filename=NULL, lty=1, 
				col=set1[1:3], fill=set1[1:3], alpha=0.3, euler.d=TRUE, fontfamily="Helvetica", cat.fontfamily="Helvetica", euler.diagram=TRUE))



n <- c(5, 10, 50, 100, 500, 1000)
s <- 10^(seq(-2,0.5,0.1))

f <- function(p,s,n=1000){
	X <- matrix(rnorm(p*n, sd = sqrt(1/p)), ncol=p)
	b <- rnorm(p, sd=sqrt(s))
	r <- X %*% b
	u <- CoxHD:::SimSurv(r)
	CoxRFX(X,u, nu=0,  sigma0 = 0.1)$sigma2
}

x <- mclapply(s, function(x){sapply(n, function(y) f(y,x))}, mc.cores=13)

for(f in 1:10){
	print(f)
}


set.seed(42)
cvIdx <- sample(1:5, nrow(os), replace=TRUE)
ggOrNotGg <- mclapply(1:10, function(x){
			i <- rep(1:5, each=2)[x]
			if(x%%2)
				fit <-  CoxRFX(dataFrame[cvIdx!=i,intersect(whichRFXOs, which(groups!="GeneGene"))], os[cvIdx!=i], groups=groups[intersect(whichRFXOs, which(groups!="GeneGene"))])
			else
				fit <-  CoxRFX(dataFrame[cvIdx!=i,whichRFXOs], os[cvIdx!=i], groups=groups[whichRFXOs])
			survConcordance(os[cvIdx==i] ~ PredictRiskMissing(fit, dataFrame[cvIdx==i, names(coef(fit))])[,1])$concordance
		}, mc.cores=10
)
matrix(unlist(ggOrNotGg), ncol=2, byrow=TRUE, dimnames=list(NULL, c("w/o int","w int")))

i <- intersect(whichRFXOs, which(groups %in% mainGroups))
fit <- CoxRFX(dataFrame[i], os, groups[i])

survConcordance(tcgaSurvival ~ PredictRiskMissing(coxRFXFitOs, tcgaData[names(coef(coxRFXFitOs))])[,1])

whichMainOs <- intersect(whichRFXOs, which(groups %in% mainGroups))
coxRFXFitOsMain <- CoxRFX(dataFrame[whichMainOs], os, groups[whichMainOs])
survConcordance(tcgaSurvival ~ PredictRiskMissing(coxRFXFitOsMain, tcgaData[names(coef(coxRFXFitOsMain))])[,1])


whichMainEfs <- intersect(whichRFXEfs, which(groups %in% mainGroups))
coxRFXFitEfsMain <- CoxRFX(dataFrameEfsTD[whichMainEfs], efsTD, groups[whichMainEfs])

PlotVarianceComponents(coxRFXFitOsMain, col=col1)
PlotVarianceComponents(coxRFXFitEfsMain, col=col1)

whichRFXOsHF <- which((colSums(dataFrame)>=16 | mainIdx) & osIdx) # ie, > 1%
coxRFXFitOsHF <- CoxRFX(dataFrame[whichRFXOsHF], os, groups[whichRFXOsHF])
survConcordance(tcgaSurvival ~ PredictRiskMissing(coxRFXFitOsHF, tcgaData[names(coef(coxRFXFitOsHF))])[,1])
PlotVarianceComponents(coxRFXFitOsHF, col=col1)


layout(matrix(1:6, nrow=2), width=c(6,6,2),height=c(3,12))
par(bty="n", mar=c(0,4,4,2), mgp=c(2,.5,0), tcl=-.25)
v <- VarianceComponents(coxRFXFitOsMain, type="rowSums")
w <- TRUE#coxRFXFitOsTD$groups %in% mainGroups
f <- factor(as.character(coxRFXFitOsMain$groups[w]), levels=levels(coxRFXFitOsMain$groups[w])[order(abs(v), decreasing=TRUE)])
o <- order(f,coxRFXFitOsMain$coefficients[w])
c <- coef(coxRFXFitOsMain)[w][o]
plot(c, type='p',  ylab="log hazard/variable", xaxt='n', col=col1[as.character(f)[o]], xaxs="i", lwd=2, pch=NA)
y <- seq(-2,2,l=100)
par(xpd=FALSE)
abline(h=seq(-.5,.5,.25), col="grey", lty=c(1,3))
colRamp <- sapply(col1[levels(f)], function(x) c(colorRampPalette(c("white",x))(11)[-1]))
#for(l in levels(f)){
#	x <- dnorm(y,coxRFXFitOs$mu[l], sqrt(coxRFXFitOs$sigma2[l]))
#	x <- x/max(x)*15
#	#polygon(x + (cumsum(table(f)) - table(f))[l], y, col=colRamp[1,l], lty=0)
#	#lines(x + (cumsum(table(f)) - table(f))[l], y, col=colRamp[1,l], lwd=1)
#	#segments((cumsum(table(f)) - table(f))[l],coxRFXFitOs$mu[l],(cumsum(table(f)) - table(f))[l]+max(x),coxRFXFitOs$mu[l], col=colRamp[8,l])
#	rect((cumsum(table(f)) - table(f))[l], coxRFXFitOs$mu[l] - 2*sqrt(coxRFXFitOs$sigma2[l]) , cumsum(table(f))[l] ,coxRFXFitOs$mu[l]+ 2*sqrt(coxRFXFitOs$sigma2[l]), col=colRamp[1,l], lty=0)
#	rect((cumsum(table(f)) - table(f))[l], coxRFXFitOs$mu[l] - 1*sqrt(coxRFXFitOs$sigma2[l]) , cumsum(table(f))[l] ,coxRFXFitOs$mu[l]+ 1*sqrt(coxRFXFitOs$sigma2[l]), col=colRamp[4,l], lty=0)	
#	segments((cumsum(table(f)) - table(f))[l], coxRFXFitOs$mu[l]  , cumsum(table(f))[l] ,coxRFXFitOs$mu[l], col=colRamp[8,l])	
#}
par(xpd=NA)
segments(seq_along(c), coxRFXFitOsMain$mu[as.character(f)[o]], seq_along(c), c, col=col1[as.character(f)[o]])
#points(c, col=col1[as.character(f)[o]], pch=16, cex=.66)
#mtext(side=3, levels(f), at= table(f)/2+ cumsum(c(0,table(f)[-nlevels(f)])), cex=.66)
rotatedLabel(table(f)/2+ cumsum(c(0,table(f)[-nlevels(f)])), y = rep(par("usr")[4]*.8, nlevels(f)), pos=3, labels=levels(f))
par(bty="L", mar=c(4,4,1,2))
X <- as.matrix(dataFrame[whichRFXOs][w][o])
p <- order(X %*% c)
X <- X[p,]
x <- t(X)/apply(X,2,max)
#h <- hclust(dist(t(x)))
image(x=1:nrow(x),y=1:ncol(x),x*.9 + as.numeric(f[o])-1 + 1e-5, useRaster=TRUE, 
		col=colRamp, 
		breaks=seq(0,nlevels(f), 0.1),  ylab="Patients",xlab="Variable", xlim=c(0,nrow(x)), ylim=c(0,ncol(x)))

par(mar=c(3,.5,1,0.5), xpd=NA)
PlotVarianceComponents(coxRFXFitOsMain, col = col1)
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
	segments(h + i*xo, seq_along(h), + i*xo, seq_along(h), col=col1[l])
	#mtext(side=3, l, cex=.66, line=.5)
	axis(side=1, at=-1:1 + xo*i, labels=-1:1)
	x <- seq(-1.5,1.5,l=101)
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



p <- 10^(-abs(interactions)) 
w <- which(matrix(p.adjust(p, "BH")<0.1, ncol=ncol(p)), arr.ind=TRUE)
n <- paste(rownames(interactions)[w[,1]],colnames(interactions)[w[,2]], sep=":")
n[!n%in% colnames(dataFrame)] <- sapply(strsplit(n[!n%in% colnames(dataFrame)], ":"), function(x) paste(rev(x), collapse=":"))
n <- n[n %in% colnames(dataFrame)]
#apply(w, 1, function(x) sum(dataFrame[rownames(interactions)[x[1]]]*dataFrame[rownames(interactions)[x[2]]]) )

m <- which(colnames(dataFrame) %in% n)

fit <- CoxRFX(dataFrame[c(whichMainOs,m)], os, groups[c(whichMainOs, m)])

## variance of C with errors
sd(sapply(1:100, function(i){
			r <- tcgaRiskRFXOs[,1] + rnorm(nrow( tcgaRiskRFXOs)) *  sqrt(tcgaRiskRFXOs[,2])
	survConcordance(tcgaSurvival ~ r)$concordance
}))

## Post relapse survival
prs <- Surv(time=as.numeric(clinicalData$Date_LF - clinicalData$Recurrence_date), event=clinicalData$Status)
t <- clinicalData$TPL_date - clinicalData$Recurrence_date
t[t < 0] <- NA
t[! clinicalData$TPL_type %in% c("ALLO","FREMD")] <- NA
index <- mainIdx &! grepl("TPL", names(dataFrame))
prsTD <- splitSurv(prs, t)
prsSplit <- splitIndex(prs, t)
prsTPL <- c(rep(0, nrow(dataFrame)), rep(1, sum(t < prs[,1], na.rm=TRUE)))
prsData <- cbind(dataFrame[prsSplit ,index],TPL=duplicated(prsSplit)+0)
plot(survfit(prsTD ~ prsTPL + clinicalData$M_Risk[prsSplit]), col=set1[1:4], lty=rep(c(1,3), each=4))
coxRFXPrsTD <-  CoxRFX(prsData, prsTD, groups=c(as.character(groups[index]), "Treatment"), sigma0 = 0.1, nu=0.1)
coxCPSSPrs <-  CoxCPSS(prsData[which(prs[,1]>0 & !is.na(prs)),], prs[prs[,1]>0 &! is.na(prs)], control="BH")


## rfs
rfs <- Surv(clinicalData$rfs, clinicalData$rfsstat)
t <- clinicalData$TPL_date - clinicalData$ERDate
t[t < 0] <- NA
rfsTD <- splitSurv(rfs, t)
rfsSplit <- splitIndex(rfs, t)
rfsTPL <- c(rep(0, nrow(dataFrame)), rep(1, sum(t < rfs[,1], na.rm=TRUE)))
rfsData <- cbind(dataFrame[splitIndex(rfs, t) ,index],TPL=rfsTPL)
coxRFXRfsTD <-  CoxRFX(rfsData, rfsTD, groups=c(as.character(groups[index]), "Treatment"))


c <- coef(coxRFXRfsTD)
plot(c, coef(coxRFXPrsTD) [names(c)], xlab="logHR RFS", ylab="logHR PRS")
text(c, coef(coxRFXPrsTD)[names(c)], names(c))
cor(c, coef(coxRFXPrsTD) [names(c)])

survConcordance(prsTD ~ PredictRiskMissing(coxRFXRfsTD, prsData)[,1])
survConcordance(rfsTD ~ PredictRiskMissing(coxRFXPrsTD, rfsData)[,1])

d <- as.data.frame(rbind(rfsTD, prsTD))
colnames(d) <- c("time","time2", "event")
s <- do.call("Surv",d)

d <- cbind(rbind(rfsData, prsData), REC=c(rep(0,nrow(rfsData)), rep(1,nrow(prsData))))
fit <- CoxRFX(d, s, groups=c(as.character(groups[index]), "Treatment","Nuisance"))

s2 <- s
s2[which(!is.na(clinicalData$Recurrence_date)[rfsSplit]),3] <- 0
fit2 <- CoxRFX(d, s2, groups=c(as.character(groups[index]), "Treatment","Nuisance"))

## CIR
c <- clinicalData$Recurrence_date - clinicalData$CR_date
d <- clinicalData$Date_LF - clinicalData$CR_date
c[is.na(c)] <- d[is.na(c)]
t <- clinicalData$TPL_date - clinicalData$CR_date
t[t < 0] <- NA
cirData <- makeTimeDependent(dataFrame[index], timeTpl=t, timeSurv=c, event=!is.na(clinicalData$Recurrence_date)+0)
coxph(Surv(time1, time2, event) ~ transplant, cirData)
plot(survfit(Surv(time1, time2, event) ~ transplant, data=cirData), col=1:2, conf.int=TRUE)
coxRFXCirTD <-  CoxRFX(cirData[c(1:101,105)], Surv(cirData$time1, cirData$time2, cirData$event), groups=c(as.character(groups[index]), "Treatment"))
summary(coxRFXCirTD)

plot(survfit(Surv(time1, time2, event) ~ clinicalData$M_Risk[cirData$index], data=cirData), col=set1[1:4])
riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)
r <- coxRFXCirTD$X %*% coef(coxRFXCirTD) - cirData$transplant * coef(coxRFXCirTD)["transplant"]
g <- "Inter-2"
w <- which(clinicalData$M_Risk[cirData$index]==g)
q <- cut(r[w], quantile(r[w], seq(0,1,.33)))

plot(survfit(Surv(time1, time2, event) ~ clinicalData$M_Risk[cirData$index], data=cirData), col=riskCol, ylab="Recurrence-free fraction", xlab="Time after CR", main="Molecular risk groups, all cases")
legend("topright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)
plot(survfit(Surv(time1, time2, event) ~ q + transplant, data=cirData[w,]), col=riskCol[g], lty=c(1,3), ylab="Recurrence-free fraction", main=paste(g, "terciles"), xlab="Time after CR")
legend("topright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[g])

plot(survfit(prsTD ~ clinicalData$M_Risk[prsSplit]), col=set1[1:4], ylab="Post-relapse survival", xlab="Time after relapse", main="Molecular risk groups, all cases")
plot(survfit(prs[w[w<=1540]] ~ q[w<=1540]), col=set1[2], lty=1:3, ylab="Post-relapse survival", xlab="Time after relapse", main="Favourable; terciles")

## NRM
t <- clinicalData$Time_1CR_TPL
t[t > (clinicalData$Recurrence_date - clinicalData$CR_date)] <- NA
s <- clinicalData$Date_LF - clinicalData$CR_date
s[s<0] <- NA
nrmData <- makeTimeDependent(dataFrame[index], timeTpl=t, timeSurv=s, event=clinicalData$Status &  is.na(clinicalData$Recurrence_date))
plot(survfit(Surv(time1, time2, event) ~ transplant, data=nrmData), col=set1[1:2], ylab="Non-relapse mortality", xlab="Time after CR", fun = function(x) 1- x)
legend("topleft", lty=1, bty='n', c("no TPL", "TPL"), col=set1[1:2])

preRecAlloTime <- clinicalData$TPL_date
preRecAlloTime[preRecAlloTime > clinicalData$Recurrence_date ] <- NA
s <- splitSurv(os, preRecAlloTime - clinicalData$ERDate)
plot(survfit(s ~ c(rep(0, nrow(os)), rep(1, nrow(s)-nrow(os)))), col=1:2)


pdf("survival.pdf", 3,3, pointsize=8)
par(mar=c(3,3,2,1), bty="n", mgp=c(2,.5,0))
riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)
r <- coxRFXCirTD$X %*% coef(coxRFXCirTD) - cirData$transplant * coef(coxRFXCirTD)["transplant"]

f <- function(x) 1-x
plot(survfit(Surv(time1, time2, event) ~ clinicalData$M_Risk[cirData$index], data=cirData), col=riskCol, ylab="CIR", xlab="Time after CR", main="Molecular risk groups, all cases", fun=f , ylim=c(0,1))
legend("bottomright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)
for(g in levels(clinicalData$M_Risk)){
	w <- which(clinicalData$M_Risk[cirData$index]==g)
	q <- cut(r[w], quantile(r[w], seq(0,1,.33)))
	
	plot(survfit(Surv(time1, time2, event) ~ q + transplant, data=cirData[w,]), col=riskCol[g], lty=c(1,3), ylab="CIR", main=paste(g, "terciles"), xlab="Time after CR", fun=f, ylim=c(0,1))
	legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[g])
}
dev.off()

d <- cbind(dataFrameOsTD[mainIdxOs], TPL = dataFrameOsTD$TPL_os, REC=(!is.na(clinicalData$Recurrence_date))[tplSplitOs])
fit <- CoxRFX(d, osTD, groups=c(as.character(groups[mainIdxOs]), "Treatment","Nuisance"))

preRecAllo <- (clinicalData$TPL_date < clinicalData$Recurrence_date | (is.na(clinicalData$Recurrence_date)  & !is.na(clinicalData$TPL_date))) & clinicalData$TPL_type %in% c("ALLO","FREMD")
preRecAllo[is.na(preRecAllo)] <- 0
postRecAllo <- (clinicalData$TPL_date > clinicalData$Recurrence_date) & clinicalData$TPL_type %in% c("ALLO","FREMD")
postRecAllo[is.na(postRecAllo)] <- 0

coxph(Surv(time1, time2, event) ~ CR + RD + REC + I(TPL*preRecAllo[tdData$i]) + I(TPL*postRecAllo[tdData$i]) , data=tdData[,-8], subset=clinicalData$M_Risk[tdData$i]=="Adverse")

s <- I(tdData$TPL*preRecAllo[tdData$i])
t <- I(tdData$TPL*postRecAllo[tdData$i])
plot(survfit(Surv(time1, time2, event) ~ REC + s + t, data=tdData), col=set1[1:5], conf.int=TRUE)


##
w <- mainIdx & !grepl("TPL", names(dataFrame))
cr <- efs
cr[!clinicalData$ereignart %in% c("RD_after_Ind1","RD_or_PR_after_Ind2"),2] <- 0
coxCPSSCr <- CoxCPSS(cbind(dataFrame[w & !names(dataFrame) %in% whichSplice], SF=anySplice)[!is.na(cr),], cr[!is.na(cr)], control="BH") 
coxRFXCr <- CoxRFX(dataFrame[w], cr, groups[w]) 


## Competing risk comparison
library(cmprsk)
r <- clinicalData$Recurrence_date - clinicalData$CR_date
d <- clinicalData$Date_LF - clinicalData$CR_date
s <- pmin(r,d)
s[is.na(r)] <- d[is.na(r)]
e <- clinicalData$Status
e[r < d] <- 2
e <- factor(e, labels=c("alive","dead","recurrence"))
#crr(ftime = as.numeric(pmin(d, t)), fstatus=s, cov1 = dataFrame[names(which(fitLrm$Pi > 0.8))])
compRisk <- cuminc(ftime = as.numeric(s), fstatus=e, cencode="alive")
plot(compRisk, xaxs="i", yaxs="i")

f <- survfit(cir ~ 1)
lines(f, lty=2, col=2, fun=function(x) 1-x, mark=NA, conf.int=FALSE)
g <- survfit(nrm ~ 1)
lines(g, lty=1, col=2, fun=function(x) 1-x, mark=NA, conf.int=FALSE)

cf <- cumsum(c(1,diff(f$surv)) * splinefun(g$time, g$surv, method="monoH.FC")(g$time))
nrsCR <- cumsum(c(1,diff(g$surv)) * splinefun(f$time, f$surv, method="monoH.FC")(f$time))
#for(i in 1){
#	cg <-  1-(1-g$surv)* splinefun(f$time, cf, method="monoH.FC")(g$time)
#	cf <-  1-(1-f$surv)* splinefun(g$time, cg, method="monoH.FC")(f$time)
#}

lines(f$time, 1-cf, col="blue", lty=2)
lines(g$time, 1-nrsCR, col="blue")


## data for web

alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD")
alloTime1CR <- clinicalData$Time_1CR_TPL + .5
alloTime1CR[!alloIdx] <- NA
alloTimeRel <- clinicalData$TPL_date - clinicalData$Recurrence_date + .5
alloTimeRel[!alloIdx] <- NA
alloTimeRel[alloTimeRel <0] <- NA


w <- mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Nuisance"

t <- clinicalData$Recurrence_date
t[is.na(t)] <- as.Date(1e6, origin="2000-01-01")
cirData <- makeTimeDependent(dataFrame[w], timeTpl=alloTime1CR, timeSurv=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), event=!is.na(clinicalData$Recurrence_date)+0)
nrmData <- makeTimeDependent(dataFrame[w], timeTpl=alloTime1CR, timeSurv=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), event=is.na(clinicalData$Recurrence_date) & clinicalData$Status)
#prsData <- makeTimeDependent(dataFrame[w], timeTpl=alloTimeRel, timeSurv=as.numeric(clinicalData$Date_LF - clinicalData$Recurrence_date)+1, event=clinicalData$Status)
prsData <- makeTimeDependent(dataFrame[w], timeTpl=alloTime1CR, timeSurv=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), time0Surv = as.numeric(clinicalData$Recurrence_date- clinicalData$CR_date), event=clinicalData$Status)
cirData$transplant1CR <- cirData$transplant
cirData$transplant <- NULL
cirData$transplantRel <- 0
nrmData$transplant1CR <- nrmData$transplant
nrmData$transplant <- NULL
nrmData$transplantRel <- 0
prsData$transplant1CR <- nrmData$transplant1CR[1:nrow(dataFrame)][prsData$index]
prsData$transplantRel <- prsData$transplant
prsData$transplant <- NULL

g <- c(as.character(groups[w]), "Treatment","Treatment")
names(g) <- c(names(dataFrame)[w],"transplant1CR","transplantRel")
coxRFXNrmTD <- CoxRFX(nrmData[names(g)], Surv(nrmData$time1, nrmData$time2, nrmData$event), groups=g)
coxRFXNrmTD$coefficients["transplantRel"] <- 0
coxRFXPrsTD <-  CoxRFX(prsData[names(g)], Surv(prsData$time1, prsData$time2, prsData$event), groups=g, sigma0 = 0.1, nu=0.1)
coxRFXCirTD <-  CoxRFX(cirData[names(g)], Surv(cirData$time1, cirData$time2, cirData$event), groups=g)
coxRFXCirTD$coefficients["transplantRel"] <- 0

osCR <- Surv(time=as.numeric(clinicalData$Date_LF - clinicalData$CR_date), event=clinicalData$Status)

c <- CoxRFX(nrmData[c("AOD_10","transplant1CR","gender")], Surv(nrmData$time1, nrmData$time2, nrmData$event))
coxRFXNrmTD$coefficients[] <-0
coxRFXNrmTD$coefficients[names(coef(c))] <- coef(c)
coxRFXNrmTD$var2[] <- 0
colnames(coxRFXNrmTD$var2) <- rownames(coxRFXNrmTD$var2) <- names(coxRFXNrmTD$coefficients)
coxRFXNrmTD$var2[names(coef(c)),names(coef(c))] <- c$var2


save(coxRFXCirTD, coxRFXNrmTD, coxRFXPrsTD, nrmData, cirData, prsData, osCR, g,crGroups, file="../../../Projects/sandbox/relapse/predict.RData")


pdf("survival.pdf", 3,3, pointsize=8)
par(mar=c(3,3,2,1), bty="n", mgp=c(2,.5,0))
riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)
r <- coxRFXCirTD$X %*% coef(coxRFXCirTD) - coxRFXCirTD$X[,"transplant1CR"] * coef(coxRFXCirTD)["transplant1CR"]
f <- function(x) 1-x
plot(survfit(Surv(time1, time2, event) ~ clinicalData$M_Risk[cirData$index], data=cirData), col=riskCol, ylab="CIR", xlab="Time after CR", main="Molecular risk groups, all cases", fun=f , ylim=c(0,1))
legend("bottomright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)
for(g in levels(clinicalData$M_Risk)){
	w <- which(clinicalData$M_Risk[cirData$index]==g)
	q <- cut(r[w], quantile(r[w], seq(0,1,.33)))
	
	plot(survfit(Surv(time1, time2, event) ~ q + transplant1CR, data=cirData[w,]), col=riskCol[g], lty=c(1,3), ylab="CIR", main=paste(g, "terciles"), xlab="Time after CR", fun=f, ylim=c(0,1))
	legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[g])
}
plot(survfit(Surv(time1, time2, event) ~ transplant1CR, data=nrmData), col="black", lty=c(1,3), ylab="NRM", xlab="Time after CR", fun=f, ylim=c(0,1))
legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col="black")

plot(survfit(Surv(time1, time2, event)   ~ clinicalData$M_Risk[prsData$index] + transplantRel, data=prsData), col=riskCol, ylab="Survival", xlab="Time after relapse", main="Molecular risk groups, all cases", ylim=c(0,1), lty=rep(c(3,1), each=4))
legend("topright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)

dev.off()



f <- "os ~ FLT3_ITD + FLT3_TKD"
s <- survfit(as.formula(f), data=dataList$Genetics)
c <- coxph(as.formula(sub("\\+","\\*",f)), data=dataList$Genetics)
summary(c)
p <- 2*pnorm(coef(c)[3],sd=sqrt(diag(c$var)[3]), lower.tail=FALSE)
plot(s, col=set1, mark=NA)
legend("topright", bty="", rownames(summary(s)$table), col=set1, lty=1)
title(main=paste("P =",signif(p,2)), font=1)



##MDS
VAF <- matrix(0,nrow = length(IDs), ncol = length(levels(mds_gen$Gene)))
rownames(VAF) <- IDs
colnames(VAF) <- levels(mds_gen$Gene)
VAF <- VAF[rownames(VAF)!="", colnames(VAF)!=""]
for(index in seq_along(mds_gen$Gene)){
	#if(mds_gen$Decision[i] %in% c("ONCOGENIC", "POSSIBLE ONCOGENIC"))
	if(mds_gen$SAMPLE.NAME[index] %in% IDs)
		if(mds_gen$SAMPLE.NAME[index] != "" & mds_gen$Decision[index] %in% c("ONCOGENIC")){
			VAF[as.character(mds_gen$SAMPLE.NAME[index]), as.character(mds_gen$Gene[index])] <- max(mds_gen$X._MUT_IN_TUM[index],VAF[as.character(mds_gen$SAMPLE.NAME[index]), as.character(mds_gen$Gene[index])]) 
			#if(!is.na(match(mds_gen$CHR[i], colnames(copy_numbers))))
			#	clonality[as.character(mds_gen$SAMPLE.NAME[i]), as.character(mds_gen$Gene[i])] <- clonality[as.character(mds_gen$SAMPLE.NAME[i]), as.character(mds_gen$Gene[i])] * copy_numbers[match(mds_gen$SAMPLE.NAME[i], mds_clin$PDID),match(mds_gen$CHR[i], colnames(copy_numbers))]/2
		}
}
VAF = VAF[rownames(VAF) != "PD6785c",]
VAF = VAF[,colSums(VAF, na.rm=TRUE)>0]
VAF = VAF[,order(colSums(VAF>0, na.rm=TRUE), decreasing = TRUE)]


#Jackknife
set.seed(42)
k <- lapply(1:100, function(b){
			w <- sample(1:nrow(dataFrame)%%2 +1 )!=1
			m <- CoxRFX(dataFrame[w,mainIdxOs], os[w], groups=groups[mainIdxOs])
			p <- as.matrix(dataFrame[!w,mainIdxOs]) %*% coef(m)
			list(coef=coef(m), conc=survConcordance(os[!w] ~ p)$concordance)
})
c <- sapply(k, `[[`, 1)
d <- sapply(k, `[[`,2) 
coxRFXKnife <- coxRFXFitOsMain
coxRFXKnife$coefficients <- rowMeans(c, na.rm=TRUE)
coxRFXKnife$var <- cov(t(c), use="c")

PlotVarianceComponents(coxRFXKnife, col=col1)
tcgaRiskRFXOs <- PredictRiskMissing(coxRFXKnife, tcgaData[which(mainIdxOs)])
survConcordance(tcgaSurvival ~ tcgaRiskRFXOs[,1])

boxplot(t(c))
points(coef(coxRFXFitOsMain), col='red',pch=19)
segments(seq_along(coef(coxRFXFitOsMain)),coef(coxRFXFitOsMain) - sqrt(diag(coxRFXFitOsMain$var2)),seq_along(coef(coxRFXFitOsMain)),coef(coxRFXFitOsMain) + sqrt(diag(coxRFXFitOsMain$var2)), col="red" )
abline(h=0)

coxCPSSKnife <- coxCPSSIntOs$coxph
set.seed(42)
h <- lapply(1:100, function(b){
			w <- sample(1:nrow(dataFrame)%%3 +1 )!=1
			f <- os[w] ~ NPM1 + RUNX1 + TP53 + CEBPA_bi + FLT3_ITD + minus7 + complex + Date_1000 + ATRA + AOD_10 + Performance_ECOG + wbc_100 + LDH_1000 + ATRA:AOD_10 + Date_1000:wbc_100
			m <- coxph(f, data=dataFrame[w,])
			p <- predict(m, newdata=dataFrame[!w,])
			list(coef=coef(m), conc=survConcordance(os[!w] ~ p)$concordance)
		})
c <- sapply(h, `[[`, "coef")
e <- sapply(h, `[[`, "conc")
boxplot(t(c))
points(coef(coxCPSSKnife), col='red',pch=19)
segments(1:15,coef(coxCPSSKnife) - sqrt(diag(coxCPSSKnife$var)),1:15,coef(coxCPSSKnife) + sqrt(diag(coxCPSSKnife$var)), col="red" )

index <- lapply(1:100, function(b){
			w <- sample(1:nrow(dataFrame)%%3 +1 )!=1
			m <- coxph(coxBICOs$formula, data=dataFrame, subset=w)
			p <- predict(m, newdata=dataFrame[!w,])
			list(coef=coef(m), conc=survConcordance(os[!w] ~ p)$concordance)
		})
c <- sapply(h, `[[`, "coef")
f <- sapply(index, `[[`, "conc")

boxplot(d,e,f)

survConcordance(tcgaSurvival ~ predict(coxBICOs, newdata=tcgaDataImputed))
survConcordance(tcgaSurvival ~ predict(coxCPSSIntOs, newdata=tcgaDataImputed))

sapply(1:100, function(i){
			w <- sample(1:nrow(tcgaSurvival), 100)
			survConcordance(tcgaSurvival[w] ~ predict(coxRFXFitOsMain, newdata=tcgaDataImputed)[w])$concordance
		})

glmPrediction <- glm$coefficients[,1:31] %*% t(tcgaDesign[,1:31])
pcPrediction <- t(glmPrediction - tcgaPca$center) %*% tcgaPca$rotation
set.seed(42)
concPcPred <- sapply(1:100, function(i){
			v <- sample(1:nrow(dataFrame) %% 5 + 1) == 1
			fit <- CoxRFX(pcPrediction[!v, 1:20], tcgaSurvival[!v], which.mu=NULL)
			p <- as.matrix(pcPrediction[v, 1:20]) %*% coef(fit)
			pcp <- survConcordance(tcgaSurvival[v]~p)$concordance[1]
			fit <- CoxRFX(tcgaPca$x[!v, 1:20], tcgaSurvival[!v], which.mu=NULL)
			p <- as.matrix(tcgaPca$x[v, 1:20]) %*% coef(fit)
			pc <- survConcordance(tcgaSurvival[v]~p)$concordance[1]
			fit <- CoxRFX(dataFrame[!v,survivalGroups %in% c("Genetics","Cytogenetics")], tcgaSurvival[!v], which.mu=NULL)
			p <- as.matrix(dataFrame[v, survivalGroups %in% c("Genetics","Cytogenetics")]) %*% coef(fit)
			gen <- survConcordance(tcgaSurvival[v]~p)$concordance[1]
			c(PC=pc, Genomics=gen, `PC (predicted)`=pcp)
		})

boxplot(t(concPcPred), notch=TRUE, ylab="C")
title(main="AML")

##MDS
pcPrediction <-  as.matrix(design[,2:18]) %*% t(glm$coefficients[,2:18])  %*% pca$rotation
set.seed(42)
mdsConcPcPred <- sapply(1:100, function(i){
			v <- sample(1:nrow(amlFreeSurvival) %% 5 + 1) == 1
			fit <- CoxRFX(pcPrediction[as.character(mdsSamples), 1:20][!v,], amlFreeSurvival[!v], which.mu=NULL)
			p <- as.matrix(pcPrediction[as.character(mdsSamples), 1:20][v,]) %*% coef(fit)
			pcp <- survConcordance(amlFreeSurvival[v]~p)$concordance[1]
			fit <- CoxRFX(Z$expression[!v, 1:20], amlFreeSurvival[!v], which.mu=NULL)
			p <- as.matrix(Z$expression[v, 1:20]) %*% coef(fit)
			pc <- survConcordance(amlFreeSurvival[v]~p)$concordance[1]
			fit <- CoxRFX(data.frame(Z$geneticsCytogenetics[!v,]), amlFreeSurvival[!v], which.mu=NULL)
			p <- as.matrix(Z$geneticsCytogenetics[v,]) %*% coef(fit)
			gen <- survConcordance(amlFreeSurvival[v]~p)$concordance[1]
			c(PC=pc, Genomics=gen, `PC (predicted)`=pcp)
		})

#+ TCGA-, fig.width=1.25, fig.height=2.5, cache=TRUE

pdf("predictedPC.pdf", width=1.25, height=2.5, pointsize=8)
par(mar=c(10,4,1,1), bty="n", mgp=c(2,.5,0))
c <- set1[c(4,2,4)]
boxplot(t(mdsConcPcPred), notch=TRUE, ylab="Concordance", names=NA, lty=1, staplewex=0, pch=16, cex=.5, xaxt="n", at=1:3, border=c, col=c(NA,NA, brewer.pal(4,"Paired")[1]))
u <- par("usr")
rotatedLabel(1:3, rep(u[3] ,3), c("Expression","Genomics","Predicted expression"))
dev.off()



f <- character(nrow(clinicalData))
for(l in levels(clinicalData$M_Risk)){
	w <- which(clinicalData$M_Risk==l)
	q <- cut(r[w], quantile(r[w], seq(0,1,.33), include.lowest=TRUE), labels=c("T1","T2","T3"))
	f[w] <- as.character(q)
}
f <- factor(f, levels=c("T1","T2","T3"))

getGT <- function(X){
	n <- colnames(X)
	apply(X,1, function(x){
				w <- x!=0
				l <- sub("^1$","",paste(signif(x[w],2)))
				paste(sub(":$","",paste(n[w], l, sep=":")), collapse=";")
			})
}

whichRFXOsGG <- which((colSums(dataFrame)>=8 | mainIdx) & osIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%
ggc <- mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,whichRFXOsGG], os[trainIdx], groups[whichRFXOsGG])
			return(c(
							RFX=survConcordance(os[!trainIdx]~as.matrix(dataFrame[!trainIdx,whichRFXOsGG]) %*% coef(coxRFXOsTrain))$concordance
					))
		}, mc.cores=20)
coxRFXOsGG<-CoxRFX(dataFrame[whichRFXOsGG], os, groups[whichRFXOsGG])
coxRFXOsTDGG<-CoxRFX(dataFrameOsTD[whichRFXOsGG], osTD, groups[whichRFXOsGG])

allIntC <- mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,whichRFXOs], os[trainIdx], groups[whichRFXOs])
			return(c(
							RFX=survConcordance(os[!trainIdx]~as.matrix(dataFrame[!trainIdx,whichRFXOs]) %*% coef(coxRFXOsTrain))$concordance
					))
		}, mc.cores=12)

confIntVar <- function(x, ci=c(0.025,0.0975)){
	v <- var(x)
	s <- sum((x-mean(x))^2)
	n <- length(x)
	s/qchisq(ci,n-1)
}


X <- cbind(dataFrame[whichRFXOsGG], `NPM1:FLT3_ITD:DNMT3A` = (rowSums(dataFrame[c('NPM1',"FLT3_ITD","DNMT3A")])==3)+0)
p <- ncol(X)
g <- c(as.character(groups[whichRFXOsGG]),"GeneGene")
gtc2 <- mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			fitGG <- CoxRFX(X[trainIdx,-p], os[trainIdx], g[-p])
			fitGT <- CoxRFX(X[trainIdx,], os[trainIdx], g)
			fitMain <- CoxRFX(dataFrame[trainIdx,mainIdxOs], os[trainIdx], groups=groups[mainIdxOs])
			
			return(c(
							GG=survConcordance(os[!trainIdx]~as.matrix(X[!trainIdx,-p]) %*% coef(fitGG))$concordance,
							GT=survConcordance(os[!trainIdx]~as.matrix(X[!trainIdx,]) %*% coef(fitGT))$concordance,
							Main=survConcordance(os[!trainIdx]~as.matrix(dataFrame[!trainIdx,mainIdxOs]) %*% coef(fitMain))$concordance
			
					))
		}, mc.cores=20)


plotSurv <- function(f, col=1:(length(terms(f))-1)^2 , ...){
	s <- survfit(f, ...)
	c <- coxph(f, ...)
	summary(c)
	#p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
	p <- pchisq(c$coefficients[3]^2/diag(c$var)[3], 1,lower.tail=FALSE)
	plot(s, col=col, mark=NA)
	legend("topright", bty="n", rownames(summary(s)$table), col=col, lty=1)
	title(main=paste("P =",signif(p,2)), font=1)
}

PredictOS <- function(coxRFXNrmTD, coxRFXCirTD, coxRFXPrsTD, data, x =365){
	getS <- function(coxRFX, data) {
		coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
		r <- PredictRiskMissing(coxRFX, matrix(data, ncol=ncol(coxRFX$Z)), var="var2")
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

w <- mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Nuisance"

osCR <- Surv(as.numeric(clinicalData$Date_LF- clinicalData$CR_date), event = clinicalData$Status)
osData <- MakeTimeDependent(dataFrame[w], timeEvent=alloTime1CR, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
osData$transplant1CR <- osData$event
osData$transplantRel <- osData$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date)  
osData$transplant1CR[osData$index %in% w] <- 0
osData$transplantRel[!osData$index %in% w] <- 0

coxRFXOsCR <- CoxRFX(osData[names(crGroups)], Surv(osData$time1, osData$time2, osData$status), groups=crGroups)

allRisk <- PredictOS(coxRFXNrmTD = coxRFXNrmTD, coxRFXPrsTD = coxRFXPrsTD, coxRFXCirTD = coxRFXCirTD, coxRFXOsCR$Z, 365)

survConcordance(Surv(osData$time1, osData$time2, osData$status) ~ -allRisk[,1])


bimp <- function(x){
	 if(any(is.na(x))) x[is.na(x)] <- rbinom(sum(is.na(x)),1, mean(x[!is.na(x)]))
	 x 
}

library(pcalg)
set.seed(42)
X <- sapply(cbind(dataList$Genetics, dataList$Cytogenetics), bimp)
fit2 <- pc(list(dm=X, adaptDF = FALSE), binCItest, labels=colnames(X), alpha=0.05)
library(Rgraphviz)
plot(fit2)

dtrbeta <- function(x, shape1,shape2, xmin = 0)
	dbeta(x,shape1,shape2, log=TRUE) - pbeta(xmin, shape1, shape2, lower.tail=FALSE, log.p=TRUE) 
p <- optim(c(0.2,8), function(y)  -sum(dtrbeta(x, y[1],y[2], xmin=1/1540)))

rtrbeta <- function(n, shape1,shape2, xmin = 0){
	x <- numeric(n)
	i <- 1
	while(i<=n){
		r <- rbeta(1,shape1,shape2)
		if(r > xmin){
			x[i] <- r
			i <- i+1
		}
	}
	return(x)
}


p1 <- c(.4, .6)
p2 <- c(.5,.5)

X <- c(50,100)

L <- apply(P,2,prod)
pi <- L/sum(L)
pi <- c(.5,.5)

for(i in 1:10000){
	P <- matrix(c(rep(pi[1]*p1[1], X[1]), rep(pi[1]*p1[2], X[2]),rep(pi[2]*p2[1], X[1]), rep(pi[2]*p2[2], X[2])), ncol=2)
	pi <- colMeans( P / rowSums(P))
}

whichRFXOsGG <- which((colSums(dataFrame)>=8 | mainIdxOs) & osIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%


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
		}, mc.cores=20)


tcgaConcordanceCV <- mclapply(1:100, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,mainIdxOs], os[trainIdx], groups=groups[mainIdxOs])
			p <- PredictRiskMissing(coxRFXOsTrain, newZ=tcgaData[mainIdxOs])
			return(c(
							RFX=survConcordance(tcgaSurvival~ p[,1])$concordance
					))
		}, mc.cores=10)


allModelsTcgaConcordance <- sapply(allModelsCV, function(x){
			sapply(x, function(y){
						p <- predictAllModels(y, newdata=tcgaDataImputed)
						survConcordance(tcgaSurvival ~ p)$concordance
					})
		})

apply(sapply(allModelsCV, function(x){
			sapply(x[1:4], function(y){
						length(attr(terms(y), "term.labels"))
					})
		}),1,quantile)


apply(sapply(allModelsCV, function(x){
			sapply(x, function(y){
						p <- predictAllModels(y, newdata=tcgaDataImputed)
						survConcordance(tcgaSurvival ~ p)$concordance
					})
		}), 1, quantile)


t <- sort(unique(tcgaSurvival[,1]))
allModelsTcgaAUC <- sapply(allModelsCV, function(x){
			sapply(x[1:4], function(y){
						p <- predict(y, newdata=tcgaDataImputed)
						a <- AUC.uno(na.omit(os), tcgaSurvival[!is.na(p) & !is.na(tcgaSurvival)], scale(p)[!is.na(tcgaSurvival) &! is.na(p)], t)$auc
					})
		})
allModelsTcgaAUC <- array(allModelsTcgaAUC, dim=c(length(t), 4,100))
for(i in 1:4)
lines(t,rowMeans(allModelsTcgaAUC, dims=2)[,i], type='l', new=i==1, col=i)

t <- sort(unique(tcgaSurvival[,1]))
tcgaAUC <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], t)$auc)
tcgaAUC


t <- sort(unique(tcgaSurvival[,1]))
tcgaAucCV <- sapply(tcgaRisk, function(x) {
			sapply(1:100, function(foo){
						ix <-  sample(1:nrow(tcgaSurvival)%%5 +1 )!=1 ## sample 1/5
			AUC.uno(na.omit(os), tcgaSurvival[ix & !is.na(x) & !is.na(tcgaSurvival)], scale(x)[ix & !is.na(tcgaSurvival) &! is.na(x)], t)$auc})})
tcgaAucCV <- aperm(array(tcgaAucCV, dim=c(length(t),100,ncol(tcgaRisk))),c(2,1,3))

tcgaAucCVm <- colMeans(tcgaAucCV)
colModels <- c("#888888", set1[c(2,1,4,3,5)])
w <- c(1,6,7,8,5,3)
plot(t,tcgaAucCVm[,1], type="l", ylim=c(.5,1), lty=0, xlab="Years", ylab="AUC", xlim=c(0,6))
for(i in seq_along(w)){
	lines(t,tcgaAucCVm[,w[i]],  ylim=c(.5,1), lty=1, col=colModels[i], type='p', pch=16,cex=.5)
	s <- smooth.spline(t, tcgaAucCVm[,w[i]], df=10)
	#sigma <- sqrt(var(s$yin - s$y)/(1-s$lev))
	sigma <- apply(tcgaAucCV[,,w[i]],2,sd)
	lines(predict(s), col=colModels[i])
	polygon(c(t, rev(t)), c(s$y + 2.0*sigma*sqrt(s$lev) , rev(s$y - 2*sigma*sqrt(s$lev) )), col=paste(colModels[i],"00", sep=""), lty=1, border=colModels[i], lwd=0.5)
}
legend("bottomright", colnames(tcgaRisk)[w], col=colModels, lty=1)



allModelsCvTcgaConcordance <- sapply(allModelsCV, function(x){
			sapply(x, function(y) {
						p <- predictAllModels(y, newdata=tcgaDataImputed)
						survConcordance(tcgaSurvival ~ p)$concordance
					})
		})
boxplot(t(allModelsCvTcgaConcordance))
apply(allModelsCvTcgaConcordance, 1, quantile)


t <- sort(unique(tcgaSurvival[,1]))
tcgaAucCVm <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], t)$auc)
tcgaAucCVm

ResidualNonp <- function(fit){
	H0 <- basehaz(fit, centered = FALSE)
	FHaz <- splinefun(c(0,H0$time), c(0,H0$hazard), method="monoH.FC")
	s <- fit$surv
	-log(FHaz(s[,ncol(s)-1]))
}




apply(allModelsCvC,1,quantile)
boxplot(t(allModelsCvC), ylim=c(0,0.5), notch=TRUE, ylab="OXS R2", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)



replicates <- 100 ## number of replicates
tmp <- mclapply(1:replicates, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			coxCPSSOsTrain <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx, mainIdxOs], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
		}, mc.cores=10)

foo <- 1
tmp2 <- sapply(tmp, function(x){
			set.seed(foo)
			foo <<- foo+1
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			p <- predictAllModels(x, newdata=dataFrame[!trainIdx,])
			survConcordance(os[!trainIdx]~p)$concordance
		})
quantile(tmp2)

foo <- 1
tmp3 <- sapply(allModelsCV, function(x){
			set.seed(foo)
			foo <<- foo+1
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			p <- predictAllModels(x[[4]], newdata=dataFrame[!trainIdx,])
			survConcordance(os[!trainIdx]~p)$concordance
		})

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


 
v <- 10^seq(-3,3, 0.01)
pdf("CvVar.pdf", 2.5,2, pointsize=8)
par(mar=c(3,3,1,1)+.1, mgp=c(2,0.5,0), bty="n")
plot(v, CoxHD:::ConcordanceFromVariance(v), type="l", log='x', xlab="Var[h]", ylab="Concordance", ylim=c(0.5,1))
abline(h=seq(0.5,1,0.1), col="grey")
abline(v=10^seq(-3,3,1), col="grey")
lines(v, CoxHD:::ConcordanceFromVariance(v), type="l", lwd=2)
dev.off()


cppFunction('int add(int x, int y, int z) {
				int sum = x + y + z;
				return sum;
				}')


rfxglm <- function(x,y, family, tol=1e-3, ...){
	lambda <- 1
	lambda0 <- 0
	lambdaSeq <- function(x) {
		xl <- log10(x)
		10^seq(xl + 2, xl, length.out=100)
	}
	i <- 0
	while(abs(lambda - lambda0) > tol){
		fit <- glmnet(x,y,family, alpha=0, standardize=FALSE, lambda=lambdaSeq(lambda))
		d <- dim(coef(fit))
		lambda0 <- lambda
		lambda <- 1/var(coef(fit)[,d[2]])
		i <- i+1
		cat(i,"\t",lambda,"\n")
	}
	fit
}

library(flexmix)
n <- names(which(colSums(dataList$Genetics) > 50))
#n <- c("NPM1","DNMT3A","FLT3_ITD","FLT3_TKD","IDH1")
dl <- dataFrame[,n]
dl <- reshape(dl, direction='long', varying=colnames(dl), v.names="presence", timevar="gene" )
dl$gene <- factor(dl$gene, labels=n)
fit1<-stepFlexmix(cbind(presence,1-presence)~ -1 + gene|id, data=dl, k=5, model=FLXglm(family="binomial"), nrep=1)
1/(1+exp(-parameters(fit1)))
fit2 <- flexmix(presence~gene|id, data=dl, k=5, model=FLXMRmultinom())
1/(1+exp(-parameters(fit2)))


y<-c(rep(0,30-2),rep(1,2),
		rep(0,30-8),rep(1,8),
		rep(0,30-15),rep(1,15),
		rep(0,30-23),rep(1,23),
		rep(0,30-27),rep(1,27))

x<-c(rep(0,30),
		rep(1,30),
		rep(2,30),
		rep(3,30),
		rep(4,30))

# Initial state
state <- NULL

# MCMC parameters
nburn<-5000
nsave<-10000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
		tune=1.1)

prior <- list(a0=2,b0=1,beta0=rep(0,2), Sbeta0=diag(10000,2))

# Fit the model
fit1 <- DPbinary(y~x,prior=prior,mcmc=mcmc,state=state,status=TRUE) 
fit1

# Summary with HPD and Credibility intervals
summary(fit1)
summary(fit1,hpd=FALSE)

# Plot model parameters (to see the plots gradually set ask=TRUE)
plot(fit1)
plot(fit1,nfigr=2,nfigc=2)	

# Plot an specific model parameter (to see the plots gradually 
# set ask=TRUE)
plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="x")	
plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="ncluster")	
plot(fit1,ask=FALSE,param="link",nfigc=1,nfigr=1)

# Table of Pseudo Contour Probabilities
anova(fit1)

# Predictive Distribution
npred<-40  
xnew<-cbind(rep(1,npred),seq(0,4,length=npred))

pp<-predict(fit1,xnew)       

plot(seq(0,4,length=npred),pp$pmean,type='l',ylim=c(0,1),
		xlab="log2(concentration)",ylab="Probability")

# Adding MLE estimates
points(c(0,1,2,3,4),c(0.067,0.267,0.500,0.767,0.900),col="red")


## Lazy load all knitr cache
for(f in dir("cache", pattern=".rdx", full.names = TRUE))
	lazyLoad(sub(".rdx","",f))


d <- cirData[c("time1","time2","status","transplantCR1","transplantRel","index")]
d$status <- 0
d <- rbind(d, prsData[c("time1","time2","status","transplantCR1","transplantRel","index")])
d$relapse <- 0
d$relapse[1:nrow(cirData)] <- cirData$status
survRel <- survfit(Surv(time1, time2, status)~1, data=d)
plot(survRel, mark=NA)
relFree <- survfit(Surv(time1, time2, relapse)~1, data=d, subset=1:nrow(cirData))
survPrs <- survfit(Surv(time1, time2, status, type="counting")~1, data=prsData)
lines(relFree, col=2, mark=NA)
lines(survPrs, col=3, mark=NA)

crAdjust <- function(surv1, surv2){
	surv2 <- cumsum(c(1,diff(surv1$surv) * splinefun(surv2$time, surv2$surv, method="monoH.FC")(surv1$time[-1])))
	data.frame(time=surv1$time, surv=surv2)
}

s2 <- cumsum(c(1,diff(relFree$surv) * splinefun(survPrs$time, 1-survPrs$surv, method="monoH.FC")(relFree$time[-1])))
lines(relFree$time, s2, col=4)

survPredict <- function(surv){
	s <- survfit(surv~1)
	splinefun(s$time, s$surv, method="monoH.FC")
}

### Different approach
prsP <- survPredict(Surv(prsData$time2-prsData$time1, prsData$status))(0:5000)
cirP <- survPredict(Surv(cirData$time1, cirData$time2, cirData$status))(0:5000)

m <- matrix(0,5001,5001)
for(j in 1:5000)
	m[j,j:5001 ] <- (cirP[j]-cirP[j+1]) * (1-prsP[1:(5001-j+1)])

n <- rep(1,5001)
df <- diff(cirP)
for(j in 1:5000)
	n[j:5001 ] <- df[j] * (1-prsP[1:(5001-j+1)]) + n[j:5001] 

plot(survRel)
lines(n, col="brown")

coxphPrs <- coxph(Surv(prsData$time2-prsData$time1, prsData$status)~ pspline(prsData$time1, df=10))
timeDepPrs <- splinefun(prsData$time1[-coxphPrs$na.action], predict(coxphPrs))

survPrsTdAdj <- rep(1,5001)
df <- diff(cirP)
td <- exp(timeDepPrs(1:5000))
for(j in 1:5000)
	survPrsTdAdj[j:5001 ] <- df[j] * (1-prsP[1:(5001-j+1)]^td[j]) + survPrsTdAdj[j:5001] 

lines(survPrsTdAdj, col="blue")


#lines(1-colSums(m, na.rm=TRUE), col='brown')

#lines(s$time, 1-(1-s$surv)*(1-splinefun(t$time, t$surv, method="monoH.FC")(s$time)), col='brown')

## double check
s3 <- Surv(as.numeric(clinicalData$Date_LF-clinicalData$CR_date), clinicalData$Status & !is.na(clinicalData$Recurrence_date))
lines(survfit(s3 ~ 1), col="orange", mark=NA, conf.int=FALSE)

##
predictAbsCox <- function(fit, data, surv){
	H0 <- basehaz(fit, centered = FALSE)
	hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
	max.x <- max(surv[,2])
	x <- surv[,1]
	S <- exp(-hazardDist(x))
	S^exp(as.matrix(data[,names(coef(fit))])%*%coef(fit))
}

predictNpCox <- function(fit, risk = predict(fit), surv=fit$surv){
	#H0 <- basehaz(fit, centered = FALSE)
	#hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
	#x <- seq(0,max(surv[,1]), l=1000)
	#S <- exp(-hazardDist(x))
	S <- survfit(fit)
	
	f <- surv
	f[,2] <- 1-f[,2]
	F <- survfit(f~1)
	FCens <- splinefun( F$time,F$surv, method="monoH.FC")
	
	t <- S$time
	SS <- sapply(seq_along(risk), function(i){
				dS <- diff(S$surv^exp(risk[i]))
				w <- which.min(abs(t - surv[i,1]))
				#cat(w,"\n")
				sum(c(0,-t[-1] * dS * pmin(pmax(ifelse(surv[i,2], FCens(t[-1]), 1-FCens(t[-1])),0),1)))
			})
	return(SS)
}


p <- PredictOS(coxRFXNrmTD, coxRFXCirTD, coxRFXPrsTD, allData, x =365)
s <- survfit(coxRFXOsCR)
q <- s$surv[which.min(abs(s$time-365))] ^ exp(predict(coxRFXOsCR, newdata=allData))

osCR <- Surv(osData$time1, osData$time2, osData$status)
survConcordance(osCR ~ q)
survConcordance(osCR ~ p$os)

p <- sapply(1:10*365/2, function(i)  PredictOS(coxRFXNrmTD, coxRFXCirTD, coxRFXPrsTD, d, x =round(i))$os)
plot(1:10*365/2,cor(p,q), xlab="time", ylab='cor')
plot(1:10*365/2,1-apply(p,2, function(x) survConcordance(osCR ~ x)$concordance))

apply(apply(-sapply(concordanceCIRcv[[1]], `[[` , "C")[4:5,],2,rank),1,function(x) table(factor(x, levels=1:6)))


a <- sapply(concordanceCIRcv[[2]], `[[` , "coef")
a <- array(a, dim = c(dim(concordanceCIRcv[[2]][[1]]$coef), ncol(a)), dimnames=c(dimnames(concordanceCIRcv[[2]][[1]]$coef),NULL))
m <- rowMeans(a,dim=2)
m <- apply(a, 1:2, median)

boxplot(t(a[,1,]))


n <- coxRFXNrmTD
n$coefficients <- m[,"NRMrfx"]

r <- coxRFXPrsTD
r$coefficients <- m[,"PRSrfx"]

c <- coxRFXCirTD
c$coefficients <- m[,"CIRrfx"]

p <- PredictOS(coxRFXNrmTD, coxRFXCirTD, coxRFXPrsTD, allData, x =365)
q <- PredictOS(n, c, r, allData, x =365)

plot(p$os, q$os)
cor(p$os, q$os)

survConcordance(osCR ~ p$os)
survConcordance(osCR ~ q$os)

## more bagging
concordanceCIRbag <- mclapply(1:replicates, function(foo){
			g <- crGroups
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%2 +1 )!=1 ## sample 1/2
			dNrm <- nrmData[nrmData$index %in% which(trainIdx),names(g)]
			sNrm <- Surv(nrmData$time1, nrmData$time2, nrmData$status)[nrmData$index %in% which(trainIdx)]
			coxRFXNrmTD <- CoxRFX(dNrm, sNrm, groups=g, nu=1, which.mu = mainGroups)
			coxRFXNrmTD$coefficients["transplantRel"] <- 0
			dPrs <- prsData[prsData$index %in% which(trainIdx), names(g)]
			sPrs <- Surv(prsData$time2 - prsData$time1, prsData$status)[prsData$index %in% which(trainIdx)]
			coxRFXPrsTD <-  CoxRFX(dPrs, sPrs, groups=g, nu=1, which.mu = mainGroups)
			dCir <- cirData[cirData$index %in% which(trainIdx), names(g)]
			sCir <- Surv(cirData$time1, cirData$time2, cirData$status)[cirData$index %in% which(trainIdx)]
			coxRFXCirTD <-  CoxRFX(dCir, sCir, groups=g, which.mu = mainGroups)
			coxRFXCirTD$coefficients["transplantRel"] <- 0
			dOs <- osData[osData$index %in% which(trainIdx), names(g)]
			sOs <- Surv(osData$time1, osData$time2, osData$status)[osData$index %in% which(trainIdx)]
			coxRFXOsCR <- CoxRFX(dOs, sOs, groups=g, which.mu = mainGroups)
			
			allRisk365 <- PredictOS(coxRFXNrmTD = coxRFXNrmTD, coxRFXPrsTD = coxRFXPrsTD, coxRFXCirTD = coxRFXCirTD, allData, 365)
			allRisk1000 <- PredictOS(coxRFXNrmTD = coxRFXNrmTD, coxRFXPrsTD = coxRFXPrsTD, coxRFXCirTD = coxRFXCirTD, allData, 1000)
			
			p365 <- -allRisk365[,1]
			p1000 <-  -allRisk1000[,1]
			pCIR <- as.matrix(cirData[names(g)]) %*% coef(coxRFXCirTD)
			pPRS <- as.matrix(prsData[names(g)]) %*% coef(coxRFXPrsTD)
			pNRM <- as.matrix(nrmData[names(g)]) %*% coef(coxRFXNrmTD)
			pOS <- as.matrix(osData[names(g)]) %*% coef(coxRFXOsCR)
			
			Ctest <- c(
					CIRrfx = survConcordance(Surv(time1, time2, status)~ pCIR, data=cirData, subset = cirData$index %in% which(!trainIdx) )$concordance,
					PRSrfx = survConcordance(Surv(time2 - time1, status) ~ pPRS, data=prsData, subset=prsData$index %in% which(!trainIdx) )$concordance,
					NRMrfx = survConcordance(Surv(time1, time2, status)~  pNRM, data=nrmData, subset=nrmData$index %in% which(!trainIdx) )$concordance,
					OSrfx = survConcordance(Surv(time1, time2, status) ~ pOS, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
					OS365 = survConcordance(Surv(time1, time2, status) ~ p365, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
					OS1000 = survConcordance(Surv(time1,time2, status) ~ p1000, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance
			)
			
			Ctrain <- c(
					CIRrfx = survConcordance(Surv(time1, time2, status)~ pCIR, data=cirData, subset = cirData$index %in% which(trainIdx) )$concordance,
					PRSrfx = survConcordance(Surv(time2 - time1, status) ~ pPRS, data=prsData, subset=prsData$index %in% which(trainIdx) )$concordance,
					NRMrfx = survConcordance(Surv(time1, time2, status)~  pNRM, data=nrmData, subset=nrmData$index %in% which(trainIdx) )$concordance,
					OSrfx = survConcordance(Surv(time1, time2, status) ~ pOS, data=osData, subset=osData$index %in% which(trainIdx) )$concordance,
					OS365 = survConcordance(Surv(time1, time2, status) ~ p365, data=osData, subset=osData$index %in% which(trainIdx) )$concordance,
					OS1000 = survConcordance(Surv(time1,time2, status) ~ p1000, data=osData, subset=osData$index %in% which(trainIdx) )$concordance
			)
			
			coef <- cbind(CIRrfx=coef(coxRFXCirTD), PRSrfx=coef(coxRFXPrsTD), NRMrfx=coef(coxRFXNrmTD),  OSrfx=coef(coxRFXOsCR))
			
			return(list(Ctest=Ctest, Ctrain=Ctrain, coef=coef, allRisk365=allRisk365, allRisk1000=allRisk1000))
		}, mc.cores=20)

		
w <- sapply(concordanceCIRbag, class) != 'try-error'
x <- sapply(concordanceCIRbag[w], function(x) c(x$Ctest, x$Ctrain))
plot(x[,5], x[,11])
		
a <- sapply(concordanceCIRbag[w], `[[` , "coef")
a <- array(a, dim = c(dim(concordanceCIRbag[w][[1]]$coef), ncol(a)), dimnames=c(dimnames(concordanceCIRbag[w][[1]]$coef),NULL))
m <- rowMeans(a,dim=2)
m <- apply(a, 1:2, median)

n <- coxRFXNrmTD
n$coefficients <- m[,"NRMrfx"]

r <- coxRFXPrsTD
r$coefficients <- m[,"PRSrfx"]

c <- coxRFXCirTD
c$coefficients <- m[,"CIRrfx"]

r <- PredictOS(n, c, r, allData, x =365)

survConcordance(osCR ~ r$os)

dtf <- as.data.frame(sapply(rev(h$height)[1:20], function(x){
			c <- cut(d, x)
			o <- order(unlist(c$lower))
			Reduce("c",sapply(seq_along(c$lower), function(i) rep(i, length(unlist(c$lower[[i]])))))[o]
		}))
dtf$all <- 1:nrow(dtf)
dtf$size <- 1
treemap(dtf, index=c(paste0("V",1:5),"V20","all"), vSize='size', vCol="V1")


vcTest <- function(x){
	r <- t(sapply(levels(x$groups), function(l) {
				n <- c(names(coef(x))[which(x$groups==l)],l)
				w <- intersect(colnames(x$Hinv),n)
				df <- sum(diag(solve(x$Hinv[w,w]) %*% x$V[w,w]))
				c <- c(coef(x) - x$mu[l], x$mu[l])[w]
				#z <- c %*% solve(x$Hinv[w,w]) %*% c
				z <- sum(c^2/diag(x$Hinv[w,w]))
				c(z=z, df=df)
			}))
	p <- pchisq(r[,"z"], r[,"df"], lower.tail=FALSE)
	data.frame(r, p.value=p, sig=mg14:::sig2star(p))
}


## Clinical and splines
clinicalSpline <- as.data.frame(sapply(dataFrame[groups=="Clinical"], function(x){
			if(all(x[1:5] %in% 0:10)) return(x)
			y <- log(x+min(x)+1e-3)
			fit <- coxph(os ~ pspline(y, df=3), subset=trainIdx)
			predict(fit, newdata=data.frame(y=y))
		}))
for(n in names(clinicalSpline)) if(!all(dataFrame[1:5,n] %in% 0:10))
	plot(dataFrame[,n], clinicalSpline[,n], log='x', xlab=paste(n, '[observed]'), ylab = paste(n, '[spline]'))

summary(coxph(os ~ ., data=clinicalSpline, subset=!trainIdx))$concordance
summary(coxph(os ~ ., data=dataFrame[groups=="Clinical"]), subset=!trainIdx)$concordance

vcSimple <- function(fit){
	groups <- fit$groups
	Z <- fit$Z
	sigma2 <- fit$sigma2
	sapply(levels(groups), function(x) {
				ix <- groups == x
				sum(cov(as.matrix(Z[,ix, drop=FALSE]))) #* (sigma2[x] + fit$mu[x]^2)  
			})
}

EvalAbsolutePred <- function(prediction, surv, time, bins=seq(0,1,0.05)){
	c <- cut(prediction, bins)
	f <- survfit(surv ~ c)
	e <- summary(f, time)
	x <- sapply(strsplit(gsub("[a-z\\=\\(]|]","",e$strata),","), function(x) mean(as.numeric(x))); 
	#w <- 1/(e$std.err+.Machine$double.eps)^2
	w <- e$n[e$strata]
	std.err = 1/sum(w, na.rm=TRUE)
	mean.error = sum((e$surv-x)^2*w, na.rm=TRUE)*std.err
	return(list(mean.error=mean.error, std.err=std.err, survfit=e, x=x))
}

absPredError <- EvalAbsolutePred(allPredict$os, Surv(allData$time1, allData$time2, allData$status), time=365)

plot(absPredError$x, absPredError$survfit$surv, xlim=c(0,1), ylim=c(0,1))
segments(absPredError$x, absPredError$survfit$lower,absPredError$x, absPredError$survfit$upper)
abline(0,1)

PredictAbsoluteCoxph <- function(coxRFXOsCR, allData, time) {
	s <- survfit(coxRFXOsCR)
	q <- s$surv[which.min(abs(s$time-time))] ^ exp(predict(coxRFXOsCR, newdata=allData))
}
q <- PredictAbsoluteCoxph(coxRFXOsCR = coxRFXOsCR, allData = allData, time=365)

absPredErrorOs <- EvalAbsolutePred(q, Surv(allData$time1, allData$time2, allData$status), time=365)
plot(absPredErrorOs$x, absPredErrorOs$survfit$surv, xlim=c(0,1), ylim=c(0,1))
segments(absPredErrorOs$x, absPredErrorOs$survfit$lower,absPredErrorOs$x, absPredErrorOs$survfit$upper)
abline(0,1)

i <- 0
absoluteErrorsCIRcv <- lapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			i <- i+1
			sapply(1:replicates, function(foo){
						set.seed(foo)
						time <- 365
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 4/5
						coef <- concordanceCIRcv[[i]][[foo]][["coef"]]

						lpCIR <- as.matrix(cirData[names(coef[,"CIRrfx"])]) %*% coef[,"CIRrfx"]
						s <- survfit(Surv(time1, time2, status)~1, data=cirData, subset=cirData$index %in% which(trainIdx))
						pCIR <- s$surv[which.min(abs(s$time-time))] ^ exp(lpCIR-mean(lpCIR[cirData$index %in%which(trainIdx)]))
												
						lpPRS <- as.matrix(prsData[names(coef[,"PRSrfx"])]) %*% coef[,"PRSrfx"] 
						s <- survfit(Surv(time2- time1, status)~1, data=prsData, subset=prsData$index %in% which(trainIdx))
						pPRS <- s$surv[which.min(abs(s$time-time))] ^ exp(lpPRS-mean(lpPRS[prsData$index %in% which(trainIdx)]))
												
						lpNRM <- as.matrix(nrmData[names(coef[,"NRMrfx"])]) %*% coef[,"NRMrfx"]
						s <- survfit(Surv(time1, time2, status)~1, data=nrmData, subset=nrmData$index %in% which(trainIdx))
						pNRM <- s$surv[which.min(abs(s$time-time))] ^ exp(lpNRM-mean(lpNRM[nrmData$index %in% which(trainIdx)]))
						
						lpOS <- as.matrix(osData[names(coef[,"OSrfx"])]) %*% coef[,"OSrfx"]
						s <- survfit(Surv(time1, time2, status)~1, data=osData, subset=osData$index %in% which(trainIdx))
						pOS <- s$surv[which.min(abs(s$time-time))] ^ exp(lpOS-mean(lpOS[osData$index %in% which(trainIdx)]))
						
						p365 <- concordanceCIRcv[[i]][[foo]][["allRisk365"]]$os
						p1000 <- concordanceCIRcv[[i]][[foo]][["allRisk1000"]]$os
						err <- sapply(list(train=which(trainIdx), test=which(!trainIdx)), function(w)
									c(
											CIRrfx = EvalAbsolutePred(pCIR[cirData$index %in% w ], Surv(cirData$time1, cirData$time2, cirData$status)[cirData$index %in% w ], time=365)$mean.error,
											PRSrfx =  EvalAbsolutePred(pPRS[prsData$index %in% w ], Surv(prsData$time2- prsData$time1, prsData$status)[prsData$index %in% w ], time=365)$mean.error,
											NRMrfx =  EvalAbsolutePred(pNRM[nrmData$index %in% w ], Surv(nrmData$time1, nrmData$time2, nrmData$status)[nrmData$index %in% w ], time=365)$mean.error,
											OSrfx =  EvalAbsolutePred(pOS[osData$index %in% w ], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% w ], time=365)$mean.error,
											OS365 =  EvalAbsolutePred(p365[osData$index %in% w ], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% w ], time=365)$mean.error,
											OS1000 =  EvalAbsolutePred(p1000[osData$index %in% w ], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% w ], time=1000)$mean.error
									))
						return(err)
					}, simplify='array')
		})


ed <- survfit(Surv(c,is.na(clinicalData$CR_date))~1)
rem <- survfit(Surv(c,!is.na(clinicalData$CR_date))~1)

crAdjust <- function(fit1, fit2){
	int2 <- splinefun(fit2$time, fit2$surv,  method="monoH.FC")
	fit1$surv <- cumsum(c(1,diff(fit1$surv)) * int2(fit1$time))
	fit1
}

coxRFXCrTD <- CoxRFX(osData[1:1540, names(crGroups)], Surv(cr[,1], cr[,2]==2), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXEsTD <- CoxRFX(osData[1:1540, names(crGroups)], Surv(cr[,1], cr[,2]==1), groups=crGroups, which.mu = NULL)
save(coxRFXCirTD, coxRFXNrmTD, coxRFXPrsTD, coxRFXOsCR, coxRFXEsTD, coxRFXCrTD, cr, nrmData, cirData, prsData, osData, crGroups, data, file="../../code/predict/predictTest.RData")


plot(survfit(c ~ 1))