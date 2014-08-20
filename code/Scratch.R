# TODO: Add comment
# 
# Author: mg14
###############################################################################

data <- cbind(dataList$Genetics, dataList$Cytogenetics)
t <- table(rowSums(data))
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

pdf("nOnc.pdf",4,3,pointsize=8)
par(mar=c(3,3,0,0)+.1, mgp=c(2,0.5,0), bty="n")
barplot(t,col=rev(brewer.pal(11,"Spectral")), xlab="Oncogenic mutations", ylab="# Cases")
dev.off()

pdf("nOncSurv.pdf",3,2.5, pointsize=8)
par(mar=c(3,3,0,0)+.1, mgp=c(2,0.5,0), bty="n")
spec <- rainbow(11, 1, 0.75) #rev(brewer.pal(11,"Spectral"))
plot(survfit(survival ~ rowSums(data, na.rm=TRUE)), col=spec, xlab="Days", ylab="EFS")
legend("topright", legend=0:10, lty=1, bty="n",  col=spec, cex=0.8, title="Oncogenic mut")
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
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
x <- 10^seq(1,log10(2500), 0.1)
plot(x, exp(-hazardDist(x)), log="x", ylim=c(0,1), type="l")
xx <- seq(250, 2000, 250)
o <- order(coxRFXFit$mu)
boxplot(exp(-exp(coef(coxRFXFit))) / exp(hazardDist(xx[order(o)[coxRFXFit$groups]]))~ factor(coxRFXFit$groups, levels=levels(groups)[o]), at=xx,border=col1[o], horizontal=FALSE, las=1, lty=1, pch=NA, cex=.66, staplewex=0, add=TRUE, xaxt="n", yaxt="n", boxwex=200, col="NA")

boxplot(-coef(coxRFXFit) + ~ factor(coxRFXFit$groups, levels=levels(groups)[o]), border=col1[o], horizontal=FALSE, las=1, lty=1, pch=NA, cex=.66, staplewex=0, ylab="", ylim=c(-1,1), yaxt="n", labels="", xaxt="n")
abline(h=0)
x <- seq(-1,1,0.01)
p <- pretty(0:2000)
axis(side=1, at = p/250+.5,labels=p ,las=1)
dev.off()

par(mfrow=c(3,3))
i <- 1
for(g in levels(groups)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(coxRFXFit$coef[groups[whichRFX] == g]), each=length(x)) + mean(r) )), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(sqrt(coxRFXFit$sigma2[g]) * c(-1,1), each=length(x)) + mean(r) + coxRFXFit$mu[g])), col=paste(col1[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( coxRFXFit$mu[g] + mean(r))), col=col1[i], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) * exp(mean(r))), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}

par(mfrow=c(3,3))
i <- 1
for(g in levels(groups)){
	plot( x, exp(-hazardDist(x)), col="white", type="l", ylim=c(0,1), ylab="Expected survival", xlab="Days")
	m <- mean(rowSums(p[,colnames(p)!=g]))
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(range(partRiskTD[,g]), each=length(x)) +m)), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskTD[,g], c(0.25,0.75)), each=length(x))+m)), col=paste(col1[i],"44",sep=""), border=NA)
	polygon( c(x, rev(x)), exp(-hazardDist(c(x,rev(x)))*exp( rep(quantile(partRiskTD[,g], c(0.05,0.95)), each=length(x))+m)), col=paste(col1[i],"44",sep=""), border=NA)
	lines( x, exp(-hazardDist(x)*exp( median(partRiskTD[,g])+m)), col=col1[i], type="l", lwd=2)	
	lines( x, exp(-hazardDist(x) *exp(+m)), col="black", lwd=2)
	mtext(side=3, g)
	i <- i+1
}


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
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
#plot(survfit(osTD ~1))
#lines(0:5000, exp(-hazardDist(0:5000)),  col='red')
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
for(i in seq_along(l))
lines(x, invHazardDist(-log(l[i]) /exp(x) )/10000, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000)
mtext(side=4, line=2.5, "Time")
mtext(side=3, at = log(-log(l)/hazardDist(par("usr")[4]*10000)), text=paste(100*l, "% survive", sep=""))


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
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
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
mtext(side=3, at = log(-log(l)/hazardDist(par("usr")[4])), text=paste(100*l, "%", sep=""), cex=.66)
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