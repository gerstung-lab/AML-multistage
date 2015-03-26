# TODO: Add comment
# 
# Author: mg14
###############################################################################


library(shiny)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)
load("predictTest.RData", envir=globalenv())
cr <<- cr
set1 <- brewer.pal(8, "Set1")
VARIABLES <- names(crGroups)[!crGroups %in% c("Nuisance","GeneGene")] 
VARIABLES <- VARIABLES[order(apply(data[VARIABLES],2,var)*c(coef(coxRFXCirTD)[VARIABLES]^2+coef(coxRFXPrsTD)[VARIABLES]^2+coef(coxRFXNrmTD)[VARIABLES]^2), decreasing=TRUE)]
INTERACTIONS <- names(crGroups)[crGroups %in% "GeneGene"] 
NUISANCE <- names(crGroups)[crGroups %in% "Nuisance"] 

scaleFactors <- rep(1, length(VARIABLES))
names(scaleFactors) <- VARIABLES
w <- crGroups[VARIABLES] %in% c("Demographics","Clinical")
r <- regexpr("(?<=_)[0-9]+$", VARIABLES[w], perl=TRUE)
scaleFactors[w][r!=-1] <- as.numeric(regmatches(VARIABLES[w],r))

computeTotalPrs <- function(x, diffCir, prsP, tdPrmBaseline, risk) {
	xLen <- length(x)
	rs <- rep(1,xLen)
	for(j in x[-1])
		rs[j:xLen ] <- diffCir[j] * (1-prsP[1:(xLen-j+1)]^(tdPrmBaseline[j] * exp(risk))) + rs[j:xLen]
	return(rs)
}

cppFunction('NumericVector computeTotalPrsC(NumericVector x, NumericVector diffCir, NumericVector prsP, NumericVector tdPrmBaseline, double risk) {
				int xLen = x.size();
				double h;
                double r = exp(risk);
				NumericVector rs(xLen);
				for(int i = 0; i < xLen; ++i) rs[i] = 1;
				for(int j = 1; j < xLen; ++j){
                 h = tdPrmBaseline[j-1] * r;
				 for(int i = j; i < xLen; ++i){
				  rs[i] += diffCir[j-1] * (1-pow(prsP[i-j], h));
				 }
                }
				return rs;
				}')

cppFunction('NumericVector computeHierarchicalSurvival(NumericVector x, NumericVector diffS0, NumericVector S1Static, NumericVector haz1TimeDep) {
				int xLen = x.size();
				double h;
				NumericVector overallSurvival(xLen);
				for(int i = 0; i < xLen; ++i) overallSurvival[i] = 1;
				for(int j = 1; j < xLen; ++j){
				h = haz1TimeDep[j-1];
				for(int i = j; i < xLen; ++i){
				overallSurvival[i] += diffS0[j-1] * (1-pow(S1Static[i-j], h));
				}
				}
				return overallSurvival;
				}')

#data <- coxRFXCirTD$Z
# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
			getData <- reactive({
						l <- list()
						for(n in VARIABLES){
							l[[n]] <- ifelse(input[[n]]=="NA",NA,as.numeric(input[[n]]))
							if(is.null(input[[n]])) l[[n]] <- NA
						}
						for(n in INTERACTIONS){
							s <- strsplit(n, ":")[[1]]
							l[[n]] <- l[[s[1]]] * l[[s[2]]]
						}
						for(n in NUISANCE)
							l[[n]] <- NA
						out <- do.call("data.frame",l)
						names(out) <- names(l)
						out[VARIABLES] <- out[VARIABLES]/scaleFactors
						return(out)
						
					})
			output$ui <- renderUI({
						pdid <- input[["pdid"]]
						if(is.null(pdid)) pdid <- "reset"
						if( pdid=="reset"){
							#cat("reset\n")
							defaults <- data[1,]
							defaults[] <- NA
						}else
							defaults <- data[pdid,]
						
						defaults <- as.numeric(defaults)
						names(defaults) <- colnames(data)
						defaults[VARIABLES] <- defaults[VARIABLES] * scaleFactors
						#cat(defaults,"\n")
						
						lapply(VARIABLES, 
								function(x) {
									d <- defaults[x]
									if(crGroups[x] %in% c("Genetics","CNA","Fusions","Treatment")){
										if(!d %in% c(0,1)) d <- NA
										d <- paste(d)
										radioButtons(x, x, choices=c("present"= "1", "absent"="0", "NA"="NA"), selected=d)
									}else{
										r <- range(data[,x]*scaleFactors[x], na.rm=TRUE)
										numericInput(inputId=x, label=paste0(sub(paste0("_",scaleFactors[x],"$"),"",x), " [",r[1],"-",r[2],"]"), value=d, min=r[1], max=r[2] )
									}}
						
						)
					})
			x <- 0:2000
			
			computeIncidence <- function(coxRFX, r, x) {
				#r=PredictRiskMissing(coxRFX, data, var="var2")
				if(!is.null(coxRFX$na.action))
					coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
				#r <- PredictRiskMissing(coxRFX, data, var="var2")
				H0 <- basehaz(coxRFX, centered = FALSE)
				hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
				lambda0 <- hazardDist(x)
				r0 <- coxRFX$means %*% coef(coxRFX)
				inc <- exp(-lambda0* exp(r[,1]))
				ciup2 <- exp(-lambda0*exp( rep(r[,1] + 2*sqrt(r[,2]) * c(1), each=length(x))))
				cilo2 <- exp(-lambda0*exp( rep(r[,1] + 2*sqrt(r[,2]) * c(-1), each=length(x))))
				ciup <- exp(-lambda0*exp( rep(r[,1] + sqrt(r[,2]) * c(1), each=length(x))))
				cilo <- exp(-lambda0*exp( rep(r[,1] + sqrt(r[,2]) * c(-1), each=length(x))))
				#p <- PartialRisk(coxRFX, dataImputed)
				return(list(inc=inc, r=r, x=x, hazardDist=hazardDist, r0 = r0, ciup=ciup, cilo=cilo, ciup2=ciup2, cilo2=cilo2))
			}
			
			plotRisk <- function(coxRFX, incidence, p, xlab="Days after diagnosis", ylab="Incidence", col="#FF0000",mark=NA) {
				plot(survfit(coxRFX), xlab=xlab, ylab=ylab, mark=mark, conf.int=FALSE, fun=function(x) 1-x, ylim=c(0,1), xlim=c(0,2000), lty=2)
				#lines(survfit(coxRFX$surv ~ 1), lty=3, mark=NA, fun=function(x) 1-x)
				polygon( c(incidence$x, rev(incidence$x)), 1-c(incidence$ciup2, rev(incidence$cilo2)), col=paste0(col,"44"), border=NA)
				abline(h=seq(0,1,.2), lty=3)
				abline(v=seq(0,2000,365), lty=3)
				#polygon( c(risk$x, rev(risk$x)), 1-c(ciup, rev(cilo)), col=paste0(col,"44"), border=NA)
				lines( incidence$x, 1-incidence$inc, col=col, lwd=2)
				legend(ifelse((1-incidence$inc[length(incidence$inc)])>.5, "bottomright","topright"), c("Population avg","Predicted","95% CI"),bty="n", lty=c(2,1,NA), fill=c(NA,NA,paste0(col,"44")), border=c(NA,NA,NA), col=c(1,col,NA))
				u <- par("usr")
				par(new=T)
				m <- colMeans(PartialRisk(coxRFX))
				rds <- .05
				p <- (matrix(p, nrow=1) - m)/10 + rds
				colnames(p) <- names(m)
				l <- cbind(x=0.25,y=ifelse((1-incidence$inc[length(incidence$inc)])>.8,0.3,.85))
				c <- cut(incidence$r[,1]-incidence$r0, quantile(coxRFX$linear.predictor,seq(0,1,l=12)))
				stars((p[,c("Demographics","Treatment","Fusions","CNA","Genetics","GeneGene","Clinical"), drop=FALSE]), scale=FALSE, locations = l, xlim=c(0,1), ylim=c(0,1), lwd=1, col.stars=rev(brewer.pal(11,"RdBu"))[c])
				symbols(l, circles=rds, inches=FALSE, add=TRUE)
				text(l[1]+cos(2*pi*0:6/7)*2*rds,l[2]+sin(2*pi*0:6/7)*2*rds,substr(c("Demographics","Treatment","Fusions","CNA","Genetics","GeneGene","Clinical"),1,5), cex=.66)
				par(usr=u)
			}
			
			dataImputed <- reactive({ImputeMissing(data[1:1540,], getData()[,colnames(data)])})
			models <- c("Es","Cr","Cir","Nrm","Prs")
			riskMissing <- reactive({sapply(models, function(m){
									fit <- get(paste("coxRFX",m,"TD", sep=""))
									if(!is.null(fit$na.action))
										fit$Z <- fit$Z[-fit$na.action,]
									PredictRiskMissing(fit, getData(),  var="var2")}, simplify = FALSE)})
			partialRiskMissing <- reactive({sapply(models, function(m){
									fit <- get(paste("coxRFX",m,"TD", sep=""))
									if(!is.null(fit$na.action))
										fit$Z <- fit$Z[-fit$na.action,]
									PartialRisk(fit, dataImputed())}, simplify = FALSE)})
			
			output$Tab <- renderDataTable({
						x <- dataImputed()
						data.frame(Covariate=colnames(x),signif(data.frame(Input=as.numeric(getData()[,colnames(data)]), Imputed=as.numeric(x), 
												`Coef CIR`=coef(coxRFXCirTD), `Value CIR`= as.numeric(x)*coef(coxRFXCirTD),
												`Coef NRM`=coef(coxRFXNrmTD), `Value NRM`= as.numeric(x)*coef(coxRFXNrmTD),
												`Coef PRS`=coef(coxRFXPrsTD), `Value PRS`= as.numeric(x)*coef(coxRFXPrsTD)),2))
					})
			output$Risk <- renderDataTable({
						data.frame(Value=c(levels(coxRFXCirTD$groups),"total","s.d"),sapply(c("Cir","Nrm","Prs"), function(m){
											r <- riskMissing()[[m]]
											x <- get(paste("coxRFX",m,"TD", sep=""))
											p <- PartialRisk(x, newZ= rbind(dataImputed(),colMeans(data[1:1540,])))
											p <- p[1,]-p[2,]
											#p <- p[-length(p)]
											c(round(p,3), `total`=round(r[1,1] - mean(x$Z %*% coef(x)),3),
													`sd`=round(sqrt(r[1,2]),3))
										}))
					})
			## Convolution approach to PRM
			survPredict <- function(surv){
				s <- survfit(surv~1)
				splinefun(s$time, s$surv, method="monoH.FC")
			}
			prsP <- survPredict(Surv(prsData$time2-prsData$time1, prsData$status))(x) # Baseline Prs (measured from relapse)
			coxphPrs <- coxph(Surv(time2-time1, status)~ pspline(time1, df=10), data=prsData) # PRS baseline with spline-based dep on CR length)
			#timeDepPrs <- splinefun(prsData$time1[-coxphPrs$na.action], predict(coxphPrs)) # spline interpolation
			tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time1=x[-1])))	
			
			# CR adjustments to obtain absolute probabilities
			crAdjust <- function(x, y, time) {
				nrs <- cumsum(c(1,diff(x$inc) * splinefun(time, y$inc)(time[-1])))
			}
			
			output$KM <- renderPlot({
						par(mfrow=c(3,2), cex=1, bty="n")
						kmEs <- computeIncidence(coxRFX = coxRFXEsTD, r = riskMissing()[["Es"]], x=x)
						kmCr <- computeIncidence(coxRFX = coxRFXCrTD, r = riskMissing()[["Cr"]], x=x)
						es <- crAdjust(x= kmEs, time=x, y=kmCr) ## Correct KM estimate for competing risk
						cr <- crAdjust(x= kmCr, time=x, y=kmEs) ## Correct KM estimate for competing risk
						
						plotRisk(coxRFX = coxRFXEsTD, kmEs,  p=partialRiskMissing()[["Es"]], ylab="Incidence", xlab="Days after diagnosis", col=set1[5])
						lines(x, 1-es, col=set1[5], lty=2, lwd=2)
						title("Early deaths")
						
						plotRisk(coxRFX = coxRFXCrTD, kmCr,  p=partialRiskMissing()[["Cr"]], ylab="Incidence", xlab="Days after diagnosis", col=set1[4])
						lines(x, 1-cr, col=set1[4], lty=2, lwd=2)
						title("Remission")
						
						kmCir <- computeIncidence(coxRFX = coxRFXCirTD, r = riskMissing()[["Cir"]], x=x)
						kmNrm <-  computeIncidence(coxRFX = coxRFXNrmTD, r = riskMissing()[["Nrm"]], x=x)

						plotRisk(coxRFX = coxRFXCirTD, kmCir,  p=partialRiskMissing()[["Cir"]], ylab="Incidence", xlab="Days after remission", col=set1[3])
						cir <- crAdjust(x= kmCir, time=x, y=kmNrm) ## Correct KM estimate for competing risk
						lines(x, 1-cir, col=set1[3], lty=2, lwd=2)
						title("Relapse")

						plotRisk(coxRFX = coxRFXNrmTD, kmNrm, p=partialRiskMissing()[["Nrm"]], ylab="Incidence", xlab="Days after remission", col=set1[2])
						nrs <- crAdjust(x = kmNrm, time = x, y = kmCir)
						lines(x, 1-nrs, col=set1[2], lty=2, lwd=2)
						title("Non-relapse mortality")
						
						kmPrs <-  computeIncidence(coxRFX = coxRFXPrsTD, r = riskMissing()[["Prs"]], x=x)
						plotRisk(coxRFX = coxRFXPrsTD, kmPrs, p=partialRiskMissing()[["Prs"]], ylab="Incidence", xlab="Days after remission", col=set1[1])
						title("Post-relapse mortality")
						l <- length(kmCir$inc)
						

## Correct KM estimate for competing risk
						
						# Survival after Relapse
						rs <- computeHierarchicalSurvival(x = x, diffS0 = diff(cir), S1Static = prsP, haz1TimeDep = tdPrmBaseline * exp(kmPrs$r[,1]-kmPrs$r0))
						# Sum up to OS
						os <- 1-(1-nrs)-(1-rs)#nrs*rs
						
						xLen <- length(x)						
						plot(survfit(coxRFXOsCR), xlab="Days from remission", ylab="Survival", mark=NA, conf.int=FALSE,  xlim=c(0,2000), ylim=c(0,1), lty=2, xaxs='r')
						polygon(c(x, x[xLen]), c(nrs,1)  , border=NA, col=paste0(set1[2],"44"))
						polygon(c(x, rev(x)), c(nrs, rev(1-(1-nrs)-(1-rs))),  border=NA, col=paste0(set1[3],"44"))
						abline(h=seq(0,1,.2), lty=3)
						abline(v=seq(0,2000,365), lty=3)
			
						lines(x, os, col=set1[1], lwd=3)

						z <- c(365,3*365)
						y <- (os)[z+1]
						points(z,y, pch=16, col=set1[1])
						if("analytical" %in% input$ciType){
							## Confidence intervals
							PlogP2 <- function(x) {(x * log(x))^2}
							errOs <- kmNrm$r[,2] * PlogP2(kmNrm$inc) * (1-(1-kmCir$inc) * (1-kmPrs$inc))^2 + kmCir$r[,2] * PlogP2(kmCir$inc) * (1-kmPrs$inc)^2* kmNrm$inc^2 +  kmPrs$r[,2] * PlogP2(kmPrs$inc) * (1-kmCir$inc)^2* kmNrm$inc^2 
							errOs <- sqrt(errOs / PlogP2(os))
							osUp <- os ^ exp(2* errOs)
							osLo <- os ^ exp(-2*errOs)
							lines(x, osUp, col=set1[1], lty=3)
							lines(x, osLo, col=set1[1], lty=3)
							segments(z, osLo[z+1] ,z,osUp[z+1], col=set1[1], lwd=2)
						}
						if("simulated" %in% input$ciType){
							## Simulate CI
							nSim <- 200
							osCiMc <- sapply(1:nSim, function(i){
										r <- exp(rnorm(3,0,sqrt(c(kmCir$r[,2],kmNrm$r[,2],kmPrs$r[,2]))))
										nrs <- cumsum(c(1,diff(kmNrm$inc^r[2]) * kmCir$inc[-1]^r[1])) ## Correct KM estimate for competing risk
										diffCir <- diff(kmCir$inc^r[1]) * kmNrm$inc[-1]^r[2] ## Correct KM estimate for competing risk							
										rs <- computeTotalPrsC(x = x, diffCir = diffCir, prsP = prsP, tdPrmBaseline = tdPrmBaseline, risk = kmPrs$r[,1]-kmPrs$r0+log(r[3]))
										return(1-(1-nrs)-(1-rs))
									})
							osCiMcQ <- apply(osCiMc,1,quantile, c(0.025,0.975))
							lines(x, osCiMcQ[1,], col=set1[1], lty=2)
							lines(x, osCiMcQ[2,], col=set1[1], lty=2)
						}
						
						#polygon(c(x, rev(x)), c(nrslo2, rev(nrsup2))*c(rslo2, rev(rsup2)), col=paste("#FF000044",sep=""), border=NA)
						#polygon(c(x, rev(x)), c(nrslo, rev(nrsup))*c(rslo, rev(rsup)), col=paste("#FF000044",sep=""), border=NA)
						text(z, y, labels=round(y,2), pos=1)
						legend(ifelse(os[2000] > .5,"bottomright","topright"), col=set1[c(2,3,1)], lty=c(NA,NA,1), fill=c(set1[c(2,3)],NA), border=NA, lwd=2 , c("Non-relapse","Relapse","Total"), box.lwd = 0, title="Death by", bg="#FFFFFF88")
						title("Overall survival after remission")
						
						### Overall survival from enrollment
						#osEr <- computeHierarchicalSurvival(x = x, diffS0 = diff(1-kmEs$inc), S1Static = os, haz1TimeDep = rep(1,length(os)))
						# Sum up to OS
						#plot(x, 1-(1-kmEs$inc) * os, type="l", xlab="Days from diagnosis", ylab="Survival", main="Overall survival after enrollment", ylim=c(0,1)) 
						#lines(x, osEr, type="l", lty=2) 
						#abline(h=seq(0,1,.2), lty=3)
						#abline(v=seq(0,2000,365), lty=3)
						
					})
			
		})
