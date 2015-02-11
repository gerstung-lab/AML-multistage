# TODO: Add comment
# 
# Author: mg14
###############################################################################


library(shiny)
library(RColorBrewer)
library(CoxHD)
load("predict2.RData")
set1 <- brewer.pal(8, "Set1")
#data <- coxRFXCirTD$Z
# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
			getData <- reactive({
						l <- list()
						for(n in colnames(data)){
							l[[n]] <- ifelse(input[[n]]=="NA",NA,as.numeric(input[[n]]))
							if(is.null(input[[n]])) l[[n]] <- NA
						}
						out <- do.call("data.frame",l)
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
						#cat(defaults,"\n")
						
						variables <- c("transplantRel",names(sort(apply(data,2,var)*coef(coxRFXCirTD)^2, decreasing=TRUE)[-100]))
						lapply(variables, 
								function(x) {
									d <- defaults[x]
									if(crGroups[x] %in% c("Genetics","CNA","BT","Treatment")){
										if(!d %in% c(0,1)) d <- NA
										d <- paste(d)
										radioButtons(x, x, choices=c("present"= "1", "absent"="0", "NA"="NA"), selected=d)
									}else{
										numericInput(inputId=x, label=x, value=d, min=min(data[,x], na.rm=TRUE), max=max(data[,x],na.rm=TRUE) , step=1e-3)
									}}
										
						)
					})
			plotRisk <- function(coxRFX, data,xlab="Days after diagnosis", ylab="Incidence", mark=NA) {
				plot(survfit(coxRFX), xlab=xlab, ylab=ylab, mark=mark, conf.int=FALSE, fun=function(x) 1-x, ylim=c(0,1), xlim=c(0,2000))
				#lines(survfit(coxRFX$surv ~ 1), lty=3, mark=NA, fun=function(x) 1-x)
				abline(h=seq(0,1,.2), lty=3)
				abline(v=seq(0,2000,365), lty=3)
				coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
				r <- PredictRiskMissing(coxRFX, data, var="var2")
				H0 <- basehaz(coxRFX, centered = FALSE)
				hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
				x <- 0:2000
				ciup2 <- exp(-hazardDist(x)*exp( rep(r[,1] + 2*sqrt(r[,2]) * c(1), each=length(x))))
				cilo2 <- exp(-hazardDist(x)*exp( rep(r[,1] + 2*sqrt(r[,2]) * c(-1), each=length(x))))
				ciup <- exp(-hazardDist(x)*exp( rep(r[,1] + sqrt(r[,2]) * c(1), each=length(x))))
				cilo <- exp(-hazardDist(x)*exp( rep(r[,1] + sqrt(r[,2]) * c(-1), each=length(x))))
				polygon( c(x, rev(x)), 1-c(ciup2, rev(cilo2)), col=paste("#FF000044",sep=""), border=NA)
				polygon( c(x, rev(x)), 1-c(ciup, rev(cilo)), col=paste("#FF000044",sep=""), border=NA)
				inc <- exp(-hazardDist(x)* exp(r[,1]))
				lines( x, 1-inc, col="red", lwd=2)
				legend(ifelse((1-inc[length(inc)])>.5, "bottomright","topright"), c("Population avg","Predicted","95% CI"),bty="n", lty=c(1,1,NA), fill=c(NA,NA,paste("#FF000044",sep="")), border=c(NA,NA,NA), col=c(1,2,NA))
				return(list(inc=inc, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX), ciup=ciup, cilo=cilo, ciup2=ciup2, cilo2=cilo2))
			}
			output$Tab <- renderDataTable({
						d <- getData()
						x <- t(ImputeMissing(data[1:1540,], getData()))
						data.frame(Covariate=colnames(d),signif(data.frame(Input=as.numeric(t(d)), Imputed=x, `Coef CIR`=coef(coxRFXCirTD), `Value CIR`= x*coef(coxRFXCirTD),
								`Coef NRM`=coef(coxRFXNrmTD), `Value NRM`= x*coef(coxRFXNrmTD),
								`Coef PRS`=coef(coxRFXPrsTD), `Value PRS`= x*coef(coxRFXPrsTD)),2))
					})
			output$Risk <- renderDataTable({
						data.frame(Value=c("log hazard","s.d"),sapply(c("Cir","Nrm","Prs"), function(m){
									fit <- get(paste("coxRFX",m,"TD", sep=""))
									r <- PredictRiskMissing(fit, getData(),  var="var2")
									c(`log hazard`=round(r[1,1] - mean(get(paste("coxRFX",m,"TD", sep=""))$Z %*% coef(get(paste("coxRFX",m,"TD", sep="")))),3),
											`sd`=round(sqrt(r[1,2]),3))
								}))
					})
			output$KM <- renderPlot({
						par(mfrow=c(2,2), cex=1)
						hazCir <- plotRisk(coxRFX = coxRFXCirTD, data = getData(), ylab="Incidence")
						#lines(compRisk$`1 recurrence`$time, compRisk$`1 recurrence`$est, lty=2)
						title("Relapse")
						hazNrm <- plotRisk(coxRFX = coxRFXNrmTD, data = getData(), ylab="Mortality")
						#lines(compRisk$`1 dead`$time, compRisk$`1 dead`$est, lty=2)
						title("Non-relapse mortality")
						hazPrs <- plotRisk(coxRFX = coxRFXPrsTD, data = getData(), ylab="Mortality")
						title("Post-relapse mortality")
						l <- length(hazCir$inc)
						#M <- matrix(0, nrow=l-1, ncol=l)
						#for(i in seq_along(hazCir$inc)[-l])
						#	M[i,i:l] <- 1-hazPrs$inc[1:(l-i+1)]
						#lrs <- rowMeans(1-(1-hazCir$inc) * sapply(seq_along(hazCir$x)[-length(hazCir$x)], function(i) c(rep(1, i), 1-hazPrs$inc[1:(length(hazPrs$inc)-i)])))
						#rs <- 1-(1-hazCir$inc) * (1-hazPrs$inc[1000])
						#rs <- 1- diff(1-hazCir$inc) %*% M
						#nrm <- 1-(1-hazNrm$inc) * hazCir$inc
						#cir <- 1-(1-hazCir$inc) * hazNrm$inc
						
						cirKM <- survfit(coxRFXCirTD$surv ~ 1)
						nrmKM <- survfit(coxRFXNrmTD$surv ~ 1)							
						
						#nrs <- cumsum(c(1,diff(hazNrm$inc) * hazCir$inc[-1]*.9)) ## Correct KM estimate for competing risk
						nrs <- cumsum(c(1,diff(hazNrm$inc) * splinefun(cirKM$time, cirKM$surv)(hazNrm$x[-1]))) ## Correct KM estimate for competing risk
					    nrslo <- cumsum(c(1,diff(hazNrm$cilo)) * hazCir$inc)
						nrsup <- cumsum(c(1,diff(hazNrm$ciup)) * hazCir$inc)
						nrslo2 <- cumsum(c(1,diff(hazNrm$cilo2)) * hazCir$inc)
						nrsup2 <- cumsum(c(1,diff(hazNrm$ciup2)) * hazCir$inc)
						
						
						#cir <- cumsum(c(1,diff(hazCir$inc) * hazNrm$inc[-1] )) ## Correct KM estimate for competing risk
						cir <- cumsum(c(1,diff(hazCir$inc)* splinefun(nrmKM$time, nrmKM$surv)(hazCir$x[-1])) ) ## Correct KM estimate for competing risk
						cirlo <- cumsum(c(1,diff(hazCir$cilo))  * hazNrm$inc) 
						cirup <- cumsum(c(1,diff(hazCir$ciup))  * hazNrm$inc) 
						cirlo2 <- cumsum(c(1,diff(hazCir$cilo2))  * hazNrm$inc) 
						cirup2 <- cumsum(c(1,diff(hazCir$ciup2))  * hazNrm$inc) 
						
						
						## Adjust cumulative distributions for competing risks
						##cirAdj <- cumsum(c(1,diff(hazCir$inc)) * hazNrm$inc)
						rs <- 1- (1-cir) * (1-hazPrs$inc) 
						#p <- cumsum(c(1,diff(hazPrs$inc) * splinefun(nrmKM$time, nrmKM$surv)(hazPrs$x[-1])))
						#rs <- cumsum(c(1,diff(cir)*(1-hazPrs$inc[-1]))) ### TODO: double check
						rslo <- 1 - (1-hazCir$cilo) * (1-hazPrs$cilo)
						rsup <- 1 - (1-hazCir$ciup) * (1-hazPrs$ciup)
						rslo2 <- 1 - (1-hazCir$cilo2) * (1-hazPrs$cilo2)
						rsup2 <- 1 - (1-hazCir$ciup2) * (1-hazPrs$ciup2)
						#rsAdj <- cumsum(c(1,diff(rs)) * hazNrm$inc)
						##rsAdj <- 1- (1-cirAdj) * (1-hazPrs$inc) 
			            ##nrsAdj <- cumsum(c(1,diff(hazNrm$inc)) * rsAdj)
						
						## Prob of relapse and death
						## Hazard
						hrs = -log(rs)
						
						## Total hazard
						hnrs = -log(nrs)
						
						plot(survfit(coxRFXOsCR$surv ~ 1), xlab="Days", ylab="Fraction", mark=NA, conf.int=FALSE,  xlim=c(0,2000))
						
						lines(hazCir$x, rs, xlab="Time", ylab="Survival", ylim=c(0,1), type='l', col=set1[3])
						abline(h=seq(0,1,.2), lty=3)
						abline(v=seq(0,2000,365), lty=3)
						#lines(hazNrm$x, cirAdj, col=set1[4])
						#lines(hazNrm$x, hazCir$x, col=set1[4], lty=2)
						#
						lines(hazNrm$x, nrs, col=set1[2])
						#lines(hazCir$x, exp(-hrs -hnrs), col=set1[1])
						lines(hazNrm$x,nrs*rs, col=set1[1])
						polygon(c(hazNrm$x, rev(hazNrm$x)), c(nrslo2, rev(nrsup2))*c(rslo2, rev(rsup2)), col=paste("#FF000044",sep=""), border=NA)
						polygon(c(hazNrm$x, rev(hazNrm$x)), c(nrslo, rev(nrsup))*c(rslo, rev(rsup)), col=paste("#FF000044",sep=""), border=NA)
						#lines(survfit(osCR ~ 1), mark=NA, conf.int=FALSE)
						x <- c(365,3*365)
						y <- (nrs*rs)[x+1]
						points(x,y, pch=16, col=set1[1])
						text(x, y, labels=round(y,2), pos=1)
						legend("bottomright", col=set1[3:1], lty=1, c("Relapse","Non-relapse","Total"), bty="n")
						title("Overall survival")
					})

		})
