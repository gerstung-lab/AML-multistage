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
pastel1 <- brewer.pal(8, "Pastel1")
VARIABLES <- names(crGroups)[!crGroups %in% c("Nuisance","GeneGene")] 
VARIABLES <- VARIABLES[order(apply(data[VARIABLES],2,var)*c(coef(coxRFXRelTD)[VARIABLES]^2+coef(coxRFXPrdTD)[VARIABLES]^2+coef(coxRFXNrdTD)[VARIABLES]^2), decreasing=TRUE)]
INTERACTIONS <- names(crGroups)[crGroups %in% "GeneGene"] 
NUISANCE <- names(crGroups)[crGroups %in% "Nuisance"] 

SCALEFACTORS<- rep(1, length(VARIABLES))
names(SCALEFACTORS) <- VARIABLES
w <- crGroups[VARIABLES] %in% c("Demographics","Clinical")
r <- regexpr("(?<=_)[0-9]+$", VARIABLES[w], perl=TRUE)
SCALEFACTORS[w][r!=-1] <- as.numeric(regmatches(VARIABLES[w],r))

LABELS <- sapply(VARIABLES, function(x){
			r <- range(data[,x]*SCALEFACTORS[x], na.rm=TRUE)
			i <- paste0(" [",r[1],"-",r[2],"]")
			paste0(sub(paste0("_",SCALEFACTORS[x],"$"),"",x), ifelse(crGroups[x] %in% c("Genetics","Fusions","CNA","Treatment"),"",i))
		})

LABELS["AOD_10"] <- sub("AOD", "Age at diagnosis (yr)", LABELS["AOD_10"])
LABELS["LDH_1000"] <- sub("LDH", "Lactic Acid Dehydrogenase (units/l)", LABELS["LDH_1000"])
LABELS["wbc_100"] <- sub("wbc", "White cell count (1e-9/l)", LABELS["wbc_100"])
LABELS["HB_10"] <- sub("HB", "Hemoglobin (g/l)", LABELS["HB_10"])
LABELS["BM_Blasts_100"] <- sub("BM_Blasts", "Bone marrow blasts (%)", LABELS["BM_Blasts_100"])
LABELS["PB_Blasts_100"] <- sub("PB_Blasts", "Peripheral blood blasts (%)", LABELS["PB_Blasts_100"])
LABELS["platelet_100"] <- sub("platelet", "Platelet count (1e-9/l)", LABELS["platelet_100"])
LABELS["VPA"] <- "VPA (Valproic acid)"
LABELS["transplantCR1"] <- "Allograft in CR1"
LABELS["transplantRel"] <- "Allograft after Relapse"
LABELS["gender"] <- "Gender: 1=male, 2=female"
LABELS <- sub("t_*([a-z,0-9]+)_([a-z,0-9]+)", "t(\\1;\\2)", LABELS)

#* AOD: Age on diagnosis 
#* LDH: Lactic Acid Dehydrogenase (units/l)
#* WBC: White cell count (1e-9/l), 
#* HB: Hemoglobin (g/l), 
#* BM_Blasts: Bone marrow blasts (%)
#* PB_Blasts: Peripheral blood blasts (%)


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

addGrid <- function() {
	abline(h=seq(0,1,.2), lty=3)
	abline(v=seq(0,2000,365), lty=3)
}

#data <- coxRFXRelTD$Z
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
						out[VARIABLES] <- out[VARIABLES]/SCALEFACTORS
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
						defaults[VARIABLES] <- defaults[VARIABLES] * SCALEFACTORS
						#cat(defaults,"\n")
						
						lapply(VARIABLES, 
								function(x) {
									d <- defaults[x]
									if(crGroups[x] %in% c("Genetics","CNA","Fusions","Treatment")){
										if(!d %in% c(0,1)) d <- NA
										d <- paste(d)
										radioButtons(x, label=if(crGroups[x]=="Genetics") tags$em(LABELS[x]) else LABELS[x], choices=c("present"= "1", "absent"="0", "NA"="NA"), selected=d)
									}else{
										r <- range(data[,x]*SCALEFACTORS[x], na.rm=TRUE)
										numericInput(inputId=x, label=LABELS[x], value=d, min=r[1], max=r[2] )
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
			
			plotRisk <- function(coxRFX, incidence, p, xlab="Days after diagnosis", ylab="Incidence", col="#FF0000",mark=NA, lty=2) {
				plot(survfit(coxRFX), xlab=xlab, ylab=ylab, mark=mark, conf.int=FALSE, fun=function(x) 1-x, ylim=c(0,1), xlim=c(0,2000), lty=3)
				#lines(survfit(coxRFX$surv ~ 1), lty=3, mark=NA, fun=function(x) 1-x)
				polygon( c(incidence$x, rev(incidence$x)), 1-c(incidence$ciup2, rev(incidence$cilo2)), col=paste0(col,"44"), border=NA)
				addGrid()
				#polygon( c(risk$x, rev(risk$x)), 1-c(ciup, rev(cilo)), col=paste0(col,"44"), border=NA)
				lines( incidence$x, 1-incidence$inc, col=col, lwd=2, lty=lty)
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
			models <- c("Ncd","Cr","Rel","Nrd","Prd")
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
												`Coef NCD`=coef(coxRFXNcdTD), `Value NCD`= as.numeric(x)*coef(coxRFXNcdTD),
												`Coef CR`=coef(coxRFXCrTD), `Value CR`= as.numeric(x)*coef(coxRFXCrTD),
												`Coef NRD`=coef(coxRFXNrdTD), `Value NRD`= as.numeric(x)*coef(coxRFXNrdTD),
												`Coef Rel`=coef(coxRFXRelTD), `Value Rel`= as.numeric(x)*coef(coxRFXRelTD),
												`Coef PRD`=coef(coxRFXPrdTD), `Value PRD`= as.numeric(x)*coef(coxRFXPrdTD)),2))
					})
			output$Risk <- renderDataTable({
						t <- sapply(c("Ncd","Cr","Nrd","Rel","Prd"), function(m){
									r <- riskMissing()[[m]]
									x <- get(paste("coxRFX",m,"TD", sep=""))
									p <- PartialRisk(x, newZ= rbind(dataImputed(),colMeans(data[1:1540,])))
									p <- p[1,]-p[2,]
									#p <- p[-length(p)]
									c(round(p,3), `total`=round(r[1,1] - mean(x$Z %*% coef(x)),3),
											`sd`=round(sqrt(r[1,2]),3))
								})
						colnames(t) <- c("Death without CR (NCD)", "Complete remission (CR)", "Death without relapse (NRD)", "Relapse","Death after relapse (PRD)")
						data.frame(Value=c(levels(coxRFXRelTD$groups),"total","s.d"), t, check.names=FALSE)
					})
			## Convolution approach to PRM
			survPredict <- function(surv){
				s <- survfit(surv~1)
				splinefun(s$time, s$surv, method="monoH.FC")
			}
			prsP <- survPredict(Surv(prdData$time2-prdData$time1, prdData$status))(x) # Baseline Prs (measured from relapse)
			coxphPrs <- coxph(Surv(time2-time1, status)~ pspline(time1, df=10), data=prdData) # PRS baseline with spline-based dep on CR length)
			tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time1=x[-1])))	

			coxphOs <- coxph(Surv(time2-time1, status)~ pspline(cr[osData$index,1], df=10), data=osData) # PRS baseline with spline-based dep on CR length)
			tdOsBaseline <- exp(predict(coxphOs, newdata=data.frame(time1=x[-1])))	
			
			# CR adjustments to obtain absolute probabilities
			crAdjust <- function(x, y, time=x$x) {
				nrs <- cumsum(c(1,diff(x$inc) * splinefun(time, y$inc)(time[-1])))
			}
			
			crAdjustIncidence <- function(inc1, inc2, time=inc1$x) {
				incOut <- incidence
				incOut$inc <- crAdjust(inc1, inc2)
				incOut$ciUp <- crAdjust(inc1, inc2)
			}
			
			output$KM <- renderPlot({
						par(bty="n", mar=c(3,3,2,1), mgp=c(2,0.5,0), tcl=-.25, xaxs="i", yaxs="i")
						layout(matrix(1:3, ncol=3), widths=c(1,1,0.5))
						par(cex=1)
						
						## KM incidence of NCD and CR
						kmNcd <- computeIncidence(coxRFX = coxRFXNcdTD, r = riskMissing()[["Ncd"]], x=x)
						kmCr <- computeIncidence(coxRFX = coxRFXCrTD, r = riskMissing()[["Cr"]], x=x)
						
						## Correct KM estimate for competing risk
						ncd <- crAdjust(x= kmNcd, time=x, y=kmCr) ## Correct KM estimate for competing risk
						cr <- crAdjust(x= kmCr, time=x, y=kmNcd) ## Correct KM estimate for competing risk
						
						## KM incidence of Relapse and NRD
						kmRel <- computeIncidence(coxRFX = coxRFXRelTD, r = riskMissing()[["Rel"]], x=x)
						kmNrd <-  computeIncidence(coxRFX = coxRFXNrdTD, r = riskMissing()[["Nrd"]], x=x)

						## Correct KM estimate for competing risk
						relCr <- crAdjust(x= kmRel, time=x, y=kmNrd) ## Correct KM estimate for competing risk
						nrsCr <- crAdjust(x = kmNrd, time = x, y = kmRel)
						
						## KM incidence of PRS
						kmPrs <-  computeIncidence(coxRFX = coxRFXPrdTD, r = riskMissing()[["Prd"]], x=x)
						
						## Outcome after Remission
						rsCr <- computeHierarchicalSurvival(x = x, diffS0 = diff(relCr), S1Static = prsP, haz1TimeDep = tdPrmBaseline * exp(kmPrs$r[,1]-kmPrs$r0))
						osCr <- 1-(1-nrsCr)-(1-rsCr)
						
						## Outcome from diagnosis
						osDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = osCr, haz1TimeDep = tdOsBaseline)
						nrsDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = nrsCr, haz1TimeDep = tdOsBaseline)
						rsDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = rsCr, haz1TimeDep = tdOsBaseline)
						relDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = relCr, haz1TimeDep = tdOsBaseline)
						
						xLen <- length(x)						
						
						# Sum up to OS
						plot(x, 1-(1-ncd)-(1-osDiag), type="l", xlab="Days from diagnosis", ylab="Probability", main="Outcome after diagnosis", ylim=c(0,1), lwd=3, lty=0) 
						y0 <- 1
						y <- ncd
						polygon(c(x, x[xLen]), c(y,y0), border=NA, col=pastel1[1])
						y0 <- y0 - (1-ncd)
						y <- y - (1-nrsDiag)
						polygon(c(x, rev(x)), c(y, rev(y0))  , border=NA, col=pastel1[2])
						y0 <- y0 - (1-nrsDiag)
						y <- y - (1-rsDiag)
						osDiag <- y
						polygon(c(x, rev(x)), c(y, rev(y0)),  border=NA, col=pastel1[3])
						y0 <- y0 - (1-rsDiag)
						y <- y - (1-relDiag) + (1-rsDiag) 
						polygon(c(x, rev(x)), c(y, rev(y0)),  border=NA, col=pastel1[5])
						polygon(c(x, rev(x)), c(cr - (1-ncd), rev(y)),  border=NA, col=pastel1[4])
						polygon(c(x, rev(x)), c(cr - (1-ncd), rev(rep(0, length(x)))),  border=NA, col="#DDDDDD")
						lines(x, osDiag, lwd=3)

						addGrid()
						

						### Plot outcome after remission
						plot(NA,NA, xlab="Days from remission", ylab="Probability",  xlim=c(0,2000), ylim=c(0,1), lty=2)
						polygon(c(x, x[xLen]), c(nrsCr,1)  , border=NA, col=pastel1[2])
						polygon(c(x, rev(x)), c(nrsCr, rev(osCr)),  border=NA, col=pastel1[3])
						polygon(c(x, rev(x)), c(osCr, rev(1-(1-nrsCr)-(1-relCr))),  border=NA, col=pastel1[5])
						polygon(c(x, rev(x)), c(1-(1-nrsCr)-(1-relCr), rep(0,length(x))),  border=NA, col=pastel1[4])
						abline(h=seq(0,1,.2), lty=3)
						abline(v=seq(0,2000,365), lty=3)
						
						lines(x, osCr, col=1, lwd=3)
						
						z <- c(365,3*365)
						y <- (osCr)[z+1]
						points(z,y, pch=16, col=1)
						if("analytical" %in% input$ciType){
							## Confidence intervals
							PlogP2 <- function(x) {(x * log(x))^2}
							errOs <- kmNrd$r[,2] * PlogP2(kmNrd$inc) * (1-(1-kmRel$inc) * (1-kmPrs$inc))^2 + kmRel$r[,2] * PlogP2(kmRel$inc) * (1-kmPrs$inc)^2* kmNrd$inc^2 +  kmPrs$r[,2] * PlogP2(kmPrs$inc) * (1-kmRel$inc)^2* kmNrd$inc^2 
							errOs <- sqrt(errOs / PlogP2(osCr))
							osUp <- osCr ^ exp(2* errOs)
							osLo <- osCr ^ exp(-2*errOs)
							lines(x, osUp, col=1, lty=2)
							lines(x, osLo, col=1, lty=2)
							#segments(z, osLo[z+1] ,z,osUp[z+1], col=1, lwd=2)
						}
						if("simulated" %in% input$ciType){
							## Simulate CI
							nSim <- 200
							osCiMc <- sapply(1:nSim, function(i){
										r <- exp(rnorm(3,0,sqrt(c(kmRel$r[,2],kmNrd$r[,2],kmPrs$r[,2]))))
										nrs <- cumsum(c(1,diff(kmNrd$inc^r[2]) * kmRel$inc[-1]^r[1])) ## Correct KM estimate for competing risk
										diffCir <- diff(kmRel$inc^r[1]) * kmNrd$inc[-1]^r[2] ## Correct KM estimate for competing risk							
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
						title("Outcome after remission")
						
						par(mar=c(0,0,0,0))
						plot(NA,NA, xlab="",ylab="", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1))
						legend(x=0,y=1, col=c(NA,NA,NA,NA,NA,NA,"black"), lty=c(NA,NA,NA,NA,NA,NA,1), fill=c(pastel1[c(1,2,3,5,4)],"#DDDDDD",NA), border=c(1,1,1,1,1,1,NA), lwd=2 , y.intersp = 1.5, c("Death without \nremission","Death without \nrelapse","Death after \nrelapse","Alive after \nrelapse","Alive in CR1", "Alive in \ninduction", "Overall survival"), box.lwd = 0,  bg="#FFFFFF88", seg.len=1)

					})
			
		})
