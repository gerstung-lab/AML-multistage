library(shiny)
library(xtable)
library(RColorBrewer)
library(CoxHD)
library(Rcpp)
load("multistage.RData", envir=globalenv())
cr <<- cr
set1 <- brewer.pal(8, "Set1")
pastel1 <- brewer.pal(8, "Pastel1")
s <- !crGroups %in% c("Nuisance","GeneGene")  & ! names(crGroups) %in% c("ATRA","VPA")
VARIABLES <- names(crGroups)[s] 
rg <- c("Fusions"=5, "CNA"=4,"Genetics"=3, "Clinical"=7, "Demographics"=8, "Treatment"=6)
o <- order(rg[crGroups[s]],((coef(coxRFXPrdTD)^2/diag(coxRFXPrdTD$var2) + coef(coxRFXNrdTD)^2/diag(coxRFXNrdTD$var2) + coef(coxRFXRelTD)^2/diag(coxRFXRelTD$var2)) * apply(data[names(crGroups)], 2, var))[VARIABLES], decreasing=TRUE)
VARIABLES <- VARIABLES[o]
NEWGRP <- c(0,diff(as.numeric(as.factor(crGroups))[s][o])) != 0
names(NEWGRP) <- VARIABLES
INTERACTIONS <- names(crGroups)[crGroups %in% "GeneGene"] 
NUISANCE <- names(crGroups)[crGroups %in% "Nuisance" | names(crGroups) %in%  c("ATRA","VPA")] 

SCALEFACTORS<- rep(1, length(VARIABLES))
names(SCALEFACTORS) <- VARIABLES
w <- crGroups[VARIABLES] %in% c("Demographics","Clinical")
r <- regexpr("(?<=_)[0-9]+$", VARIABLES[w], perl=TRUE)
SCALEFACTORS[w][r!=-1] <- as.numeric(regmatches(VARIABLES[w],r))

CATEGORIES <- sapply(VARIABLES, function(x){
			if(length(unique(data[,x])) <= 10){
				c <- min(data[,x]):max(data[,x])
				if(all(c %in% 0:1))
					names(c) <- c("absent","present")
				else if(x =="gender")
					names(c) <- c("male","female")
				return(c)
			}
			else
				NULL
		})

LABELS <- sapply(VARIABLES, function(x){
			r <- round(range(data[,x]*SCALEFACTORS[x], na.rm=TRUE),1)
			i <- paste0(" [",r[1],"-",r[2],"]")
			paste0(sub(paste0("_",SCALEFACTORS[x],"$"),"",x), ifelse(is.null(CATEGORIES[[x]]),i,""))
		})

LIMITS <- sapply(VARIABLES, function(x){
			r <- round(range(data[,x], na.rm=TRUE),1)})

LABELS["AOD_10"] <- sub("AOD", "Age at diagnosis (yr)", LABELS["AOD_10"])
LABELS["LDH_1000"] <- sub("LDH", "Lactic Acid Dehydrogenase (units/l)", LABELS["LDH_1000"])
LABELS["wbc_100"] <- sub("wbc", "White cell count (1e-9/l)", LABELS["wbc_100"])
LABELS["HB_10"] <- sub("HB", "Hemoglobin (g/dl)", LABELS["HB_10"])
LABELS["BM_Blasts_100"] <- sub("BM_Blasts", "Bone marrow blasts (%)", LABELS["BM_Blasts_100"])
LABELS["PB_Blasts_100"] <- sub("PB_Blasts", "Peripheral blood blasts (%)", LABELS["PB_Blasts_100"])
LABELS["platelet_100"] <- sub("platelet", "Platelet count (1e-9/l)", LABELS["platelet_100"])
LABELS["VPA"] <- "VPA (Valproic acid)"
LABELS["transplantCR1"] <- "Allograft in CR1"
LABELS["transplantRel"] <- "Allograft after Relapse"
LABELS["gender"] <- "Gender"
LABELS <- sub("t_*([a-z,0-9]+)_([a-z,0-9]+)", "t(\\1;\\2)", LABELS)
LABELS[crGroups[VARIABLES] %in% c("Fusions","CNA")] <- gsub("_","/",LABELS[crGroups[VARIABLES] %in% c("Fusions","CNA")])
LABELS <- sub("plus","+",LABELS)
LABELS <- sub("minus|^mono","-",LABELS)
LABELS <- sub("(_|/)*other"," (other)", LABELS)
LABELS <- sub("_([0-9a-zA-Z]+)"," (\\1)", LABELS)


COMPVAR <- list(`Allogeneic HSCT`=c(none="none", `in first CR`="transplantCR1", `after relapse`="transplantRel"), `AML type`=c(primary='AML', secondary='sAML',tertiary='tAML',other='oAML')) ## Compound variables (factors)
COMPIDX <- numeric(length(VARIABLES))
names(COMPIDX) <- VARIABLES
COMPIDX[c("transplantRel","oAML")] <- 1 ## Index of last elements for display
VAR2COMP <- unlist(sapply(names(COMPVAR), function(n) rep(n, length(COMPVAR[[n]])))) 
names(VAR2COMP) <- unlist(COMPVAR)

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
					if(diffS0[j-1] != 0){
						h = haz1TimeDep[j-1];
						for(int i = j; i < xLen; ++i){
							overallSurvival[i] += diffS0[j-1] * (1-pow(S1Static[i-j], h));
						}
					}
				}
				return overallSurvival;
				}')

addGrid <- function(scale=1) {
	abline(h=seq(0,1,.2), lty=3)
	abline(v=seq(0,2000,365.25)/scale, lty=3)
}

wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px;"


# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
			getData <- reactive({
						input$compute
						isolate({
									l <- list()
									for(n in VARIABLES){
										if(!n %in% unlist(COMPVAR)){
											l[[n]] <- ifelse(input[[n]]=="NA",NA,as.numeric(input[[n]]))
											if(is.null(input[[n]])) l[[n]] <- NA
										}else{
											l[[n]] <- ifelse(input[[VAR2COMP[n]]]=="NA", NA, input[[VAR2COMP[n]]]==n) + 0
											if(is.null(input[[VAR2COMP[n]]])) l[[n]] <- NA
										}
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
					})
			
			makeWarning <- function(title="Warning", message=HTML(""), id="warningModal") {
				list(tags$div(HTML(paste0('<div class="modal fade" id="',id,'" tabindex="-1" role="dialog" aria-labelledby="myModalLabel">
												<div class="modal-dialog" role="document">
												<div class="modal-content">
												<div class="modal-header" style="background-color:rgb(242,222,222);color:rgb(169,68,66);border-top-left-radius:6px; border-top-right-radius:6px">
												<h4 class="modal-title" id="myModalLabel">',title,'</h4>
												</div>
												<div class="modal-body">')), message,
								HTML(paste0('</div>
												<div class="modal-footer">
												<button type="button" class="btn btn-default" data-dismiss="modal">Dismiss</button>
												</div>
												</div>
												</div>
												</div>
												<script>$(\'#',id,'\').modal(\'show\')</script>'))))
			}
			
			output$inputCheck <- renderUI({
						g <- as.numeric(getData()[VARIABLES])
						cond <- g < LIMITS[1,] | g > LIMITS[2,]
						if(any(na.omit(cond))){
							makeWarning(message=list(HTML('<h5>The following values are out of range:</h5>'),
											renderTable(data.frame(Variable=LABELS[VARIABLES[which(cond)]], `Entered value`=(g * SCALEFACTORS[VARIABLES])[which(cond)], check.names=FALSE), include.rownames=FALSE),
											HTML('This is likely to lead to uncontrolled behaviour of the predictions.')
									), id='warningInput')
						}
					})
			
			output$ui <- renderUI({
						pdid <- input[["pdid"]]
						if(is.null(pdid)) pdid <- "reset all variables"
						if( pdid=="reset all variables"){
							#cat("reset\n")
							defaults <- data[1,]
							defaults[] <- NA
						}else
							defaults <- data[pdid,]
						
						defaults <- as.numeric(defaults)
						
						## Obfuscation
						defaults <- signif(defaults * 20,1)/20
						
						names(defaults) <- colnames(data)
						defaults[VARIABLES] <- defaults[VARIABLES] * SCALEFACTORS
						#cat(defaults,"\n")
						
						makeMenu <- function(x) {
							d <- defaults[x]
							f <- if(x %in% unlist(COMPVAR)){
										if(!COMPIDX[x]) return(NULL)
										s <- defaults[COMPVAR[[VAR2COMP[x]]][-1]]
										w <- if(any(is.na(s))) 'N/A' else if(all(s==0)) 1 else if(any(!s %in% c(0,1))) 'N/A' else which(s==1)+1
										c <- c(COMPVAR[[VAR2COMP[x]]], "N/A"="NA")
										radioButtons(VAR2COMP[x], label=VAR2COMP[x], choices=c, selected=c[w], inline=FALSE)
									}else if(crGroups[x] %in% c("Genetics","CNA","Fusions","Treatment")){
										if(!d %in% c(0,1)) d <- NA
										d <- paste(d)
										radioButtons(x, label=if(crGroups[x]=="Genetics") tags$em(LABELS[x]) else LABELS[x], choices=c("present"= "1", "absent"="0", "N/A"="NA"), selected=d, inline=TRUE)
									}else{
										r <- round(quantile(data[,x]*SCALEFACTORS[x], c(0.05,0.95), na.rm=TRUE),1)
										if(is.null(CATEGORIES[[x]]))
											numericInput(inputId=x, label=LABELS[x], value=d, min=r[1], max=r[2], step = if(round(min(data[,x]*SCALEFACTORS[x], na.rm=TRUE),1) %% 1 ==0) 1 else 0.1)
										else{
											if(!d %in% 0:10) d <- NA
											d <- paste(d)
											radioButtons(x, label=LABELS[x], choices=c(CATEGORIES[[x]],"N/A"="NA"), selected=d, inline=TRUE)
										}
									}
							h <- if(NEWGRP[x]) list(tags$em(tags$b(crGroups[x]))) else NULL
							list(h,f)}
						
						list(wellPanel(	list(	tags$b("2. or enter/amend variables"),	
												#tags$em(tags$small("Click to see a list of variables")),
												#tags$hr(),
												wellPanel(
														actionLink("showClinical", 					
																tags$div("Clinical variables", HTML("&#9662;")),
																style = "color:rgb(0,0,0);"
														),
														style = wellStyle),
												conditionalPanel(condition = 'input.showClinical % 2',
														wellPanel(
																tags$em(tags$b(crGroups[VARIABLES[1]])),
																lapply(VARIABLES[crGroups[VARIABLES] %in% c("Clinical","Demographics")], makeMenu),
																style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
														)
												),
												#tags$hr(),
#												conditionalPanel(condition = 'input.showClinical % 2', tags$hr()),
#										)
#								),
#								wellPanel(	list(
												wellPanel(
														actionLink("showDrivers", 					
																tags$div("Driver mutations", HTML("&#9662;")),
																style = "color:rgb(0,0,0);"
														),
														style = wellStyle),
												
												conditionalPanel(condition = 'input.showDrivers % 2',
														wellPanel(
																lapply(VARIABLES[crGroups[VARIABLES] %in% c("Genetics","Fusions","CNA")], makeMenu),
																tags$hr(),
																style = paste(wellStyle,"margin-top:-20px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
														)
												),
												wellPanel(
														actionLink("showTreatment", 					
																tags$div("Treatment", HTML("&#9662;")),
																style = "color:rgb(0,0,0);"
														),
														style = paste(wellStyle, "margin-bottom: 0px")),
												conditionalPanel(condition = 'input.showTreatment % 2',
														wellPanel(
																lapply(VARIABLES[crGroups[VARIABLES] == "Treatment"], makeMenu),
																style = paste(wellStyle,"margin-bottom:0px; overflow-y:scroll; max-height: 400px; position:relative; 2px 1px 1px rgba(0, 0, 0, 0.05) inset")
														)
												)
										)
								)
						)
					})
			
			x <- seq(0,2000,1)#0:2000
			
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

	
			dataImputed <- reactive({
						ImputeMissing(data[1:1540,], getData()[,colnames(data)])
					})
			models <- c("Ncd","Cr","Rel","Nrd","Prd")
			riskMissing <- reactive({
						sapply(models, function(m){
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
									Z <- if(!is.null(x$na.action)) x$Z[-x$na.action,] else x$Z 
									p <- PartialRisk(x, newZ= rbind(dataImputed(),colMeans(data[1:1540,])))
									p <- p[1,]-p[2,]
									#p <- p[-length(p)]
									c(round(p,3), `total`=round(r[1,1] - mean(Z %*% coef(x)),3),
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
			prsP <- survPredict(Surv(prdData$time1, prdData$time2, prdData$status))(x) # Baseline Prs (measured from relapse)

			coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])) 
			tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=x[-1]))) ## Hazard (function of CR length)	
			
			coxphOs <- coxph(Surv(time1,time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))) 
			tdOsBaseline <- exp(predict(coxphOs, newdata=data.frame(time0=x[-1])))	 ## Hazard (function of induction length), only for OS (could do CIR,NRM,PRS seperately)
						
			# CR adjustments to obtain absolute probabilities
			crAdjust <- function(x, y, time=x$x) {
				xadj <- .crAdjust(x$inc, y$inc, time)
			}
			
			.crAdjust <- function(inc1, inc2, time) {
				cumsum(c(1,diff(inc1) * splinefun(time, inc2)(time[-1])))
			}
			
			computeAbsoluteProbabilities <- reactive({
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
						osDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = osCr, haz1TimeDep = tdOsBaseline) - (1-ncd)
						nrsDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = nrsCr, haz1TimeDep = tdOsBaseline)
						rsDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = rsCr, haz1TimeDep = tdOsBaseline)
						relDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = relCr, haz1TimeDep = tdOsBaseline)
						
						
						## Confidence intervals
						osLoDiag <- osUpDiag <- rep(NA, length(osDiag))
						
						if("analytical" %in% input$ciType){
							PlogP2 <- function(x) {(x * log(x))^2}
							errOsCr <- kmNrd$r[,2] * PlogP2(kmNrd$inc) * (1-(1-kmRel$inc) * (1-kmPrs$inc))^2 + kmRel$r[,2] * PlogP2(kmRel$inc) * (1-kmPrs$inc)^2* kmNrd$inc^2 +  kmPrs$r[,2] * PlogP2(kmPrs$inc) * (1-kmRel$inc)^2* kmNrd$inc^2 
							errOsCr <- sqrt(errOsCr / PlogP2(osCr))
							osUpCr <- osCr ^ exp(2* errOsCr)
							osLoCr <- osCr ^ exp(-2*errOsCr)
							#segments(z, osLo[z+1] ,z,osUp[z+1], col=1, lwd=2)
						}
						if("simulated" %in% input$ciType){
							## Simulate CI
							nSim <- 200
							osCrMc <- sapply(1:nSim, function(i){
										r <- exp(rnorm(5,0,sqrt(c(kmRel$r[,2],kmNrd$r[,2],kmPrs$r[,2], kmNcd$r[,2], kmCr$r[,2]))))
										nrsCr <- .crAdjust(kmNrd$inc^r[2],  kmRel$inc^r[1], time=x) ## Correct KM estimate for competing risk
										diffCir <- diff(kmRel$inc^r[1]) * kmNrd$inc[-1]^r[2] ## Correct KM estimate for competing risk							
										rsCr <- computeHierarchicalSurvival(x = x, diffS0 = diffCir, S1Static = prsP, haz1TimeDep = tdPrmBaseline * exp(kmPrs$r[,1]-kmPrs$r0+log(r[3])))
										osCr <- 1-(1-nrsCr)-(1-rsCr)
										
										cr <- .crAdjust(kmCr$inc^r[5], kmNcd$inc^r[4], time=x)
										ncd <- .crAdjust(kmNcd$inc^r[4], kmCr$inc^r[5], time=x)
										osDiag <- computeHierarchicalSurvival(x = x, diffS0 = diff(cr), S1Static = osCr, haz1TimeDep = tdOsBaseline)
										
										return(cbind(osCr, osDiag - (1-ncd)))
									}, simplify="array")
							osCrMcQ <- apply(osCrMc,1:2,quantile, c(0.025,0.975))
							osLoCr <- osCrMcQ[1,,1]
							osUpCr <- osCrMcQ[2,,1]
							osLoDiag <- osCrMcQ[1,,2]
							osUpDiag <- osCrMcQ[2,,2]
						}
						
						absolutePredictions=data.frame(x=x, cr=cr, ncd=ncd, osDiag=osDiag, nrsDiag=nrsDiag, rsDiag=rsDiag, relDiag=relDiag, osCr=osCr, nrsCr=nrsCr, relCr=relCr, rsCr=rsCr, osUpCr=osUpCr, osLoCr=osLoCr, osLoDiag=osLoDiag, osUpDiag=osUpDiag)
						
						return(absolutePredictions)
						
						
					})
			
			output$multistageCheck <- renderUI({
						if(any(is.na(computeAbsoluteProbabilities()[,1:11])))
							makeWarning(title="Error", message=HTML('An error has occurred in calculating the multistage probabilities. Please check your input values.'), 
									id = 'warningMultistage')
					})
					
			output$multistageDiag <- renderPlot({
						par(bty="n", mar=c(3,3,1,1), mgp=c(2,0.5,0), tcl=-.25, xaxs="i", yaxs="i")
						layout(matrix(1:2, ncol=2), widths=c(1,0.33))
						par(cex=1)
						
						with(computeAbsoluteProbabilities(),{
						
						xLen <- length(x)
						
						scale <- 365.25/12
						xScaled <- x/scale
						
						## Plot probabilities
						plot(xScaled, 1-(1-ncd)-(1-osDiag), type="l", xlab="Months from diagnosis", ylab="Probability", ylim=c(0,1), lwd=3, lty=0) 
						y0 <- 1
						y <- ncd
						polygon(c(xScaled, xScaled[xLen]), c(y,y0), border=NA, col=pastel1[1])
						y0 <- y0 - (1-ncd)
						y <- y - (1-nrsDiag)
						polygon(c(xScaled, rev(xScaled)), c(y, rev(y0))  , border=NA, col=pastel1[2])
						y0 <- y0 - (1-nrsDiag)
						y <- y - (1-rsDiag)
						osDiag <- y
						polygon(c(xScaled, rev(xScaled)), c(y, rev(y0)),  border=NA, col=pastel1[3])
						y0 <- y0 - (1-rsDiag)
						y <- y - (1-relDiag) + (1-rsDiag) 
						polygon(c(xScaled, rev(xScaled)), c(y, rev(y0)),  border=NA, col=pastel1[5])
						polygon(c(xScaled, rev(xScaled)), c(cr - (1-ncd), rev(y)),  border=NA, col=pastel1[4])
						polygon(c(xScaled, rev(xScaled)), c(cr - (1-ncd), rev(rep(0, length(xScaled)))),  border=NA, col="#DDDDDD")
						lines(xScaled, osDiag, lwd=3)

						z <- round(c(365.25,3*365.25))
						y <- (osDiag)[z+1]
						points(z/scale,y, pch=16, col=1)
						text(z/scale, y, labels=round(y,2), pos=1)
						addGrid(scale)
						
						
						## CI
						lines(xScaled, osUpDiag, col=1, lty=2)
						lines(xScaled, osLoDiag, col=1, lty=2)

						## Legend
						par(mar=c(0,0,0,0))
						plot(NA,NA, xlab="",ylab="", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1))
						legend(x=0,y=1, col=c(NA,NA,NA,NA,NA,NA,"black","black"), lty=c(NA,NA,NA,NA,NA,NA,1,4), fill=c(pastel1[c(1,2,3,5,4)],"#DDDDDD",NA,NA), border=c(1,1,1,1,1,1,NA,NA), lwd=c(NA,NA,NA,NA,NA,NA,3,1), y.intersp = 1.5, c("Death without \nremission","Death without \nrelapse","Death after \nrelapse","Alive after \nrelapse","Alive in CR1", "Alive in \ninduction", "Overall survival", "95% C.I."), box.lwd = 0,  bg="#FFFFFF88", seg.len=1)
						
						})
					})
					
					output$multistageCR <- renderPlot({
								par(bty="n", mar=c(3,3,1,1), mgp=c(2,0.5,0), tcl=-.25, xaxs="i", yaxs="i")
								layout(matrix(1:2, ncol=2), widths=c(1,0.33))
								par(cex=1)
								
								with(computeAbsoluteProbabilities(),{
											
											xLen <- length(x)
											
											scale <- 365.25/12
											xScaled <- x/scale
																			
											z <- round(c(365.25,3*365.25))
											
											## Plot outcome after remission
											plot(NA,NA, xlab="Months from remission", ylab="Probability",  xlim=c(0,2000)/scale, ylim=c(0,1), lty=2)
											polygon(c(xScaled, xScaled[xLen]), c(nrsCr,1)  , border=NA, col=pastel1[2])
											polygon(c(xScaled, rev(xScaled)), c(nrsCr, rev(osCr)),  border=NA, col=pastel1[3])
											polygon(c(xScaled, rev(xScaled)), c(osCr, rev(1-(1-nrsCr)-(1-relCr))),  border=NA, col=pastel1[5])
											polygon(c(xScaled, rev(xScaled)), c(1-(1-nrsCr)-(1-relCr), rep(0,length(xScaled))),  border=NA, col=pastel1[4])
											addGrid(scale)
											lines(xScaled, osCr, col=1, lwd=3)
											#title("Outcome after remission")
											
											y <- (osCr)[z+1]
											points(z/scale,y, pch=16, col=1)
											text(z/scale, y, labels=round(y,2), pos=1)
											
											## CI
											lines(xScaled, osUpCr, col=1, lty=2)
											lines(xScaled, osLoCr, col=1, lty=2)
											
											## Legend
											par(mar=c(0,0,0,0))
											plot(NA,NA, xlab="",ylab="", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1))
											legend(x=0,y=1, col=c(NA,NA,NA,NA,"black","black"), lty=c(NA,NA,NA,NA,1,4), fill=c(pastel1[c(2,3,5,4)],NA,NA), border=c(1,1,1,1,NA,NA), lwd=c(NA,NA,NA,NA,3,1), y.intersp = 1.5, c("Death without \nrelapse","Death after \nrelapse","Alive after \nrelapse","Alive in CR1", "Overall survival", "95% C.I."), box.lwd = 0,  bg="#FFFFFF88", seg.len=1)
											
										})
							})
					
					
					printMutations <- function(data) {res <- paste(LABELS[colnames(data)[which(data[1,]==1)]], collapse=", "); if(res=="") NULL else res}
					
					output$patientSummary <- renderText({
								d <- getData()[, colnames(data)]
								#x <- dataImputed()
								bloodVariables <- c("BM_Blasts_100","PB_Blasts_100","wbc_100","LDH_1000","HB_10","platelet_100")
								paste0( "<table><tr><td style='width:20%'><b>Patient:</b></td><td> ", paste(c(na.omit(d[["AOD_10"]]*10), if(is.na(d[["AOD_10"]])) "" else "yr old ", c(`1`="male",`2`='female',`NA`="")[paste(d[["gender"]])]), collapse=""), "</td></tr>",
										"<tr><td><b>Driver mutations:</b></td><td> ", paste(c(printMutations(d[,crGroups %in% "Genetics", drop=FALSE]), 
												printMutations(d[,crGroups %in% "Fusions", drop=FALSE]),
												printMutations(d[,crGroups %in% "CNA", drop=FALSE])), collapse="; "),"<br>",
										"<tr><td><b>Blood counts:</b></td><td> ", paste((d[, bloodVariables] * SCALEFACTORS[bloodVariables])[!is.na(d[, bloodVariables])], sub("(.+) \\((.+)\\).+", "\\2 \\1",LABELS[bloodVariables][!is.na(d[, bloodVariables])], perl=TRUE), collapse=", "), "</td></tr>",
										"<tr><td><b>Treatment:</b></td><td> ", if(!is.na(d[,'transplantRel'])) if(d[,"transplantRel"]) "HSCT after relapse" else if(d[,"transplantCR1"]) "HSCT in CR1" else if(d[,"transplantCR1"]==0 & d[,"transplantRel"]==0) "No HSCT", "</td></tr></table>")
							})
					
					round100 <- function(x){
						y <- floor(x)
						d <- x-y
						o <- order(d, decreasing=TRUE)
						i <- 1
						while(sum(y) < 100){
							y[o[i]] <- y[o[i]] + 1
							i <- i+1
						}
						return(y)
					}
					
					printRisk <- function(x, img="human.svg") paste0(paste0(c("",rep(paste0("<img src=\"",img,"\" alt=\"%\" width=\"8\"></img>"), x)), collapse=""), " ",x, "%")
					#printRisk <- function(x, img="human.svg") paste0(paste0(c("<style type=\"text/css\">.st0{fill:#A9C1D9;}</style>",rep(paste0("<object type=\"image/svg+xml\" data=\"",img,"\" alt=\"%\" width=\"8\"></object>"), x)), collapse=""), " ",x, "%")

					
					output$absoluteRiskDiag <- renderText({
								r <- computeAbsoluteProbabilities()
								r <- r[which.min(abs(round(3*365.25)-r$x)),,drop=FALSE]
								p <- c((1-r[,"ncd"]),(1-r[,"nrsDiag"]),(1-r[,"rsDiag"]),(1-r[,"relDiag"] - (1-r[,"rsDiag"])),(r[,"osDiag"] -(1-r[,"relDiag"] - (1-r[,"rsDiag"])) - (r[,"cr"] - (1-r[,"ncd"]))  ),(r[,"cr"] - (1-r[,"ncd"]) ))
								p <- round100(p*100)
								paste0( '<div class="table-responsive"><table style="width:100%"><tr><td style="width:20%">Death without remission</td><td style="width:80%">', printRisk(p[1], "human-rose.svg"), "</td></tr>",
										"<tr><td>Death without relapse</td><td>", printRisk(p[2], "human-blue.svg"),"</td></tr>",
										"<tr><td>Death after relapse</td><td>", printRisk(p[3], "human-green.svg"),"</td></tr>",
										"<tr><td>Alive after relapse</td><td>", printRisk(p[4], "human-yellow.svg") ,"</td></tr>",
										"<tr><td>Alive in CR1</td><td>", printRisk(p[5], "human-violet.svg"),"</td></tr>",
										"<tr><td>Alive without CR</td><td>", printRisk(p[6], "human-grey.svg"),"</td></tr></table></div>")
#								data.frame(paste("<b>",c("Death without remission","Death without relapse","Death after relapse","Alive after relapse","Alive in CR1","Alive without CR"),"</b>"),
#										sapply(p, printRisk))
							})#, include.colnames=FALSE, include.rownames=FALSE)
					
					output$absoluteRiskCr <- renderText({
								r <- computeAbsoluteProbabilities()
								r <- r[which.min(abs(round(3*365.25)-r$x)),,drop=FALSE]
								p <- c((1-r[,"nrsCr"]),(1-r[,"rsCr"]),(1-r[,"relCr"] - (1-r[,"rsCr"])),(r[,"osCr"] -(1-r[,"relCr"] - (1-r[,"rsCr"])) ))
								p <- round100(p*100)
								paste0( "<table  style='width:100%'><tr><td style='width:20%'>Death without relapse</td><td style='width:80%'>", printRisk(p[1], "human-blue.svg"),"</td></tr>",
										"<tr><td>Death after relapse</td><td>", printRisk(p[2], "human-green.svg"),"</td></tr>",
										"<tr><td>Alive after relapse</td><td>", printRisk(p[3], "human-yellow.svg") ,"</td></tr>",
										"<tr><td>Alive in CR1</td><td>", printRisk(p[4], "human-violet.svg"),"</td></tr></table>")
							})
				
		})
