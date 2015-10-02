#save(genes, file="genesImputation.RData")

load("loo.RData")
load("genesImputation.RData")
library(mg14)
library(CoxHD)
library(Rcpp)

i <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

cvIdx <- 1:nrow(dataFrame)
whichTrain <- which(cvIdx != i)

e <- new.env()
t <- try(load(paste0("loo/",i,".RData"), env=e))
if(class(t)=="try-error"){
	stop() 
}else{
	whichTrain <<- (1:nrow(data))[-i]
	dMiss <- do.call("rbind", lapply(c(0,seq_along(genes)), function(g){
				na.genes <- if(g==0) genes else genes[-(1:g)]
				if(length(na.genes)==0) na.genes <- "FOO42"
				d <- data[i,, drop=FALSE]
				d[grepl(paste(na.genes, collapse="|"), colnames(d))] <- NA
				d
			}))
	
	xx <- 0:2000
	coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])[prdData$index %in% whichTrain,]) 
	tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
	
	coxphOs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))[osData$index %in% whichTrain,]) 
	tdOsBaseline <- exp(pmin(predict(coxphOs, newdata=data.frame(time0=500)),predict(coxphOs, newdata=data.frame(time0=xx[-1])))) ## cap predictions at induction length 500 days.
	
	multiRfx5Imputed <- MultiRFX5(e$rfxEs, e$rfxCr, e$rfxNrs, e$rfxRel, e$rfxPrs, dMiss, tdPrmBaseline = tdPrmBaseline, tdOsBaseline = tdOsBaseline, x=2000)
	save(multiRfx5Imputed, file=paste0("imputed/",i,".RData"))
}	