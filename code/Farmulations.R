load("simulationData.RData")
library(survival)
library(parallel)
interactionEffects <- seq(-2,2,0.25)
nData <- c(100, 1000, 10000)
treatmentEffects <- seq(0,-2,-0.25)
nSim <- 100
nJobs <- as.numeric(Sys.getenv("NJOBS"))

jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

X <- expand.grid(nData=nData, gene=colnames(simDataFrame)[which(groups=="Genetics")], treatmentEffects=treatmentEffects, interactionEffects=interactionEffects, nSim=1:nSim)

#+ eval=FALSE
#set.seed(seed)
#simResults <- lapply(nData, 
#		function(n) mclapply(colnames(simDataFrame)[which(groups=="Genetics")], 
#					function(g) lapply(treatmentEffects,
#								function(t) sapply(interactionEffects, 
#											function(e) sapply(1:nSim, 
#														function(i) simulate(simDataFrame[1:n,whichRFX], coxRFXFit, g, treatEffect = t, geneTreatEffect = e)[1,1]))), mc.cores=1))
indeces <- which(1:nrow(X) %% nJobs == jobIndex %% nJobs)

warn <- pWald <- pLR <- coefHat <- seeds <-  numeric(length(indeces))
for(i in seq_along(indeces)){
	seed <- round(runif(1) * 1e8)
	set.seed(seed)
	x <- X[indeces[i],]
	seeds[i] <- seed
	whichMain <- colnames(simDataFrame)[groups %in% c("Clinical","Genetics","Cytogenetics","Treatment")]
	sim <- SimAndTest(simDataFrame[1:x$nData,whichRFX], coxRFXFit, whichGene = as.character(x$gene), treatEffect = x$treatmentEffects, geneTreatEffect = x$interactionEffects, whichMain=whichMain)
	pWald[i] <- sim[1,1]
	pLR[i] <- sim[1,2]
	coefHat[i] <- sim[1,3]
	warn[i] <- sim[1,4]
}
simResults <- data.frame(pWald, pLR, coefHat, warn, indeces, seeds)
save(simResults, file=paste("output/",Sys.getenv("LSB_JOBNAME"),"_",Sys.getenv("LSB_JOBINDEX"),".RData", sep=""))