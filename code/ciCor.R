load("ciCor.RData")

library(mg14)
library(CoxHD)
library(Rcpp)

#save(allDataTpl, coxRFXNrdTD, coxRFXPrdTD, coxRFXRelTD, MultiRFX3, prdData, relData, nrdData, crGroups, nSim, file="../code/ciCor.RData")

cppFunction('NumericVector computeTotalPrsC(NumericVector x, NumericVector diffCir, NumericVector prsP, NumericVector tdPrmBaseline, double risk) {
				int xLen = x.size();
				double hj;
				double r = exp(risk);
				NumericVector rs(xLen);
				for(int i = 0; i < xLen; ++i) rs[i] = 1;
				for(int j = 1; j < xLen; ++j){ 
				hj = tdPrmBaseline[j-1] * r;
				for(int i = j; i < xLen; ++i){
				rs[i] += diffCir[j-1] * (1-pow(prsP[i-j], hj));
				}
				}
				return rs;
				}', rebuild=TRUE)

jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

load(paste0("loo/",jobIndex,".RData"))

cvIdx <- 1:nrow(dataFrame)
whichTrain <- which(cvIdx != jobIndex)

multiRFX3TplCiCorLoo <- sapply(1:nSim, function(foo){
			set.seed(foo)
			cNrd <- rfxNrs
			cNrd$coefficients <- mvtnorm::rmvnorm(1, mean=cNrd$coefficients, sigma=coxRFXNrdTD$var2, method="chol")[1,]
			cRel <- rfxRel
			cRel$coefficients <- mvtnorm::rmvnorm(1, mean=cRel$coefficients, sigma=coxRFXRelTD$var2, method="chol")[1,]
			cPrd <- rfxPrs
			cPrd$coefficients <- mvtnorm::rmvnorm(1, mean=cPrd$coefficients, sigma=coxRFXPrdTD$var2, method="chol")[1,]
			multiRFX3Tpl <- MultiRFX3(cNrd, cRel, cPrd, data=allDataTpl[3*jobIndex + (-2:0),], x=3*365, prdData=prdData[prdData$index %in% whichTrain,])
			multiRFX3Tpl <- matrix(multiRFX3Tpl$os, ncol=3, byrow=TRUE, dimnames=list(NULL, c("None","CR1","Relapse")))
		})

save(multiRFX3TplCiCorLoo, file=paste0("ciCorLoo/",jobIndex,".RData"))

