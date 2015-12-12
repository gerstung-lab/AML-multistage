load("ciCor.RData")
library(mg14)
library(CoxHD)
library(Rcpp)

#save(allDataTpl, coxRFXNrdTD, coxRFXPrdTD, coxRFXRelTD, MultiRFX3, prdData, relData, nrdData, crGroups, file="../code/ciCor.RData")

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

set.seed(jobIndex)

cNrd <- coxRFXNrdTD
cNrd$coefficients <- mvtnorm::rmvnorm(1, mean=cNrd$coefficients, sigma=cNrd$var2, method="chol")[1,]
cRel <- coxRFXRelTD
cRel$coefficients <- mvtnorm::rmvnorm(1, mean=cRel$coefficients, sigma=cRel$var2, method="chol")[1,]
cPrd <- coxRFXPrdTD
cPrd$coefficients <- mvtnorm::rmvnorm(1, mean=cPrd$coefficients, sigma=cPrd$var2, method="chol")[1,]
multiRFX3Tpl <- MultiRFX3(cNrd, cRel, cPrd, data=allDataTpl, x=3*365, prdData=prdData)
multiRFX3Tpl <- matrix(multiRFX3Tpl$os, ncol=3, byrow=TRUE, dimnames=list(NULL, c("None","CR1","Relapse")))


save(multiRFX3Tpl, file=paste0("ciCor/",jobIndex,".RData"))

