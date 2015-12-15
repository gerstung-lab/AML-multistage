#' # SF/Chromatin AML
#' ### Load data
library(knitr)
library(RColorBrewer)
col <- c(brewer.pal(9,"Set1")[c(9,1:8)], brewer.pal(8,"Dark2"))
load("/Volumes/mg14/Git/AML/code/AML.RData")
load("/Volumes/mg14/Git/AML/code/dpClass.RData")

#' ### Classes
kable(posteriorMeans)
table(dpClass)

#' Poor man's pairwise interactions within 8 (do properly in future): Some structure
lesions <- names(which(posteriorMeans[,9]>10)) # > 10 or any other cutoff
lesions
heatmap(cor(dataFrame[dpClass==8, lesions]), breaks=seq(-1,1,l=12), col=brewer.pal(11,"RdBu"), scale="none")

#' ### Survival
#' Class 8 has different outcome than other groups
library(survival)
summary(coxph(os ~ dpClass))

#' Heterogenetiy within class 8
summary(coxph(os ~ ., data=dataFrame[lesions]), subset=dpClass=="8")

#' ### Phenotype
#' #### Overall differences among classes
#' x35% less whites than avg ..
boxplot(clinicalData$wbc ~ dpClass, log='y', xlab="Class",ylab="wbc", col=col)
summary(lm(log(clinicalData$wbc) ~ dpClass==8))

#' .. -4.7% Less blasts ..
boxplot(clinicalData$BM_Blasts ~ dpClass, xlab="Class",ylab="Blast %", col=col)
summary(lm(clinicalData$BM_Blasts ~ dpClass==8))

#' .. and 6.16 yrs older than average
boxplot(clinicalData$AOD ~ dpClass, xlab="Class",ylab="Age", col=col)
summary(lm(clinicalData$AOD ~ dpClass==8))

#' Heterogeneity with class 8 - 
summary(lm(log(clinicalData$wbc) ~ ., data=dataFrame[lesions], subset=dpClass==8))
summary(lm(clinicalData$BM_Blasts ~ ., data=dataFrame[lesions], subset=dpClass==8))
summary(lm(clinicalData$AOD ~ ., data=dataFrame[lesions], subset=dpClass==8))


