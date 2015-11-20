#' ---
#' title: Supplementary Methods
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     number_sections: true
#'     auto_identifiers: true
#'     table_captions: true
#'     includes:
#'       in_header: mathjax.html
#' author: Moritz Gerstung
#' bibliography: lit.bib
#' 
#' ---
#' 
#' # Data
#' 
#' 
#' ## Variables
#' 
#' We use the data from N=1,540 AML cases as described in our companion paper [@PapaemmanuilSM2015]. 
#' 
#' These can be summarised as follows:
#' 
#' 
#' |Group        |   p|Variables                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
#' |:------------|---:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
#' |Fusion genes           |   8|t_MLL, inv3_t3_3, t_9_22, t_15_17, t_8_21, inv16_t16_16, t_6_9, t_9_11, t_v_11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
#' |CNA          |  18|minus5_5q, minus7, minus7q, abn7other, plus8_8q, minus9q, mono12_12p_abn12p, plus13, mono17_17p_abn17p, minus18_18q, minus20_20q, plus21, plus22, minusY, abn3q_other, plus11_11q, mono4_4q_abn4q, complex                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
#' |Genetics     |  58|ASXL1, ATRX, BCOR, BRAF, CBL, CBLB, CDKN2A, CREBBP, CUX1, DNMT3A, EP300, ETV6, EZH2, FBXW7, GATA2, GNAS, IDH1, IKZF1, JAK2, KDM5A, KDM6A, KIT, KRAS, MLL, MLL2, MLL3, MLL5, MPL, MYC, NF1, NPM1, NRAS, PHF6, PRPF40B, PTEN, PTPN11, RAD21, RB1, RUNX1, SF1, SF3A1, SF3B1, SFRS2, SH2B3, STAG2, TET2, TP53, U2AF1, U2AF2, WT1, ZRSR2, CEBPA_mono, CEBPA_bi, FLT3_ITD, FLT3_TKD, FLT3_other, IDH2_p172, IDH2_p140                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
#' |Gene:Gene product terms     | 126|BCOR:DNMT3A, ASXL1:EZH2, DNMT3A:IDH1, DNMT3A:KRAS, DNMT3A:MLL, IDH1:MLL, DNMT3A:MYC, DNMT3A:NF1, CBL:NPM1, DNMT3A:NPM1, GATA2:NPM1, IDH1:NPM1, KIT:NPM1, KRAS:NPM1, MYC:NPM1, NF1:NPM1, ASXL1:NRAS, BCOR:NRAS, DNMT3A:NRAS, EZH2:NRAS, GATA2:NRAS, IDH1:NRAS, KIT:NRAS, KRAS:NRAS, MLL:NRAS, NPM1:NRAS, NPM1:PHF6, DNMT3A:PTPN11, IDH1:PTPN11, KRAS:PTPN11, NPM1:PTPN11, NRAS:PTPN11, DNMT3A:RAD21, NPM1:RAD21, NRAS:RAD21, PTPN11:RAD21, ASXL1:RUNX1, BCOR:RUNX1, DNMT3A:RUNX1, EZH2:RUNX1, IDH1:RUNX1, MLL:RUNX1, NRAS:RUNX1, PHF6:RUNX1, NRAS:SF3B1, ASXL1:SFRS2, DNMT3A:SFRS2, IDH1:SFRS2, NPM1:SFRS2, NRAS:SFRS2, RUNX1:SFRS2, ASXL1:STAG2, DNMT3A:STAG2, EZH2:STAG2, MLL:STAG2, NPM1:STAG2, NRAS:STAG2, RUNX1:STAG2, SFRS2:STAG2, ASXL1:TET2, DNMT3A:TET2, KIT:TET2, MLL:TET2, NPM1:TET2, NRAS:TET2, PTPN11:TET2, RUNX1:TET2, SFRS2:TET2, STAG2:TET2, DNMT3A:TP53, NRAS:TP53, NRAS:U2AF1, NPM1:WT1, NRAS:WT1, DNMT3A:CEBPA_mono, NPM1:CEBPA_mono, TET2:CEBPA_mono, GATA2:CEBPA_bi, NRAS:CEBPA_bi, WT1:CEBPA_bi, DNMT3A:FLT3_ITD, EZH2:FLT3_ITD, IDH1:FLT3_ITD, MLL:FLT3_ITD, MYC:FLT3_ITD, NPM1:FLT3_ITD, NRAS:FLT3_ITD, PHF6:FLT3_ITD, PTPN11:FLT3_ITD, RAD21:FLT3_ITD, RUNX1:FLT3_ITD, STAG2:FLT3_ITD, TET2:FLT3_ITD, WT1:FLT3_ITD, CEBPA_mono:FLT3_ITD, CEBPA_bi:FLT3_ITD, DNMT3A:FLT3_TKD, IDH1:FLT3_TKD, MLL:FLT3_TKD, NPM1:FLT3_TKD, NRAS:FLT3_TKD, RAD21:FLT3_TKD, RUNX1:FLT3_TKD, TET2:FLT3_TKD, WT1:FLT3_TKD, FLT3_ITD:FLT3_TKD, DNMT3A:FLT3_other, NPM1:FLT3_other, NRAS:FLT3_other, PTPN11:FLT3_other, RUNX1:FLT3_other, TET2:FLT3_other, FLT3_ITD:FLT3_other, DNMT3A:IDH2_p172, ASXL1:IDH2_p140, DNMT3A:IDH2_p140, MLL:IDH2_p140, NPM1:IDH2_p140, NRAS:IDH2_p140, PTPN11:IDH2_p140, RUNX1:IDH2_p140, SFRS2:IDH2_p140, STAG2:IDH2_p140, FLT3_ITD:IDH2_p140, FLT3_TKD:IDH2_p140, NPM1:FLT3_ITD:DNMT3A |
#' |Clinical     |  11|Performance_ECOG, BM_Blasts_100, PB_Blasts_100, wbc_100, LDH_1000, HB_10, platelet_100, Splenomegaly, oAML, sAML, tAML                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
#' |Demographics |   2|AOD_10, gender                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
#' |Treatment    |   4|ATRA, VPA, TPL_os, TPL_rel                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
#' |Nuisance     |   4|AMLHD98A, AMLHD98B, Date_1000, MissingCyto                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
#' 
#' Table: Table 1. Variables 
#' 
#' ## Preprocessing
#' 
#' The following preprocessing steps were applied:
#'  
#' * Agnostic imputation of missing variables by mean.
#' * Quantitative variables linearly rescaled by a power of 10 to a magnitude of 1. This is necessary, as we are working with a penalty on the coefficient size. 
#' 
#' ### Genomic variables
#' Fusion genes, Copy number alterations, Genetics and Gene:Gene interactions were encoded as 0 (absent) and 1 (present) based on the same annotation as in [@PapaemmanuilSM2015]. 
#' A gene was considered mutated and encoded as 1 if it contained at least one oncogenic mutation, and 0 otherwise. In addition, we used the following rules:
#' 
#' * For _CEBPA_ we additionally differentiated mono- and bi-allelic lesions. 
#' * For _FLT3_ we distinguished between ITD, TKD and other mutations. 
#' * For _IDH2_ we separately encoded  p172 and p140 point mutations.
#' 
#' ### Gene:gene product terms
#' 
#' Gene:Gene product terms were computed indicating whether a combination of two genes was present. This allows to
#' account for non-additive genetic interaction. To limit the number of variables, product terms included if there were at least 8 occurrences. 
#' 
#' 
#' ### Clinical and demographic variables
#' 
#' Quantitative clinical variables were rescaled to a magnitude of 1 as described above.  To assess the validity of out log-linear risk model
#' we computed spline fits, that allow for a non-linear dependence between log-hazard and each variable. We did not observe a measurable improvement 
#' of our model fits in cross-validation. 
#' 
#' ### Treatment
#' 
#' The following variables were included in the model:
#' 
#' 1.	Allograft (MRD, MUD) in CR1 as a time-dependent covariate. 
#' Not considered for RFS, CPSS (no time-dependence allowed)
#' 2. Allograft (MRD, MUD) after relapse as a time-dependent covariate. 
#' Multi-stage model only.
#' 3.	Extra cycles of ATRA encoded as 0/1.
#' 4.	VPA (AMLSG 07/04 only).
#' 
#' ### Nuisance
#' 
#' In addition to the aforementioned explanatory variables we used the following multiplicative strata to
#' account for potential confounding factors:
#' 
#' 1.	Missing cytogenetic information (0/1). We observed that cytogenetic data was missing more frequently for patients having died early. To avoid a negative
#' bias we included this as an additional factor 
#' 2.	Trial (2/3), AMLSG07/04, AMLHD98A, AMLHD98B. A factor was included to account for systematic differences between trials.
#' 3.	Date. The date of diagnosis was included to account for and improvement of patient care over time.
#' 
#' 
#' ## Code

#+ Preliminaries, echo=FALSE
options(width=120)
pdf.options(pointsize=8)
library(knitr)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' ### Libraries
#' Load a few libraries - see end of document for a full list of libraries and their versions.
library(CoxHD)
library(mg14)
set1 <- brewer.pal(8, "Set1")

#' ### Raw data
#' Load clinical data
#+ clinicalData, cache=TRUE
clinicalData <- read.table("../../data/Ulm1.17_MG_Clinical.txt", sep="\t", header=TRUE, na.strings = "na", comment.char = "", quote="\"")
clinicalData <- clinicalData[order(clinicalData$PDID),]
clinicalData$ERDate <- as.Date(as.character(clinicalData$ERDate), "%d-%b-%y")
clinicalData$CR_date <- as.Date(as.character(clinicalData$CR_date), "%d-%b-%y")
clinicalData$TPL_date <- as.Date(as.character(clinicalData$TPL_date), "%d-%b-%y")
clinicalData$Date_LF <- as.Date(as.character(clinicalData$Date_LF), "%d-%b-%y")
clinicalData$Recurrence_date <- as.Date(as.character(clinicalData$Recurrence_date), "%d-%b-%y")
levels(clinicalData$Study) <- c(`_07-04`="AMLSG0704" ,   `98A`="AMLHD98A" ,  `98B`="AMLHD98B")[levels(clinicalData$Study)]
clinicalData$VPA[is.na(clinicalData$VPA)] <- 0
clinicalData$ATRA_arm[is.na(clinicalData$ATRA_arm)] <- 0
colnames(clinicalData) <- gsub('\\.',"",colnames(clinicalData))
clinicalData <- clinicalData[!is.na(clinicalData$TypeAML),] ## remove unknown patients
clinicalData$PDID <- factor(as.character(clinicalData$PDID))
t <- read.table("../../data/Ulm1.9_sm_Clinical.txt", header=T, sep="\t", na.strings = "na",comment.char = "", quote="\"")
karyotypes <- t$karyotype[match(clinicalData$PDID,t$PDID)]
rm(t)
clinicalData$t_9_11 <- grepl("t\\(9;11\\)\\(p22;q23\\)", karyotypes) + 0 # t(9;11)
clinicalData$t_v_11 <- clinicalData$t_MLL &! clinicalData$t_9_11
clinicalData$t_MLL <- NULL

dim(clinicalData)

#' Load mutation data
mutationData = read.table("../../data/Ulm1.14_MG_Genetic.txt", sep="\t", header=TRUE, strip.white = TRUE)
mutationData$SAMPLE_NAME <- factor(as.character(mutationData$SAMPLE_NAME), levels = levels(clinicalData$PDID)) ## Refactor
mutationTable <- (table(mutationData[mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & mutationData$FINAL_CALL == "OK" ,c("SAMPLE_NAME","GENE")]) > 0)+0
dim(mutationTable)

all(rownames(mutationTable)==clinicalData$PDID)

#' ### Survival data
#+ survival, cache=TRUE
os <- Surv(clinicalData$OS, clinicalData$Status) #OS
t <- clinicalData$Time_Diag_TPL
t[is.na(t) | !clinicalData$TPL_Phase %in% "CR1" | !clinicalData$TPL_type %in% c("ALLO","FREMD") ] <- Inf ## Only allografts in CR1
o <- clinicalData$OS
tplIndexOs <-  t < o
osTD <-  Surv(time = rep(0, nrow(clinicalData)), time2=pmin(o, t), event=ifelse(tplIndexOs, 0, clinicalData$Status) )
osTD <- rbind(osTD, 
		Surv(time=t[which(tplIndexOs)],
				time2=o[which(tplIndexOs)], 
				event=clinicalData$Status[which(tplIndexOs)])
)
osTD = Surv(osTD[,1],osTD[,2],osTD[,3])
rm(o,t)
tplSplitOs <- c(1:nrow(clinicalData), which(tplIndexOs))
osYr <- os
osYr[,1] <- osYr[,1]/365
osYrTD <- osTD
osYrTD[,1] <- osYrTD[,1]/365

#' ### Covariates
#' All data as list
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
		Cytogenetics = clinicalData[,grep("^(t_)|(inv)|(abn)|(plus)|(minus)|(mono)|(complex)",colnames(clinicalData))],
		Nuisance = data.frame( MakeInteger(clinicalData$Study)[,1:2], Date=scale(as.numeric(clinicalData$ERDate), scale=FALSE), MissingCyto=is.na(clinicalData$t_15_17)+0),
		Treatment = data.frame(ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL_os=tplIndexOs),
		Demographics = clinicalData[,c("AOD","gender")],
		Clinical = cbind(clinicalData[, c("Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH","HB","platelet","Splenomegaly")], MakeInteger(clinicalData$TypeAML)[,-1]))#,
#MolRisk = makeInteger(clinicalData$M_Risk))
#dataList$Genetics$CEBPA <-  clinicalData$CEBPA # encoded as 0,1,2
dataList$Genetics$CEBPA_mono <-  clinicalData$CEBPA == 1 # encoded as 0,1,2
dataList$Genetics$CEBPA_bi <-  clinicalData$CEBPA == 2 # encoded as 0,1,2
dataList$Genetics$CEBPA <- NULL
dataList$Genetics$FLT3 <- NULL
dataList$Genetics$FLT3_ITD <- clinicalData$FLT3_ITD != "0"
dataList$Genetics$FLT3_TKD <- clinicalData$FLT3_TKD != "0"
dataList$Genetics$FLT3_other <- clinicalData$FLT3_other != "0"
dataList$Genetics$IDH2_p172 <- table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("172", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2_p140 <-  table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("140", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2 <- NULL
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Cytogenetics$MLL_PTD <- NULL
dataList$Genetics = dataList$Genetics + 0
dataList$GeneGene <- MakeInteractions(data.frame(dataList$Genetics), data.frame(dataList$Genetics))[,as.vector(upper.tri(matrix(0,ncol=ncol(dataList$Genetics), nrow=ncol(dataList$Genetics))))]
dataList$GeneGene <- dataList$GeneGene[,colSums(dataList$GeneGene, na.rm=TRUE)>0] 
dataList$GeneGene$`NPM1:FLT3_ITD:DNMT3A` <- (rowSums(dataList$Genetics[c('NPM1',"FLT3_ITD","DNMT3A")])==3)+0 ## Add NPM1:FLT3_ITD:DNMT3A product term as well
dataList$CytoCyto <- MakeInteractions(dataList$Cytogenetics, dataList$Cytogenetics)[,sapply(1:ncol(dataList$Cytogenetics), `<`, 1:ncol(dataList$Cytogenetics))]
dataList$CytoCyto <- dataList$CytoCyto[, colSums(dataList$CytoCyto, na.rm=TRUE) > 0]
dataList$GeneCyto <- MakeInteractions(dataList$Genetics, dataList$Cytogenetics)
dataList$GeneCyto <- dataList$GeneCyto[,colSums(dataList$GeneCyto, na.rm=TRUE) > 0]
dataList$GeneTreat <- MakeInteractions(dataList$Genetics, dataList$Treatment)
dataList$GeneTreat <- dataList$GeneTreat[,colSums(dataList$GeneTreat, na.rm=TRUE) > 0]
dataList$CytoTreat <- MakeInteractions(dataList$Cytogenetics, dataList$Treatment)
dataList$CytoTreat <- dataList$CytoTreat[,colSums(dataList$CytoTreat, na.rm=TRUE) > 0]

#' Condensing to a data.frame
dataRaw <- do.call(cbind,dataList)
names(dataRaw) <- unlist(sapply(dataList, names))
dataFrame <- StandardizeMagnitude(dataRaw)
dim(dataFrame)

#+ groups, cache=TRUE
groups <- unlist(sapply(names(dataList), function(x) rep(x, ncol(dataList[[x]]))))
groups[grepl("^(t_)|(inv)", colnames(dataFrame)) &! grepl(":", colnames(dataFrame))] <- "Fusions"
groups[groups=="Cytogenetics"] <- "CNA"
groups <- factor(groups)
names(groups) <- colnames(dataFrame)
table(groups)

#' Poor man's imputation by column means
#+ dataFrame, cache=TRUE
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
dataFrame <- as.data.frame(sapply(dataFrame, poorMansImpute))
rownames(dataFrame) <- clinicalData$PDID



#' 
#' 
#' 
#' # Models for overall survival
#' 
#' 
#' We use overall survival, measured from date of diagnosis, as the endpoint.
#' 
#' ## Random effects modelling
#' 
#' We implemented sparse random effects for the Cox proportional hazards model in the `CoxHD` R package available at [http://github.com/mg14/CoxHD].
#' This implementation can handle constant covariate and time-dependent models. The latter is important to quantify the effects of allografts,
#' which are typically administered well after diagnosis.
#' `CoxHD::CoxRFX()`
#' 
#' Let the hazard be:
#' $$ \lambda = \lambda_0(t) \exp(u^T Z)$$
#' 
#' Define $h=u^T Z$ as the log hazard. $\lambda_0(t)$ is the normal baseline hazard in a `coxph` model.
#' 
#' The random effects model used here is an example of a hierarchical model with an additional assumption about the
#' distribution of the parameters $u$. We assume that these follow a normal distributions. This additional assumption
#' leads to a ridge-type regularisation of the log-likelihood.
#' 
#' Let there be $p$ covariates and $\{g\}$ be a partitioning of the $p$ variables into $|g|$ groups. For each group assume
#' that the parameters $u_j$ are iid Normally distributed in each group:
#' 
#' $$\forall j \in g: u_j \sim \operatorname{N}(\mu_g;\sigma^2_g) \qquad iid. $$
#' 
#' The shared means $\mu_g$ are motivated by the observation that the effect of oncogenic lesions is, on average, deleterious.
#' 
#' We use the convention that variables without indexes refer to the set of variables. In particular $u = \{u_j: j=1,...,p\}$, $u_g=\{u_j: j\in g\}$.
#' 
#' The full logarithmic likelihood reads:
#' 
#' $$\begin{align} \ell(u,\sigma^2,\mu;Z) &= \ell_0(u;Z) - \sum_g \frac{\sum_{j\in g}(u_j-\mu_g)^2}{\sigma^2_g} \cr &=  \ell_0(u;Z) + \ell_2(u,\mu,\sigma^2). \end{align}$$
#' 
#' The term $\ell_0(u)$ is the likelihood of an unpenalised coxph model. The second term is a sum of ridge penalties resulting
#' from the constraints imposed by the normal distribution of $u$, which penalises large values of $u_j-\mu_g$ with strength $1/\sigma_g$.
#' 
#' Note that the likelihood can be reparametrised by introducing the auxiliary variables $z_g = \sum_{j\in g} Z_{.j}$ and the centred
#' effects $u_j = u_j - \mu_g$:
#' 
#' $$\ell(u,\sigma^2,\mu;Z) = \ell_0(u,\mu;Z,z) + \ell_2(u,\sigma^2) =: \ell(u, \mu,\sigma^2;Z)$$
#' 
#' ### Implementation
#' All of the following steps are implemented in the `CoxHD` `R` package, available
#' at [http://www.github.com/mg14/CoxHD]
#' It can be installed using the `devtools::install_github("mg14/CoxHD/CoxHD")`. 
#' The implementation makes heavy use of the `survival` package [@Therneau2014]. The implementation is about 100x faster
#' than the `coxme` `R` package for mixed effects Cox models by @Therneau2012, as it exploits that $u$ are iid.
#' 
#' ### Parameter estimation
#' We use an EM algorithm as suggested by @PerperoglouSM2014 for Cox models, based on the work by @SchallB1991. 
#' The algorithm iteratively estimates the following quantities:
#' 
#' 1.  Given $\hat\sigma^2$, jointly estimate 
#' 
#' 	1.1. the **shared means** $\hat\mu_g$ as the effect of the auxiliary variables $z_g$.
#' 	
#' 	1.2. the **centred variables** $\hat u$ as a ridge estimate,
#' 	$$\hat\mu, \hat u = \arg \max \ell(u, \mu,\hat\sigma^2;Z)$$
#' 
#' 2.  Given $\hat\mu$ and $\hat u$ the **variances** are estimated as:
#' 	$$\hat\sigma_g^2 = \sum_{j\in g}\hat u_j^2/df_g, \qquad df_g = \operatorname{tr} [\mathcal{I_gg} H_{gg}^{-1}],$$
#'     where $H$ is the Hessian matrix of the penalised model 
#' 	and $\mathcal{I}$ the observed Fisher information of the unpenalised model (each evaluated for variables of group $g$).
#' 
#' Iterate until convergence of parameters and penalised likelihood.
#' 
#' The final parameters are given by uncentering $\hat{u}_j = \hat u_j + \hat \mu_g$.
#' 
#' **Note**: There estimates $\hat{u}$ are maximum a posteriori (MAP) from a Bayesian interpretation with
#'  $\hat\sigma^2$ and $\hat\mu$ being empirical Bayes estimates. 
#' 
#' 
#' ### Semiparametric bootstrap
#' 
#' To assess the sampling distributions of our estimates, e.g., to assess the their variances, we use the following semi-parametric bootstrap approach:
#' 
#' For i=1:100 simulate $n$ semiparametric survial times `Y` (see [#survival]):
#' 	
#' * Using MAP estimates $\hat{u}$
#' * Using full covariate set $Z$
#' 	
#' This allows to assess the distribution of all estimates in a semi-parametric way.
#' 
#' ### Analytical confidence intervals of individual parameters
#' Two estimates exists for the covariance matrices of the parameters $\hat u$ and $\mu$[@TherneauJCAGS2003]:
#' 
#' 1. $\hat V_1 = H^{-1}$ 
#' 
#' 2.  $\hat V_2 = H^{-1} \mathcal{I} H^{-1}$, where $H$ is the Hessian matrix of the penalised model 
#' 	and $\mathcal{I}$ the observed Fisher information of the unpenalised model. Semi-parametric bootstrap simulations show that $V_2$ is more accurate in our context.
#' 
#' The estimates $\hat V_\cdot$ have dimension $(p+|g|)\times(p+|g|)$.
#' 
#' The uncentered variance estimates of the parameter $u_j = u_j + \mu_g$ are given by 
#' 
#' $$\hat V_\cdot[u_j] = \hat V_\cdot[u_j,u_j] + \hat V_\cdot[\mu_g,\mu_g] + 2 \hat V_\cdot[u_j,\mu_g],$$
#' 
#' thus accounting for the correlation of $u_j$ and $\mu_g$.
#' 
#' 
#' ### Wald test of individual parameters
#' 
#' Using variance estimate $\hat V_2$, allows for computing a Wald-type test with one degree of freedom. 
#' 
#' $$
#' \begin{align}
#' z &= \hat{u}^2 / \hat V_2[u] \\\\
#' Z &\sim \chi^2_1
#' \end{align} 
#' $$
#' 
#' This is implemented as `CoxHD::WaldTest()`
#' 
#' P-values of each test are corrected for multiple testing. Due to dependence imposed 
#' by the shared distribution we use the Benjamini-Yektuieli method for controlling
#' the false discovery rate (Q < FDR), implemented as `p.adjust(x,method="BY")`.
#' 
#' **Note**: There exists a lively debate about how, and if at all, random effects
#' shall be tested or not, see [http://glmm.wikidot.com/faq] or [https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q2/000743.html]. 
#' Here we use an approach outlined by @GrayJTASA1992, @TherneauJCAGS2003 and @WoodB2013.
#' However, it is important to check that the variances are correctly specified using
#' a parametric bootstrap approach. 
#' 
#' ### Variance components
#' 
#' #### Partial log hazard
#' In an additive model the linear predictor of the log hazard $h$ is given by:
#' 
#' $$h = u^T Z = \sum_g \sum_{j\in g} u_j^T Z_{.j} = \sum_g h_g$$
#' 
#' Where the set of ${g}$ is partitioning of all covariates. We define $h_g$ as the partial logarithmic hazard contributed by group $g$.
#' 
#' #### Variance components
#' The variance of the logarithmic hazard is given by:
#' 
#' $$Var[h] = \sum_{g,h} Cov(h_g,h_h)$$
#' 
#' Taking just the diagonal elements of $Cov(h_g,h_h)$ guarantees positive values, which do not necessarily add to the total variance.
#' Using $V_g = \sum_h Cov(h_h,h_g)$ yields additive variance components, albeit at the cost of being negative in cases with strong 
#' collinearity of the components.
#' 
#' Variance components are implemented as `CoxHD::VarianceComponents()`.
#' 
#' **Note**: Unlike a classical mixed model $V$ is not computed by marginalising the random effects, but by the MAP estimates. This
#' can be seen as a first order approximation.
#' 
#' The standard deviation $\sqrt{Var[h]}$ determines the average difference between any two patients in logarithmic hazard.
#' 
#' #### Relation to concordance
#' For a normally distributed hazard, the variance $\sigma^2_h$ of the log hazard is related to the concordance metric [@GonenB2005]
#' 
#' $$C = \int \frac{1}{1+\exp(-|x|)}f(x;0,\sigma^2_h) dx$$  
#' 
#' where $f(x,\mu,\sigma_2)$ is the density of normal distribution. There exists no analytical 
#' solution to the above equation, but it may be computed numerically. For a variance of 1, the concordance is 72.5%.
#' 
#' **Note**: For a Cox proportional hazards model even perfect knowledge of the hazard does not guaratuee perfect concordance (i.e. C=1)
#' due to the sampling of the survival times. The limit $Var[h] \rightarrow \infty$, in which the hazard ratio between any two patients is infinite, yields a 
#' deterministic behaviour with $C=1$.
#' 
#' ### Prediction error
#' 
#' The prediction error of a the log hazard for patient $i$ is given by
#' 
#' $$\hat V[h_i] = V[\hat{u}^T Z_{i\cdot}] = Z_{i\cdot}^T \hat V[u] Z_{i\cdot}$$
#' 
#' where $V[\hat{u}]$ is the covariance matrix of the parameters defined in [Analytical confidence intervals](#analytical-confidence-intervals-of-individual-parameters).
#' 
#' **Note**:  In a linear model, the lhs corresponds to the the residual $r_i$ of observation $i$ and
#' the identity $\hat V = Z^T Z \times RSS/n$ holds. In our case $V$ is derived from the Fisher information,
#' but it can be intuitive to think about the average prediction error $\sum_i \hat V[h_i]/n$ as a
#' pseudo residual variance.
#' 
#' ### Covariance-based imputation
#' 
#' To predict the log-hazard in the presence of missing variables, we can use the following
#' imputation, leveraging the covariance in the training set:
#' 
#' Suppose that $Z = (Z_o,Z_m)$, where $Z_o$ are observed and $Z_m$ missing parts of 
#' the data set. Suppose we know the means $\mu$ and covariance $\Sigma$. Then
#' 
#' $$E[Z_m] = \mu_m + \Sigma_{m,o} \Sigma_{o,o}^{-1} (Z_o - \mu_o)$$
#' 
#' $$V[Z_m] = \Sigma_{mm} - \Sigma_{mo} \Sigma_{oo}^{-1} \Sigma_{om} $$
#' 
#' #### Prediction error
#' 
#' The uncertainty in $Z_m$ adds another term to the prediction error:
#' 
#' $$\hat V[h_i] = Z_{io}^T \hat V[u]_{oo} Z_{io} + u_m^T V[Z_{im}] u_m$$ 
#' 
#' ## Other survival models
#' 
#' ### Stepwise variable selection
#' 
#' Coxph + AIC or BIC forward and backward selection beginning from empty model.
#' The implementation in the `survival` `R` package [@Therneau2014] handles constant covariate and time-dependent models.
#' 
#' 
#' ### Complementary pairs stability selection
#' 
#' Complementary pairs stability selection (CPSS) is an extension of the stability selection protocol, which combines
#' subsampling and LASSO-regularised regression to obtain a robust subset of predictor variables [@MeinshausenJOTRSSSBSM2010]. 
#' Using complimentary pairs subsamples @ShahJTRSSSSM2013 derived 
#' a tighter bound for error control.
#' 
#' We have recently used CPSS to analyse the association of genomic predictors and outcome in Myelodysplastic syndromes [@PapaemmanuilB2013].
#' To this end, we have implemented CPSS in the `CoxHD` `R` package. Our implementation fits the CPSS model using the `glmnet` algorithm [@FriedmanSS2010, @SimonJSS2011]. 
#' 
#' The algorithm `CoxHD::CoxCPSS()` uses the following parameters:
#' 
#' * Parameters
#'     - 50 pairs
#'     - Selection probability 80%
#'     - Penalty range chosen to conform FDR < 10%
#' * Refit `coxph()` with selected variables for predictions
#' 
#' Note that the `glmnet` algorithm cannot handle time-dependent covariates.
#' 
#' ### Random survival forests
#' 
#' Random survival forests are an intrinsically non-linear alternative to Cox proportional hazards based regression [@IshwaranTAAS2008]. The idea is to
#' fit an ensemble of regression trees based on subsampling of patients and/or covariates. The resulting predictions are averaged across
#' the forest of regression trees. We used version 1.6 of the `randomForestSRC` package and default options for `randomForestSRC::rfsrc()`. 
#' Note that the model can only handle constant covariates.
#' 
#' ## Code
#' ### Number of oncogenic mutations

#' Construct data.frame for OS, replicating patients (rows) before and after allograft.
#+ dataFrameOsTD, cache=TRUE
dataFrameOsTD <- dataFrame[tplSplitOs,]
dataFrameOsTD[which(tplIndexOs), grep("TPL", colnames(dataFrameOsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero 

#' Define some indexes relating to subsets of variables used by the random effects model.
#+ indeces, cache=TRUE
mainGroups <- grep("[A-Z][a-z]+[A-Z]",levels(groups), invert=TRUE, value=TRUE)
mainGroups
mainIdx <- groups %in% mainGroups
osIdx <- !grepl("TPL", colnames(dataFrame)) ## Exclude TPL from OS analyses..
whichRFXOs <- which((colSums(dataFrame)>=8 | mainIdx) & osIdx) # ie, > 0.5%
mainIdxOs <- mainIdx & osIdx
osTDIdx <- !grepl("TPL_efs", colnames(dataFrame))
whichRFXOsTD <- which((colSums(dataFrame)>=8 | mainIdx) & osTDIdx) # ie, > 0.5%
mainIdxOsTD <- mainIdx & osTDIdx
whichRFXOsGG <- which((colSums(dataFrame)>=8 | mainIdxOs) & osIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%

#' Compute the number of oncogenics, excluding complex karyotype.
NONC <- rowSums(cbind(dataList$Cytogenetics[names(dataList$Cytogenetics)!="complex"], dataList$Genetics), na.rm=TRUE)

#' #### Number of oncogenics
#+ NONCs, fig.width=2.5, fig.height=2
NONCs <- factor(ceiling(pmin(NONC,7)/2), labels=c("0","1-2","3-4","5-6","7+"))
c <- set1[c(3,2,4,1,5)]
f <- survfit(osYr ~ NONCs)
s <- summary(f)
plot(f, col=c, xlim=c(0,10),xlab="Years", ylab="Survival",mark="|", cex=.5)
legend('topright', bty='n', col=c, legend=paste0(levels(NONCs)," (n=",table(NONCs),")"), lty=1)

#' ##### Linearity of continuous variables
#' Fit a spline through continuous covariates
#+ clinicalSpline, fig.width=6, fig.height=6
set.seed(42)
trainIdx <- sample(c(TRUE,FALSE), nrow(dataFrame), replace=TRUE, prob=c(0.66,0.34))
trainIdxOsTD <- trainIdx[tplSplitOs]
par(mfrow=c(3,3))
clinicalSpline <- as.data.frame(sapply(dataFrame[groups %in% c("Clinical","Demographics")], function(x){
					if(all(x[1:5] %in% 0:10)) return(x)
					y <- log(x+min(x)+1e-3)
					fit <- coxph(os ~ pspline(y, df=3), subset=trainIdx)
					predict(fit, newdata=data.frame(y=y))
				}))
for(n in names(clinicalSpline)) if(!all(dataFrame[1:5,n] %in% 0:10))
		plot(dataFrame[,n], clinicalSpline[,n], log='x', xlab=paste(n, '[observed]'), ylab = paste(n, '[spline]'))

summary(coxph(os ~ ., data=clinicalSpline, subset=!trainIdx))$concordance
summary(coxph(os ~ ., data=dataFrame[groups %in% c("Clinical","Demographics")]), subset=!trainIdx)$concordance

#' No measurable improvement over (scaled) linear terms thus. 
#' 
#' ### Random effects models
#' Here we fit the random effects model using our implementation in the `CoxHD` package. First for main effects only.
 
#+ coxRFXFitOsTDMain, cache=TRUE
coxRFXFitOsTDMain <- CoxRFX(dataFrameOsTD[,mainIdxOsTD], osTD, groups[mainIdxOsTD])

#' Now including gene:gene interaction terms (min. recurrence = 8)
 
#+ coxRFXFitOsTDGG, cache=TRUE
whichRFXOsTDGG <- which((colSums(dataFrame)>=8 | mainIdxOsTD) & osTDIdx & groups %in% c(mainGroups,"GeneGene")) # ie, > 0.5%
coxRFXFitOsTDGGc <- CoxRFX(dataFrameOsTD[,whichRFXOsTDGG], osTD, groups[whichRFXOsTDGG], which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 

#' Compute Harrel's concordance index
survConcordance(osTD~coxRFXFitOsTDGGc$linear.predictors)

#' #### Variance components
#' Here we compute the variance components.
#+ varianceComponents
colGroups <- c(brewer.pal(12, "Paired")[c(10)],brewer.pal(12, "Paired")[c(6,4,3,5,12,9,1,2,7)],"#999999", brewer.pal(12, "Paired")[c(8)])
colGroups <- colGroups[c(2:6,1,7:14)]
names(colGroups) <- levels(groups)[order(toupper(levels(groups)))]
PlotVarianceComponents(coxRFXFitOsTDGGc, col=colGroups)
title("Risk contributions OS (time-dep)")

#' ### Confidence intervals and significance tests
#' Estimate confidence intervals by parametric bootstrap and compare with Wald Test. Note that the usual sample with replacement 
#' yields inconsistencies for the interaction terms due to the overdispersed correlations.
#' The theoretical description of the survival time simulation is given in section [survival](#survival).
#+ parBoot, cache=TRUE, eval=TRUE
set.seed(42)
risk <- as.matrix(dataFrame[whichRFXOsTDGG]) %*% coxRFXFitOsTDGGc$coefficients
risk <- risk - mean(risk)
parBoot <- mclapply(1:100, function(i) {
			s <- SimSurvNonp(risk, os)
			c <- try(CoxRFX(dataFrame[whichRFXOsTDGG], s, groups=groups[whichRFXOsTDGG], sigma0=0.1, nu=0))
			if(class(c)=="try-error")
				return(s)
			c$Z <- NULL # set X to zero to save mem
			return(c)
		}, mc.cores=10)

#' Distributions of mean, sigma and df
#+ parBootPlots, fig.width=3, fig.height=2,  eval=TRUE
boxplot(t(sapply(parBoot, `[[`, "sigma2")), border=colGroups[names(parBoot[[1]]$sigma2)], lty=1, pch=16, staplewex=0, ylab="sigma2", las=2, log="y", ylim=c(1e-3,1))
abline(h=0, lty=3)
points(coxRFXFitOsTDGGc$sigma2,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "mu")), border=colGroups[names(parBoot[[1]]$mu)], lty=1, pch=16, staplewex=0, ylab="mu", las=2)
abline(h=0, lty=3)
points(coxRFXFitOsTDGGc$mu,  pch=19)

boxplot(t(sapply(parBoot, `[[`, "df")), border=colGroups[names(parBoot[[1]]$mu)], lty=1, pch=16, staplewex=0, ylab="df", las=2)
abline(h=0, lty=3)
points(coxRFXFitOsTDGGc$df,  pch=19)

#' Coefficients
#+ parBootSignif, fig.width=2.5, fig.height=2.5
v <- apply(sapply(parBoot, `[[`, "coefficients"), 1, var, na.rm=TRUE)
w <- diag(coxRFXFitOsTDGGc$var) ## H^{-1}
w2 <- diag(coxRFXFitOsTDGGc$var2) ## H^{-1} I H^{-1}
c <- coef(coxRFXFitOsTDGGc)
plot(c^2/v, c^2/w, log="xy", xlab="Chi2 (bootstrap)", ylab="Chi2 (analyt.)", cex=.66)
par(xpd=NA)
points(c^2/v, c^2/w2, pch=16, cex=.7)
arrows(c^2/v, c^2/w, c^2/v,c^2/w2, length=0.05)
abline(0,1)
abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))
#' The plot indicates a good agreement of the variance estimate `var2`, see [section 2.1.4.](#Analytical-confidence-intervals-of-individual-parameters). Knowing the distribution
#' of the variance allows us to compute a Wald test of the coefficients.

#' #### Supplementary Table S1
#' Table with significance
#+ parBootTable, results='asis'
library(DT)
library(htmlwidgets)
pBoot <- pchisq(c^2/v,1, lower.tail=FALSE)
pVar2 <- pchisq(c^2/w2,1, lower.tail=FALSE)
pVar <- pchisq(c^2/w,1, lower.tail=FALSE)
waldOut <- data.frame(group = groups[whichRFXOsTDGG], 
		`beta (log-hazard)`= c, 
		`hazard exp(beta)` = exp(c),
		n = ifelse(groups[whichRFXOsTDGG] %in% c("CNA","Fusions","Genetics","GeneGene"), colSums(dataRaw[sub("_10*$","",names(whichRFXOsTDGG))], na.rm=TRUE), NA),
		sd = sqrt(w2), 
		`sd (bootstrap)` = sqrt(v), 
		`sd (var)`= sqrt(w),
		`P-value`=pVar2,
		`Q (Benjamini-Yekutieli)` = p.adjust(pVar2, "BY"),
		`Q (Benjamini-Hochberg)` = p.adjust(pVar2, "BH"),
		check.names=FALSE
)
datatable(as.data.frame(lapply(waldOut, function(x) if(class(x)=="numeric") round(x,4) else x), check.names=FALSE, row.names=row.names(waldOut)))
library(xlsx)
wb <- createWorkbook("xlsx")
sheet  <- createSheet(wb, sheetName="Overall survival")
addDataFrame(waldOut, 
		sheet,
		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
)

#' Volcano plot
#+ volcanoGGc, fig.width=3, fig.height=3
par(mar=c(3,3,1,1)+.1,  bty="n", mgp=c(2,.5,0))
i <- coxRFXFitOsTDGGc$groups %in% c("Genetics", "CNA","Fusions","GeneGene","Treatment")#apply(coxRFXFitOsTDGGc$Z,2,min) == 0 & apply(coxRFXFitOsTDGGc$Z,2,max) == 1
p <- pVar2 ## pvalues coxRFX
plot(c, 1/p, log='y', col=paste(colGroups[as.character(coxRFXFitOsTDGGc$groups)],"BB", sep=""), pch=ifelse(i,16,16), ylab="P-value",xlab="log hazard", cex=ifelse(i, sqrt(colMeans(coxRFXFitOsTDGGc$Z[!rev(duplicated(rev(tplSplitOs))),])*50),1), xlim=range(c*1.2))
#abline(h=qchisq(c(0.95,0.99,0.999), 1, lower.tail=TRUE), lty=c(1,2,3))
w <- which(p.adjust(p,"BY") < 0.1)
points(c[w], 1/p[w],  pch=1, cex=ifelse(i[w], sqrt(colMeans(coxRFXFitOsTDGGc$Z[!rev(duplicated(rev(tplSplitOs))),w])*50),1))
w <- which(p.adjust(p,"bonf") < 0.05)
par(xpd=NA)
text(c[w], 1/p[w], names(c[w]), pos=3)
u <- par("usr")
f <- c(0.01,0.05,0.1,0.2,0.5)
s <- sqrt(f*50)
legend("topright",legend=f, pch=16, pt.cex=s, bty='n', col=paste("#88888888"))
par(xpd=FALSE)
abline(h=1/0.05, lty=2)
abline(h=1/max(p[which(p.adjust(p,"BY") < 0.1)]), lty=3)

#' P-values and random model
#+ pValuesMain, fig.width=2.5, fig.height=2.5, cache=TRUE
set.seed(42)
Z <- apply(coxRFXFitOsTDGGc$Z, 2,sample)[1:nrow(dataFrame),] ## random covariates
coxRFXFitOsRain <- CoxRFX(Z, os, groups=coxRFXFitOsTDGGc$groups, nu=1) ## model
w2 <- diag(coxRFXFitOsRain$var2) 
c <- coef(coxRFXFitOsRain)
p2 <- pVar2
plot(seq(0,1,l=length(p2)+1)[-1],sort(p2), xlab="P-value (expected)", ylab="P-value (observed)", pch=16, col="grey")
abline(0,1)
points(seq(0,1,l=length(p)+1)[-1],sort(p), pch=16)
legend("topleft",bty="n", c("observed","randomised"), pch=16, col=c("black","grey"))

#' Distribution of the variance components
#+ parBootVarianceComp, fig.width=3, fig.height=2, cache=TRUE, eval=TRUE
v <- t(sapply(parBoot, function(x) {t <- try(VarianceComponents(x, newZ=dataFrame[whichRFXOsTDGG])); if(class(t)=="try-error") rep(NA, nlevels(x$groups)+1) else t}))
boxplot(v, border=colGroups[colnames(v)], lty=1, pch=16, staplewex=0, ylab="variance comp.", las=2)
abline(h=0, lty=3)
points(VarianceComponents(coxRFXFitOsTDGGc),  pch=19)

rm(parBoot)

#' ### Risk constellation plots
#' #### Supplementary Figure S3A
#' Plot of log hazard v outcome
#+ logHazOutcome, fig.width=3, fig.height=2
par(mar=c(3,3,3,1), mgp=c(2,.5,0))
t <- os
s <- survfit(os~1)
q <- quantile(t[,1], seq(0,1,.1))# q <- splinefun( s$surv, s$time,"monoH.FC")(seq(1,min(s$surv),l=10))
c <- cut(t[,1], q, na.rm=TRUE)
h <- coxRFXFitOsTDGGc$linear.predictors[rev(!duplicated(rev(tplSplitOs)))][order(tplSplitOs[rev(!duplicated(rev(tplSplitOs)))])]
o <- order(h)
plot(h[o], col= (brewer.pal(10,'RdBu'))[c[o]], type='h', xaxt="n", xlab='Patient', las=2, ylab="log hazard")
u <- par("usr")
q <- pmin(q,365*12)
image(x=q/max(q)*500, y=c(u[4]-(u[4]-u[3])/20, u[4]), matrix(1:10), col= (brewer.pal(10,'RdBu')), add=TRUE)
#axis(side=3, at=seq(1,500,l=11), labels=seq(0,1,.1))
axis(side=3, at=pretty(q/365)/max(q)*365*500, labels=pretty(q/365))
lines(ksmooth(seq_along(o),t[o,2]==0, bandwidth=50))

#' #### Supplementary Figure S3C
#' Risk constellation plots using the `stars()` function
#+ stars, fig.width=12, fig.height=12
set.seed(42)
library(HilbertVis)
nStars <- 32
s <- sample(nrow(dataFrame),nStars^2) #1:(nStars^2)
l <- "coxRFXFitOsTDGGc"
t <- os#get(l)$surv
p <- PartialRisk(get(l),  newZ=dataFrame[, colnames(get(l)$Z)])
p <- p[,colnames(p)!="Nuisance"]
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
h <- hclust(dist(p[s,]))
x <- p - rep(colMeans(p), each=nrow(p))
x <- x/(2*sd(x)) + 1
c <- cut(t[s,1][h$order], quantile(t[,1], seq(0,1,0.1), na.rm=TRUE))
if(l=="coxRFXFitOsTDGGc")
	x <- x[,c("Demographics","Treatment","Fusions","CNA","Genetics","GeneGene","Clinical")]
mg14:::stars(x[s,][h$order,]/2, scale=FALSE, locations=locations, key.loc=c(0,-3), col.lines=ifelse(t[s,2][h$order],1,NA), col.stars = (brewer.pal(11,'RdBu'))[c], density=ifelse(t[s,2][h$order],NA,NA))
symbols(locations[,1], locations[,2], circles=rep(.5,(nStars^2)), inches=FALSE, fg="grey", add=TRUE, lty=1)
title(main=l)


#' #### Supplementary Figure S3B
#' Randomly select 5 patients according to their genotype/outcome
#+ patientStars, fig.width=4, fig.height=3
patients <- c(
		which(dataFrame$`TP53`==1 & dataFrame$complex==1 & os[,1] < 300 & os[,2]==1)[1],
		which(dataFrame$`NPM1:FLT3_ITD:DNMT3A`==1 & os[,1] < 300 & os[,2]==1)[1],
		which(dataFrame$SFRS2==1 & clinicalData$WHOcat=='no' & os[,2]==1)[1],
		which(dataFrame$NPM1==1 & dataFrame$FLT3_ITD==0 & dataFrame$DNMT3A==0 & os[,1] > 2000)[1],
		which(dataFrame$t_15_17==1 & os[,1] > 2000)[1]
)
genotype <- apply(dataFrame[groups %in% c("Fusions","CNA","Genetics")]==1, 1,function(x) paste(names(which(x)), collapse=";"))

t <- os
p <- PartialRisk(coxRFXFitOsTDGGc, newZ=dataFrame[, whichRFXOsTDGG])
p <- p[,colnames(p)!="Nuisance"]
locations <- 1.5*hilbertCurve(log2(nStars)) #2*expand.grid(1:nStars,1:nStars)
x <- p - rep(colMeans(p), each=nrow(p))
x <- x/(2*sd(x)) + 1
c <- cut(t[patients,1], quantile(t[,1], seq(0,1,0.1), na.rm=TRUE))
x <- x[patients,c("Demographics","Treatment","Fusions","CNA","Genetics","GeneGene","Clinical")]
locations <- expand.grid(seq_along(patients)* 1.5, 1)
mg14:::stars(x/2, scale=FALSE, locations=locations, key.loc=NA, col.lines=ifelse(t[patients,2],1,NA), col.stars = (brewer.pal(11,'RdBu'))[c])
symbols(locations[,1], locations[,2], circles=rep(.5,length(patients)), inches=FALSE, fg="grey", add=TRUE, lty=1)
text(locations[,1], locations[,2]-1,labels=clinicalData$PDID[patients], pos=1)
l <- apply(dataFrame[patients,c("gender","AOD_10","TPL_os","wbc_100")], 1,paste, collapse=";")
par(xpd=NA)
text(locations[,1], locations[,2]+1,labels=paste(gsub(";","\n",genotype[patients]),l, paste(round(os[patients,1],2), osYr[patients,2]), sep="\n"), pos=3)



 
#' # Multistage modelling
#' 
#' ## Definitions
#' 
#' ### Nomenclature
#' We use the following nomenclature: $f(T=t)=f(t)$ denotes a probability density, $F(T=t)=P(T<t)=F(t)$ the corresponding cumulative distribution function. $S(T=t)=1-F(t)$ is the survivor function, 
#' the name being motivated by the situation that $t$ is a death time. In cases where it is clear to which variable a (cumulative) density refers to, we may drop the stochastic variable and 
#' simply use its value as the argument, $f(t)=f(T=t)$. We use the convention of lower case variables $t$ to denote the values of the corresponding upper case stochastic 
#' variable $T=t$, $U=u$ and so on. For a categorical stochastic process $X_t$, $t\in \mathbb{R}^+$ we use the symbol $P(X_t)$ to denote the 
#' probability distribution at time $t$. The symbol $Z$ denotes the covariates.
#' 
#' ### Transitions
#' 
#' We use a hierarchical multistage model to quantify the rates at which a patient progresses from one disease/treatment stage to another (Figure 3a).
#' After learning the marginal time-dependent transition probabilities for each event, we can combine these into a time-dependent joint probability. 
#' 
#' In particular, we model the following transition times:
#' 
#' * Time of complete remission (complete remission CR, $T_{CR}$)
#' * Time of non-remission death (Non-complete remission death NCD, $T_{NCD}$)
#' * Time from CR to relapse (Relapse R, $T_{R}$)
#' * Time from CR to non-relapse death (Non-relapse death NRD, $T_{NRD}$)
#' * Time from relapse to post-relapse death (Post-relapse mortality PRD, $T_{PRD}$)
#' 
#' CR and CIR and midpoints, allowing for further events, NCD, NCD and PRD are endpoints. 
#' Due to the hierarchical nature of the model, only one endpoint can ever occur and midpoints are transient. 
#' 
#' ### States
#' 
#' The probability to be in a given state is then given by the combination of
#' event times, such that 
#' 
#' * the transition to a given stage has happened before the other competing transitions
#' * no subsequent transition has occurred yet
#' * no other endpoint has been reached yet
#'  
#' To be alive in CR at time t, for example, requires that CR occurred before t, CR was achieved before NCD, and neither relapse nor NCD have occurred yet.
#' Overall, a patient can only be in one of the following six states at time t, each corresponding to
#' a particular ordering of event times:
#' 
#' 
#' | Stage | Abbreviation | Ordering of times | Symbol |
#' |-------|--------------|-------------------|--------| 
#' | Alive in induction | AI| $t < T_{CR}, T_{NCD}$ | $\mathcal{I}_{AI}(t)$ |
#' | Death without complete remission | NCM | $T_{NCD} < t$; $T_{NCD} < T_{CR}$ |  $\mathcal{I}_{NCM}(t)$ |
#' | Alive in complete remission | ACR | $T_{CR} < T_{NCD}$; $t < T_R, T_{NCD}$ |  $\mathcal{I}_{ACR}(t)$ |
#' | Death without relapse | NRM | $T_{CR} < T_{NCD}$; $T_{NRD} < t < T_R$| $\mathcal{I}_{NRM}(t)$ |
#' | Alive after relapse | AR | $T_{CR} < T_{NCD}$; $T_R < T_{NRD}$; $t < T_{PRD}$ |  $\mathcal{I}_{AAR}(t)$ |
#' | Death after relapse | PRM | $T_{CR} < T_{NCD}$; $T_R < T_{NRD}$ ; $T_{PRD} < t$| $\mathcal{I}_{PRM}(t)$ |
#' 
#' Table: Table 2. Stages 
#' 
#' This defines a stochastic process $X_t$ on the given set of six states, $X_t \in \{AI, ACR, AR, NCM, NRM, PRM\}$. Initially, all patients will be alive in induction, $X_0=AI$.
#' 
#' ### Factorisation of the joint probability
#' 
#' The hierarchical nature of the model implies that the joint probability of event time factorises
#' 
#' $$ f(T_{CR},T_{NCD}, T_{R}, T_{NRD}, T_{PRD}) = f(T_{CR}) \times f(T_{NCD}) \times f(T_R \mid  T_{CR})\times f(T_{NRD}\mid T_{CR})\times f(T_{PRD} \mid  T_R). \label{eq:joint-density}$$
#' 
#' The above factorisation lays out a strategy in which each of the 5 factors may be estimated separately. 
#' The probability of each state $P(X_t)$, defined in section [states](#states), are then computed by integrating the joint density $f$, Eq.$\eqref{eq:joint-density}$ over the
#' simplexes $\mathcal{I}_{\cdot}(t)$ defining a particular ordering of transitions detailed in table 2:
#' $$ P(X_t = x) =  \iiiint\!\!\!\!\!\int_{\mathcal{I}_x(t)} f(t_{CR},t_{NCD}, t_{R}, t_{NRD}, t_{PRD})\ dt_{CR}\ dt_{NCD}\ dt_{R}\ dt_{NRD}\ dt_{PRD}. \label{eq:mult-prob}$$
#' 
#' The integral can be successively evaluated as described below.
#' 
#' 
#' ## Static multistage models
#' 
#' To estimate the population average transition probabilities and absolute incidence of each each individual stage we use the `msSurv` R package [@FergusonJOSS2012]. 
#' The resulting time-dependent joint distribution $P(X_t)$ is shown in Figure 3b.
#' 
#' ## Multistage random effects modelling
#' 
#' To estimate how each transition $T$ depends on the set of variables $Z$ introduced in section [variables](transitions), we use a random effects model for each transition to
#' obtain $f(T\mid Z)$. Competing events are considered to be censored. We apply a separate random effects model to estimate all five terms in Eq.$\eqref{eq:joint-density}$.
#' 
#' ### Unconditional densities
#' The estimation of the unconditional densities $f(T=t\mid Z)$ is straightforward.
#' The random effects model yield the marginal survivor function $S(t \mid  Z) = S_0(t) ^{\exp(u Z)}$, quantifying 
#' the hypothetical scenario that there were no competing events $T'$. From $S(t\mid Z)$ we can derive
#' the marginal densities $f(t \mid  Z) = -dS(t)/dt = -\exp(u Z) S_0(t)^{\exp(u Z) -1} dS_0(t)/dt$ for each transition $T$ given the covariates $Z$. 
#' The Kaplan-Meier estimate of $S_0(t)$ is a step function, so we may compute $f(t\mid Z)$ via a numerical differentiation.
#' 
#' ### Conditional densities
#' To estimate the conditional densities of the type $f(U=u\mid Z, T)$, we use the following approach: 
#' 
#' $$S(u \mid  Z, T=t) =  S_0(u-t \mid  Z) ^{\exp(g(t))} = S_0(u-t) ^{\exp(u Z + g(t))}. \label{eq:cond-dens}$$
#' 
#' This allows us to estimate the incidence of each event from the beginning of each stage $S_0(u-t \mid  Z)$ and express the time-dependence as
#' a smooth function $g(t)$. For example, the duration of CR1 is a prognostic factor for post-relapse mortality, e.g. [@BurnettJCO2013].
#' 
#' The above corresponds to a Cox proportional hazards model with a time-dependent smooth covariate g(t). 
#' 
#' Here we estimate $g(t)$ with a spline term with 10 degrees of freedom. We 
#' estimate $S_0(u-t)$ and $u$ using a random effects model and subsequently estimate $g(t)$ using $u Z$ as an
#' intercept: 
#' 
#' ```{r, eval=FALSE}
#' fit_u_minus_t_given_Z <- CoxRFX(Z, Surv(U-T, event))
#' beta <- coef(fit_u_minus_t_given_Z)
#' fit_u_given_Z_t <- coxph(Surv(U-T, event) ~ I(beta %*% Z) + spline(T, df=10))
#' ```
#' 
#' We may thus obtain $S(u\mid T=t, Z) = \left(S_0(u-t)^{\exp(u Z)}\right)^{\exp(g(t))}$ by offsetting the baseline hazard $S_0(u-t)$ and exponentiating
#' for the effect of covariates $u Z$ and exponentiating for the effect of time-dependence.
#' 
#' The absolute probability to be in state $U$ is given by integrating over the conditional probabilities $S(u\mid T=t, Z)$, weighted by the probability of 
#' the preceding event probabilities $f(T=t\mid Z)$:
#' 
#' $$P(U < u \mid  Z) = \int_0^t f(T=t \mid  Z) \int_t^u f(U=v \mid  T=t, Z) dt dv = \int_0^u f(T=t \mid  Z) F(u \mid  T=t, Z) dt \label{eq:cond-prob}$$ 
#' 
#' With $F(u\mid T=t, Z) = 1-S(u\mid T=t, Z) $, we can use the above definition to numerically solve the above integral.
#' 
#' The pseudo code for this is given by:
#' 
#' ```{r, eval=FALSE}
#' S0_given_Z <- S0 ^ exp(beta %*% Z)
#' gt <- predict(fit_u_given_t, data=data.frame(T=1:length(u)))
#' ft <- -diff(St)
#' for(t in 1:length(u)){
#'    Fu_given_Zt <- 1-S0_given_Z[-(1:t)] ^ exp(gt)
#'    Pu <- cumsum(ft * Fu_given_Zt)
#' }
#' ```
#' 
#' ### Competing risk adjustment
#' In cases of competing events (CR and NCD; NCD and CIR), we use a competing risk adjustment between two event times $T$, $U$, to obtain 
#' $$S(T=t \mid  Z, T < U) = \int_t^v \int_v^\infty f(T=t'\mid  Z) f(U=u' \mid  Z) dt' du' = \int_0^t f(T=t' \mid  Z) S(U=t'\mid Z) dv.$$
#' 
#' In practical terms, $S(t\mid Z) = S_0(t) ^{\exp(u Z)}$ denotes the survivor function estimated by the Kaplan-Meyer estimate $S_0(t)$, exponentiated by the hazard $\exp(u Z)$. The differential $f(t\mid Z)$ is
#' obtained by evaluating the difference of $S(t+1\mid Z) - S(t\mid Z)$ at intervals of length 1 day, pseudo code
#' 
#' `S_t_cr <- cumsum(diff(S_t) * S_u)`
#' 
#' 
#' ### Encoding of events
#' 
#' | Endpoint | Censored | Model | Competing | Interval | Time-dependency |
#' |----------|----------|-------|-----------|----------|-----------------|
#' | Complete remission | Non-remission death | CoxRFX | Non-remission death | From ER | |
#' | Non-remission death | Complete remission | CoxRFX | Complete remission | From ER | |
#' | Relapse | Non-relapse death | CoxRFX | Non-relapse death | From CR1 | Time to CR1 |
#' | Non-relapse death | Relapse | CoxRFX | Relapse | From CR1 | Time to CR1 |
#' | Post-relapse death | Last follow up | CoxRFX | - | From relapse | Duration of CR1 | 
#' 
#' Table: Table 3.
#' 
#' ### Probabilities of each state
#' 
#' As the density $\eqref{eq:joint-density}$ factorises we can successively evaluate each term, beginning with the first transition. 
#' The probability to be in a given state is then computed according to the rules outlined in the previous two subsections. 
#' 
#' 
#' #### Death without complete remission
#' 
#' The probability is given by a simple competing risk adjustment between T_{CR} and T_{NCD}:
#' 
#' $$ P(X_t = NCM) \mid  Z) = P(T_{NCD} < t, T_{NCD} < T_{CR} \mid  Z) =  1 - \int_0^t f(T_{NCD} = u \mid  Z) F(T_{CR} = u\mid Z) du.$$
#' 
#' #### Complete remission
#' 
#' We first compute the probability that CR is achieved, irrespective of the subsequent events, using a competing risk adjustment with T_{NCD}:
#' 
#' $$P(X_t = CR \mid  Z) = P(T_{CR} < T_{NCD}, T_{CR} < t\mid Z) = 1 - \int_0^t f(T_{CR} = u \mid  Z) F(T_{NCD} = u\mid Z) du ,$$
#' 
#' where $CR = \{ACR \cup NRD \cup AAR \cup PRD\}$, which can then be further subdivided according to the possible subsequent events. 
#' To be alive, neither relapse nor non-relapse death may have occurred.
#' 
#' #### Alive in induction
#' The probability to be alive in induction is given by 
#' $$P(X_t = AAI\mid Z) = 1-P(X_t = CR\mid Z) -P(X_t=NRD\mid Z).$$
#' 
#' #### Non-relapse deaths
#' 
#' The probability of non-relapse deaths $P(X_t = NRM \mid Z)$ is computed in the following way. We first estimate transition rates
#' for non-relapse deaths and relapses,  $f(T_{NRM} \mid  Z, T_{CR})$  
#' and  $f(T_{NRM} \mid  Z, T_{CR})$, as outlined in 
#' [conditional-densities](#conditional densities). Instead of the marginal density $f(T_{CR} \mid  Z)$ we use the differential
#' of $dP(X_t = CR \mid  Z) / dt$ in Eq.$\eqref{eq:cond-prob}$.
#' 
#' As only one of the two events can ever
#' occur we then use a [competing risk adjustment](#competing-risk-adjustment) to obtain the absolute probability
#' probability of $P(X_t = NRD \mid Z)$ and $P(X_t = R \mid Z)$, respectively, where $R = \{ AAR \cup PRM \}$ denotes a relapse.
#' 
#' #### Alive in complete remission
#' The probability to be alive in first complete remission equals the probability of neither dying nor relapsing:
#' $$P(X_t=ACR) = 1 - P(X_t = NRD) - P(X_t = R). $$
#' 
#' #### Post-relapse death
#' The probability of post-relapse deaths $P(X_t = PRM \mid Z)$ is computed as described in 
#' [conditional-densities](#conditional densities). We first estimate the rate
#' for post-relapse deaths and relapses,  $f(T_{PRM} \mid  Z, T_{R})$, with the derivative
#' of $dP(X_t = R \mid  Z) / dt$ in Eq.$\eqref{eq:cond-prob}$. 
#' 
#' #### Alive after relapse
#' Finally, the probability to be alive after relapse is given by
#' $$ P(X_t = AAR \mid  Z) = P(R\mid Z) - P(PRM\mid Z).$$
#' 
#' ### Comments
#' In the absence of an established estimator of the joint density Eq.$\eqref{eq:joint-density}$, we assumed that each factor of
#' the density may be separately estimated using a random effects model. We note that the interdependence of observed events
#' could in general introduce a bias as the censoring is not independent. The precise magnitude of this effect still needs to be investigated.
#' 
#' We observed, however, a good consistency of the average
#' predictions and static multistage probabilities, indicating that those biases, on average, tend to cancel. Moreover cross-validation of our methodology ascertained 
#' a very good predictive performance despite all potential shortcomings.
#' 
#'     
#' ## Confidence intervals
#' 
#' ### Marginal probabilities
#' 	
#' For each predicted variable we can derive 95% confidence intervals from the prediction error of
#' the log hazard,  $(h_{0.025},h_{0.975})\approx  h + (-2,2) \times \hat V[h \mid  Z]$, with $V[h\mid Z]$. 
#' This translates to
#'  the survival function as follows using the log-log approach:
#' 	$$S_{0.025}(t \mid  Z) = S_0(t)^{\exp(h_{0.025})}$$
#' 	$$S_{0.975}(t \mid  Z) = S_0(t)^{\exp(h_{0.975})}.$$
#' 
#' Note that this does not model the error of the baseline survival estimate $S_0(t)$.
#' 
#' ### Survival after remission
#' 
#' Let the symbol PCS denote post remission survival. In the following sections all quantities are conditional on the data $Z$.
#' 
#' #### Analytical confidence intervals
#' 	
#' Analytical confidence intervals can be calculated using a the propagation of errors based on a Taylor expansion of the PCS probability:
#' $$\begin{align}
#' V[h_{PCS}]  &\approx \sum_i \left(\frac{\partial h_{PCS}}{\partial h_i}\right)^2 V[h_i] \\
#' h_{PCS} &=  \log\log P_{PCS} +  \log\log P_0(t)\\
#' \frac{\partial h_{PCS}}{\partial h_i} &= \frac{\partial \log\log P_{PCS}}{\partial h_i} \\
#' &= \frac{1}{P_{PCS}\log(P_{PCS})} \frac{\partial P_{PCS}}{\partial h_i}
#' \end{align}
#' $$
#' 
#' To facilitate an efficient computation of the derivatives of $P_{PCS}$, which is given by the integrals above, we use the pointwise approximation: 
#' $$\begin{align}
#' P_{PCS} &\approx S_{NRD}(1 -  (1-S_{R})(1- S_{PRD})) = S_{PCS}\\
#' \end{align}$$ 
#' where $S_\cdot$ denote the Kaplan-Meyer estimates of the survival probabilities.
#' 
#' The partial derivative of the loglog is given by
#' $$\frac{\partial \log\log S}{\partial x} = \frac{1}{\log(S)} \frac{1}{S} \frac{\partial S}{\partial x}.$$
#' 
#' So the variance of the loglog overall survival reads: 
#' $$V[h_{PCS}]\approx \frac{1}{(S_{PCS}\log S_{PCS})^2} \left( V[S_{NRD}](1-(1-S_{R})(1- S_{PRD}))^2 + V[S_{R}] S_{NRD}^2 (1-S_{PRD})^2  + V[S_{PRD}] S_{NRD}^2 (1-S_{R})^2 \right )$$
#' with 
#' $$V[S_i] = \left(\frac{\partial S_0(t)^{\exp(h)}}{\partial h}\right)^2 V[h] = (S_i \log S_i)^2 V[h]. $$
#' 
#' This allows to define 95% confidence intervals of $h_{PCS}$ as:
#' $$h_{PCS}^{0.975} = h_{PCS} + 2V[h_{PCS}].$$
#' 
#' This translates to an overall survival:
#' $$P_{PCS}^{0.975} = \exp \exp(h_{PCS}^{0.975}) = \exp (\exp(h_{PCS}) \exp( 2V[h_{PCS}]) ) = P_{PCS} ^{\exp(2 V[h_{PCS}])}$$
#' 
#' Note that in the last step uses the competing risk and time-adjusted estimate $P_{PCS}$ again
#' 
#' #### Simulated
#' A more accurate account comes from simulations of errors in the predicted log hazard. The cumulative
#' survival functions are given by
#' $$S{_\cdot}(t) = S_{0\cdot}(t)^{\exp(h_\cdot + \epsilon_\cdot)}$$
#' So drawing 
#' $$\epsilon_\cdot \sim N(0,\hat V[h_\cdot])$$ 
#' for each event type and repeating the computations outlined in [combined-os] yields
#' an empirical distribution of the survival distribution of $S_{PCS}(t)$.
#' 
#' We use `i=200` simulations to compute the empirical confidence intervals.
#' 
#' **Note**: In all cases the prediction errors are assumed to be independent.
#' 
#' ### Differential survival
#' 
#' Confidence intervals for differential survival, e.g. with and without allograft are computed as in the previous section. A complication 
#' arises as the errors are correlated. We hence sample errors for all common variables and then sample those variable that differ. This approach
#' allows to assess the uncertainty resulting from a subset of variables, and on the background of the joint variation in the set of common features.
#' 
#' ### Overall survival from diagnosis
#' 
#' We use a numerical approach similar to the one outlined above to compute the confidence intervals for overall survival measured from diagnosis.
#' Note that it is in principle also possible to derive analytical confidence intervals analogous to [section 5.4.2.1](#analytical-confidence-intervals).
#' 

#' ## Code
#' 
#' ### Static multistage model
#' #### Figure 3b 
#' Multi-state using msSurv  [@FergusonJOSS2012].
#+ mstate, fig.width=3,fig.height=2.5
library(msSurv)
d <- sapply(1:nrow(clinicalData), function(i){
			i <<- i
			t <- c(as.numeric(clinicalData[i,c("CR_date","Recurrence_date","Date_LF")]) - as.numeric(clinicalData$ERDate[i]))
			o <- order(t, na.last=NA)
			stages <- c(1:3,0)
			r <- stages[c(1, o+1)]
			if(clinicalData$Status[i])
				r[length(r)] <- r[length(r)-1] +3
			tt <- c(0,t[o])
			if(length(o)==0)
				return(c(rep(NA,7),i))
			s <- cbind(id=i, stop=tt[-1], start.stage=r[-length(r)], end.stage=r[-1])[diff(tt)!=0,]
			#s <- cbind(time1 = tt[-length(tt)], time2=tt[-1], death=c(rep(0, length(o)-1), clinicalData$Status[i]), outer(0:(length(o)-1), r[-3], `>=`)+0, i=i)[diff(tt)!=0,]
			return(s)
		})
d <- as.data.frame(do.call("rbind",d))
nodes <- as.character(1:6)
edges <- list(`1`=list(edges=c("2","4")), `2`=list(edges=c("3","5")), `3`=list(edges="6"), `4`=list(edges=NULL), `5`=list(edges=NULL),`6`=list(edges=NULL))
struct <-  new("graphNEL", nodes = nodes, edgeL = edges, edgemode = "directed")
msurv <- msSurv(d, struct, bs = FALSE)
y <- t(apply(cbind(1,-msurv@ps[,c(4:6, 3:1)]),1,cumsum))
par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0), las=1) 
plot(msurv@et/365.25, y[,1], ylim=c(0,1), type="s",lty=0, xlab="Time after diagnosis", ylab="Fraction of patients", xlim=c(0,10), xaxs="i", yaxs="i")
steps <- function(x, type="s") rep(x, each=2)[if(type=="s") -1 else -2*length(x)]
x <- steps(msurv@et/365.25, type="S")
for(i in 1:6)
	polygon(c(x, rev(x)), c(steps(y[,i]), rev(steps(y[,i+1])) ), col=c(brewer.pal(5,"Pastel1")[c(1:3,5,4)],"#DDDDDD")[i], border=NA)
abline(h=seq(0,1,.2), col='white', lty=3)
abline(v=seq(0,10,1), col='white', lty=3)
lines(x, steps(y[,4]), lwd=2)
w <- which.min(abs(msurv@et/365.25-10))
text(x=par("usr")[2], y= y[w,-7]+diff(y[w,])/2, labels=c("early death","death in CR","death after relapse","alive with relapse","alive in remission","induction/LOF"), pos=2)

#' ### Prepare covariates
#' Times for allografts pre and post relapse, after 1CR only
#+ alloIdx
alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD") # only allografts
alloTimeCR1 <- clinicalData$Time_1CR_TPL + .5 # +.5 to make > 0
alloTimeCR1[!alloIdx | (clinicalData$TPL_date < clinicalData$Recurrence_date & !clinicalData$TPL_Phase %in% c("CR1","RD"))] <- NA

#' Create data frames for each phase
#+ postCR1Data, cache=TRUE
whichRFXRel <- whichRFXOsTDGG[grep("TPL",names(whichRFXOsTDGG), invert=TRUE)] #mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Nuisance"
t <- clinicalData$Recurrence_date
t[is.na(t)] <- as.Date(1e6, origin="2000-01-01")
relData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=!is.na(clinicalData$Recurrence_date)+0)
relData$transplantCR1 <- relData$event
relData$event <- NULL
relData$transplantRel <- 0
nrdData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=is.na(clinicalData$Recurrence_date) & clinicalData$Status)
nrdData$transplantCR1 <- nrdData$event
nrdData$event <- NULL
nrdData$transplantRel <- 0
alloTimeRel <- clinicalData$TPL_date - clinicalData$Recurrence_date + .5 # +.5 to make > 0
alloTimeRel[!alloIdx | (clinicalData$TPL_date < clinicalData$Recurrence_date & !clinicalData$TPL_Phase %in% c("CR1","RD"))] <- NA
i <- !is.na(clinicalData$Recurrence_date)
prdData <- MakeTimeDependent(dataFrame[i,whichRFXRel], timeEvent=alloTimeRel[i], timeStop=as.numeric(clinicalData$Date_LF- clinicalData$Recurrence_date)[i], status=clinicalData$Status[i])
prdData$transplantCR1 <- rep(0,nrow(prdData))
w <- sub("\\.1","",rownames(relData))[relData$status==1 & relData$transplantCR1==1]
prdData$transplantCR1[sub("\\.1","",rownames(prdData)) %in% w] <- 1
prdData$transplantRel <- prdData$event
prdData$event <- NULL
w <- which(prdData$time1 == prdData$time2) ## 5 cases with LF=Rec
prdData$time2[w] <- prdData$time2[w] + .5
prdData$time0 <- as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index]

#' ### RFX fit of transitions
#+ postCR1Fits, cache=TRUE
crGroups <- c(as.character(groups[whichRFXRel]), "Treatment","Treatment")
names(crGroups) <- c(names(dataFrame)[whichRFXRel],"transplantCR1","transplantRel")
coxRFXNrdTD <- CoxRFX(nrdData[names(crGroups)], Surv(nrdData$time1, nrdData$time2, nrdData$status), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXNrdTD$coefficients["transplantRel"] <- 0
#prsData$time1[!is.na(prsData$time1)] <- 0
coxRFXPrdTD <-  CoxRFX(prdData[names(crGroups)], Surv(prdData$time1, prdData$time2, prdData$status), groups=crGroups, nu=1, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXRelTD <-  CoxRFX(relData[names(crGroups)], Surv(relData$time1, relData$time2, relData$status), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXRelTD$coefficients["transplantRel"] <- 0

#' #### OS
#+ coxRFXOs, cache=TRUE
osData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
osData$transplantCR1 <- osData$event
osData$transplantRel <- osData$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date | clinicalData$TPL_Phase != "CR1")  
osData$transplantCR1[osData$index %in% w] <- 0
osData$transplantRel[!osData$index %in% w] <- 0
data <- osData[rev(!duplicated(rev(osData$index))),colnames(coxRFXRelTD$Z)]
osData$transplantRel <- 0 # Note: confounded by relapse
rownames(data) <- sub("\\.1$","", rownames(data))
data <- data[rownames(dataFrame),]

coxRFXOsCR <- CoxRFX(osData[names(crGroups)], Surv(osData$time1, osData$time2, osData$status), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))

#' #### Early deaths
#+ coxRFXCr, cache=TRUE
table(CR=!is.na(clinicalData$CR_date), os[,2])

c <- as.numeric(clinicalData$CR_date - clinicalData$ERDate)
c[is.na(c)] <- clinicalData$OS[is.na(c)]
cr <- Surv(c, factor(pmin(2 * (!is.na(clinicalData$CR_date))+os[,2],2), levels=0:2, labels=c("cens","ED","CR")), type="mstate")

coxRFXCrTD <- CoxRFX(osData[1:1540, names(crGroups)], Surv(cr[,1], cr[,2]==2), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXNcdTD <- CoxRFX(osData[1:1540, names(crGroups)], Surv(cr[,1], cr[,2]==1), groups=crGroups, which.mu = NULL)

#save(coxRFXRelTD, coxRFXNrdTD, coxRFXPrdTD, coxRFXOsCR, coxRFXNcdTD, coxRFXCrTD, cr, nrdData, relData, prdData, osData, crGroups, data, clinicalData, file="../../code/multistage/multistage.RData")

#' ### Variance components
#+ allVarComp, fig.width=6, fig.height=4
par(mfrow=c(3,2), xpd=FALSE)
o <- c(1,4,6,5,2,3,7,8)
PlotVarianceComponents(coxRFXNcdTD, col=colGroups, order=o)
title(main="Early deaths")
PlotVarianceComponents(coxRFXCrTD, col=colGroups, order=o)
title(main="Remission")
PlotVarianceComponents(coxRFXRelTD, col=colGroups, order=o)
title(main="Relapse")
PlotVarianceComponents(coxRFXNrdTD, col=colGroups, order=o)
title(main="Non-relapse deaths")
PlotVarianceComponents(coxRFXPrdTD, col=colGroups, order=o)
title(main="Post-relapse deaths")

#' #### Figure 3c
#' As barplot
#+ allVarCompBar, fig.width=2, fig.height=2
par(mar=c(4,3,1,5))
allVarComp <- sapply(c("NcdTD","CrTD","NrdTD","RelTD","PrdTD"), function(x){
			m <- get(paste0("coxRFX",x))
			Z <- get(sub("\\[.+","",as.character(m$call["data"])))
			i <- if(x%in%c("CrTD","EsTD")) 1:1540 else Z$index
			VarianceComponents(m, newZ=Z[!rev(duplicated(rev(i))),colnames(m$Z)])})
colnames(allVarComp) <- c("Early deaths","Remission","Non-relapse d.","Relapse","Post-relapse d.")
w <- c("CNA","Fusions","Genetics","GeneGene","Clinical","Demographics","Treatment","Nuisance")
z <- allVarComp[w,]#/rep(colSums(allVarComp[-9,]), each=8)
b <- barplot(z, col=colGroups[w], ylab="Variance [log hazard]", names.arg=rep("",ncol(z)))
rotatedLabel(x0=b, labels=colnames(z))
Z <- rbind(0,apply(z,2,cumsum))
n <- ncol(z)
segments(b[-n]+.5,t(Z[,-n]),b[-1]-.5 ,t(Z[,-1]))

z <- allVarComp[w,]/rep(colSums(allVarComp[-9,]), each=8)
b <- barplot(z, col=colGroups[w], ylab="Relative importance", names.arg=rep("",ncol(z)))
rotatedLabel(x0=b, labels=colnames(z))
Z <- rbind(0,apply(z,2,cumsum))
n <- ncol(z)
segments(b[-n]+.5,t(Z[,-n]),b[-1]-.5 ,t(Z[,-1]))
mtext(side=4, at=Z[-1,n] - diff(Z[,n])/2, text=rownames(Z)[-1], las=2)

v <- c(1,3,5,4,2)
z <- allVarComp[w,v]/rep(colSums(allVarComp[-9,v]), each=8)
b <- barplot(z, col=colGroups[w], ylab="Relative importance", names.arg=rep("",ncol(z)))
rotatedLabel(x0=b, labels=colnames(z))
Z <- rbind(0,apply(z,2,cumsum))
n <- ncol(z)
segments(b[-n]+.5,t(Z[,-n]),b[-1]-.5 ,t(Z[,-1]))
mtext(side=4, at=Z[-1,n] - diff(Z[,n])/2, text=rownames(Z)[-1], las=2)


#' Pairwise scatter plots of the log hazard for each transition
#+ allStagesRisk, fig.width=4,fig.height=4
allStagesRisk <- as.data.frame(sapply(c("NcdTD","CrTD","NrdTD","RelTD","PrdTD"), function(x){
					m <- get(paste0("coxRFX",x))
					#Z <- get(sub("\\[.+","",as.character(m$call["data"])))
					#i <- if(x=="Cr") 1:1540 else Z$index
					Z <- if(x=="Cr") dataFrame else data[rownames(dataFrame),]
					predict(m, newdata=as.data.frame(Z))}))
f <- function(x,y,...) {points(x,y, col=densCols(x,y),...); lines(lowess(x,y), col='red')}
pairs(allStagesRisk, panel=f, pch=19)

#' #### Supplementary Tables 2-6 
#' Non-complete remission deaths
w <- WaldTest(coxRFXNcdTD)
w$Q.BH <- p.adjust(w$p.value, "BH")
w$Q.BY <- p.adjust(w$p.value, "BH")
datatable(w)
sheet  <- createSheet(wb, sheetName="Non-complete remission deaths")
addDataFrame(w, 
		sheet,
		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
)

#' Complete remission
w <- WaldTest(coxRFXCrTD)
w$Q.BH <- p.adjust(w$p.value, "BH")
w$Q.BY <- p.adjust(w$p.value, "BH")
datatable(w)
sheet  <- createSheet(wb, sheetName="Complete remission")
addDataFrame(w, 
		sheet,
		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
)


#' Relapses
w <- WaldTest(coxRFXRelTD)
w$Q.BH <- p.adjust(w$p.value, "BH")
w$Q.BY <- p.adjust(w$p.value, "BY")
datatable(w)
sheet  <- createSheet(wb, sheetName="Relapse")
addDataFrame(w, 
		sheet,
		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
)


#' Post-relapse survival
w <- WaldTest(coxRFXPrdTD)
w$Q.BH <- p.adjust(w$p.value, "BH")
w$Q.BY <- p.adjust(w$p.value, "BY")
datatable(w)
sheet  <- createSheet(wb, sheetName="Post-relapse deaths")
addDataFrame(w, 
		sheet,
		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
)


#' Non-relapse deaths
w <- WaldTest(coxRFXNrdTD)
w$Q.BH <- p.adjust(w$p.value, "BH")
w$Q.BY <- p.adjust(w$p.value, "BY")
datatable(w)
sheet  <- createSheet(wb, sheetName="Non-relapse deaths")
addDataFrame(w, 
		sheet,
		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
)


saveWorkbook(wb, file="SupplementaryTables.xlsx") 

#' ### Predicting outcome after CR
#' We use the following function to compute the hierarchical adjustment for two subsequent stages. It is implemented in C++ for efficiency using the `Rcpp` package [@EddelbuettelJOSS2011]. 
library(Rcpp)
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

#' Function to predict OS from Relapse, PRS and NRM, as described in [Section 4.3.5](#probabilities-of-each-state).
MultiRFX3 <- function(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, x =365, ciType="analytical", prdData){
	## Step 1: Compute KM survival curves and log hazard
	getS <- function(coxRFX, data, max.x=5000) {		
		if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
		data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data))])
		r <- PredictRiskMissing(coxRFX, data, var="var2")
		H0 <- basehaz(coxRFX, centered = FALSE)
		hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
		x <- c(0:max.x,max.x)
		S <- exp(-hazardDist(x))
		return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
	}
	kmRel <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
	kmNrd <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
	kmPrd <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
	
	## Step 2: Adjust CIR and NRM curve for competing risks, accounting for hazard
	kmRel$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmRel$S^exp(kmRel$r[i,1]))) * kmNrd$S ^ exp(kmNrd$r[i,1])))
	kmNrd$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmNrd$S^exp(kmNrd$r[i,1]))) * kmRel$S ^ exp(kmRel$r[i,1]))) ## array times x nrow(data)
	
	stopifnot(length(x)==1 | length(x) == nrow(data))
	if(length(x)==nrow(data))
		w <- match(x,kmRel$x)
	else if(length(x)==1)
		w <- rep(match(x, kmRel$x), nrow(data))
	y <- mapply(function(i,j) kmNrd$Sadj[i,j], w,1:length(w) ) # select time for each sample
	nrs <- y
	nrsUp <- y^exp(2*sqrt(kmNrd$r[,2]))
	nrsLo <- y^exp(- 2*sqrt(kmNrd$r[,2]))
	
	y <- mapply(function(i,j) kmRel$Sadj[i,j], w,1:length(w) ) # select time for each sample
	cir <- y
	cirLo <- y^exp( 2*sqrt(kmRel$r[,2]))
	cirUp <- y^exp( - 2*sqrt(kmRel$r[,2]))
	
	## Step 3: Compute post-relapse survival
	survPredict <- function(surv){
		s <- survfit(surv~1)
		splinefun(s$time, s$surv, method="monoH.FC")
	}
	xx <- 0:max(x)
	# Baseline Prs (measured from relapse)
	kmPrs0 <- survPredict(Surv(prdData$time1, prdData$time2, prdData$status))(xx) 
	# PRS baseline with spline-based dep on CR length)

	coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=prdData ) 
	tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
	rs <- sapply(1:nrow(data), function(i){
				### Different approach				
				xLen <- 1+floor(x)
				cir <- kmRel$Sadj[1:xLen,i]
				rs <- computeTotalPrsC(x = xx, diffCir = diff(cir), prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = kmPrd$r[i,1]-kmPrd$r0)
				rs[xLen]
			})
	
	## Step 4: Combine into overall survival
	if(any(1-(1-rs)-(1-nrs)<0)) warning("OS < 0 occured.")	
	os <- pmax(pmin(1-(1-rs)-(1-nrs),1),0)
	
	## Step 5: Confidence intervals for OS
	osCi <- sapply(1:nrow(data), function(i){
				if("analytical" == ciType){
					## Confidence intervals
					PlogP2 <- function(x) {(x * log(x))^2}
					errOs <- kmNrd$r[i,2] * PlogP2(kmNrd$S[w[i]]) * (1-kmRel$S[w[i]] * kmPrd$S[w[i]])^2 + kmRel$r[i,2]  * (1-kmNrd$S[w[i]])^2* kmPrd$S[w[i]]^2 * PlogP2(kmRel$S[w[i]]) +  kmPrd$r[i,2]  * (1-kmNrd$S[w[i]])^2* kmRel$S[w[i]]^2 * PlogP2(kmPrd$S[w[i]])
					errOs <- errOs / PlogP2(1-(1-kmNrd$S[w[i]])*(1-kmRel$S[w[i]]*kmPrd$S[w[i]]))
					return(c(osUp=os[i] ^ exp(-2* errOs), osLo= os[i] ^ exp(+2*errOs)))
				} else if("simulated" == ciType){
					## Simulate CI
					nSim <- 200
					osCiMc <- sapply(1:nSim, function(foo){
								H <- exp(rnorm(3,c(kmRel$r[i,1],kmNrd$r[i,1],kmPrd$r[i,1]),sqrt(c(kmRel$r[i,2],kmNrd$r[i,2],kmPrd$r[i,2]))))
								nrs <- cumsum(c(1,diff(kmNrd$S^H[2]) * kmRel$S[-1]^H[1])) ## Correct KM estimate for competing risk
								diffCir <- diff(kmRel$inc^H[1]) * kmNrd$inc[-1]^H[2] ## Correct KM estimate for competing risk							
								rs <- computeTotalPrsC(x = x, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrd$r0+log(H[3]))
								return((1-(1-nrs)-(1-rs))[w[i]])
							})
					osCiMcQ <- quantile(osCiMc, c(0.025,0.975))
					return(c(osUp = osCiMcQ[2], osLo = osCiMcQ[1]))
				}
			})
	
	return(data.frame(os=os, osLo = osCi[2,], osUp = osCi[1,],  cir=cir, cirLo=cirLo, cirUp=cirUp, nrs=nrs, nrsLo=nrsLo, nrsUp=nrsUp, rs=rs ))
}

#' Create a data.frame with all data in cr
#+ allPredict, cache=TRUE
allData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
allData$transplantCR1 <- allData$event
allData$transplantRel <- allData$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date)  
allData$transplantCR1[allData$index %in% w] <- 0
allData$transplantRel[!allData$index %in% w] <- 0

multiRFX3 <-  MultiRFX3(coxRFXNrdTD = coxRFXNrdTD, coxRFXPrdTD = coxRFXPrdTD, coxRFXRelTD = coxRFXRelTD, data=allData, x=3*365, prdData=prdData)


#' #### Model assessment
#' ##### Random cross-validation
#+concordanceCIRcv, cache=TRUE
replicates <- 100 ## number of replicates
concordanceCIRcv <- lapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			mclapply(1:replicates, function(foo){
						set.seed(foo)
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 4/5
						dNrm <- nrdData[nrdData$index %in% which(trainIdx),names(g)]
						sNrm <- Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% which(trainIdx)]
						coxRFXNrdTD <- CoxRFX(dNrm, sNrm, groups=g, nu=1, which.mu = mainGroups)
						coxRFXNrdTD$coefficients["transplantRel"] <- 0
						dPrs <- prdData[prdData$index %in% which(trainIdx), c(names(g),"time0","time1","time2","status")]
						sPrs <- Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% which(trainIdx)]
						coxRFXPrdTD <-  CoxRFX(dPrs, sPrs, groups=g, nu=1, which.mu = mainGroups)
						dCir <- relData[relData$index %in% which(trainIdx), names(g)]
						sCir <- Surv(relData$time1, relData$time2, relData$status)[relData$index %in% which(trainIdx)]
						coxRFXRelTD <-  CoxRFX(dCir, sCir, groups=g, which.mu = mainGroups)
						coxRFXRelTD$coefficients["transplantRel"] <- 0
						dOs <- osData[osData$index %in% which(trainIdx), names(g)]
						sOs <- Surv(osData$time1, osData$time2, osData$status)[osData$index %in% which(trainIdx)]
						coxRFXOsCR <- CoxRFX(dOs, sOs, groups=g, which.mu = mainGroups)
						
						allRisk365 <- MultiRFX3(coxRFXNrdTD = coxRFXNrdTD, coxRFXPrdTD = coxRFXPrdTD, coxRFXRelTD = coxRFXRelTD, data=allData, x=365, prdData=dPrs)
						allRisk1000 <- MultiRFX3(coxRFXNrdTD = coxRFXNrdTD, coxRFXPrdTD = coxRFXPrdTD, coxRFXRelTD = coxRFXRelTD, data=allData, x=1000, prdData=dPrs)
						
						p365 <- -allRisk365[,1]
						p1000 <-  -allRisk1000[,1]
						pCIR <- as.matrix(relData[names(g)]) %*% coef(coxRFXRelTD)
						pPRS <- as.matrix(prdData[names(g)]) %*% coef(coxRFXPrdTD)
						pNRM <- as.matrix(nrdData[names(g)]) %*% coef(coxRFXNrdTD)
						pOS <- as.matrix(osData[names(g)]) %*% coef(coxRFXOsCR)
						
						C <- c(
								CIRrfx = survConcordance(Surv(time1, time2, status)~ pCIR, data=relData, subset = relData$index %in% which(!trainIdx) )$concordance,
								PRSrfx = survConcordance(Surv(time1, time2, status) ~ pPRS, data=prdData, subset=prdData$index %in% which(!trainIdx) )$concordance,
								NRMrfx = survConcordance(Surv(time1, time2, status)~  pNRM, data=nrdData, subset=nrdData$index %in% which(!trainIdx) )$concordance,
								OSrfx = survConcordance(Surv(time1, time2, status) ~ pOS, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
								OS365 = survConcordance(Surv(time1, time2, status) ~ p365, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
								OS1000 = survConcordance(Surv(time1,time2, status) ~ p1000, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance
						)
						
						coef <- cbind(CIRrfx=coef(coxRFXRelTD), PRSrfx=coef(coxRFXPrdTD), NRMrfx=coef(coxRFXNrdTD),  OSrfx=coef(coxRFXOsCR))
						
						return(list(C=C, coef=coef, allRisk365=allRisk365, allRisk1000=allRisk1000))
					}, mc.cores=10)
		})

apply(apply(-sapply(concordanceCIRcv[[1]], `[[` , "C")[4:6,],2,rank),1,function(x) table(factor(x, levels=1:3)))
apply(apply(-sapply(concordanceCIRcv[[2]], `[[` , "C")[4:6,],2,rank),1,function(x) table(factor(x, levels=1:3)))

#' Test and train errors
i <- 0
concordanceCIRcvTrain <- lapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			i <- i+1
			sapply(1:replicates, function(foo){
						set.seed(foo)
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 4/5
						coef <- concordanceCIRcv[[i]][[foo]][["coef"]]
						pCIR <- as.matrix(relData[names(coef[,"CIRrfx"])]) %*% coef[,"CIRrfx"] 
						pPRS <- as.matrix(prdData[names(coef[,"PRSrfx"])]) %*% coef[,"PRSrfx"] 
						pNRM <- as.matrix(nrdData[names(coef[,"NRMrfx"])]) %*% coef[,"NRMrfx"]
						pOS <- as.matrix(osData[names(coef[,"OSrfx"])]) %*% coef[,"OSrfx"]
						p365 <- -concordanceCIRcv[[i]][[foo]][["allRisk365"]]$os
						p1000 <- -concordanceCIRcv[[i]][[foo]][["allRisk1000"]]$os
						C <- sapply(list(train=which(trainIdx), test=which(!trainIdx)), function(w)
									c(
											CIRrfx = survConcordance(Surv(time1, time2, status)~ pCIR, data=relData, subset = relData$index %in% w )$concordance,
											PRSrfx = survConcordance(Surv(time1, time2, status) ~ pPRS, data=prdData, subset=prdData$index %in% w )$concordance,
											NRMrfx = survConcordance(Surv(time1, time2, status)~  pNRM, data=nrdData, subset=nrdData$index %in% w )$concordance,
											OSrfx = survConcordance(Surv(time1, time2, status) ~ pOS, data=osData, subset=osData$index %in% w )$concordance,
											OS365 = survConcordance(Surv(time1, time2, status) ~ p365, data=osData, subset=osData$index %in% w )$concordance,
											OS1000 = survConcordance(Surv(time1,time2, status) ~ p1000, data=osData, subset=osData$index %in% w )$concordance
									))
						return(C)
					}, simplify='array')
		})

#' Plot test and training errors
#+ concordanceCIRcvTrainTest
for(i in 1:4){
	plot(t(concordanceCIRcvTrain[[2]][i,,] ), main=rownames(concordanceCIRcvTrain[[2]])[i])
	abline(0,1)
}

#' Plot coefficients v mean of subsampled coef
#+ concordanceCIRcvMeanCoef
r <- rowMeans(sapply(concordanceCIRcv[[2]], `[[` , "coef", simplify="array"), dim=2)
plot(r[,1], coef(coxRFXRelTD)); abline(0,1)
plot(r[,2], coef(coxRFXPrdTD)); abline(0,1)
plot(r[,3], coef(coxRFXNrdTD)); abline(0,1)
plot(r[,4], coef(coxRFXOsCR)); abline(0,1)


#' Variance-based concordance estimate
#+ concordanceCIRcvVar, cache=TRUE
i <- 0
concordanceCIRcvVar <- lapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			i <- i+1
			sapply(1:replicates, function(foo){
						set.seed(foo)
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 4/5
						coef <- concordanceCIRcv[[i]][[foo]][["coef"]]
						pCIR <- as.matrix(relData[names(coef[,"CIRrfx"])]) %*% coef[,"CIRrfx"] 
						pPRS <- as.matrix(prdData[names(coef[,"PRSrfx"])]) %*% coef[,"PRSrfx"] 
						pNRM <- as.matrix(nrdData[names(coef[,"NRMrfx"])]) %*% coef[,"NRMrfx"]
						pOS <- as.matrix(osData[names(coef[,"OSrfx"])]) %*% coef[,"OSrfx"]
						C <- sapply(list(train=which(trainIdx), test=which(!trainIdx)), function(w){
									sapply(ls(sys.frame(-3),pattern='^p[A-Z]+'), function(x)
												CoxHD:::ConcordanceFromVariance(var(get(x)[w], na.rm=TRUE)))[c(1,4,2,3)] 
								})
					}, simplify="array")})

for(i in 1:4)
{cat(rownames(concordanceCIRcvTrain[[2]])[i],"\n"); print(summary(data.frame(harrel=t(concordanceCIRcvTrain[[2]][i,1:2,]) , var=t(concordanceCIRcvVar[[2]][i,1:2,]))))}

#' ##### Cross-validation across trials
#+ concordanceCIRcvTrial, cache=TRUE
concordanceCIRcvTrial <- mclapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			mclapply(levels(clinicalData$Study), function(study){
						trainIdx <- clinicalData$Study != study
						g <- g[colSums(allData[allData$index %in% which(trainIdx), names(g)])>0]
						if(study == "AMLSG0704") g <- g[names(g) != "AMLHD98B"] # avoid collinearity
						dNrm <- nrdData[nrdData$index %in% which(trainIdx),names(g)]
						sNrm <- Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% which(trainIdx)]
						coxRFXNrdTD <- CoxRFX(dNrm, sNrm, groups=g, nu=1, which.mu = mainGroups)
						coxRFXNrdTD$coefficients["transplantRel"] <- 0
						dPrs <- prdData[prdData$index %in% which(trainIdx), c(names(g),"time0","time1","time2","status")]
						sPrs <- Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% which(trainIdx)]
						coxRFXPrdTD <-  CoxRFX(dPrs, sPrs, groups=g, nu=1, which.mu = mainGroups)
						dCir <- relData[relData$index %in% which(trainIdx), names(g)]
						sCir <- Surv(relData$time1, relData$time2, relData$status)[relData$index %in% which(trainIdx)]
						coxRFXRelTD <-  CoxRFX(dCir, sCir, groups=g, which.mu = mainGroups)
						coxRFXRelTD$coefficients["transplantRel"] <- 0
						dOs <- osData[osData$index %in% which(trainIdx), names(g)]
						sOs <- Surv(osData$time1, osData$time2, osData$status)[osData$index %in% which(trainIdx)]
						coxRFXOsCR <- CoxRFX(dOs, sOs, groups=g, which.mu = mainGroups)
						
						allRisk365 <- MultiRFX3(coxRFXNrdTD = coxRFXNrdTD, coxRFXPrdTD = coxRFXPrdTD, coxRFXRelTD = coxRFXRelTD, data=allData, x=365, prdData=dPrs)
						allRisk1000 <- MultiRFX3(coxRFXNrdTD = coxRFXNrdTD, coxRFXPrdTD = coxRFXPrdTD, coxRFXRelTD = coxRFXRelTD, data=allData, x=1000, prdData=dPrs)
						
						p365 <- -allRisk365[,1]
						p1000 <-  -allRisk1000[,1]
						pCIR <- as.matrix(relData[names(g)]) %*% coef(coxRFXRelTD)
						pPRS <- as.matrix(prdData[names(g)]) %*% coef(coxRFXPrdTD)
						pNRM <- as.matrix(nrdData[names(g)]) %*% coef(coxRFXNrdTD)
						pOS <- as.matrix(osData[names(g)]) %*% coef(coxRFXOsCR)
						
						C <- c(
								CIRrfx = survConcordance(Surv(time1, time2, status)~ pCIR, data=relData, subset = relData$index %in% which(!trainIdx) )$concordance,
								PRSrfx = survConcordance(Surv(time2 - time1, status) ~ pPRS, data=prdData, subset=prdData$index %in% which(!trainIdx) )$concordance,
								NRMrfx = survConcordance(Surv(time1, time2, status)~  pNRM, data=nrdData, subset=nrdData$index %in% which(!trainIdx) )$concordance,
								OSrfx = survConcordance(Surv(time1, time2, status) ~ pOS, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
								OS365 = survConcordance(Surv(time1, time2, status) ~ p365, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance,
								OS1000 = survConcordance(Surv(time1,time2, status) ~ p1000, data=osData, subset=osData$index %in% which(!trainIdx) )$concordance
						)
						
						coef <- cbind(RELrfx=coef(coxRFXRelTD), PRSrfx=coef(coxRFXPrdTD), NRSrfx=coef(coxRFXNrdTD),  OSrfx=coef(coxRFXOsCR))
						
						return(list(C=C, coef=coef, allRisk365=allRisk365, allRisk1000=allRisk1000))
					}, mc.cores=3)
		}, mc.cores=2)

dotplot(sapply(concordanceCIRcvTrial[[1]], `[[` , "C")[4:6,])
dotplot(sapply(concordanceCIRcvTrial[[2]], `[[` , "C")[4:6,])

#' ##### Test for TPL:Age interactions
#+tplAge
# CIR
c <- coxph(Surv(time1,time2,status)~transplantCR1*AOD_10, data=relData)
print(c)
anova(c)
#NRM
c <- coxph(Surv(time1,time2,status)~transplantCR1*AOD_10, data=nrdData)
print(c)
anova(c)
#PRS
c <- coxph(Surv(time1,time2,status)~ transplantRel*AOD_10, data=prdData)
print(c)
anova(c)

#' #### Absolute survival probabilities
#' This function computes the average accuracy of multiple absolute survival predictions at a given point in time by subdividing them
#' into equally sized bins and computing the weighted average absolute difference of the KM estimated survival probability and predicted.
EvalAbsolutePred <- function(prediction, surv, time, bins=seq(0,1,0.05)){
	c <- cut(prediction, bins)
	f <- survfit(surv ~ c)
	e <- summary(f, time)
	x <- sapply(strsplit(gsub("[a-z\\=\\(]|]","",e$strata),","), function(x) mean(as.numeric(x))); 
	#w <- 1/(e$std.err+.Machine$double.eps)^2
	w <- e$n[e$strata]
	std.err = 1/sum(w, na.rm=TRUE)
	mean.error = sum(abs(e$surv-x)*w, na.rm=TRUE)*std.err
	return(list(mean.error=mean.error, std.err=std.err, survfit=e, x=x))
}

#' #### Supplementary Figure 5
#+ absError
absPredError <- EvalAbsolutePred(multiRFX3$os, Surv(allData$time1, allData$time2, allData$status), time=3*365)

plot(absPredError$x, absPredError$survfit$surv, xlim=c(0,1), ylim=c(0,1), xlab="Predicted probability", ylab="Observed", main="Prediction tool")
segments(absPredError$x, absPredError$survfit$lower,absPredError$x, absPredError$survfit$upper)
abline(0,1)

PredictAbsoluteCoxph <- function(coxRFXOsCR, allData, time) {
	s <- survfit(coxRFXOsCR)
	q <- s$surv[which.min(abs(s$time-time))] ^ exp(predict(coxRFXOsCR, newdata=allData))
}
q <- PredictAbsoluteCoxph(coxRFXOsCR = coxRFXOsCR, allData = allData, time=365)

absPredErrorOs <- EvalAbsolutePred(q, Surv(allData$time1, allData$time2, allData$status), time=365)
plot(absPredErrorOs$x, absPredErrorOs$survfit$surv, xlim=c(0,1), ylim=c(0,1), xlab="Predicted probability", ylab="Observed", main="RFX on OS")
segments(absPredErrorOs$x, absPredErrorOs$survfit$lower,absPredErrorOs$x, absPredErrorOs$survfit$upper)
abline(0,1)

#' Eval cross-validated samples
#+ absoluteErrorsCIRcv, cache=TRUE
i <- 0
absoluteErrorsCIRcv <- lapply(list(crGroups[crGroups %in% mainGroups], crGroups), function(g){ 
			i <- i+1
			sapply(1:replicates, function(foo){
						set.seed(foo)
						time <- 365
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 4/5
						coef <- concordanceCIRcv[[i]][[foo]][["coef"]]
						
						lpCIR <- as.matrix(relData[names(coef[,"CIRrfx"])]) %*% coef[,"CIRrfx"]
						s <- survfit(Surv(time1, time2, status)~1, data=relData, subset=relData$index %in% which(trainIdx))
						pCIR <- s$surv[which.min(abs(s$time-time))] ^ exp(lpCIR-mean(lpCIR[relData$index %in%which(trainIdx)]))
						
						lpPRS <- as.matrix(prdData[names(coef[,"PRSrfx"])]) %*% coef[,"PRSrfx"] 
						s <- survfit(Surv(time2- time1, status)~1, data=prdData, subset=prdData$index %in% which(trainIdx))
						pPRS <- s$surv[which.min(abs(s$time-time))] ^ exp(lpPRS-mean(lpPRS[prdData$index %in% which(trainIdx)]))
						
						lpNRM <- as.matrix(nrdData[names(coef[,"NRMrfx"])]) %*% coef[,"NRMrfx"]
						s <- survfit(Surv(time1, time2, status)~1, data=nrdData, subset=nrdData$index %in% which(trainIdx))
						pNRM <- s$surv[which.min(abs(s$time-time))] ^ exp(lpNRM-mean(lpNRM[nrdData$index %in% which(trainIdx)]))
						
						lpOS <- as.matrix(osData[names(coef[,"OSrfx"])]) %*% coef[,"OSrfx"]
						s <- survfit(Surv(time1, time2, status)~1, data=osData, subset=osData$index %in% which(trainIdx))
						pOS <- s$surv[which.min(abs(s$time-time))] ^ exp(lpOS-mean(lpOS[osData$index %in% which(trainIdx)]))
						
						p365 <- concordanceCIRcv[[i]][[foo]][["allRisk365"]]$os
						p1000 <- concordanceCIRcv[[i]][[foo]][["allRisk1000"]]$os
						err <- sapply(list(train=which(trainIdx), test=which(!trainIdx)), function(w)
									c(
											CIRrfx = EvalAbsolutePred(pCIR[relData$index %in% w ], Surv(relData$time1, relData$time2, relData$status)[relData$index %in% w ], time=365)$mean.error,
											PRSrfx =  EvalAbsolutePred(pPRS[prdData$index %in% w ], Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% w ], time=365)$mean.error,
											NRMrfx =  EvalAbsolutePred(pNRM[nrdData$index %in% w ], Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% w ], time=365)$mean.error,
											OSrfx =  EvalAbsolutePred(pOS[osData$index %in% w ], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% w ], time=365)$mean.error,
											OS365 =  EvalAbsolutePred(p365[osData$index %in% w ], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% w ], time=365)$mean.error,
											OS1000 =  EvalAbsolutePred(p1000[osData$index %in% w ], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% w ], time=1000)$mean.error
									))
						return(err)
					}, simplify='array')
		})

summary(t(absoluteErrorsCIRcv[[2]][,1,]))
boxplot(t(absoluteErrorsCIRcv[[2]][,1,]), main="Training")

summary(t(absoluteErrorsCIRcv[[2]][,2,]))
boxplot(t(absoluteErrorsCIRcv[[2]][,2,]), main="Test")

#' ##### KM curves after remission
riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)

#+ mortality, fig.width=3, fig.height=3
i <- 1
rsStatus <- osData$status
rsStatus[osData$index %in% nrdData$index[nrdData$status==1]] <- 0
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	plot(NA,NA, ylim=c(0,1),  xlab="Years", ylab="Mortality", xlim=c(0,10), yaxs='i', xaxs='i')
	abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	#abline(v=seq(1,9), col='lightgrey')
	lines(survfit(Surv(time1/365, time2/365, status) ~ clinicalData$M_Risk[osData$index], data=osData, subset=clinicalData$M_Risk[osData$index]==l), col=riskCol[l],fun=function(x) 1-x ,mark=NA, lty=1, conf.int=FALSE)
	rsKM <- survfit(Surv(time1/365, time2/365, rsStatus) ~ 1, data=osData, subset=  clinicalData$M_Risk[osData$index]==l)
	nrsKM <- survfit(Surv(time1/365, time2/365, status) ~ 1, data=nrdData, subset=  clinicalData$M_Risk[nrdData$index]==l)
	
	rsCR <- cumsum(c(1,diff(rsKM$surv)) * splinefun(nrsKM$time, nrsKM$surv, method="monoH.FC")(rsKM$time))
	nrsCR <- cumsum(c(1,diff(nrsKM$surv)) * splinefun(rsKM$time, rsKM$surv, method="monoH.FC")(nrsKM$time))
	
	
	lines(rsKM$time, 1-rsCR, col=riskCol[l], lty=2, type='s')
	lines(nrsKM$time, 1-nrsCR, col=riskCol[l], lty=3, type='s')
	if(i ==1)
		legend(ifelse(i<=3,"topleft","bottomright"), c("total","relapse","non-rel"), lty=c(1,2,3), col="black", box.lty = 0, bg="white")
	i <- i+1
	mtext(l, side=3, font=2)
}

f <- function(x) 1-x
plot(survfit(Surv(time1/365, time2/365, status) ~ clinicalData$M_Risk[relData$index], data=relData), col=riskCol, ylab="CIR", xlab="Time after CR", main="Molecular risk groups, all cases", fun=f , ylim=c(0,1))
legend("bottomright", lty=1, bty="n", paste(levels(clinicalData$M_Risk), table(clinicalData$M_Risk[!is.na(c)])), col=riskCol)

#' Incidence of relapse v risk tercile
#+ cirSplits, fig.width=3, fig.height=3
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
riskCirTD <- coxRFXRelTD$Z %*% coef(coxRFXRelTD) - relData$transplantCR1 * coef(coxRFXRelTD)["transplantCR1"]
quantileRiskCirTD <- numeric(nrow(relData))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	w <- which(clinicalData$M_Risk[relData$index]==l)
	q <- cut(riskCirTD[w], quantile(riskCirTD[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
	quantileRiskCirTD[w] <- q
	plot(NA,NA,  ylab="CIR", main=paste(l, "terciles"),  xlab="Years after CR", ylim=c(0,1), xlim=c(0,10), xaxs="i", yaxs="i")
	#abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	fit <- survfit(Surv(time1/365, time2/365, status) ~ q + transplantCR1, data=relData[w,])
	## adjust for competing risk (NRM)
	i <- c(0,diff(fit$surv))
	s <- split(fit$surv, cumsum(i>0)) # split into strata
	u <- split(fit$upper, cumsum(i>0)) # split into strata
	v <- split(fit$lower, cumsum(i>0)) # split into strata
	
	t <- split(fit$time, cumsum(i>0))
	nrsKM <- survfit(Surv(time1/365, time2/365, status) ~ 1, data=nrdData, subset=  clinicalData$M_Risk[nrdData$index]==l)
	
	fit$surv <- unlist(sapply(seq_along(s), function(i) cumsum(c(1,diff(s[[i]])) * splinefun(nrsKM$time, nrsKM$surv, method="monoH.FC")(t[[i]])))) #adjust increments by nrs KM est
	fit$lower <- unlist(sapply(seq_along(s), function(i) cumsum(c(1,diff(v[[i]])) * splinefun(nrsKM$time, nrsKM$surv, method="monoH.FC")(t[[i]])))) #adjust increments by nrs KM est
	fit$upper <- unlist(sapply(seq_along(s), function(i) cumsum(c(1,diff(u[[i]])) * splinefun(nrsKM$time, nrsKM$surv, method="monoH.FC")(t[[i]])))) #adjust increments by nrs KM est
	lines(fit, col=rep(sapply(2:0,function(x) colTrans(riskCol[l],x)), each=2), lty=c(1,0), mark=NA, xlab="Time after CR", fun=f)
	#legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[l])
}

#' #### Figure 4a
#' We use the `survival` package to compute the following mstate fits of CIR and NRM
t <- clinicalData$Recurrence_date
t[is.na(t)] <-  clinicalData$Date_LF[is.na(t)]
time <- as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date)
status <- factor(ifelse(!is.na(clinicalData$Recurrence_date), "relapse", ifelse(clinicalData$Status==1,"dead","alive"  ))) 
status[is.na(clinicalData$CR_date)] <- NA
alloCR1 <- 1:1540 %in% osData$index[osData$transplantCR1==1]
mSurv <- Surv(time/365.25, status, type="mstate")

f <- function(x) 1-x
#+ cirSplitsCR, fig.width=3, fig.height=3
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
riskCir <- (coxRFXRelTD$Z %*% coef(coxRFXRelTD) - relData$transplantCR1 * coef(coxRFXRelTD)["transplantCR1"])[1:1540] # Risk w/o allograft
qtl <- numeric(nrow(dataFrame))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	w <- which(clinicalData$M_Risk==l)
	q <- cut(riskCir[w], quantile(riskCir[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
	qtl[w] <- q
	plot(NA,NA,  ylab="Fraction relapsed", main=paste0(l,", n=",sum(clinicalData$M_Risk[!is.na(mSurv)]==l, na.rm=TRUE)),  xlab="Years after CR", ylim=c(0,1), xlim=c(0,5), xaxs="i", yaxs="i", font.main=1)
	#abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	fit <- survfit(mSurv~ qtl, subset= clinicalData$M_Risk==l)
	
	lines(fit, col=sapply(2:0, function(x) c(colTrans(set1[2],x), colTrans(set1[5],x))), lty=c(1,1), mark=NA, xlab="Time after CR", fun=f)
	#legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[l])
}


#' Overall survival after remission v risk tercile
#+ osSplits, fig.width=3, fig.height=3
par(mfrow=c(2,2), mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0))
riskOsCR <- coxRFXOsCR$Z %*% coef(coxRFXOsCR) - osData$transplantCR1 * coef(coxRFXOsCR)["transplantCR1"]
quantileRiskOsCR <- numeric(nrow(osData))
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	w <- which(clinicalData$M_Risk[osData$index]==l)
	q <- cut(riskOsCR[w], quantile(riskOsCR[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
	quantileRiskOsCR[w] <- q
	plot(NA,NA,  ylab="OS", main=paste(l, "terciles"),  xlab="Years after CR", ylim=c(0,1), xlim=c(0,10), xaxs="i", yaxs="i")
	abline(h=seq(0.2,0.8,0.2),lty=1, col='lightgrey')
	fit <- survfit(Surv(time1/365, time2/365, status) ~ q + transplantCR1, data=osData[w,])
	lines(fit, col=rep(sapply(2:0,function(x) colTrans(riskCol[l],x)), each=2), lty=c(1,3), mark=NA, xlab="Time after CR")
	legend("bottomright", lty=c(1,3), bty="n", c("no TPL","TPL"), col=riskCol[l])
}

#' ##### Risk factors of relapse and survival
p <- lapply(levels(clinicalData$M_Risk), function(l) {
			w <- which(clinicalData$M_Risk==l)
			q <- cut(riskOsCR[w], quantile(riskOsCR[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
			sapply(split(relData[w, names(crGroups)], q), colMeans)
		})
names(p) <- levels(clinicalData$M_Risk)

#+relapseFactors, fig.width=5,fig.heigh=4
par(mfrow=c(4,1), xpd=NA)
for(l in levels(clinicalData$M_Risk)[c(2,4,3,1)]){
	t <- t((p[[l]]) * coef(coxRFXRelTD))[,crGroups != "Treatment"]
	z <- (coef(coxRFXRelTD)/sqrt(diag(coxRFXRelTD$var2)))[crGroups != "Treatment"]
	o <- order(z)
	w <- c(1:15,ncol(t)-14:0)
	b <- barplot(t[,o][,w], las=2, col=sapply(2:0,function(x) colTrans(riskCol[l],x)), beside=TRUE, ylim=c(-.5,.5), names.arg=rep("", length(w)))
	rotatedLabel(b[2,],pmin(0,apply(t,2,min)[o][w]), colnames(t)[o][w])
	s <- matrix(rep(sqrt(diag(coxRFXRelTD$var2)[crGroups != "Treatment"]), each=3) * t/rep(coef(coxRFXRelTD)[crGroups != "Treatment"], each=3), nrow=3)[,o][,w]
	segments(b[1,], (colMeans(coxRFXRelTD$Z)*coef(coxRFXRelTD))[crGroups != "Treatment"][o][w] ,b[3,], (colMeans(coxRFXRelTD$Z)*coef(coxRFXRelTD))[crGroups != "Treatment"][o][w])
	segments(b,t[,o][,w]-s, b,t[,o][,w]+s)
}

#p <- as.data.frame(PartialRisk(coxRFXRelTD)[1:nrow(clinicalData),])
partialRiskOsCR <- as.data.frame(PartialRisk(coxRFXOsCR)[1:nrow(clinicalData),])
s <- do.call("rbind",lapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(l) {
					w <- which(clinicalData$M_Risk==l)
					q <- cut(riskOsCR[w], quantile(riskOsCR[w], seq(0,1,.33)), include.lowest=TRUE, labels=c("T1","T2","T3"))
					t(sapply(split(partialRiskOsCR[w, ], q), colMeans) +.5 - colMeans(partialRiskOsCR))
				}))

#' Risk constellation for OS after remission
#+relapseStars, fig.width=3,fig.heigh=3
c <- sapply(2:0, function(t) sapply(riskCol[c(2,4,3,1)], function(c) colTrans(c,t)))
g <- expand.grid(1:3,1:4-1)*3
stars(2*s[,c("Clinical","Demographics","Genetics","GeneGene","CNA","Fusions","Treatment")], scale=FALSE, col.stars=t(c), key.loc = c(13,0), locations=g, labels=NA)
symbols(g[,1], g[,2], circles=rep(1,12), inches=FALSE, add=TRUE)
text(1, 0:3*3, names(riskCol[c(2,4,3,1)]), pos=2)
text(1:3*3, 11, c("Best","Intermediate","Worst"), pos=3)

#' Prototypical risk constellations
prototypes <- sapply(levels(clinicalData$M_Risk)[c(2,4,3,1)], function(l) sapply(1:3, function(i){
						#d <- dist(as.data.frame(coxRFXRelTD$Z[which(clinicalData$M_Risk[cirData$index]==l & quantileRiskCirTD==i &! is.na(clinicalData$CR_date[cirData$index])), ])) 
						w <- which(clinicalData$M_Risk[relData$index]==l & quantileRiskOsCR==i &! is.na(clinicalData$CR_date[relData$index]))
						d <- dist(t(t(coxRFXOsCR$Z[w, ]) ))
						osData$index[w][which.min(rowMeans(as.matrix(d), na.rm=TRUE))]
					}))

c <- sapply(2:0, function(t) sapply(riskCol[c(2,4,3,1)], function(c) colTrans(c,t)))
g <- expand.grid(1:3,1:4-1)*3
stars(2*t(t(partialRiskOsCR[prototypes,])- colMeans(partialRiskOsCR))[,c("Clinical","Demographics","Genetics","CNA","Fusions","Treatment")] +1, scale=FALSE, col.stars=t(c), key.loc = c(13,0), locations=g, labels=NA)
symbols(g[,1], g[,2], circles=rep(1,12), inches=FALSE, add=TRUE)
text(1, 0:3*3, names(riskCol[c(2,4,3,1)]), pos=2)
text(1:3*3, 11, c("Low","Intermediate","High"), pos=3)

#+ starsOS, fig.width=4, fig.height=4
s <- partialRiskOsCR - rep(colMeans(partialRiskOsCR), each=nrow(partialRiskOsCR))
w <- sapply(split(1:1540, paste(clinicalData$M_Risk, quantileRiskOsCR[1:1540])), `[`, 1:12)
w <- w[,!grepl("NA", colnames(w))][,c(4:6,10:12,7:9,1:3)]
l <- stars(s[w,c("Demographics","Treatment","Fusions","CNA","Genetics","GeneGene","Clinical")] + .5, scale=FALSE, col.stars = mapply(function(i,j) {t <- try(c[i,j]); if(class(t)=="try-error") NA else t}, as.character(clinicalData$M_Risk[w]),quantileRiskCirTD[w]), labels="")
symbols(l[,1],l[,2], circles=rep(0.5, nrow(l)), inches=FALSE,add=TRUE)

#+ starsCIR, fig.width=4, fig.height=4
layout(matrix(c(1:4), ncol=2),heights = c(10,1), widths = c(10,1))
partialRiskCirTD <- as.data.frame(PartialRisk(coxRFXRelTD))
s <- partialRiskCirTD[1:nrow(clinicalData),] - rep(colMeans(partialRiskCirTD), each=nrow(clinicalData))
u <- unique(relData$index[!is.na(relData$time2)])
w <- sapply(split(u, paste(clinicalData$M_Risk, quantileRiskCirTD[1:1540])[u]), `[`, 1:12)
w <- w[,!grepl("NA", colnames(w))][,c(4:6,10:12,7:9,1:3)]
i <- which(rev(!duplicated(rev(relData$index))))
m <- i[order(relData$index[i])]
c <- cut(relData$time2, quantile(relData$time2[m], seq(0,1,0.1), na.rm=TRUE))
l <- mg14:::stars(s[w,c("Demographics","Treatment","Fusions","CNA","Genetics","GeneGene","Clinical")] + .5, scale=FALSE, col.stars = brewer.pal(11,"RdBu")[-6][c[w]], labels="", density=ifelse(relData$status[m][w]==1,NA,48),  col.lines=rep(1,(12^2)))
symbols(l[,1],l[,2], circles=rep(0.5, nrow(l)), inches=FALSE,add=TRUE, fg='lightgrey')
par(mar=c(2,2,0,2))
barplot(matrix(diff(quantile(relData$time2[m], na.rm=T, seq(0,1,0.1))), ncol=1)/365.25, col=brewer.pal(11,"RdBu")[-6], horiz=TRUE, border=NA, xlim=c(0,20))


#' #### Allogeneic hematopoietic stem cell transplants
#' Create a data.frame with all possibilities for allografts - none, CR1, after relapse.
#+survivalTpl, cache=TRUE
w <- sort(unique(osData$index[which(quantileRiskOsCR==3 & clinicalData$M_Risk[osData$index]=="Favorable")]))
allDataTpl <- osData[rep(1:nrow(dataFrame), each=3),]
allDataTpl$transplantCR1 <- rep(c(0,1,0), nrow(dataFrame))
allDataTpl$transplantRel <- rep(c(0,0,1), nrow(dataFrame))
multiRFX3Tpl <- MultiRFX3(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data=allDataTpl, x=3*365, prdData=prdData)
multiRFX3Tpl <- as.data.frame(matrix(multiRFX3Tpl$os, ncol=3, byrow=TRUE, dimnames=list(NULL, c("None","CR1","Relapse"))), row.names=rownames(dataFrame))
survivalTpl <- data.frame(multiRFX3Tpl, os=osYr, age=clinicalData$AOD, ELN=clinicalData$M_Risk, tercile=quantileRiskOsCR[1:nrow(multiRFX3Tpl)])

#+ survivalTplOut, results='asis'
datatable(format(survivalTpl[order(survivalTpl$CR1 -survivalTpl$Relapse),], digits=4))
datatable(multiRFX3Tpl[patients,])

#' Function to predict OS from  Relapse, PRS and NRM. This one also computes confidence intervals for each type of allograft and the predicted differences in outcome between allograft types.
MultiRFX3TplCi <- function(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, x =365, prsData, ciType="simulated", nSim = 200, mc.cores=10){
	## Step 1: Compute KM survival curves and log hazard
	getS <- function(coxRFX, data, max.x=5000) {		
		if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
		data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data))])
		r <- PredictRiskMissing(coxRFX, data, var="var2")
		H0 <- basehaz(coxRFX, centered = FALSE)
		hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
		x <- c(0:max.x,max.x)
		S <- exp(-hazardDist(x))
		return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
	}
	
	data$transplantCR1 <- 0
	data$transplantRel <- 0
	
	kmCir <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
	kmNrs <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
	kmPrs <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
	
	survPredict <- function(surv){
		s <- survfit(surv~1)
		splinefun(s$time, s$surv, method="monoH.FC")
	}
	xx <- 0:max(x)
	
	# Baseline Prs (measured from relapse)
	kmPrs0 <- survPredict(Surv(prdData$time1, prdData$time2, prdData$status))(xx) 
	
	# PRS baseline with spline-based dep on CR length)
	coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=prdData) 
	tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
	
	stopifnot(length(x)==1 | length(x) == nrow(data))
	if(length(x)==nrow(data))
		w <- match(x,kmCir$x)
	else if(length(x)==1)
		w <- rep(match(x, kmCir$x), nrow(data))
	
	survival <- sapply(c("None","Rel","CR1"), function(type){
				if(type=="None"){
					data$transplantCR1 <- 0
					data$transplantRel <- 0
				}else if(type=="Rel"){
					data$transplantCR1 <- 0
					data$transplantRel <- 1					
				}else if(type=="CR1"){
					data$transplantCR1 <- 1
					data$transplantRel <- 0
				}
				
				
				kmCir <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
				kmNrm <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
				kmPrs <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
				
				## Step 2: Adjust CIR and NRM curve for competing risks, accounting for hazard
				kmCir$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmCir$S^exp(kmCir$r[i,1]))) * kmNrm$S ^ exp(kmNrm$r[i,1])))
				kmNrm$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmNrm$S^exp(kmNrm$r[i,1]))) * kmCir$S ^ exp(kmCir$r[i,1]))) ## array times x nrow(data)
				
				y <- mapply(function(i,j) kmNrm$Sadj[i,j], w,1:length(w) ) # select time for each sample
				nrs <- y
				nrsUp <- y^exp(2*sqrt(kmNrm$r[,2]))
				nrsLo <- y^exp(- 2*sqrt(kmNrm$r[,2]))
				
				y <- mapply(function(i,j) kmCir$Sadj[i,j], w,1:length(w) ) # select time for each sample
				cir <- y
				cirLo <- y^exp( 2*sqrt(kmCir$r[,2]))
				cirUp <- y^exp( - 2*sqrt(kmCir$r[,2]))
				
				## Step 3: Compute post-relapse survival
				rs <- sapply(1:nrow(data), function(i){
							### Different approach				
							xLen <- 1+floor(x)
							cir <- kmCir$Sadj[1:xLen,i]
							rs <- computeTotalPrsC(x = xx, diffCir = diff(cir), prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = kmPrs$r[i,1]-kmPrs$r0)
							rs[xLen]
						})
				
				## Step 4: Combine into overall survival
				if(any(1-(1-rs)-(1-nrs)<0)) warning("OS < 0 occured.")	
				os <- pmax(pmin(1-(1-rs)-(1-nrs),1),0)
				cbind(os, rs, nrs, aar=rs-cir)
			}, simplify='array')
	
	## Step 5: Confidence intervals for OS
	osCi <- sapply(mclapply(1:nrow(data), function(i){
						{
							## Simulate CI
							osCiMc <- sapply(1:nSim, function(foo){
										r0 <- rnorm(3,c(kmCir$r[i,1],kmNrs$r[i,1],kmPrs$r[i,1]),sqrt(c(kmCir$r[i,2],kmNrs$r[i,2],kmPrs$r[i,2])))
										H0 <- exp(r0)
										nrs0 <- cumsum(c(1,diff(kmNrs$S^H0[2])) * kmCir$S^H0[1]) ## Correct KM estimate for competing risk
										diffCir <- diff(c(1,kmCir$S^H0[1])) * kmNrs$S^H0[2] ## Correct KM estimate for competing risk			
										cir0 <- 1+cumsum(diffCir)
										rs0 <- computeTotalPrsC(x = xx, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrs$r0+log(H0[3]))
										aar0 <- rs0-cir0
										
										Hcr1 <- exp(r0 + rnorm(3,c(coxRFXRelTD$coefficients["transplantCR1"],coxRFXNrdTD$coefficients["transplantCR1"],coxRFXPrdTD$coefficients["transplantCR1"]), 
														sqrt(c(coxRFXRelTD$var2["transplantCR1","transplantCR1"],coxRFXNrdTD$var2["transplantCR1","transplantCR1"],coxRFXPrdTD$var2["transplantCR1","transplantCR1"])))) 
										nrsCr1 <- cumsum(c(1,diff(kmNrs$S^Hcr1[2])) * kmCir$S^Hcr1[1]) ## Correct KM estimate for competing risk
										diffCir <- diff(c(1,kmCir$S^Hcr1[1])) * kmNrs$S^Hcr1[2] ## Correct KM estimate for competing risk	
										cirCr1 <- 1+cumsum(diffCir)
										rsCr1 <- computeTotalPrsC(x = xx, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrs$r0+log(Hcr1[3]))
										aarCr1 <- rsCr1-cirCr1
										
										
										Hrel <- exp(r0 + rnorm(3,c(coxRFXRelTD$coefficients["transplantRel"],coxRFXNrdTD$coefficients["transplantRel"],coxRFXPrdTD$coefficients["transplantRel"]), 
														sqrt(c(coxRFXRelTD$var2["transplantRel","transplantRel"],coxRFXNrdTD$var2["transplantRel","transplantRel"],coxRFXPrdTD$var2["transplantRel","transplantRel"])))) 
										nrsRel <- cumsum(c(1,diff(kmNrs$S^Hrel[2])) * kmCir$S^Hrel[1]) ## Correct KM estimate for competing risk
										diffCir <- diff(c(1,kmCir$S^Hrel[1])) * kmNrs$S^Hrel[2] ## Correct KM estimate for competing risk							
										cirRel <- 1+cumsum(diffCir)
										rsRel <- computeTotalPrsC(x = xx, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrs$r0+log(Hrel[3]))
										aarRel <- rsRel-cirRel
										
										
										os0 <- (1-(1-nrs0[1:w[i]])-(1-rs0))[w[i]]
										osCr1 <- (1-(1-nrsCr1[1:w[i]])-(1-rsCr1))[w[i]]
										osRel <- (1-(1-nrsRel[1:w[i]])-(1-rsRel))[w[i]]
										return(cbind(os=c(none=os0, cr1=osCr1, rel=osRel, dCr1=osCr1-os0, dRel=osRel-os0, dCr1Rel=osCr1-osRel),
														rs=c(none=rs0[w[i]], cr1=rsCr1[w[i]], rel=rsRel[w[i]], dCr1=rsCr1[w[i]]-rs0[w[i]], dRel=rsRel[w[i]]-rs0[w[i]], dCr1Rel=rsCr1[w[i]]-rsRel[w[i]]),
														nrs=c(none=nrs0[w[i]], cr1=nrsCr1[w[i]], rel=nrsRel[w[i]], dCr1=nrsCr1[w[i]]-nrs0[w[i]], dRel=nrsRel[w[i]]-nrs0[w[i]], dCr1Rel=nrsCr1[w[i]]-nrsRel[w[i]]),
														aar=c(none=aar0[w[i]], cr1=aarCr1[w[i]], rel=aarRel[w[i]], dCr1=aarCr1[w[i]]-aar0[w[i]], dRel=aarRel[w[i]]-aar0[w[i]], dCr1Rel=aarCr1[w[i]]-aarRel[w[i]])))
									}, simplify='array')
							osCiMcQ <- apply(osCiMc,1:2,quantile, c(0.025,0.5,0.975))
							return(sapply(c("os","rs","nrs","aar"), function(t) 
												cbind(hat = c(survival[i,t,1], survival[i,t,3], survival[i,t,2], survival[i,t,3]-survival[i,t,1], survival[i,t,2]-survival[i,t,1], survival[i,t,3]-survival[i,t,2]), 
														median = osCiMcQ[2,,t], lower = osCiMcQ[1,,t], upper = osCiMcQ[3,,t]), simplify="array"))
						}
					}, mc.cores=mc.cores), I, simplify="array")
	#cat(os, "\n")
	return(osCi)
}

#+predictOsTplCi, cache=TRUE
set.seed(42)
d <- osData[1:nrow(dataFrame),]
d$transplantCR1 <- 0
d$transplantRel <- 0
p <- grep("PD11104a|PD8314a|PD11080a",rownames(dataFrame))
predict3 <- MultiRFX3TplCi(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data=d[p,colnames(coxRFXNrdTD$Z)], x=3*365, nSim=1000, prsData=prsData) ## selected with 1000
dimnames(predict3)[[4]] <- rownames(dataFrame)[p]
predict3
set.seed(42)
multiRFX3TplCi <- MultiRFX3TplCi(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data=d[,colnames(coxRFXNrdTD$Z)], x=3*365, nSim=200, prsData=prsData) ## others with 200
dimnames(multiRFX3TplCi)[[4]] <- rownames(dataFrame)

#' The figure shows the mortality reduction of allograft CR1 v none, allograft in Rel v none, and CR1 v Relapse.
#+mortalityReduction, fig.width=3.5, fig.height=2.5
set.seed(42)
par(mar=c(3,3,1,3), las=2, mgp=c(2,.5,0), bty="n")
patients <- grep("PD11104a|PD8314a",rownames(dataFrame))
for(t in c("dCr1","dRel","dCr1Rel")){
	s <- clinicalData$AOD < 60 #sample(1:1540,100)
	x <- 1-multiRFX3TplCi["none","hat","os",]
	y <- multiRFX3TplCi[t,"hat","os",]
	plot(x[s], y[s], pch=NA, ylab="Mortality reduction from allograft", xlab="3yr mortality with standard chemo", col=riskCol[clinicalData$M_Risk], cex=.8, las=1, ylim=range(multiRFX3Tpl$CR1-multiRFX3Tpl$None))
	abline(h=seq(-.1,.3,.1), col='grey', lty=3)
	abline(v=seq(.2,.9,0.2), col='grey', lty=3)
	points(x[s], y[s], pch=16,  col=riskCol[clinicalData$M_Risk[s]], cex=.8)
	segments(1-multiRFX3TplCi["none","lower","os",patients], y[patients],1-multiRFX3TplCi["none","upper","os",patients],y[patients])
	segments(x[patients], multiRFX3TplCi[t,"lower","os",patients],x[patients], multiRFX3TplCi[t,"upper","os",patients])
	xn <- seq(0,1,0.01)
	p <- predict(loess(y~x, data=data.frame(x=x[s], y=y[s])), newdata=data.frame(x=xn), se=TRUE)
	yn <- c(p$fit + 2*p$se.fit,rev(p$fit - 2*p$se.fit))
	polygon(c(xn, rev(xn))[!is.na(yn)],yn[!is.na(yn)], border=NA, col="#00000044", lwd=1)
	lines(xn,p$fit, col='black', lwd=2)
	legend("topleft", pch=c(16,16,16,16,NA),lty=c(NA,NA,NA,NA,1), col=c(riskCol[c(2,4,3,1)],1),fill=c(NA,NA,NA,NA,"grey"), border=NA, c(names(riskCol)[c(2,4,3,1)],"loess average"), box.lty=0)
	n <- c(100,50,20,10,5,4,3)
	axis(side=4, at=1/n, labels=n, las=1)
	mtext("Number needed to treat", side=4, at=.2, line=2, las=0)
	axis(side=4, at=-1/n, labels=n, las=1)
	mtext("Number needed to harm", side=4, at=-.1, line=2, las=0)
}

#' The following shows boxplots of the mortality reduction v the risk terciles.
#+survivalTplBoxPlot
par(mar=c(7,5,1,1))
f <- factor(clinicalData$M_Risk, levels=levels(clinicalData$M_Risk)[c(2,4,3,1)])
boxplot(multiRFX3Tpl$CR1 - multiRFX3Tpl$None ~ quantileRiskOsCR[1:1540]  + f, las=2, col=t(outer(riskCol[c(2,4,3,1)], 2:0, colTrans)), ylab="Survival gain TPL CR1 at 3yr")
boxplot(multiRFX3Tpl$Relapse - multiRFX3Tpl$None ~ quantileRiskOsCR[1:1540]  + f, las=2, col=t(outer(riskCol[c(2,4,3,1)], 2:0, colTrans)), ylab="Survival gain TPL Relapse at 3yr")
boxplot(multiRFX3Tpl$CR1 - multiRFX3Tpl$Relapse ~ quantileRiskOsCR[1:1540] + f, las=2, col=t(outer(riskCol[c(2,4,3,1)], 2:0, colTrans)), ylab="Survival gain TPL in CR1 over salvage at 3yr")
abline(h=0)

#' Mortality reduction v age
par(c(3,3,1,1))
y <- multiRFX3Tpl$CR1 - multiRFX3Tpl$None
x <- dataFrame$AOD_10*10
plot(y ~ x)
lines(lowess(x[x<60], y[x<60]), col="green")
#' Note: The jump after 60 arises from patients after 60 in AMLHD98B not having received allografts. Based on the trial stratum they are hence (incorrectly) predicted to have very low non-relapse mortality upon allograft. However,
#' this doesn't affect novel patients.

plot(multiRFX3Tpl$CR1 - multiRFX3Tpl$None ~ predict(coxRFXOsCR, newdata=osData[1:1540,]), xlab="Risk", ylab="Survival gain TPL CR1 at 1000d")
lines(lowess(predict(coxRFXOsCR, newdata=osData[1:1540,]), multiRFX3Tpl$CR1 - multiRFX3Tpl$None), col='green')


#' #### Best treatment options
#' We can explore the the hypothetical survival gain if each patient had received the optimal treatment strategy.
#' Distribution of ranks.
apply(apply(-multiRFX3Tpl,1,rank),1,function(x) table(factor(x, levels=1:3)))

#' Split by ELN risk
table(clinicalData$M_Risk, factor(apply(multiRFX3Tpl, 1, which.max), levels=1:3, labels=colnames(multiRFX3Tpl)))[c(2,4,3,1),]

#' Split by ELN risk, requiring TPL in CR1 to offer 5% advantage over salvage
table(clinicalData$M_Risk, apply(multiRFX3Tpl, 1, function(x) x[2] > x[3]+.05))[c(2,4,3,1),]

#' Observed outcome
summary(survfit(Surv(time1, time2, status) ~ 1, data=osData), time=3*365)

#' Under different scenarios
colMeans(multiRFX3Tpl[!is.na(clinicalData$CR_date),])

#' Using observed 
mean(sapply((1:nrow(data))[!is.na(clinicalData$CR_date) & clinicalData$AOD < 60], function(i) multiRFX3Tpl[i, 1+data[i,"transplantCR1"] + 2*data[i, "transplantRel"] ]))

#' Best possible - everyone had received the optimum 
mean(apply(multiRFX3Tpl[!is.na(clinicalData$CR_date),],1,max))

#' Benefit v number of allografts in CR1
#+ survNallo
par(bty="L")
o <- order(-multiRFX3TplCi["dCr1Rel","hat","os",] + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>=60,NA,0), na.last=NA)
f <- sum(prdData$transplantRel & clinicalData$AOD[ !is.na(clinicalData$Recurrence_date)][prdData$index] < 60)/sum(relData$status & clinicalData$AOD[relData$index] < 60 ) # fraction of patients that have received a salvage transplant
s <- sapply(seq_along(o), function(i) mean(c(multiRFX3TplCi["cr1","hat","os",o[1:i]], (1-f)*multiRFX3TplCi["none","hat","os",o[-(1:i)]] + f*multiRFX3TplCi["rel","hat","os",o[-(1:i)]]), na.rm=TRUE))
plot(seq_along(s)/length(s), s, type='l', xlab="Fraction of allografts in CR1", ylab="Survival of eligible patients 3yrs after CR", col=set1[1], lty=3)
s <- rowMeans(sapply(1:10, function(foo){ set.seed(foo)
					o <- order(-multiRFX3TplCi["dCr1Rel","hat","os",] + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>=60,NA,0) + rnorm(1540,sd=(multiRFX3TplCi["dCr1Rel","upper","os",]-multiRFX3TplCi["dCr1Rel","lower","os",])/4), na.last=NA)
					s <- sapply(seq_along(o), function(i) mean(c(multiRFX3TplCi["cr1","hat","os",o[1:i]], (1-f)*multiRFX3TplCi["none","hat","os",o[-(1:i)]] + f*multiRFX3TplCi["rel","hat","os",o[-(1:i)]]), na.rm=TRUE))
				}))
lines(seq_along(s)/length(s), s, type='l',col=set1[1], lty=1)
p <- order(na.zero(c(1,4,2,3)[clinicalData$M_Risk])  + dataFrame$AOD_10/20 + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>=60,NA,0), na.last=NA)
e <- sapply(seq_along(p), function(i) mean(c(multiRFX3Tpl[p[1:i],2], f*multiRFX3Tpl[p[-(1:i)],3] + (1-f)*multiRFX3Tpl[p[-(1:i)],1]), na.rm=TRUE)) 
lines(seq_along(e)/length(e), e, type='l', col=set1[2])
legend("bottomright", c("Personalised risk", "ELN and age"), lty=1, col=set1, bty="n")

#' Consider fraction of patients relapsing with salvage allograft as well
r <- 1+multiRFX3TplCi[1:2,1,"aar",] - multiRFX3TplCi[1:2,1,"rs",] ## Relapse probabilities
a <- sapply(seq_along(o), function(i) mean(c(r[2,o[1:i]], r[1,o[-(1:i)]]), na.rm=TRUE)) # Personalised
b <- sapply(seq_along(p), function(i) mean(c(r[2,p[1:i]], r[1,p[-(1:i)]]), na.rm=TRUE)) # ELN

x <- seq_along(s)/length(s)
plot(x + (1-x)*a*f,s, type='l', xlab="Total fraction of allografts", ylab="Survival of eligible patients 3yrs after CR", col=set1[1])
lines(x + (1-x)*b*f,e, lty=1, col=set1[2])
legend("bottomright", c("Personalised risk", "ELN and age"), lty=1, col=set1, bty="n")

#' Fraction of relapse
plot(x ,a*f, type='l', xlab="Fraction with allograft in CR1", ylab="Fraction relapsing", col=set1[1])
lines(x, b*f, lty=1, col=set1[2])
legend("bottomright", c("Personalised risk", "ELN and age"), lty=1, col=set1, bty="n")

#' #### Leave one out cross-validation 
#' ##### Three state model
#' We compute LOO out-of-sample predictions for the survival gain by allograft in CR1 v relapse by training 1540 models on 1539 patients each. 
#' This we compare to in-sample predictions of the model trained on all 1540 patients.
#+ allPredictLOO, cache=TRUE
cvFold <- nrow(dataFrame)
cvIdx <- 1:nrow(dataFrame)
multiRFX3LOO <- Reduce("rbind", mclapply(1:nrow(data), function(i){
					whichTrain <- which(cvIdx != i)
					rfxNrm <- CoxRFX(nrdData[nrdData$index %in% whichTrain, names(crGroups)], Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
					rfxNrm$coefficients["transplantRel"] <- 0
					rfxPrs <-  CoxRFX(prdData[prdData$index %in% whichTrain, names(crGroups)], Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% whichTrain], groups=crGroups, nu=1, which.mu = intersect(mainGroups, unique(crGroups)))
					rfxCir <-  CoxRFX(relData[relData$index %in% whichTrain, names(crGroups)], Surv(relData$time1, relData$time2, relData$status)[relData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
					rfxCir$coefficients["transplantRel"] <- 0
					allPrd <- MultiRFX3(rfxNrm, rfxCir, rfxPrs, data=allDataTpl[rep(cvIdx, each=3) == i,], x=3*365, prdData=prdData[prdData$index %in% whichTrain,])
					allPrd <- matrix(allPrd$os, ncol=3, byrow=TRUE, dimnames=list(NULL, c("None","CR1","Relapse")))
				}, mc.cores=10))

plot(multiRFX3LOO[,2]-multiRFX3LOO[,3],multiRFX3TplCi["dCr1Rel","hat","os",] )
cor(multiRFX3LOO[,2]-multiRFX3LOO[,3],multiRFX3TplCi["dCr1Rel","hat","os",] )

#' The correlation of in-sample and out-of-sample predictions is very high, but there are some differences.
#' We now assess the accuracy of our predictions by comparing the observed survival with the out-of-sample prediction. To this end,
#' we split out the quarter of patients predicted to benefit the most. In both subsets we compare the observed 3yr survial between patients with
#' and without allograft in CR1 and compute the difference. CIs by boostrapping.
#+ allPredictLOOkM
d <- multiRFX3LOO[,2]-multiRFX3LOO[,3]
w <- which(clinicalData$AOD < 60)
q <- c(min(d), 0.1, max(d)) 
c <- cut(d, breaks=q, include.lowest=TRUE)# , paste0("[",names(q)[-length(q)],",",names(q)[-1],")"))
e <- sapply(levels(c), 
		function(cc) {
			t <- try(survfit(Surv(time1, time2, status) ~ transplantCR1, data=osData, subset=c[osData$index]==cc & osData$index %in% w)); 
			if(class(t)=="try-error") 
				rep(NA,2)
			else {
				s <- summary(t, time=3*365)
				if(length(s$surv)==2) {
					ci <- sapply(1:200, function(foo){
								set.seed(foo)
								b <- sample(1:nrow(dataFrame), replace=TRUE)
								diff(summary(survfit(Surv(time1, time2, status) ~ transplantCR1, data=osData, subset=c[osData$index]==cc & osData$index %in% w & osData$index %in% b), time=3*365)$surv)
							})
					r <- c(diff(s$surv), quantile(ci, c(0.025, 0.975)))
					names(r) <- c("delta", "lower",'upper')
					r
				}
				else rep(NA,3)
			}})
x <- sapply(split(d[w],c[w]),median)
par(xpd=NA, bty="L")
plot(x,e[1,], pch=19, xlim=c(-.05,.2), ylim=c(-.05,.2), xlab = "Predicted survival benefit", ylab="Observed survival benefit (leave-one-out CV)")
h <- density(d[w]) 
y <- h$y/diff(range(h$y))*.05 + par("usr")[3]
v <- h$x <= q[2] 
par(xpd=FALSE)
polygon(c(h$x[v], h$x[which(v)[length(which(v))]]), c(y[v],par("usr")[3]), border=NA, col=set1[1])
polygon(c(h$x[which(!v)[1]], h$x[!v]), c(par("usr")[3],y[!v]), border=NA, col=set1[3])
lines(h$x, y)
segments(x,e[2,],x,e[3,])
#rug(d, col="#00000022")
abline(0,1, lty=3)



#' #### Figure 4b
#' Violins plot of the predicted survival gain
#+ benefit_hsct, fig.width=1, fig.height=2.5
par(mar=c(3,3,1,1), mgp=c(2,0.5,0), bty="n")
h <- density(d[w]) 
y <- h$y/diff(range(h$y))*.05 + par("usr")[3]
par(xpd=FALSE)
xx <- c(h$x, rev(h$x))
yy <- c(h$y, -rev(h$y))
v <- xx <= q[2] 
plot(yy,xx, pch=NA, ylab="Predicted benefit", xlab="", xaxt="n")
polygon(yy[v], xx[v], border=NA, col=set1[1])
polygon(yy[!v], xx[!v], border=NA, col=set1[2])
lines(yy, xx)

#' #### Figure 4c
#' KM plot of the high v low benefit groups
#+ survival_hsct, fig.width=3, fig.height=2.5
par(mar=c(3,3,1,1), mgp=c(2,0.5,0), bty="L")
f <- survfit(Surv(time1/365, time2/365, status) ~ group +  transplantCR1, data=cbind(osData, group=c[osData$index]), subset=osData$index %in% w)
summary(f, time=3)
plot(f, col=rep(set1[1:nlevels(c)],each=2), lty=rep(c(1,2), nlevels(c)), xlab="TIme after CR", ylab="Survival", xlim=c(0,5), cex=.5)
t <- table(w %in% osData$index[osData$transplantCR1==1],c[w],!is.na(clinicalData$CR_date[w]))[,,"TRUE"]
legend("topright", legend=as.numeric(t), col=rep(set1[1:nlevels(c)],each=2), lty=rep(c(1,2), nlevels(c)), bty="n")



#' The model predictions appear consistent, but there is still a substantial uncertainty. We are unable to further tease apart the accuracy in the lower quartiles, as
#' the predicted effect size is too small. The following is a barplot of the predicted and observed survival benefit.
#+ allPredictLOOkMbar, fig.width=1.5
barplot(rbind(x,e[1,]), beside=TRUE, names.arg=colnames(e), args.legend=list(x="topleft",bty="n",legend=c("predicted","observed")), legend=TRUE, xlab = "Predicted survival benefit", ylab="Three-year survival benefit", ylim=c(-0.05,.2)) -> b
segments(b[2,],e[2,],b[2,],e[3,])

#' The bottom line is that we are able to confidently isolate a quarter of patients with high benefit of allografts (about 12% absolute benefit). The breakdown across 
#' ELN risk groups is:
table(c, clinicalData$M_Risk)

#' ##### Leave one out cross-validation for RFX on post-CR OS 
#+ coxRFXOsCrLOO, cache=TRUE
cvFold <- nrow(dataFrame)
cvIdx <- 1:nrow(dataFrame)
p <- Reduce("rbind", mclapply(cvIdx, function(i){
					whichTrain <- which(cvIdx != i)
					rfxOS <- CoxRFX(osData[osData$index %in% whichTrain, names(crGroups)], Surv(osData$time1, osData$time2, osData$status)[osData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
					p <- as.data.frame(predict(rfxOS, newdata=osData[!osData$index %in% whichTrain, names(crGroups)], se.fit=TRUE))
					s <- summary(survfit(Surv(osData$time1, osData$time2, osData$status)[osData$index %in% whichTrain] ~ 1), time=3*365)$surv
					cbind(p, surv=s^exp(p$fit))
				}, mc.cores=10))
d <- duplicated(sub(".1$","",rownames(p)))
coxRFXOsCrLOO <- rbind(p[!d,], p[d,])
rm(p,d)

#' Compare with corresponding multistage predictions
#+ coxRFXOsCrLOOplot
m <- c(multiRFX3LOO[,3],multiRFX3LOO[osData$index[osData$transplantCR1==1],2])
r <- c(coxRFXOsCrLOO$surv[1:1540],coxRFXOsCrLOO$surv[osData$transplantCR1==1]) 
plot(m, r)
abline(0,1)
cor(m, r)

#' #### Prediction errors
#' ##### Training error
#' 3-state model
c <- Surv(as.numeric(clinicalData$Date_LF - clinicalData$CR_date), clinicalData$Status)
p <- multiRFX3TplCi[3,1,1,]
p[osData$index[osData$transplant1CR==1]] <-  multiRFX3TplCi[2,1,1,]
ape(p, c, time=3*365)

#' RFX
unduplicate <- function(index) {u <- unique(index); u[which(rev(duplicated(rev(index))))] <- seq_along(index)[duplicated(index)]; return(u)}
q <- summary(survfit(Surv(time1,time2,status) ~ 1, data=osData), time=3*365)$surv^exp(scale(predict(coxRFXOsCR, newdata=osData[unduplicate(osData$index),]), scale=FALSE))
ape(q, c, time=3*365)

#' ##### LOO test error
#' 3-state model
p <- multiRFX3LOO[,3]
p[osData$index[osData$transplantCR1==1]] <-  multiRFX3LOO[osData$index[osData$transplantCR1==1],2]
ape(p, c, time=3*365)

#' RFX
ape(coxRFXOsCrLOO$surv[unduplicate(osData$index)], c, time=3*365)


#' ### Predicting outcome from diagnosis
#' The following function fits a 5-stage model. Note that we use a single smooth function g(t) to model the association between time of CR and all subsequent events.
MultiRFX5 <- function(coxRFXNcdTD, coxRFXCrTD, coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, x =365, tdPrmBaseline = rep(1, ceiling(max(x))+1), tdOsBaseline = rep(1, ceiling(max(x))+1), ciType="analytical"){
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
	
	
	
	## Step 1: Compute KM survival curves and log hazard
	getS <- function(coxRFX, data, max.x=5000) {		
		if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
		data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data)), drop=FALSE])
		r <- PredictRiskMissing(coxRFX, data, var="var2")
		H0 <- basehaz(coxRFX, centered = FALSE)
		hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
		x <- c(0:ceiling(max.x))
		S <- exp(-hazardDist(x))
		return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
	}
	kmCr <- getS(coxRFX = coxRFXCrTD, data = data, max.x=max(x))
	kmEs <- getS(coxRFX = coxRFXNcdTD, data = data, max.x=max(x))
	kmCir <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
	kmNrm <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
	kmPrs <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
	
	xx <- 0:ceiling(max(x))
	
	sapply(1:nrow(data), function(i){
				## Step 2: Adjust curves for competing risks, accounting for hazard
				crAbs <-  cumsum(c(1,diff(kmCr$S^exp(kmCr$r[i,1]))) * kmEs$S ^ exp(kmEs$r[i,1]))
				esAbs  <- cumsum(c(1,diff(kmEs$S^exp(kmEs$r[i,1]))) * kmCr$S ^ exp(kmCr$r[i,1])) ## array times x nrow(data)
				cirCrAbs <- cumsum(c(1,diff(kmCir$S^exp(kmCir$r[i,1]))) * kmNrm$S ^ exp(kmNrm$r[i,1]))
				nrsCrAbs <- cumsum(c(1,diff(kmNrm$S^exp(kmNrm$r[i,1]))) * kmCir$S ^ exp(kmCir$r[i,1])) ## array times x nrow(data)
				
				## Step 3: Compute hierarchical survival
				### Prs			
				rsCrAbs <- computeHierarchicalSurvival(x = xx, diffS0 = diff(cirCrAbs), S1Static = kmPrs$S, haz1TimeDep = tdPrmBaseline * exp(kmPrs$r[i,1]))
				
				## Confidence intervals (loglog)
				PlogP2 <- function(x) {(x * log(x))^2}
				errOs <- kmNrm$r[i,2] * PlogP2(kmNrm$S^exp(kmNrm$r[i,1])) * (1-(1-kmCir$S ^ exp(kmCir$r[i,1]))) * (1-kmPrs$S ^ exp(kmPrs$r[i,1]))^2 + kmCir$r[i,2] * PlogP2(kmCir$S ^ exp(kmCir$r[i,1])) * (1-kmPrs$S ^ exp(kmPrs$r[i,1]))^2 * (kmNrm$S ^ exp(kmNrm$r[i,1]))^2 +  kmPrs$r[i,2] * PlogP2(kmPrs$S ^ exp(kmPrs$r[i,1])) * (1-kmCir$S ^ exp(kmCir$r[i,1]))^2 * (kmNrm$S ^ exp(kmNrm$r[i,1]))^2 
				sdOsCr <- sqrt(errOs / PlogP2(1-(1-nrsCrAbs)-(1-rsCrAbs)))
				
				
				### Overall survival from enrollment
				nrsEr <- computeHierarchicalSurvival(x = xx, diffS0 = diff(crAbs), S1Static = nrsCrAbs, haz1TimeDep = tdOsBaseline)
				rsEr <- computeHierarchicalSurvival(x = xx, diffS0 = diff(crAbs), S1Static = rsCrAbs, haz1TimeDep = tdOsBaseline)
				cirEr <- computeHierarchicalSurvival(x = xx, diffS0 = diff(crAbs), S1Static = cirCrAbs, haz1TimeDep = tdOsBaseline)
				cbind(deathInErFromEr=1-esAbs, deathInCrFromEr=1-nrsEr, deathInRelFromEr=1-rsEr, aliveInRelFromEr=1-cirEr-(1-rsEr), aliveInCrFromEr=1-crAbs - (1-cirEr) - (1-nrsEr),
						deathInCrFromCr = 1-nrsCrAbs, deathInRelapseFromCr=(1-rsCrAbs), aliveInRelapseFromCr = (1-cirCrAbs) - (1-rsCrAbs), osInCrFromCrSd=sdOsCr
				)
			}, simplify='array')
}

#' PRS baseline with spline-based dep on CR length)
#+ fiveStagePredicted, cache=TRUE
xmax <- 2000
xx <- 0:ceiling(xmax)
coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])) 
tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1]))) ## Hazard (function of CR length)	

coxphOs <- coxph(Surv(time1,time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))) 
tdOsBaseline <- exp(predict(coxphOs, newdata=data.frame(time0=xx[-1])))	 ## Hazard (function of induction length), only for OS (could do CIR,NRM,PRS seperately)

fiveStagePredicted <- MultiRFX5(coxRFXNcdTD, coxRFXCrTD, coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, tdPrmBaseline = tdPrmBaseline, tdOsBaseline = tdOsBaseline, x=xmax)

#' Function to plot stages
sedimentPlot <- function(Y, x=1:nrow(Y), y0=0, y1=NULL, col=1:ncol(Y), ...){
	Z <- cbind(t(apply(cbind(y0,Y),1,cumsum)),y1)
	plot(x,Z[,1], xlim=range(x), ylim=range(Z), lty=0, pch=NA,...)
	for(i in 2:ncol(Z))
		polygon(c(x,rev(x)), c(Z[,i-1],rev(Z[,i])), border=NA, col=col[i-1])
}

lineStage <- function(CR_date, Recurrence_date, Date_LF, ERDate, Status, y=0, col=1:5, pch.trans=19, pch.end=19, ...){
	xpd <- par("xpd")
	par(xpd=NA) 
	t <- as.numeric(c(CR_date, Recurrence_date, Date_LF) - ERDate )
	w <- !is.na(t)
	o <- order(t)
	to <- pmin(t[o], par("usr")[2])
	l <- length(to)
	segments(c(0,to[-l]), rep(y,l), to, rep(y,l), col=col, lend=1, ...)
	status <- if(Status == 1) 3 else 0
	if(is.na(Recurrence_date))
		status <- status - 1
	if(is.na(CR_date))
		status <- status - 1
	x <- ifelse(t <= par("usr")[2], t, NA)
	points(x, rep(y, length(t)), pch=c(pch.trans,pch.trans, if(Status) pch.end else NA), col=col[c(2:3,status+3)])
	par(xpd=xpd)
}


#' Average of all multistage predictions, note the precise agreement with overall survival.
#+ fiveStagePredictedAvg, fig.width=3, fig.height=2.5
pastel1 <- brewer.pal(9, "Pastel1")
par(mfrow=c(1,1), mar=c(3,3,1,1), cex=1)
sedimentPlot(-rowMeans(fiveStagePredicted[,1:5,], dims=2), y0=1, y1=0,  col=c(pastel1[c(1:3,5,4)], "#DDDDDD"))
lines(survfit(Surv(OS, Status) ~ 1, data=clinicalData))


#' Multistage predictions v overall survival
#+ HRFXvRFX
for(i in 1:5)
plot(summary(survfit(coxRFXFitOsTDGGc), i*365)$surv^ exp(coxRFXFitOsTDGGc$linear.predictors[1:1540]), 1-rowSums(aperm(fiveStagePredicted[,1:3,], c(3,1,2)), dim=2)[,365*3],
		xlab="Survival RFX OS", ylab="Survival RFX Multistage", main=paste(i, "years"))


#' #### Leave-one-out cross-validation
#' The following code is run on the cluster
read_chunk('../../code/leaveOneOut.R', labels="leaveOneOut")
#+ leaveOneOut, eval=FALSE

#' Multistage model
#+ multiRfx5Loo, cache=TRUE
times <- round(seq(0,5,0.05)*365)
multiRfx5Loo <- sapply(mclapply(1:nrow(data), function(i){
					e <- new.env()
					t <- try(load(paste0("../../code/loo/",i,".RData"), env=e))
					if(class(t)=="try-error") rep(NA, length(times))
					else e$multiRfx5[times+1,,1]
				}, mc.cores=6), I, simplify="array")

#' Error OS
survConcordance(os ~ colSums(multiRfx5Loo[times == 3*365,1:3,]))
ape(1-colSums(multiRfx5Loo[times == 3*365,1:3,]), os, 3*365)


#' #### Figure 2
#' In order of risk constellation plots
#+ fiveStagePredictedHilbert, fig.width=12, fig.height=12
set.seed(42)
s <- sample(nrow(dataFrame),nStars^2) #1:(nStars^2)
library(HilbertVis)
nStars <- 32
l <- "coxRFXFitOsTDGGc"
t <- os#get(l)$surv
p <- PartialRisk(get(l),  newZ=dataFrame[, colnames(get(l)$Z)])
p <- p[,colnames(p)!="Nuisance"]
locations <- hilbertCurve(log2(nStars))+1 
mat <- matrix(order(locations[,1], locations[,2]), ncol=nStars)
h <- hclust(dist(p[s,]))
layout(mat[nStars:1,])
par(mar=c(0,0,0,0),+.5, bty="n")
for(i in 1:nStars^2){ # Fitted predictions
	sedimentPlot(-fiveStagePredicted[seq(1,2001,200),1:5,s[h$order[i]]], x=seq(1,2001,200),y0=1, y1=0,  col=c(pastel1[c(1:3,5,4)], "#DDDDDD"), xlab="time",ylab="fraction", xaxt="n", yaxt="n")
	lines(x=seq(1,2001,200), y=1-rowSums(fiveStagePredicted[seq(1,2001,200),1:3,s[h$order[i]]]), lwd=2)
	i <- s[h$order[i]]
	lineStage(clinicalData$CR_date[i], clinicalData$Recurrence_date[i], clinicalData$Date_LF[i], clinicalData$ERDate[i], clinicalData$Status[i], col=c(brewer.pal(8,"Dark2")[8], set1[c(4:5,1:3)]), lwd=2, pch.trans=NA, y=0.05)	
}
for(i in 1:nStars^2){ # Leave-one-out predictions
	sedimentPlot(-multiRfx5Loo[seq(1,length(times),5),1:5,s[h$order[i]]], x=times[seq(1,length(times),5)],y0=1, y1=0,  col=c(pastel1[c(1:3,5,4)], "#DDDDDD"), xlab="time",ylab="fraction", xaxt="n", yaxt="n")
	lines(x=times[seq(1,length(times),5)], y=1-rowSums(multiRfx5Loo[seq(1,length(times),5),1:3,s[h$order[i]]]), lwd=2)
	i <- s[h$order[i]]
	lineStage(clinicalData$CR_date[i], clinicalData$Recurrence_date[i], clinicalData$Date_LF[i], clinicalData$ERDate[i], clinicalData$Status[i], col=c(brewer.pal(8,"Dark2")[8], set1[c(4:5,1:3)]), lwd=2, pch.trans=NA, y=0.05)	
}

#' #### Comparison with RFX
#+ rfx5Loo, cache=TRUE
rfx5Loo <- sapply(mclapply(1:nrow(data), function(i){
					e <- new.env()
					t <- try(load(paste0("../../code/loo/",i,".RData"), env=e))
					if(class(t)=="try-error") rep(NA, length(times))
					else {
						cvIdx <- 1:nrow(dataFrame)
						whichTrain <<- which(cvIdx != i)
						pNrs <- predict(e$rfxNrs, newdata=data[cvIdx==i,])
						pRel <- predict(e$rfxRel, newdata=data[cvIdx==i,])
						pPrs <- predict(e$rfxPrs, newdata=data[cvIdx==i,])
						pCr <- predict(e$rfxCr, newdata=data[cvIdx==i,])
						pEs <- predict(e$rfxEs, newdata=data[cvIdx==i,])
						pOs <- predict(e$rfxOs, newdata=dataFrame[cvIdx==i,])
						c(pCr, pEs, pNrs, pRel, pPrs, pOs)
					}
				}, mc.cores=6), I, simplify="array")

colnames(rfx5Loo) <- rownames(data)
survConcordance(Surv(nrdData$time1, nrdData$time2, nrdData$status) ~ rfx5Loo[3,nrdData$index])
survConcordance(Surv(prdData$time1, prdData$time2, prdData$status) ~ rfx5Loo[5,rownames(prdData)[prdData$index]])
survConcordance(Surv(relData$time1, relData$time2, relData$status) ~ rfx5Loo[4,relData$index])
survConcordance(Surv(cr[,1], cr[,2]==2) ~ rfx5Loo[1,])
survConcordance(Surv(cr[,1], cr[,2]==1) ~ rfx5Loo[2,])
survConcordance(os ~ rfx5Loo[6,])

#' #### Figure 1d
#' Plot of absolute risk at 3yr v outcome
#+ survival_risk, fig.width=3, fig.height=1.5
par(mar=c(3,3,2,1), mgp=c(1.5,.5,0), bty="n")
t <- os
q <- quantile(t[,1], seq(0,1,.1))# q <- splinefun( s$surv, s$time,"monoH.FC")(seq(1,min(s$surv),l=10))
c <- cut(t[,1], q, na.rm=TRUE)
h <- colSums(multiRfx5Loo[times == 3*365,1:3,])
o <- order(h)
plot(h[o], col= (brewer.pal(10,'RdBu'))[c[o]], type='h', xaxt="n", xlab='', las=2, ylab="Survival at 3 years")
mtext(side=1, line=1, "Patient")
u <- par("usr")
q <- pmin(q,365*12)
image(x=q/max(q)*500, y=c(u[4]-(u[4]-u[3])/20, u[4]), matrix(1:10), col= (brewer.pal(10,'RdBu')), add=TRUE)
#axis(side=3, at=seq(1,500,l=11), labels=seq(0,1,.1))
axis(side=3, at=pretty(q/365)/max(q)*365*500, labels=pretty(q/365))
lines(ksmooth(seq_along(o),t[o,2]==0, bandwidth=50))

#' #### Supplementary Figure 3
#' Plots of concordance and absolute prediction measures
#+ errorsMultiRfxOsLoo, fig.width=2.5, fig.height=2.5
multiRfx5C <- sapply(seq_along(times), function(i) survConcordance(os ~ colSums(multiRfx5Loo[i,1:3,]))$concordance[1])
plot(times, multiRfx5C, type='l', xlab="Time", ylab="Concordance", ylim=c(0.65, 0.73), col=set1[1])
abline(h=survConcordance(os ~ rfx5Loo[6,])$concordance, col=set1[2], lwd=2)
legend("bottomright",c("RFX OS","RFX Multistage"), col=set1[1:2], lty=1, bty="n")

a <- sapply(times, function(t) ape(1-colSums(multiRfx5Loo[times == t,1:3,]), os, t))
s <- summary(survfit(coxRFXFitOsTDGGc), times=times)
b <- sapply(times, function(t) ape(s$surv[times==t]^exp(rfx5Loo[6,]), os, t))
e <- sapply(times, function(t) ape(s[times==t], os, t))
for(i in 1:4){
	plot(times/365.25, e[i,], type='l', xlab="Time (yr)", ylab=rownames(a)[i], col=set1[9])
	lines(times/365.25, a[i,], col=set1[1])
	lines(times/365.25, b[i,], col=set1[2])
	legend("bottomright",c("Kaplan-Meier","RFX OS","RFX Multistage"), col=set1[c(9,1:2)], lty=1, bty="n")
}

#' Figure of predicted survival for 100 patients, comparing multistage and OS predictions
#+ survivalMultiRfxOsLoo, fig.width=2.5, fig.height=2.5
plot(s$surv^exp(rfx5Loo[6,1]), 1-rowSums(multiRfx5Loo[,1:3,1]), type='l', xlim=c(0,1), ylim=c(0,1), col='grey', xlab="Predicted survival RFX", ylab="Pedicted survival Multistage")
for(i in 2:100)
	lines(s$surv^exp(rfx5Loo[6,i]), 1-rowSums(multiRfx5Loo[,1:3,i]), col='grey')

#' #### Figure 3a-f, Supplementary Figure 4
#' With and without TPL
#+ threePatientsAllo, fig.width=3, fig.height=2.5
xmax=2000
patients <- c("PD11104a","PD8314a","PD9227a")
fiveStagePredictedTplLoo <- sapply(patients, function(pd){
			e <- new.env()
			i <- which(rownames(dataFrame)==pd)
			load(paste0("../../code/loo/",i,".RData"), env=e)

			cvIdx <- 1:nrow(dataFrame)
			whichTrain <<- which(cvIdx != i)
			xx <- 0:2000
			coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])[prdData$index %in% whichTrain,]) 
			tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
			
			coxphOs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))[osData$index %in% whichTrain,]) 
			tdOsBaseline <- exp(pmin(predict(coxphOs, newdata=data.frame(time0=500)),predict(coxphOs, newdata=data.frame(time0=xx[-1])))) ## cap predictions at induction length 500 days.
			m <- MultiRFX5(e$rfxEs, e$rfxCr, e$rfxNrs, e$rfxRel, e$rfxPrs, allDataTpl[grep(pd, rownames(allDataTpl)),], tdPrmBaseline = tdPrmBaseline, tdOsBaseline = tdOsBaseline, x=2000)			
		}, simplify="array")
#par(mfrow=c(2,2))
par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
w <- seq(1,2001,10)
at <- ceiling(1:5 * 365.5)
x <- (w-1)/365.25
for(pd in patients)
	for(i in c(2,3)){
	sedimentPlot(-fiveStagePredictedTplLoo[w,6:8,i,pd],x=x, y0=1, y1=0,  col=pastel1[c(2:3,5,4)], xlab="Years from CR",ylab="Probability", xaxs='i', yaxs='i')
	o <- 1-rowSums(fiveStagePredictedTplLoo[w,6:7,i,pd])
	abline(v=c(1:5), col="white", lty=3)
	abline(h=seq(0.2,0.8,0.2), col="white", lty=3)
	lines(x,o, lwd=2)
	lines(x,o ^ exp(qnorm(0.975) * fiveStagePredictedTplLoo[w,9,i,pd]))
	lines(x,o ^ exp(-qnorm(0.975) * fiveStagePredictedTplLoo[w,9,i,pd]))
	text(x=rep(0,3), c(0.1,0.2,0.3), c("rel./al.", "rel./death", "n.r./death") )
	text(x=1:5, y=rep(0.3, 5), round(fiveStagePredictedTplLoo[at,6,i,pd],2))
	text(x=1:5, y=rep(0.2, 5), round(fiveStagePredictedTplLoo[at,7,i,pd],2))
	text(x=1:5, y=rep(0.1, 5), round(fiveStagePredictedTplLoo[at,8,i,pd],2))
	#text(x=at, y=rep(0.1, 5), round(fiveStagePredictedTpl[w,6,i],2))
}

#' #### Three patients with numerical CI's and LOO
patients <- c("PD11104a","PD8314a","PD9227a")
threePatientTplCiLoo <- sapply(patients, function(pd){
			e <- new.env()
			i <- which(rownames(dataFrame)==pd)
			set.seed(i)
			load(paste0("../../code/loo/",i,".RData"), env=e)			
			multiRFX3TplCi <- MultiRFX3TplCi(e$rfxNrs, e$rfxRel, e$rfxPrs, data=data[i,colnames(coxRFXNrdTD$Z)], x=3*365, nSim=200, prsData=prsData[prsData$index!=i,]) ## others with 200
		}, simplify="array")


#' #### Figure 4a
#' The figure shows the mortality reduction of allograft CR1 v none, allograft in Rel v none, and CR1 v Relapse, for LOO predictions similar to above.
#+mortalityReductionLoo, fig.width=3.5, fig.height=2.5
par(mar=c(3,3,1,3), las=2, mgp=c(2,.5,0), bty="n")
benefit <- multiRFX3LOO[,2]-multiRFX3LOO[,3]
absrisk <- multiRFX3LOO[,1]
s <- clinicalData$AOD < 60 & !is.na(clinicalData$CR_date)#sample(1:1540,100)
x <- 1-absrisk
y <- benefit
plot(x[s], y[s], pch=NA, ylab="Mortality reduction from allograft", xlab="3yr mortality with standard chemo", col=riskCol[clinicalData$M_Risk], cex=.8, las=1, ylim=range(benefit))
abline(h=seq(-.1,.3,.1), col='grey', lty=3)
abline(v=seq(.2,.9,0.2), col='grey', lty=3)
points(x[s], y[s], pch=16,  col=riskCol[clinicalData$M_Risk[s]], cex=.8)
segments(1-threePatientTplCiLoo["none","lower","os",1,patients], threePatientTplCiLoo["dCr1Rel","hat","os",1,patients],1-threePatientTplCiLoo["none","upper","os",1,patients],threePatientTplCiLoo["dCr1Rel","hat","os",1,patients])
segments(1-threePatientTplCiLoo["none","hat","os",1,patients], threePatientTplCiLoo["dCr1Rel","lower","os",1,patients],1-threePatientTplCiLoo["none","hat","os",1,patients], threePatientTplCiLoo["dCr1Rel","upper","os",1,patients])
xn <- seq(0,1,0.01)
p <- predict(loess(y~x, data=data.frame(x=x[s], y=y[s])), newdata=data.frame(x=xn), se=TRUE)
yn <- c(p$fit + 2*p$se.fit,rev(p$fit - 2*p$se.fit))
polygon(c(xn, rev(xn))[!is.na(yn)],yn[!is.na(yn)], border=NA, col="#00000044", lwd=1)
lines(xn,p$fit, col='black', lwd=2)
legend("topleft", pch=c(16,16,16,16,NA),lty=c(NA,NA,NA,NA,1), col=c(riskCol[c(2,4,3,1)],1),fill=c(NA,NA,NA,NA,"grey"), border=NA, c(names(riskCol)[c(2,4,3,1)],"loess average"), box.lty=0)
n <- c(100,50,20,10,5,4,3)
axis(side=4, at=1/n, labels=n, las=1)
mtext("Number needed to treat", side=4, at=.2, line=2, las=0)
axis(side=4, at=-1/n, labels=n, las=1)
mtext("Number needed to harm", side=4, at=-.1, line=2, las=0)


#' #### Imputation of missing genes
#' For RFX model on OS
#+ imputationGenes, cache=TRUE
w <- WaldTest(coxRFXFitOsTDGGc)
o <- order(w$p.value[groups[whichRFXOsTDGG] %in% c("Genetics","GeneGene")])
genes <- unique(sub("_.+","",unlist(strsplit(names(whichRFXOsTDGG[groups[whichRFXOsTDGG]%in% c("Genetics","GeneGene")])[o],":"))))

cvFold <- 1540
foo <- 42
set.seed(foo)
cvIdx <- 1:cvFold #sample(1:nrow(dataFrame)%% cvFold +1 ) ## sample 1/10

m <- unlist(sapply(1:cvFold, function(i) which(tplSplitOs %in% which(cvIdx==i))))
o <- order(m)

imputedRiskCv <- do.call("abind", c(mclapply(1:cvFold, function(i){
							whichTrain <- which(cvIdx != i)
							ix <- tplSplitOs %in% whichTrain
							cRfx <- CoxRFX(dataFrameOsTD[ix,whichRFXOsTDGG], osTD[ix], groups[whichRFXOsTDGG], which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 
							
							imputedRisk <- sapply(mclapply(c(0,seq_along(genes)), function(i){
												na.genes <- if(i==0) genes else genes[-(1:i)]
												if(length(na.genes)==0) na.genes <- "FOO42"
												d <- dataFrameOsTD[,whichRFXOsTDGG]
												d[grepl(paste(na.genes, collapse="|"), colnames(d))] <- NA
												p <- PredictRiskMissing(cRfx, d[!ix,,drop=FALSE])
											}, mc.cores=1), I, simplify="array")
							dimnames(imputedRisk)[[3]] <- c("None",genes)
							return(imputedRisk)
							
						}, mc.cores=10), along=1))[o,,]



#+ imputationGenesPlot, fig.width=5, fig.height=2.5
par(mar=c(3,3,3,1))
imputedCCv <- sapply(dimnames(imputedRiskCv)[[3]], function(i) as.numeric(survConcordance(osTD ~ imputedRiskCv[,1,i])[c("concordance","std.err")]))
x <- 0:ncol(imputedCCv)-.5
plot(x, c(imputedCCv[1,], imputedCCv[1,ncol(imputedCCv)]), type="s", xaxt="n", xlab="", ylab="Concordance", ylim=range(imputedCCv[1,]) + c(-1,1)*imputedCCv[2,1])
polygon(c(rep(x,each=2)[-c(1, 2*length(x))],rep(rev(x), each=2)[-c(1, 2*length(x))]), c(rep(imputedCCv[1,]+imputedCCv[2,], each=2), rep(rev(imputedCCv[1,]-imputedCCv[2,]), each=2)), border=NA, col="#00000044")
mtext(dimnames(imputedRiskCv)[[3]], side=1, at=1:ncol(imputedCCv)-1, las=2, font=3, cex=.9)
abline(v=seq(0,50,10), lty=3)
abline(h=seq(0.68,0.73,0.01), lty=3)
axis(side=3)


#' Genetic imputation multi stage
read_chunk('../../code/imputation.R', labels="imputationMultiRfx")
#+ imputationMultiRfx, eval=FALSE

#' Collect data
#+ multiRfx5CvImputed, cache=TRUE
multiRfx5CvImputed <- sapply(mclapply(1:nrow(data), function(i){
			e <- new.env()
			t <- try(load(paste0("../../code/imputed/",i,".RData"), env=e))
			if(class(t)=="try-error") return(rep(NA, length(genes)+1))
			else colSums(e$multiRfx5Imputed[3*365,1:3,])
		}, mc.cores=10), I)

#' #### Figure S8
#' Imputed accuracy
#+ multiRfx5CvImputedPlot, cache=TRUE, fig.width=5, fig.height=2.5
par(mar=c(3,3,3,1))
multiRfx5CvImputedC <- sapply(1:nrow(multiRfx5CvImputed), function(i) as.numeric(survConcordance(os ~ multiRfx5CvImputed[i,])[c('concordance','std.err')]))					
x <- 0:ncol(multiRfx5CvImputedC)-.5
plot(x, c(multiRfx5CvImputedC[1,], multiRfx5CvImputedC[1,ncol(multiRfx5CvImputedC)]), type="s", xaxt="n", xlab="", ylab="Concordance", ylim=range(multiRfx5CvImputedC[1,]) + c(-1,1)*multiRfx5CvImputedC[2,1])
polygon(c(rep(x,each=2)[-c(1, 2*length(x))],rep(rev(x), each=2)[-c(1, 2*length(x))]), c(rep(multiRfx5CvImputedC[1,]+multiRfx5CvImputedC[2,], each=2), rep(rev(multiRfx5CvImputedC[1,]-multiRfx5CvImputedC[2,]), each=2)), border=NA, col="#00000044")
mtext(dimnames(imputedRiskCv)[[3]], side=1, at=1:ncol(multiRfx5CvImputedC)-1, las=2, font=3, cex=.9)
abline(v=seq(0,50,10), lty=3)
abline(h=seq(0.68,0.73,0.01), lty=3)
axis(side=3)

par(mar=c(3,3,3,1))
multiRfx5CvImputedApe <- sapply(1:nrow(multiRfx5CvImputed), function(i) ape(1-multiRfx5CvImputed[i,], os, 3*365))					
x <- 0:ncol(multiRfx5CvImputedApe)-.5
for(i in 1:4){
	plot(x, c(multiRfx5CvImputedApe[i,], multiRfx5CvImputedApe[i,ncol(multiRfx5CvImputedApe)]), type="s", xaxt="n", xlab="", ylab=rownames(multiRfx5CvImputedApe)[i], col=set1[i])
	mtext(dimnames(imputedRiskCv)[[3]], side=1, at=1:ncol(multiRfx5CvImputedApe)-1, las=2, font=3, cex=.9)
	abline(v=seq(0,50,10), lty=3)
	abline(h=axTicks(side=2), lty=3)
	axis(side=3)
}


#' # Model comparison
#' 
#' We used cross-validation to evaluate the performance of different modelling strategies. The idea is to split the data
#' into a training and a test set; the model is fitted on the training part and its prognostic accuracy evaluated on the test set.
#' 
#' 
#' ## Random cross validation
#' We have randomly split the data 100 times 80% training data and 20% validation data. For each split, we 
#' evaluated the following metrics:
#' 
#' * Concordance,  defined as the probability that the survival times of two individuals are concordant to the ranking of their risk [@HarrellSM1996],
#' implemented in `survival::survConcordance()` 
#' * Area under the receiver-operating characteristic curve (AUC), a measure of the correct classification into dead and alive at
#' a given time [@UnoJOTASA2007], implemented in `survAUC::auc.Uno()`.
#' * Brier score, an absolute measure of the prediction error [@GerdsBJ2006], implemented in `survAUC::predErr()`.
#' * A generalised coefficient of explained variation $R^2$  [@NagelkerkeB1991,@OQuigleySM2005] , implemented as `survAUC::OXS()`.
#' 
#'  The latter three algorithms are implemented in the `survAUC` `R` package [@Potapov2012].
#' 
#' ## Inter-trial
#' 
#' The data comprised patients from three different trials - AMLSG07/04 (n=740), AMLHD98A (n=627), and AMLHD98B (n=173).
#' In addition to randomly splitting the data into training and test partitions, we trained the model on all three combinations of 2 trials and
#' evaluated the prognostic accuracy on the third trial. This situation is more challenging as there may 
#' be some systematic differences between the trials, but it can also be expected to more closely mimic the situation of predicting a novel cohort.
#' 
#' 
#' ## TCGA
#' 
#' As an additional and independent evaluation cohort, we downloaded data from the cancer genome atlas (TCGA) [@TCGANEJM2013]. 
#' We downloaded variant calls from exome sequencing and cytogenetic data for n=200 and annotated oncogenic mutations as described in our companion paper.
#' Overall survival was available for n=186 patients. For missing prognostic variables, we use a [covariance-based imputation](#covariance-based-imputation), with
#' a covariance matrix derived from our original data set (n=1,540). We note that there was no data available for allografts.
#' 
#' ## Code

#' ### Systematic cross-validation
#' 
#' #### Static models
#+allModelsCV, cache=TRUE, cache.lazy=FALSE
library(rpart)
library(randomForestSRC)
replicates <- 100 ## number of replicates
scope <- c("Genetics","CNA","Treatment","Fusions") ## For CPSS
scopeStep <- as.formula(paste("os ~", paste(colnames(dataFrame)[mainIdxOs& osIdx], collapse="+"))) ## For AIC&BIC
allModelsCV <- mclapply(1:replicates, function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			c <- coxph(os[trainIdx] ~ 1, data=dataFrame[trainIdx,mainIdxOs])
			scopeStep <- as.formula(paste("os[trainIdx] ~", paste(colnames(dataFrame)[mainIdxOs], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxCPSSOsTrain <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx, mainIdxOs], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,mainIdxOs], os[trainIdx], groups=groups[mainIdxOs])
			coxRFXOsTrain$Z <- NULL
			coxRFXOsGGc <- CoxRFX(dataFrame[trainIdx,whichRFXOsGG], os[trainIdx], groups=groups[whichRFXOsGG], which.mu=mainGroups)
			coxRFXOsGGc$Z <- NULL
			rForestOsTrain <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs])[trainIdx,], ntree=100, importance="none")
			return(list(
							BIC=coxBICOsTrain,
							AIC=coxAICOsTrain,
							CPSS=coxCPSSOsTrain,
							RFX=coxRFXOsTrain,
							RFXgg=coxRFXOsGGc,
							rForest=rForestOsTrain
					))
		}, mc.cores=10)

#' Compute predictions for all model fits
#' 
#+ allModelsCvPredictions, cache=TRUE
predictAllModels <- function(x, newdata){
	if("rfsrc" %in% class(x)){
		predict(x, newdata, importance="none")$predicted
	}else{
		predict(x, newdata)
	}
}

allModelsCvPredictions <- mclapply(seq_along(allModelsCV), function(foo){
			set.seed(foo)
			x <- allModelsCV[[foo]]
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[!trainIdx]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrame[!trainIdx,])
							}))
		}, mc.cores=10)

colModels <- c("#888888", set1[c(2,1,4,3,5,7)])

#' Harrel's C
#+ allModelsCv-C
foo <- 1
allModelsCvC <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x, 2 , function(p){						
						survConcordance(osYr[!trainIdx,] ~ p)$concordance
					})
		})
apply(allModelsCvC,1,quantile)

#+ allModelsCvBoxplot, fig.width=2, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(t(apply(-allModelsCvC,2,rank))),factor, levels=1:7)),table)
o <- order(apply(allModelsCvC,1,median))
boxplot(t(allModelsCvC[o,]), notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:7, rep(par("usr")[3],7), rownames(allModelsCvC)[o])

#+ allModelsCvRank, fig.width=2, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
barplot(r[,o]/replicates, col=c(set1[c(3,2,4,1,5,7)],"grey"), ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],6), colnames(allModelsCvC)[o])
legend(par("usr")[1],1.5, fill=c(set1[c(3,2,4,1,5,7)],"grey"), legend=1:6, bty="n", border=NA, horiz=TRUE, title="Rank")

#' #### Supplementary Figure S1B
#' Brier scores
#+ allModelsCv-Brier, cache=TRUE
library(survAUC)
foo <- 1
allModelsCvBrier<- sapply(allModelsCV, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			sapply(x, function(y){
						p <- predictAllModels(y, newdata=dataFrame)
						a <- predErr(Surv.rsp = osYr[trainIdx,], Surv.rsp.new = osYr[!trainIdx,], lp=p[trainIdx], lpnew = p[!trainIdx], times= c(90,365,1000)/365, type="brier")$error
					})
		})
apply(allModelsCvBrier,1,quantile)
rownames(allModelsCvBrier) <- paste(rep(names(allModelsCV[[1]]), each=3), rep(c(90,365,1000), length(allModelsCV[[1]])))
boxplot(t(allModelsCvBrier)[,rep(0:5*3, 3) + rep(1:3, each=6)],notch=TRUE, ylab="Brier score", border=rep(colModels[-1],3), las=2, lty=1, pch=16, staplewex=0)


#' GHCI
#+ allModelsCv-GHCI, cache=TRUE
allModelsCvGHCI<- sapply(allModelsCvPredictions, function(x){
			apply(x[,2:6], 2 , function(p){
						p <- GHCI(lpnew = na.omit(p))
					})
		})
apply(allModelsCvGHCI,1,quantile)
boxplot(t(allModelsCvGHCI),notch=TRUE, ylab="GH C", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)


#' OXS R2 estimates
#+ allModelsCv-OXS, cache=TRUE
foo <- 1
allModelsCvOXS <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x[,2:6], 2 , function(p){						
						a <- OXS(osYr[!trainIdx,], p, rep(0,length(p)))
					})
		})
apply(allModelsCvOXS,1,quantile)
boxplot(t(allModelsCvOXS), ylim=c(0,0.5), notch=TRUE, ylab="OXS R2", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)

#' Nagelk R2 estimates
#+ allModelsCv-Nagelk, cache=TRUE
foo <- 1
allModelsCvNagelk <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x[,2:6], 2 , function(p){						
						a <- Nagelk(osYr[!trainIdx,], p, rep(0,length(p)))
					})
		})
apply(allModelsCvNagelk,1,quantile)
boxplot(t(allModelsCvNagelk), ylim=c(0,0.4), notch=TRUE, ylab="Nagelk's R2", border=colModels[2:6], las=2, lty=1, pch=16, staplewex=0)

#' UnoC
#+ allModelsCv-UnoC, cache=TRUE
foo <- 1
allModelsCvUnoC<- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			apply(x,2, function(p){
						a <- UnoC(Surv.rsp = osYr[trainIdx,], Surv.rsp.new = osYr[!trainIdx,][!is.na(p)],  lpnew = na.omit(p), time=5)
					})
		})
apply(allModelsCvUnoC,1,quantile)
boxplot(t(allModelsCvUnoC), notch=TRUE,  ylab="Uno's C", border=colModels, lty=1, pch=16, staplewex=0)

#' AUC UNO
#+ allModelsCv-AUCuno, cache=TRUE
t <- seq(0.1,5,0.1) #times
allModelsCvAuc <- sapply(seq_along(allModelsCvPredictions), function(foo){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			apply(allModelsCvPredictions[[foo]], 2, function(p){
						AUC.uno(osYr[trainIdx,], osYr[!trainIdx, ][!is.na(p)], scale(na.omit(p)), t)$auc
					})
		})
allModelsCvAuc <- array(allModelsCvAuc, dim=c(length(t),ncol(allModelsCvPredictions[[1]]),length(allModelsCvPredictions)))
plot(NA,NA, xlab="Years",ylab="AUC", xlim=range(t), ylim=c(0.5,0.8))
for(i in 1:dim(allModelsCvAuc)[2]){
	lines(t,rowMeans(allModelsCvAuc, dims=2)[,i], type='l', new=i==1, col=colModels[i])
}
legend("bottomright", colnames(allModelsCvPredictions[[1]]), bty="n", lty=1, col=colModels)


#' Wisdom of the crowds?
foo <- 1
allModelsCvCCrowd <- sapply(allModelsCvPredictions, function(x){
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
			foo <<- foo +1		
			r <- rowMeans(apply(x, 2 , rank))
			survConcordance(osYr[!trainIdx,] ~ r)$concordance
		})
quantile(allModelsCvCCrowd)
boxplot(cbind(t(allModelsCvC),allModelsCvCCrowd), notch=TRUE, ylab="Concordance", border=c(colModels,1), las=2, lty=1, pch=16, staplewex=0)

ranks <- apply(apply(-cbind(t(allModelsCvC),kraut=allModelsCvCCrowd),1,rank, ties.method="random"),1,function(x) table(factor(x, levels=1:8)))
ranks <- ranks[,order(1:8 %*% ranks)]

#' Clean up.. 
rm(allModelsCV)

#' #### Different RFX models
#' Here we assess RFX models with interaction terms for different variable categories.
#+ allModelsCvRfx, cache=TRUE
replicates <- 100 ## number of replicates
allModelsCvRfxC <- do.call("rbind",mclapply(1:100, function(foo){
					set.seed(foo)
					trainIdx <- sample(1:nrow(dataFrameOsTD)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
					coxRFXOsMain <- CoxRFX(dataFrameOsTD[trainIdx,mainIdxOsTD], osTD[trainIdx], groups=groups[mainIdxOsTD])
					coxRFXOsGG <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTDGG], osTD[trainIdx], groups=groups[whichRFXOsTDGG])
					coxRFXOsGGc <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTDGG], osTD[trainIdx], groups=groups[whichRFXOsTDGG], which.mu=mainGroups)
					coxRFXOsAll <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTD], osTD[trainIdx], groups=groups[whichRFXOsTD])
					coxRFXOsAllc <- CoxRFX(dataFrameOsTD[trainIdx,whichRFXOsTD], osTD[trainIdx], groups=groups[whichRFXOsTD], which.mu=mainGroups)
					return(c(
									Main=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,mainIdxOsTD]) %*% coef(coxRFXOsMain))$concordance,
									GeneGene=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTDGG]) %*% coef(coxRFXOsGG))$concordance,
									GeneGeneCentred=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTDGG]) %*% coef(coxRFXOsGGc))$concordance,
									AllInt=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTD]) %*% coef(coxRFXOsAll))$concordance,	
									AllIntCentred=survConcordance(osTD[!trainIdx]~as.matrix(dataFrameOsTD[!trainIdx,whichRFXOsTD]) %*% coef(coxRFXOsAllc))$concordance	
							))
				}, mc.cores=10))
colnames(allModelsCvRfxC) <- sub(".concordant","",colnames(allModelsCvRfxC))

#' #### Supplementary Figure 2

#+ allModelsCvRfxBoxplot, fig.width=2, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(round(t(apply(-allModelsCvRfxC,1,rank)))),factor, levels=1:6)),table)
o <- order(colMeans(allModelsCvRfxC))
boxplot(allModelsCvRfxC[,o], notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:ncol(allModelsCvRfxC), rep(par("usr")[3],ncol(allModelsCvRfxC)), colnames(allModelsCvRfxC)[o])

#+ allModelsCvRfxRank, fig.width=2, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
barplot(r[,o]/replicates, col=set1[c(3,2,4,1,5,7)][1:ncol(allModelsCvRfxC)], ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],ncol(allModelsCvRfxC)), colnames(allModelsCvRfxC)[o])
legend(par("usr")[1],1.5, fill=set1[c(3,2,4,1,5,7)][1:ncol(allModelsCvRfxC)], legend=1:ncol(allModelsCvRfxC), bty="n", border=NA, horiz=TRUE, title="Rank")

#' #### Time-dependent models
#' The following models allow for quantifying the effect of a time-dependent covariate, such as a bone marrow transplant, which is typically
#' administered after diagnosis.

#' The subsequent code is executed on our LSF cluster for 100 replicates
#+ allModelsCVTDCode
read_chunk('../../code/cv100.R', labels="cv100")
#+ cv100, eval=FALSE


#' Gathering results and computing predictions
#+ allModelsCvTdPredictions, cache=TRUE
replicates <- 100
allModelsCvTdPredictions <- mclapply(1:replicates, function(foo) try({
			e <- new.env()
			load(paste0("../../code/cv100/",foo,".RData"), envir=e)
			set.seed(foo)
			x <- list(
					BIC=e$coxBICOsTD,
					AIC=e$coxAICOsTD,
					RFX=e$coxRFXOsTD,
					RFXgg=e$coxRFXOsTDGGc
			)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			pred <- cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[tplSplitOs][!trainIdx]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrameOsTD[!trainIdx,])
							}))
			pred <- cbind(pred, mRFX1yr=colSums(e$multiRfx5[365,1:3,]), mRFX3yr=colSums(e$multiRfx5[3*365,1:3,]), mRFX5yr=colSums(e$multiRfx5[5*365,1:3,]))
			return(pred)
		}), mc.cores=10)

#' Harrel's C
#+ allModelsCvTd-C
allModelsCvTdC <- sapply(1:replicates, function(foo){
			x <- allModelsCvTdPredictions[[foo]] 
			set.seed(foo)
			trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )[tplSplitOs]!=1 ## sample 1/5
			apply(x, 2 , function(p){						
						survConcordance(osYrTD[!trainIdx,] ~ p)$concordance
					})
		})
apply(allModelsCvTdC,1,quantile)


#+ allModelsCvTdCBoxplot, fig.width=1.5, fig.height=1.5
par(mar=c(3,3,1,1),bty="n", mgp=c(2,.5,0), las=2)
r <- sapply(as.data.frame(lapply(as.data.frame(t(apply(-allModelsCvTdC,2,rank))),factor, levels=1:nrow(allModelsCvTdC))),table)
o <- order(apply(allModelsCvTdC,1,median))
boxplot(t(allModelsCvTdC[o,]), notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n")
rotatedLabel(1:nrow(allModelsCvTdC), rep(par("usr")[3],nrow(allModelsCvTdC)), rownames(allModelsCvTdC)[o])

#+ allModelsCvTdCRank, fig.width=1.5, fig.height=1.5
par(mar=c(3,3,3,1), xpd=NA, las=2, mgp=c(2,.5,0))
clr <- brewer.pal(nrow(allModelsCvTdC),"PiYG")#set1[c(3,2,4,1,5,7)]
barplot(r[,o]/replicates, col=clr[1:ncol(allModelsCvTdC)], ylab="Fraction", names.arg=rep("",ncol(r))) -> b
rotatedLabel(b, rep(par("usr")[3],ncol(allModelsCvTdC)), colnames(allModelsCvTdC)[o])
legend(par("usr")[1],1.5, fill=clr[1:nrow(allModelsCvTdC)], legend=1:nrow(allModelsCvTdC), bty="n", border=NA, horiz=TRUE, title="Rank")



#' ### Inter-study CV
#+ allModelsTrial, cache=TRUE
allModelsTrial <- mclapply(levels(clinicalData$Study), function(foo){
			#set.seed(foo)
			trainIdx <- clinicalData$Study != foo 
			c <- coxph(os[trainIdx] ~ 1, data=dataFrame[trainIdx,mainIdxOs])
			scopeStep <- as.formula(paste("os[trainIdx] ~", paste(colnames(dataFrame)[mainIdxOs], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxCPSSOsTrain <- CoxCPSSInteractions(dataFrame[!is.na(os) & trainIdx, mainIdxOs], na.omit(os[trainIdx]), bootstrap.samples=50, scope = which(groups %in% scope))
			w <- colnames(dataFrame[mainIdxOs])
			w <- setdiff(w, names(which(colSums(dataFrame[trainIdx,w])==0)))
			coxRFXOsTrain <- CoxRFX(dataFrame[trainIdx,w], os[trainIdx], groups=groups[w], nu = if(foo=="AMLSG0704") 1 else 0) # add prior for 0704 (just one group member)
			coxRFXOsTrain$Z <- NULL
			w <- whichRFXOsGG
			w <- setdiff(w, which(colSums(dataFrame[trainIdx,w])==0))
			coxRFXOsGGc <- CoxRFX(dataFrame[trainIdx,w], os[trainIdx], groups=groups[w], which.mu=mainGroups, nu = if(foo=="AMLSG0704") 1 else 0)
			coxRFXOsGGc$Z <- NULL
			rForestOsTrain <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs])[trainIdx,], ntree=100, importance="none")
			return(list(
							BIC=coxBICOsTrain,
							AIC=coxAICOsTrain,
							CPSS=coxCPSSOsTrain,
							RFX=coxRFXOsTrain,
							RFXgg=coxRFXOsGGc,
							rForest=rForestOsTrain
					))
		}, mc.cores=3)
names(allModelsTrial) <- levels(clinicalData$Study)

allModelsTrialPredictions <- mclapply(names(allModelsTrial), function(foo){
			x <- allModelsTrial[[foo]]
			trainIdx <- clinicalData$Study != foo
			cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[!trainIdx]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrame[!trainIdx,])
							}))
		}, mc.cores=10)
names(allModelsTrialPredictions) <- names(allModelsTrial)

allModelsTrialC <- sapply(names(allModelsTrial), function(foo){
			trainIdx <- clinicalData$Study != foo
			apply(allModelsTrialPredictions[[foo]], 2 , function(p){						
						unlist( survConcordance(osYr[!trainIdx,] ~ p)[c("concordance","std.err")])
					})
		}, simplify="array")

allModelsTrialC

#' #### Time-dependent
#+ allModelsTrialTD, cache=TRUE
allModelsTrialTD <- mclapply(levels(clinicalData$Study), function(foo){
			#set.seed(foo)
			trainIdxTD <- clinicalData$Study[tplSplitOs] != foo 
			c <- coxph(osTD[trainIdxTD] ~ 1, data=dataFrameOsTD[trainIdxTD,mainIdxOsTD])
			scopeStep <- as.formula(paste("osTD[trainIdx] ~", paste(colnames(dataFrameOsTD)[mainIdxOsTD], collapse="+")))
			coxBICOsTrain <- step(c, scope=scopeStep, k = log(sum(trainIdxTD)), trace=0)
			coxAICOsTrain <- step(coxBICOsTrain, scope=scopeStep, k = 2, trace=0)
			coxRFXOsTrain <- CoxRFX(dataFrameOsTD[trainIdxTD,mainIdxOsTD], osTD[trainIdxTD], groups=groups[mainIdxOsTD], nu = if(foo=="AMLSG0704") 1 else 0)
			coxRFXOsTrain$Z <- NULL
			coxRFXOsGGc <- CoxRFX(dataFrameOsTD[trainIdxTD,whichRFXOsTDGG], osTD[trainIdxTD], groups=groups[whichRFXOsTDGG], which.mu=mainGroups, nu = if(foo=="AMLSG0704") 1 else 0)
			coxRFXOsGGc$Z <- NULL
			
			return(list(
							BIC=coxBICOsTrain,
							AIC=coxAICOsTrain,
							RFX=coxRFXOsTrain,
							RFXgg=coxRFXOsGGc							))
		}, mc.cores=3)
names(allModelsTrialTD) <- levels(clinicalData$Study)

allModelsTrialTdPredictions <- mclapply(names(allModelsTrialTD), function(foo){
			x <- allModelsTrialTD[[foo]]
			trainIdxTD <<- clinicalData$Study[tplSplitOs] != foo
			pred <- cbind(ELN=c(4,1,3,2)[clinicalData$M_Risk[tplSplitOs][!trainIdxTD]],
					sapply(x, function(y){
								predictAllModels(y, newdata=dataFrameOsTD[!trainIdxTD,])
							}))
			whichTrain <<- which(trainIdxTD[1:nrow(dataFrame)])
			rfxNrs <- CoxRFX(nrdData[nrdData$index %in% whichTrain, names(crGroups)], Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
			rfxNrs$coefficients["transplantRel"] <- 0
			rfxPrs <-  CoxRFX(prdData[prdData$index %in% whichTrain, names(crGroups)], Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% whichTrain], groups=crGroups, nu=1, which.mu = intersect(mainGroups, unique(crGroups)))
			rfxRel <-  CoxRFX(relData[relData$index %in% whichTrain, names(crGroups)], Surv(relData$time1, relData$time2, relData$status)[relData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
			rfxRel$coefficients["transplantRel"] <- 0
			rfxCr <- CoxRFX(osData[whichTrain, names(crGroups)], Surv(cr[,1], cr[,2]==2)[whichTrain], groups=crGroups, which.mu = NULL)#intersect(mainGroups, unique(crGroups)))
			rfxEs <- CoxRFX(osData[whichTrain, names(crGroups)], Surv(cr[,1], cr[,2]==1)[whichTrain], groups=crGroups, which.mu = NULL)
			ix <- tplSplitOs %in% whichTrain
			rfxOs <- CoxRFX(dataFrameOsTD[ix,whichRFXOsTDGG], osTD[ix], groups[whichRFXOsTDGG], which.mu=mainGroups) ## allow only the main groups to have mean different from zero.. 
			xx <- 0:2000
			coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])[prdData$index %in% whichTrain,]) 
			tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
			coxphOs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))[osData$index %in% whichTrain,]) 
			tdOsBaseline <- exp(pmin(predict(coxphOs, newdata=data.frame(time0=500)),predict(coxphOs, newdata=data.frame(time0=xx[-1])))) ## cap predictions at induction length 500 days.
			
			dataTD <- data[tplSplitOs, ]
			dataTD$transplantCR1[1:nrow(data)] <- 0
			dataTD$transplantRel[1:nrow(data)] <- 0
			multiRfx5 <- MultiRFX5(rfxEs, rfxCr, rfxNrs, rfxRel, rfxPrs, dataTD[!trainIdxTD,], tdPrmBaseline = tdPrmBaseline, tdOsBaseline = tdOsBaseline, x=2000)
			
			pred <- cbind(pred, mRFX1yr=colSums(multiRfx5[365,1:3,]), mRFX3yr=colSums(multiRfx5[3*365,1:3,]), mRFX5yr=colSums(multiRfx5[5*365,1:3,]))
			return(pred)
			
		}, mc.cores=3)
names(allModelsTrialTdPredictions) <- names(allModelsTrialTD)

allModelsTrialTdC <- sapply(names(allModelsTrialTD), function(foo){
			trainIdx <- clinicalData$Study[tplSplitOs] != foo
			apply(allModelsTrialTdPredictions[[foo]], 2 , function(p){						
						unlist( survConcordance(osYrTD[!trainIdx,] ~ p)[c("concordance","std.err")])
					})
		}, simplify="array")

allModelsTrialTdC

#' ### TCGA validation
#' 
#' #### Fit models
#' Fit a single tree [@Therneau2014A] and a random forest model [@IshwaranTAAS2008].
#+ tree, fig.width=3, fig.height=3
library(rpart)
library(randomForestSRC)
tree <- rpart(os ~ ., data=dataFrame[mainIdxOs & osIdx])
plot(tree)
text(tree)
survConcordance(na.omit(os)~predict(tree))

#' Random forest
#+ rForest, cache=TRUE
rForest <- rfsrc(Surv(time, status) ~.,data= cbind(time = os[,1], status = os[,2], dataFrame[,mainIdxOs & osIdx]), ntree=100)
boxplot(rForest$importance ~ droplevels(groups[mainIdxOs & osIdx]), border= colGroups[mainGroups], staplewex=0, pch=16, cex=0.75, ylab="RSF importance", lty=1, xaxt="n")
rotatedLabel(labels=mainGroups)
rForestVimp <- sapply(mainGroups, function(g) vimp(rForest, colnames(dataFrame)[which(groups==g)]))

survConcordance(na.omit(os)~predict(rForest, importance="none")$predicted)

#' Complementary pairs stability selection with interaction terms
#+ CoxCPSSIntOs, cache=TRUE
set.seed(42)
coxCPSSIntOs <- CoxCPSSInteractions(dataFrame[!is.na(os),groups %in% mainGroups & osIdx], na.omit(os), bootstrap.samples=50, scope = which(groups %in% scope))
selectedIntOs <- names(which(coxCPSSIntOs$Pi > 0.8))
coxCPSSIntOs

#' Stepwise model selection by BIC
#+ coxBIC, cache=TRUE, warning=FALSE
c <- coxph(os ~ 1, data=dataFrame[,mainIdxOs & osIdx])
scopeStep <- as.formula(paste("os ~", paste(colnames(dataFrame)[mainIdxOs& osIdx], collapse="+")))
coxBICOs <- step(c, scope=scopeStep, k = log(sum(trainIdx)), trace=0)
summary(coxBICOs)

#' With AIC
#+ coxAIC, cache=TRUE, warning=FALSE
coxAICOs <- step(c, scope= scopeStep, k = 2, trace=0)
summary(coxAICOs)

#' Time-dep AIC and BIC, including allografts
#+ coxBICosTD, cache=TRUE, warning=FALSE
c <- coxph(osTD ~ 1, data=dataFrameOsTD[mainIdxOsTD])
scopeStep <- as.formula(paste("osTD ~", paste(colnames(dataFrameOsTD)[mainIdxOsTD], collapse="+")))
coxBICOsTD <- step(c, scope=scopeStep, k = log(nrow(dataFrame)), trace=0)
coxAICOsTD <- step(coxBICOsTD, scope=scopeStep, k = 2, trace=0)


#' #### TCGA data
#' Load data
#+ tcgaData, cache=TRUE
tcgaClinical <- read.table("../../data/TCGA_clin.txt", sep="\t", header=TRUE)
tcgaGenetic <- read.table("../../data/TCGA_gen.txt", sep="\t", header=TRUE)
tcgaGenetic$TCGA_ID <- factor(as.character(tcgaGenetic$TCGA_ID), levels = levels(tcgaClinical$TCGA_ID))
g <- as.character(tcgaGenetic$Hugo_Symbol)
g[tcgaGenetic$Hugo_Symbol=="FLT3" & tcgaGenetic$Variant_Type == 'INS'] <- "FLT3_ITD"
g[tcgaGenetic$Hugo_Symbol=="FLT3" & tcgaGenetic$Variant_Type == 'SNP'] <- "FLT3_TKD"
tcgaMutation <- (table(tcgaGenetic$TCGA_ID,g)) + 0
t <- data.frame(tcgaMutation[,]>0, CEBPA_mono = tcgaMutation[,"CEBPA"]==1,CEBPA_bi = tcgaMutation[,"CEBPA"]>1,tcgaClinical[,-c(1,2,4,5,6,13,25)], MakeInteger(tcgaClinical$TypeAML)) + 0
w <- grep("_10+$",colnames(dataFrame), value=TRUE)
f <- as.numeric(sub(".+_","",w))
n <- sub("_10+","",w)
f <- f[n %in% colnames(tcgaClinical)]
n <- n[n %in% colnames(tcgaClinical)]
t[n] <- t[n] / rep(f, each=nrow(t))
colnames(t)[match(n,colnames(t))] <- paste(n,f,sep="_")
rm(w,n,f,g)

tcgaData <- dataFrame[1:nrow(t),]
tcgaData[,] <- NA
w <- intersect(names(t), names(tcgaData))
tcgaData[w] <- t[w]
tcgaData$TPL_os <- NA
tcgaData[groups=="Genetics"][is.na(tcgaData[groups=="Genetics"])] <- 0
tcgaData$MissingCyto <- (tcgaClinical$karyotype == '[Not Available]' )+0
rm(t,w)
tcgaSurvival <- Surv(tcgaClinical$OS/365, tcgaClinical$Status)

tb <- read.xlsx("../../data/TCGA_SupplementalTable01.xlsx", 1, colIndex=1:29)
tb <- tb[order(tb$TCGA.Patient.ID),]

tt <- strsplit(as.character(tb$Trnsplt), ", ")
tp <-  strsplit(as.character(tb$Dz.Stat....trnsplt),", ")
tcgaTpl <- tb(sapply(1:nrow(tb) , function(i){
					transplantCR1=0; transplantRel=0
					if(tt[i] != "0") {
						a <- tt[i]%in%c("MUD","sib Allo") & !tp[i] %in% c("Refr dz","refr dz","refr dz post induction","xxxxx","aplastic post chemo","0") 
						if(any(a)){
							if(any(a & tp[i] %in% c("CR1","CR 1"))) transplantCR1 <- 1
							if(any(a & !tp[i] %in% c("CR1","CR 1"))) transplantRel <- 1
						}
					}
					return(c(transplantCR1=transplantCR1, transplantRel=transplantRel))
				}) )

tcgaData$TPL_os <- tcgaTpl[,"transplantCR1"]

#' Plot mutation frequencies
#+ tcgaFreq, fig.width=2.5, fig.height=2.5
plot(colMeans(dataFrame[groups=="Genetics"])[colnames(tcgaMutation)], colMeans(tcgaMutation))
text(colMeans(dataFrame[groups=="Genetics"])[colnames(tcgaMutation)], colMeans(tcgaMutation), colnames(tcgaMutation))
cor(colMeans(dataFrame[groups=="Genetics"])[colnames(tcgaMutation)], colMeans(tcgaMutation), use='c')

#' NPM1 survival
plot(survfit(tcgaSurvival ~ NPM1, data=tcgaData), col=set1[1:2])
lines(survfit(osYr ~ NPM1, data=dataFrame), col=set1, lty=3,mark=NA)
legend("topright", col=c(set1[1:2],"black","black"), c("NPM1 wt", "NPM1 mut","TCGA","AML"), lty=c(1,1,1,3), bty='n')

#' #### Analyse risk
#' CoxRFX model and covariance-based imputation
tcgaRiskRFXOs <- PredictRiskMissing(coxRFXFitOsTDGGc, tcgaData[whichRFXOsTDGG])
survConcordance(tcgaSurvival ~ tcgaRiskRFXOs[,1])

#' CPSS model
tcgaDataImputed <- as.data.frame(ImputeMissing(dataFrame[mainIdxOs], newX=tcgaData[mainIdxOs]))
tcgaRiskCPSSOs <- predict(coxCPSSIntOs, newdata=tcgaDataImputed)
survConcordance(tcgaSurvival ~ tcgaRiskCPSSOs)

#' Blind imputation (mean only)
f <- function(X) {X <- sapply(X, poorMansImpute);X[is.na(X)] <- 0; X}
survConcordance(tcgaSurvival ~ predict(coxCPSSIntOs, newdata=as.data.frame(f(tcgaData[mainIdxOs]))))

#' Cytogenetic risk
survConcordance(tcgaSurvival ~ c(3,1,2)[tcgaClinical$C_Risk])

#' PINA score [@PastoreJCO2014] for NK AML.
PINAOs <- function(X){
	coef <- c( NPM1=-1.2,
			FLT3_ITD=-.26,
			`NPM1:FLT3_ITD`=.89,
			CEBPA_bi=-1.3,
			wbc_log10=.57,
			age=0.044,
			ecog24=.4)
	x <- cbind(X[,colnames(X) %in% names(coef)], wbc_log10 = log10(100*1e3*pmax(X[,"wbc_100"], 0.001)), age = X[,"AOD_10"]*10, ecog24 = X[,"Performance_ECOG"]>=2)
	risk <- as.matrix(x[,names(coef)]) %*% coef
	group <- cut(risk, c(min(risk), 4,5.4, max(risk)), labels = c("low","int","high"))
	return(data.frame(risk, group))
}
pinaOs <- PINAOs(dataFrame)

#+ PINAos, fig.width=3, fig.height=2.5
nkIdx <- clinicalData$NK == 1
plot(survfit(os[nkIdx] ~ pinaOs[nkIdx,2]), col=rev(set1[1:3]))
survConcordance(os[nkIdx] ~ pinaOs[nkIdx,1])

#' Compared to CPSS (AML data)
survConcordance(os[nkIdx] ~ predict(coxCPSSIntOs, newdata=dataFrame)[nkIdx])

#' And on TCGA data
tcgaPinaOs <- PINAOs(cbind(tcgaDataImputed, `NPM1:FLT3_ITD` = tcgaDataImputed[,"NPM1"]*tcgaDataImputed[,"FLT3_ITD"]))
tcgaNkIdx <- tcgaClinical$karyotype == "Normal"
survConcordance(tcgaSurvival[tcgaNkIdx] ~ tcgaPinaOs[tcgaNkIdx,1])
survConcordance(tcgaSurvival[tcgaNkIdx] ~ tcgaRiskCPSSOs[tcgaNkIdx])


#' ELN score [@DohnerB2010]
ELN <- function(X, nkIdx){
	factor(ifelse(X$inv3_t3_3==1 | X$t_6_9==1 | X$minus5_5q==1 | X$mono17_17p_abn17p==1 | X$minus7==1 | X$complex==1 | X$t_v_11==1,
					"Adverse",
					ifelse(X$t_15_17==1 | X$t_8_21==1 | X$inv16_t16_16==1 | ((X$CEBPA_bi==1 |  X$CEBPA_mono==1 | (X$NPM1==1 & X$FLT3_ITD==0)) & nkIdx),
							"Favorable",
							ifelse(nkIdx & (X$FLT3_ITD==1 | X$NPM1==0 & X$FLT3_ITD==0), 
									"Inter-1", "Inter-2"))), levels=rev(c("Adverse","Inter-1","Inter-2","Favorable")))
}

table(clinicalData$M_Risk, ELN(dataFrame, nkIdx))

#' Other models
tcgaRisk <- data.frame(
		#stdRisk = c(3,1,2)[tcgaClinical$C_Risk],
		ELN = as.numeric(ELN(tcgaDataImputed, tcgaNkIdx)),
		tree = predict(tree, newdata=tcgaDataImputed),
		rForest = predict(rForest, newdata = tcgaDataImputed, importance="none")$predicted,
		PINAos = tcgaPinaOs[,1],
		coxRFX = tcgaRiskRFXOs[,1],
		coxBIC = predict(coxBICOs, newdata=tcgaDataImputed),
		coxAIC = predict(coxAICOs, newdata=tcgaDataImputed),
		coxCPSS = tcgaRiskCPSSOs
)

#' Concordance of all models
#+ concordanceTCGA, fig.width=3, fig.height=2.5
tcgaConcordance <- sapply(tcgaRisk, function(x) {c <- survConcordance(tcgaSurvival ~ x); c(c$concordance, c$std.err)})
tcgaConcordance
o <- order(tcgaConcordance[1,])
barplot(tcgaConcordance[1,o], border=NA, col= set1[-6], las=2, xaxt="n", ylab="Concordance", ylim=c(0.5,0.75), xpd=FALSE) -> b
segments(b,tcgaConcordance[1,o]-tcgaConcordance[2,o],b,tcgaConcordance[1,o]+tcgaConcordance[2,o])
rotatedLabel(b, rep(0.49,length(b)), colnames(tcgaConcordance)[o], srt=45)

#' AUC of all models
#+ aucTCGA, fig.width=3, fig.height=2.5
library(survAUC)
library(survivalROC)
tcgaAUC <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], c(90,365,1000)/365)$auc)
tcgaAUCi <- sapply(tcgaRisk, function(x) AUC.uno(na.omit(os), tcgaSurvival[!is.na(x) & !is.na(tcgaSurvival)], scale(x)[!is.na(tcgaSurvival) &! is.na(x)], sort(na.omit(tcgaSurvival[,1])))$iauc)
o <- order(colMeans(tcgaAUC))
barplot(tcgaAUC[,o], border=1, col= rep(c("grey",set1[-6]),each=3), las=2, xaxt="n", ylab="AUC", beside=TRUE, density=c(NA, 48,24), ylim=c(0.5,0.85), xpd=FALSE) -> b
legend("topleft", bty="n", c("3mo","1yr","3yr"), fill='black', density=c(NA, 48,24))
rotatedLabel(b[seq(3, length(b), 3)], rep(0.49,length(tcgaRisk)), names(tcgaRisk)[o], srt=45)

#' KM curves for four risk categories (quartiles)
#+ kmTCGA, fig.width=3, fig.height=2.5
risk <- cut(tcgaRiskRFXOs[,1], quantile(tcgaRiskRFXOs[,1]), labels=c("1st Q","2nd Q","3rd Q","4th Q"))
s <- survfit(tcgaSurvival ~ risk)
plot(s, col=set1[c(3,2,4,1)], mark=NA, xlab="Years", ylab="Survival")
legend("topright", bty="n", rownames(summary(s)$table), col=set1[c(3,2,4,1)], lty=1)

#' Distribution of risk v cytogenic categories
#+ riskTCGA, fig.width=3, fig.height=2.5
risk <- tcgaRiskRFXOs[,1] - mean(tcgaRiskRFXOs[,1])
x <- seq(from=-4,to=4, l=512)
r <- sapply(levels(tcgaClinical$C_Risk)[c(2,3,1)], function(r){
			i <- tcgaClinical$C_Risk==r
			d <- density(na.omit(risk[i]), from=-4,to=4)$y * mean(i, na.rm=TRUE)
		})
par(mar=c(4,4,3,4)+.1, bty="n")
plot(exp(x),rowSums(r), type='l', lty=0,xlab="Hazard", ylab="Prop. patients", log='x', ylim=c(0,.55))
for(i in 1:3)
	polygon(exp(c(x, rev(x))), c(rowSums(r[,1:i, drop=FALSE]), rev(rowSums(cbind(0,r)[,1:i, drop=FALSE]))), col=set1[c(3,2,1)][i], border=NA)

H0 <- basehaz(coxph(tcgaSurvival ~ risk), centered=TRUE)
hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
invHazardDist <- splinefun(H0$hazard, H0$time, method="monoH.FC")
l <- c(0.1,.5,.9)#c(0.1,0.25,.5,.75,.9)
for(i in seq_along(l))
	lines(exp(x), pmax(0,invHazardDist(-log(l[i]) /exp(x) ))/10000*365, col='black', lty=c(2,1,2)[i])
axis(side=4, at=seq(0,.5,0.1), labels=seq(0,.5,.1)*10000*365)
mtext(side=4, "Time", line=2.5)
mtext(side=3, at = -log(l)/hazardDist(par("usr")[4]*10000*365), text=paste(100*l, "% survive", sep=""))
legend("topright", levels(tcgaClinical$C_Risk)[c(2,3,1)], fill=set1[c(3,2,1)], bty="n", title="M risk")

#' ### Multistage models
d <- tcgaData
d$transplantRel <- tcgaTpl[,"transplantRel"]
d$transplantCR1 <- tcgaTpl[,"transplantCR1"]
d$MissingCyto <- (tcgaClinical$karyotype == '[Not Available]' )+0
multiRfx5Tcga <- MultiRFX5(coxRFXNcdTD, coxRFXCrTD, coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, d, tdPrmBaseline = tdPrmBaseline, tdOsBaseline = tdOsBaseline, x=xmax)

s <- rowMeans(colSums(aperm(multiRfx5Tcga[,1:3,],c(2,1,3))))
plot(survfit(tcgaSurvival ~ 1))
lines(seq(0,2000)/365,1-s)

c <- sapply(seq(1,2000,10), function(i) survConcordance(tcgaSurvival ~  colSums(multiRfx5Tcga[i,1:3,]))$concordance)
plot(c)
plot(seq(1,2000,10),c)

r <- colSums(multiRfx5Tcga[365,1:3,])
survConcordance(tcgaSurvival ~ r)

#' TCGA concordance time-dependent models
tcgaDataTdImputed <- as.data.frame(ImputeMissing(dataFrame[mainIdxOsTD], newX=tcgaData[mainIdxOsTD]))
tcgaRiskTD <- data.frame(
		coxBICTD = predict(coxBICOsTD, newdata=tcgaDataTdImputed),
		coxAICTD = predict(coxAICOsTD, newdata=tcgaDataTdImputed),
		coxRFXTD = PredictRiskMissing(coxRFXFitOsTDGGc, tcgaData)[,1],
		mRFX365 = r
)
tcgaConcordanceTD <- sapply(tcgaRiskTD, function(x) unlist(survConcordance(tcgaSurvival ~ x)[c("concordance","std.err")]))





#' #### Figure 1A
#' Here we generate the overview shown in Figure 1b.
#+ concordanceCvTcga, fig.width=3.5, fig.height=2.5
library(abind)
par(mar=c(3,3.5,.5,.5),bty="n", mgp=c(2.5,.5,0), las=2,  lend=1, xpd=FALSE)
o <- c(1,7,2,3,4,6)
x <- rbind(allModelsCvC[o,], allModelsCvTdC[c("BIC","AIC","RFXgg","mRFX3yr"),])
col <- brewer.pal(4,"Pastel1")
#boxplot(t(x[o,]), notch=TRUE, ylab="Concordance", staplewex=0, lty=1, pch=16, xaxt="n", border="white", ylim=c(0.5,0.75), boxwex=.5)
bplot <- function(x, at=1:ncol(x),..., ylim=range(x), xlab="", col="black", col.lines="grey"){
	y <- apply(x,2,fivenum)
	plot(at,y[3,], pch=NA, ..., ylim=ylim, xlab="", xaxt="n")
	segments(at,y[1,],at,y[5,], col=col.lines, lwd=2)
	segments(at,y[2,],at,y[4,], col=col.lines, lwd=4)
	points(at,y[3,], pch=15, col=col)
}
s <- .2 #space
a <- c(1:6, 7:10+.5)
bplot(t(x), at=a-1.5*s,ylab="Concordance", ylim=c(0.5,0.75), xlim=range(a)+c(-.5,.5))
abline(h=seq(.5,.75,.05), col="lightgrey")
par(xpd=NA)
t <- tcgaConcordance[,c(1,3,6,7,8,5)]
z <- abind(abind(TCGA=t, allModelsTrialC[,o,]), abind(TCGA=tcgaConcordanceTD, allModelsTrialTdC[,-c(1,4),]), along=2)
m <- sapply(1:ncol(z),function(i){
			err <- 1 / sum(1/z[2,i,]^2)
			avg <- sum(z[1,i,] /z[2,i,]^2) * err
			c(avg,sqrt(err))})
#segments(1:6+s/2,m[1,]-m[2,],1:6+s/2,m[1,]+m[2,], lwd=2, col="#00000044")
#points(m[1,], pch=19, cex=1.5)
#segments(a-s/2,t[1,]-t[2,],a-s/2,t[1,]+t[2,], col=paste0(col[1],"FF"), lwd=2)
#points(a-s/2,t[1,], col=col[1], pch=16, cex=1)
i <- 0; for(n in dimnames(z)[[3]]) { i<-i+1;
	segments(a -s +s/2*i, z[1,,n] -  z[2,,n],a -s +s/2*i, z[1,,n]+ z[2,,n], col=paste0(col[i],"FF"), lwd=2)
	points(a -s +s/2*i, z[1,,n], col=col[i], pch=16, cex=1)
}
segments(a -3/4*s, m[1,],a+s*5/4,m[1,], lwd=3)
rotatedLabel(a, labels= rownames(x))
legend("bottomright", 
		c(
				"random CV 4/5 x100", 
				paste0("TCGA, (n=",nrow(na.omit(tcgaSurvival)),")"),
				paste0(dimnames(allModelsTrialC)[[3]]," (n=",table(clinicalData$Study),")"),
				"average"), 
		lty=c(1,1), bg="white", col=c("grey",col[1:4], "black"), pch=c(15,16,16,16,16,16,NA))



#' # Simulations
#' 
#' We use simulations to assess different properties of our risk modelling approach.
#' 
#' ## Survival
#' 
#' Simulating survival times is useful, for example, to verify the consistency of our estimators and obtain 
#' empirical confidence intervals.
#' 
#' In the Cox proportional hazards model, the hazard is given by:
#' 
#' $$ \lambda(t) = \lambda_0(t) \exp(u Z) = -\frac{dS(t)}{dt}\frac{1}{S(t)}.$$
#' 
#' On the transformed time-scale $\tau(t) = \int_0^t \lambda_0(t') dt'$, the 
#' hazard is constant and survival times are distributed exponentially. 
#' A strategy to model survival times according to the Cox proportional
#' hazards model is therefore to draw unit survival times $\tau ~ \operatorname{Exp}(u Z)$ and 
#' to scale those according to $\tau^{-1}$. 
#' 
#' The observed survival times $T_o$ are subject to censoring. The generative process can
#' be thought of as $T_o = \min\{T, T_c\}$, where $T_c$ is a censoring time and $T$ the actual survival. 
#' This process may be simulated by estimating
#' the cumulative distribution of censoring times $\hat F(T_c)$ using the Kaplan-Meier estimator and subsequently
#' simulating censoring times $T_c = \hat F^{-1}(U); U \sim \operatorname{Unif}(0,1)$.
#' 
#' The simulated times and events are then $T_o = \min\{T, T_c\}$ and the status is 1 if $T < T_c$ and 0 otherwise.
#' 
#' Hence our algorithm to simulate survival times can be summarised as follows:
#' 
#' 1. KM estimate of baseline hazard `L_0(t)`
#' 2. Monotonous spline interpolation and inversion `Linv_0(x)`
#' 3. KM estimate of cumulative censoring distribution `Fcens(t)`
#' 4. Monotonous spline interpolation of inverse `Finv`
#' 5. Sample standardised exponential survival times with linear predictor `h`, transform using `T=Linv_0(rexp(h))`
#' 6. Sample follow-up times  `T_c=Finv(runif())`
#' 7. Times: `pmin(T_c,T)`. 
#' 8. Status: `T < S`
#' 
#' This algorithm is implemented in `CoxHD::SimSurvNonp()`.
#' 
#' ## Interpolations
#' 
#' We use interpolations subsampling patients and genes to assess the influence of cohort size and breadth of genomic sequencing
#' on our predictive performance. 
#' 
#' ### Subsampling of genes
#' 
#' We take the following approach:
#' 
#' * Randomly take subsamples of $p' < p_{genes} = 58$ genetic covariates 
#' * Extract gene:gene interaction terms with $\ge 8$ occurrences
#' * Re-estimate the model using all patients.
#' 
#' We observed that the key determinant for the variance of the genetic log hazard is the
#' average number of (genetic) drivers/patient
#' 
#' ### Subsampling of patients
#' 
#' For a given size `n` repeat `r` times:
#' 
#' * Draw a random subsample of patients of size `n` using `sample()`
#' * Train `CoxRFX` model
#' * Predict on the remaining patients and compute concordance
#' 
#' Compute the average concordance across all repetitions `r`. The number of repetitions `r` was chosen for each `n` such that the test size `1540 - n` 
#' was constant in order to achieve a similar error of the average concordance.
#' 
#' ## Extrapolations
#' 
#' Here we use a non-parametric approach to simulate data sets of larger cohorts to extrapolate influence of cohort size on prognostic accuracy.
#' We also use a parametric approach to quantify the relation between the number of genes sequenced and model performance.
#' 
#' ### Patients 
#' To extrapolate to larger cohort size we need to simulate new patients, distributed according to the empirical distribution. 
#' We observed that a simple resampling exaggerates the effect of interaction terms as particular constellations will be overrepresented.
#' We therefore resampled patients and variables and used a multiple imputation package to impute the missing variables, noting that this will
#' more likely generate non-duplicate data points that still satisfy the empirical distribution of the original data.
#' 
#' So we used the following steps:
#' 
#' * Set 10% of variables to NA; impute with `mice` [@van-BuurenJOSS2011] using 10 chains. 
#' * Sample from chains
#' * Sample effect sizes from RFX mean and variance parameters
#' * Keep covariates and interaction terms fixed
#' 
#' This protocol is implemented as `CoxHD::SimData()`
#' 
#' ### Genes
#' 
#' One observation made during [subsampling of genes][#subsampling-of-genes] was that the predicted variation of risk
#' was a linear function of the average number of drivers/patient. Here we derive the theoretical groundwork supporting 
#' this observation.
#' 
#' Let $Z$ be the set of genetic predictors and $u \sim N(\mu;\sigma)$ the distribution of effect sizes. Then the 
#' variation in log hazard is given by
#' 
#' $$
#' \begin{align} Var[h] = Var[u^T Z] &= E[Var[u^T Z|Z]] + Var[E[u^T Z | Z]] \cr
#' &= E[ Z^T Var[u] Z ] + \mu^2 Var[\textstyle\sum_i Z_i]  \cr
#' &= \sigma^2 E[ \textstyle\sum_i Z^2_i] + \mu^2 Var[\textstyle\sum_i Z_i] \cr 
#' &= \sigma^2 E[ \textstyle\sum_i Z_i] + \mu^2 Var[\textstyle\sum_i Z_i] \qquad Z_i^2 = Z_i \in \{0,1\} \cr
#' &= \sigma^2 E[D] + \mu^2 Var[D]
#' \end{align} 
#' $$
#' 
#' Where $D=\sum_i Z_i$ denotes the total number of drivers per patient. The latter term $\mu^2 Var[D]$ can be ignored as long as $D \approx 1$ and $|\mu| < 1$.
#' Hence the variation in the log hazard increases proportionally to the mean number of drivers.
#' 
#' **Note**: These derivations hold for an additive model. 
#' 
#' For interactions
#' $$E[Var[Z^T B Z | Z]] = \sigma^2 E[\sum Z Z^T] =\sigma^2 E[I]$$ 
#' I being the number of interaction terms. $I < D(D-1)/2$.
#' 
#' #### TCGA
#' 
#' On TCGA data we can estimate the number of drivers by means of the sum of the excess of non-synonymous over synonymous mutations at each gene [@Martincorena2015].
#' We use the total number of indels as an upper bound for the number of driver indels.
#' 
#' Using this approach we detect an average of 2.3 point mutations and 1.4 indels adding to 3.7 drivers per AML case when considering the entire exome. This compares to 
#' an average of 1.55 driver substitutions and 0.94 driver indels, with a total of 2.3 mutations (excluding multiple mutations in the same gene) in our cohort, 
#' having sequenced the 111 most prevalent driver genes.
#' 
#' It therefore appears that the variance explained due to the number of genes considered could be increased by approximately 50%.
#' 
#' ## Code

#' ### Interpolations
#' #### Figure 5b
#' Subsampling patients
#+ subsetPatients, cache=TRUE
library(survivalROC)
set.seed(42)
subsets <- seq(100,1500,100)
subsetPatients <- lapply(subsets, function(s){
			mclapply(1:ceiling(50000/(1540-s)), function(foo){
						set.seed(s*foo)
						trn <- 1:nrow(dataFrame) %in% sample(nrow(dataFrame), s)
						tst <-  !trn 
						fit <- CoxRFX(dataFrameOsTD[tplSplitOs[trn], whichRFXOsTDGG], osTD[tplSplitOs[trn]], groups[whichRFXOsTDGG], which.mu=mainGroups, nu = 0.1)
						C <- survConcordance(osTD[tplSplitOs[tst]]~predict(fit, newdata=dataFrameOsTD[tplSplitOs[tst], whichRFXOsTDGG]))
						ROC <- survivalROC(Stime=os[!is.na(os) & tst,1], status=os[!is.na(os) & tst,2], marker = predict(fit, newdata=dataFrame[tst, whichRFXOsTDGG]), predict.time = 850, method="KM", cut.values=seq(-5,5,0.1))
						list(C, ROC, trn, tst, coef(fit))}, mc.cores=10)
		})

#+ subsetPatientsPlot, fig.width=2, fig.height=2
#pdf("subsetConcordance.pdf", 2.5,2.5, pointsize=8)
col1 <- colorRampPalette(set1[c(3,2,4,1,5)])(length(subsets))
plot(NA,NA, xlim=c(0,1),ylim=c(0,1), xlab="FPR",ylab="TPR")
abline(0,1, lty=3)
for(i in seq_along(subsets)){
	x <- sapply(subsetPatients[[i]], function(x) x[[2]]$FP)
	y <- sapply(subsetPatients[[i]], function(x) x[[2]]$TP)
	lines(rowMeans(x),rowMeans(y), col=col1[i], type="l")
}
#legend("bottomright", legend=rev(subsets), lty=1, col=col1[5:1], bty="n")


rangeplot2 <- function(x, y, col = 1, pch = 19, lty = 1, ylim=range(unlist(y)),...){
	plot(x, sapply(y, mean), col = col, pch=pch, ylim = ylim, ..., xaxt="n")
	points(jitter(unlist(sapply(seq_along(y), function(i) rep(x[i], length(y[[i]]))))),unlist(y),pch=1, col=unlist(sapply(seq_along(y), function(i) rep(col[i], length(y[[i]])))), cex=.2) 
	lines(x, sapply(y, mean), lwd=2)
	lines(x, sapply(y, mean) + 2*sapply(y,  sd)/sqrt(sapply(y,length)))
	lines(x, sapply(y, mean) - 2*sapply(y,  sd)/sqrt(sapply(y,length)))
	axis(at = x, labels=x, side=1)
	#segments(x,apply(y,2,min),x,apply(y,2,max), col=col, lty = lty)
}

rangeplot3 <- function(x, y, col = 1, pch = 19, lty = 1, ylim=range(unlist(y)),...){
	plot(x, sapply(y, mean), col = col, pch=pch, ylim = ylim, ...)
	#points(jitter(unlist(sapply(seq_along(y), function(i) rep(x[i], length(y[[i]]))))),unlist(y),pch=1, col=unlist(sapply(seq_along(y), function(i) rep(col[i], length(y[[i]])))), cex=.2) 
	#lines(x, sapply(y, mean), lwd=2)
	s <- sapply(y,  sd)/sqrt(sapply(y,length)) 
	m <- sapply(y, mean)
	segments(x, m+s*2, x, m-s*2, col=col)
	#axis(at = x, labels=x, side=1)
	#segments(x,apply(y,2,min),x,apply(y,2,max), col=col, lty = lty)
}

rangeplot2(x=subsets, y = sapply(subsetPatients, function(x) sapply(x, function(y) y[[2]]$AUC)) , col=col1, xlab="Cohort", ylab="AUC", lty=1, ylim=c(0.7,0.85))
rangeplot2(x=subsets, y = sapply(subsetPatients, function(x) sapply(x, function(y) y[[1]]$concordance)) , col=col1, xlab="Cohort", ylab="Concordance", lty=1, ylim=c(0.65,.75), log='')

rangeplot3(x=subsets, y = sapply(subsetPatients, function(x) sapply(x, function(y) y[[1]]$concordance)) , col=col1, xlab="Cohort", ylab="Concordance", lty=1, ylim=c(0.67,.73), log='')


#lines(x=subsets, y = concordanceFromVariance(sapply(subsetPatients, function(x) {
#					mean(sapply(x, function(y) {
#										h <-  var(as.matrix(dataFrameOsTD[tplSplitOs[y[[3]]],whichRFXOsTDGG]) %*% y[[5]])
#									}))
#				})) , col=1, xlab="Cohort", ylab="Concordance", ylim=c(0.65,.75))
#

#' #### Figure 5a
#' Subsampling genes
#+ subsetGenes, cache=TRUE
set.seed(42)
subsets <- seq(5,55,5)
genes <- names(whichRFXOsTDGG[groups=="Genetics"])
subsetGenes <- lapply(subsets, function(s){
			mclapply(1:100, function(foo){
						g <- sample(genes, s)
						ix <- !grepl(paste(g,collapse="|"), names(whichRFXOsTDGG))
						trainIdx <- sample(1:nrow(dataFrame)%%5 +1 )!=1 ## sample 1/5
						testIdx <-  !trainIdx 
						fit <- CoxRFX(dataFrameOsTD[tplSplitOs[trainIdx], whichRFXOsTDGG[ix]], osTD[tplSplitOs[trainIdx]], groups[whichRFXOsTDGG[ix]], which.mu=mainGroups, nu = 0.1)
						C <- survConcordance(osTD[tplSplitOs[testIdx]]~predict(fit, newdata=dataFrameOsTD[tplSplitOs[testIdx], whichRFXOsTDGG[ix]]))
						ROC <- survivalROC(Stime=os[!is.na(os) & testIdx,1], status=os[!is.na(os) & testIdx,2], marker = predict(fit, newdata=dataFrame[testIdx, whichRFXOsTDGG[ix]]), predict.time = 850, method="KM", cut.values=seq(-5,5,0.1))
						fit <- CoxRFX(dataFrameOsTD[, whichRFXOsTDGG[ix]], osTD, groups[whichRFXOsTDGG[ix]], which.mu=mainGroups, nu = 0.1)
						S <- cov(PartialRisk(fit))
						list(C, ROC, S, trainIdx, testIdx, ix, mean(rowMeans(dataFrame[setdiff(genes,g)])))
					}, mc.cores=10)
		})

#+ subsetGenesPlotTCGA, fig.width=2.5, fig.height=2.5
plot(sapply(subsetGenes, function(x) sapply(x, function(y) y[[7]]*sum(y[[6]][1:58]))), sapply(subsetGenes, function(x) sapply(x, function(y) {t <- try(sum(y[[3]][c("Genetics","GeneGene"),c("Genetics","GeneGene")])); ifelse(class(t)=="try-error",NA,t)})), xlab="Mean no. of drivers", ylab=expression(paste(Var,"[",h[g],"]")), xlim=c(0,3.8), ylim=c(0,.35), pch=16, col=c("#00000044"))
x <- c(0,3.7)
s <- coxRFXFitOsTDGGc$sigma2["Genetics"]
segments(c(2.3, 3.7), rep(par("usr")[3],2), c(2.3, 3.7), c(2.3, 3.7) * s, col="grey")
segments( rep(par("usr")[1],2),  c(2.3, 3.7) * s, c(2.3, 3.7), c(2.3, 3.7) * s, col="grey")
lines(x, x*s, col="red")
par(xpd=NA)
axis(at=c(2.3, 3.7), labels=c("111 genes", "TCGA (exome)"), tcl=0.5, side=1, mgp=c(-2.5,-2,0))

#' ### Extrapolations
#' ##### Generate new data
#' Simulate data using multiple imputation.
#+ simData, cache=TRUE
set.seed(42)
SimDataNonp
d <- as.matrix(dataFrame[mainIdxOsTD])
w <- groups[mainIdxOsTD] %in% c("Genetics","Fusions","CNA")
d[,w][! as.matrix(d[,w]) %in% c(0,1)] <- NA # remove those imputed ones
simData <- SimDataNonp(d, nData = 10000, m=10)
names(simData) <- names(dataFrame[mainIdxOsTD])

#' Merge into data.frame
#+ simDataFrame, cache=TRUE
set.seed(42)
g <- groups[mainIdxOsTD]
for(w in which(colSums(simData,na.rm=TRUE) == 0))
	simData[[w]] <- rbinom(nrow(simData),1,mean(dataFrame[mainIdxOsTD][,w]))
all(colSums(simData,na.rm=TRUE) != 0)
simDataFrame <- cbind(simData,
		MakeInteractions(simData[,g=="Genetics"], simData[,g=="Genetics"])[,as.vector(upper.tri(matrix(0,ncol=sum(g=="Genetics"), nrow=sum(g=="Genetics"))))])
for(n in unique(which(is.na(simDataFrame), arr.ind = TRUE)[,2]))
	simDataFrame[[n]] <- poorMansImpute(simDataFrame[[n]])
simDataFrame <- StandardizeMagnitude(simDataFrame)
simDataFrame <- simDataFrame[,colnames(simDataFrame)  %in% names(whichRFXOsTDGG) | colSums(simDataFrame)>=8]
simDataFrame$`NPM1:FLT3_ITD:DNMT3A` <- simDataFrame$NPM1 * simDataFrame$FLT3_ITD * simDataFrame$DNMT3A
dim(simDataFrame)

#' ##### Basic simulations
set.seed(42)
simGroups <- factor(c(as.character(g), rep("GeneGene", ncol(simDataFrame)-length(g))))
names(simGroups) <- colnames(simDataFrame)
simCoef <- CoxHD:::SimCoef(coxRFXFitOsTDGGc, groups = simGroups)

simRisk <- as.matrix(simDataFrame[names(whichRFXOsTDGG)]) %*% simCoef[names(whichRFXOsTDGG)]
simSurv <- SimSurvNonp(simRisk, os)

survConcordance(simSurv ~ simRisk)

#' ##### Save output
save(coxRFXFitOsTDGGc, whichRFXOsTDGG, simDataFrame, simGroups, os, mainGroups, file="sim2Data.RData")

#' ##### Simulation code
#' The following code is run on the farm
#+ farmulations, cache=FALSE
read_chunk('../../code/Farmulations2.R', labels="farmulationsCode")
#+ farmulationsCode, eval=FALSE

#' ##### Analysis
#' Read files
files <- dir("../../code/simRFX", pattern="Farmulations\\[1-1000\\]*", full.names = TRUE)
tmp <- new.env()
load(files[1], envir = tmp)

#' #### Supplementary Figure 7
#' ##### P-values
#' Plot the P-values as a function of Npu^2.
#+ pVarSchoenfeld, fig.width=2, fig.height=2, cache=TRUE
w <- groups[whichRFXOsTDGG] %in% c("Genetics","Fusions","CNA", "GeneGene") ## Which groups
psi <- mean(os[,2]) ## Fraction of uncensored observations
plot(colSums(simDataFrame[names(whichRFXOsTDGG[w])]) * tmp$simCoef[whichRFXOsTDGG[w]]^2 , CoxHD:::WaldTest( tmp$fit10000)$p[w], log="yx", pch=NA, xlab=expression(psi *N *p *beta^2), ylab="P-value", ylim=c(1e-50,1))
for(f in files[1:50]){
	load(f, envir = tmp)
	points(psi*colSums(simDataFrame[names(whichRFXOsTDGG[w])]) * tmp$simCoef[names(whichRFXOsTDGG[w])]^2 , CoxHD:::WaldTest( tmp$fit10000)$p[w],  col=colGroups[as.character(groups)[whichRFXOsTDGG[w]]], pch=1, cex=.5)
	points(psi*colSums(simDataFrame[tmp$w1000, names(whichRFXOsTDGG[w])]) * tmp$simCoef[names(whichRFXOsTDGG[w])]^2 , CoxHD:::WaldTest( tmp$fit1000)$p[w], col=colGroups[as.character(groups)[whichRFXOsTDGG[w]]], pch=2, cex=.5)
	if(tmp$fit100$iter[1] < 50) ## Exclude simulations without convergence
		points(psi*colSums(simDataFrame[tmp$w100, names(whichRFXOsTDGG[w])]) * tmp$simCoef[names(whichRFXOsTDGG[w])]^2 ,CoxHD:::WaldTest( tmp$fit100)$p[w],  col=colGroups[as.character(groups)[whichRFXOsTDGG[w]]], pch=3, cex=.5)
}
legend("bottomleft", lty=c(0,1),pch=c(1,NA), c("Simulations","Schoenfeld"), bty="n")
x <- 10^seq(-4,4,0.1)
lines(x, pnorm(sqrt(x), lower.tail = FALSE))


#' ##### Power
#' The theoretical power according to Schoenfeld/Schmoor is given by [@SchmoorSM2000]:
power <- function(beta, N, p, psi=0.5, alpha=0.05){
	pnorm(sqrt(N*psi*beta^2*p*(1-p))-qnorm(1-alpha/2))
}

#' #### Figure 5c
#' Plot for observed cases and overlay a few usual suspects
#+ power1540, fig.width=3, fig.height=3
x <- seq(-2,2,0.01)
y <- 10^seq(-4,0,0.01)
colLevels <- colorRampPalette(brewer.pal(9, "Reds")[-(1:2)])(11)
g <- c("Fusions","CNA","Genetics","GeneGene")
xObs <- matrix(exp(rep(coxRFXFitOsTDGGc$mu[g], each=2) + c(-1,1) * rep(sqrt(coxRFXFitOsTDGGc$sigma2[g]),each=2)), nrow=2) ## Mean log haz +/- sd
yObsQ <- sapply(split(colMeans(dataFrameOsTD[whichRFXOsTDGG]), groups[whichRFXOsTDGG]),quantile, c(0.05,0.5,0.95))[,g] ## 5,50,95% frequency quantiles

contour(outer(x,y,function(x,y) power(x,1540,y)), x=exp(x),y=y, log='xy', xlab="Hazard ratio", ylab="Mutation frequency", main="N=1540", col=colLevels)
rect(xObs[1,],yObsQ[1,],xObs[2,],yObsQ[3,], border = colGroups[c("Fusions","CNA","Genetics","GeneGene")])
#segments(exp(coxRFXFitOsTDGGc$mu[g]),yObsQ[1,],exp(coxRFXFitOsTDGGc$mu[g]),yObsQ[3,], col = colGroups[g])
#segments(xObs[1,],yObsQ[2,],xObs[2,],yObsQ[2,], col = colGroups[g])

effects <- c("NPM1","TP53","inv3_t3_3","t_15_17","inv16_t16_16","CEBPA_bi","FLT3_ITD","complex","NPM1:FLT3_ITD:DNMT3A") ## A few interesting variables
points(exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), col=colGroups[as.character(groups[effects])], pch=19)
text(labels=effects,exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), pos=ifelse(sign(coef(coxRFXFitOsTDGGc)[effects])==1,4,2))
legend("bottom", lty=c(1,NA,NA,NA,NA,NA),pch=c(NA,19,22,22,22,22), c("Power","Selected variables", paste("Dist.", g)), col=c(colLevels[10], "black", colGroups[g]), bty="n", ncol=2)

#' Compared to other cohort sizes
#+ power100-10000, fig.width=1.5, fig.height=1.5
for(N in c(100,1000,10000)){
	contour(outer(x,y,function(x,y) power(x,N,y)), x=exp(x),y=y, log='xy', xlab="Hazard ratio", ylab="Mutation frequency", main=paste("N=",N,sep=""), col=colLevels, drawlabels=FALSE)
	rect(xObs[1,],yObsQ[1,],xObs[2,],yObsQ[3,], border = colGroups[g])
}

#' ##### Concordance
#+ concordance100-10000, fig.width=2, fig.height=2, cache=TRUE
C <- sapply(files[1:500], function(f){
			load(f)
			r <- c(sapply(tmp$nData, function(n){
								survConcordance(SimSurvNonp(simRisk[get(paste0('w',n))], os)~get(paste0("fit",n))$linear.predictors)$concordance
							}),survConcordance(SimSurvNonp(simRisk, os)~simRisk)$concordance)
			names(r) <- c(nData,"Truth")
			return(r)
		})
boxplot(t(C), staplewex=0, pch=16, lty=1, ylab="", ylab="Concordance", xaxt="n")
rotatedLabel(labels=(sub(".concordant","", rownames(C))))
abline(h=CoxHD:::ConcordanceFromVariance(var(simRisk)))

#' #### Figure 5e
#' ##### Mean prediction error
#+ predError100-10000, cache=TRUE
load("../../code/sim2Data.RData")
R <- sapply(files[1:100], function(f){
			load(f, envir=.GlobalEnv)
			r <- c(sapply(tmp$nData, function(n){
								f <- get(paste0("fit",n))
								assign("s", get(paste0("w",n)), envir=.GlobalEnv)
								x <- as.matrix(simDataFrame[s, names(coef(f))])
								h <- x %*% coef(f)
								#z <- t(t(x)-colMeans(x))
								#e <- rowSums(z %*% f$var2 * z) 
								#return(mean(e))
								S <- try(survfit(f, newdata = as.data.frame(t(colMeans(x)))))
								if(class(S)[1]=="try-error") return(NA)
								hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
								sapply(seq(0,365*5, by=365/4), function(w){
											w <- which.min(abs(S$time-w))
											p <- S$surv[w]^exp(h - mean(h))
											q <- S$surv[w]^exp(simRisk[s])
											abs(p-q)})
							}))
			names(r) <- c(nData)
			return(r)
		})

#+ predError100-10000Plot, fig.width=2, fig.height=2
q <- sapply(1:nrow(R),function(i) apply(Reduce("rbind", R[i,]),2,quantile, c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE), simplify="array")
contour(q[3,,],x=seq(0,365*5, by=365/4)/365, y=nData , log='y', xlab="Time (years)", ylab="Cohort size", las=1, xlim=c(0,3))
plot(nData,q[3,5,], log='xy', ylim=c(1e-2,1), type='l', ylab="Prediction error", lwd=2, xlab="Cohort size")
polygon(c(nData,rev(nData)), c(q[1,5,], rev(q[5,5,])), border=NA, col="#88888844")
polygon(c(nData,rev(nData)), c(q[2,5,], rev(q[4,5,])), border=NA, col="#88888844")
#rotatedLabel(labels=(sub(".concordant","", rownames(q))))


#' ##### Cohort size
#+ cohort, fig.width=2.5, fig.height=2.5
par(mar=c(3,3,1,1), bty='n', mgp=c(2,0.5,0))
cohort <- function(beta, p, psi=0.5, alpha=0.05, power=0.5){
	(qnorm(1-alpha/2) + qnorm(1-power) )^2 / (beta^2 * psi * p * (1-p))
}
x <- seq(-2,2, 0.01)
y <- 10^seq(-3,0, 0.01)
contour(outer(x,y,function(x,y) cohort(x,y, alpha=0.05/100)), x=exp(x),y=y, log='xy', xlab="Hazard ratio", ylab="Mutation frequency",  col=colLevels, levels=c(10,20,50,100,200,500,1000,2000,5000,10000,20000))
rect(xObs[1,],yObsQ[1,],xObs[2,],yObsQ[3,], border = colGroups[c("Fusions","CNA","Genetics","GeneGene")])
effects <- c("NPM1","TP53","inv3_t3_3","t_15_17","inv16_t16_16","CEBPA_bi","FLT3_ITD","complex","NPM1:FLT3_ITD:DNMT3A") ## A few interesting variables
points(exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), col=colGroups[as.character(groups[effects])], pch=19)
text(labels=effects,exp(coef(coxRFXFitOsTDGGc)[effects]), colMeans(dataFrame[effects]), pos=ifelse(sign(coef(coxRFXFitOsTDGGc)[effects])==1,4,2))
#legend("bottom", lty=c(1,NA,NA,NA,NA,NA),pch=c(NA,19,22,22,22,22), c("Power","Selected variables", paste("Dist.", g)), col=c(colLevels[10], "black", colGroups[g]), bty="n", ncol=2)

#' #### Figure 5d
#' Number of cases needed
#+ cases, fig.width=2.2, fig.height=2
par(mar=c(3,5,1,1), bty='n', mgp=c(2.5,0.5,0))
cases <- function(beta, alpha=0.05, power=0.5, p = 1e-2, psi=0.5){
	(qnorm(1-alpha/2) + qnorm(1-power) )^2 / (beta^2 * (1-p) * psi) 
}
x <- seq(-1,1,0.01)

x0 <- log(c(0.01,0.02,0.05,0.1,0.2,0.5,1)+1)
plot(exp(x), cases(x, alpha=5e-2), log='yx', type='l', xlab="Hazard ratio", ylab="Minimal number of cases", las=1)
#lines(exp(x), cases(x, alpha=1e-2),  type='l', lty=2)
lines(exp(x), cases(x, alpha=1e-3),  type='l', lty=3)
segments(exp(x0), par("usr")[3],exp(x0),cases(x0, alpha=5e-2), col='grey')
segments(exp(x[1]), cases(x0, alpha=5e-2),exp(x0),cases(x0, alpha=5e-2), col='grey')
axis(side=2, at=cases(x0, alpha=5e-2), labels=exp(x0), tcl=.5, line=0, las=2, mgp=c(-2.5,-.5,0), hadj=0)
axis(side=1, at=c(seq(0.1,3,0.1)), labels=rep("",30), tcl=-.2, line=0, las=2)
axis(side=2, at=rep(c(1:10), 4) * 10^rep(1:4, each=10), labels=rep("",40), tcl=-.2, line=0, las=2)
legend("topright", legend=c("P < 0.05 *","P < 0.001 ***"), lty=c(1,3), bty="n")

#' ### Multistage simulations 

#' #### Simulation function
#' The following function simulates data from the 5-stage multistage RFX model
SimSurv5 <- function(coxRFXNcdTD, coxRFXCrTD, coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, coxphOs, coxphPrs, censInd, censCr, censRel){	
	
	## Step 1: Compute KM survival curves and log hazard
	getS <- function(coxRFX, data, max.x=5000) {		
		if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
		data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data))])
		r <- PredictRiskMissing(coxRFX, data, var="var2")
		H0 <- basehaz(coxRFX, centered = FALSE)
		hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
		invHazardDist <- splinefun(c(0,H0$hazard), c(0,H0$time), method="monoH.FC")
		x <- c(0:ceiling(max.x))
		S <- exp(-hazardDist(x))
		return(list(S=S, r=r, x=x, hazardDist=hazardDist, invHazardDist=invHazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
	}
	
	x <- 15000
	kmCr <- getS(coxRFX = coxRFXCrTD, data = data, max.x=max(x))
	kmEs <- getS(coxRFX = coxRFXNcdTD, data = data, max.x=max(x))
	kmCir <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
	kmNrm <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
	kmPrs <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
	
	getCens <- function(surv, n){
		F <- survfit(surv~1)
		FCensInv <- splinefun(F$surv, F$time)
		censTimes <- FCensInv(runif(n,0,1)) ## Simulate censoring times
	}
	
	censIndTimes <- getCens(censInd, nrow(data))
	censCrTimes <- getCens(censCr, nrow(data))
	censRelTimes <- getCens(censRel, nrow(data))
	
	as.data.frame(t(sapply(1:nrow(data), function(i){
								crTime <- edTime <- relTime <- nrdTime <- prdTime <- NA
								status <- 1
								crTime <- kmCr$invHazardDist(rexp(1, exp(kmCr$r[i,1])))
								edTime <- kmEs$invHazardDist(rexp(1, exp(kmEs$r[i,1])))
								firstTime <- pmin(edTime, crTime, censIndTimes[i])
								if(firstTime==censIndTimes[i]){
									edTime <- firstTime
									status <- 0
									crTime <- NA
								}
								if(firstTime==edTime){
									crTime <- NA
								}else{
									edTime <- NA
									rInd <- predict(coxphOs, newdata=data.frame(time0=crTime))
									relTime <- kmCir$invHazardDist(rexp(1, exp(kmCir$r[i,1] + rInd)))
									nrdTime <- kmNrm$invHazardDist(rexp(1, exp(kmNrm$r[i,1] + rInd)))
									secondTime <- pmin(relTime, nrdTime, censCrTimes[i])
									if(secondTime==censCrTimes[i]){
										nrdTime <- secondTime
										relTime <- NA
										status <- 0
									}
									if(secondTime==nrdTime){
										relTime <- NA
									}else{
										nrdTime <- NA
										rCr <- predict(coxphPrs, newdata=data.frame(time0=relTime))
										prdTime <- kmPrs$invHazardDist(rexp(1, exp(kmPrs$r[i,1] + rCr)))
										if(prdTime > censRelTimes[i]){
											prdTime <- min(prdTime, censRelTimes[i])
											status <- 0
										}
									}
								}
								times <- c(crTime=crTime, edTime=edTime, relTime=relTime+crTime, nrdTime=nrdTime+crTime, prdTime=prdTime+crTime+relTime, status=status)
								return(times)
							}, simplify='array')))
}

#' #### Simulate outcomes
#' First prepare the data. Allograft indices:
alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD") # only allografts
alloTimeRel <- clinicalData$TPL_date - clinicalData$Recurrence_date + .5 # +.5 to make > 0
alloTimeRel[!alloIdx | (clinicalData$TPL_date < clinicalData$Recurrence_date & !clinicalData$TPL_Phase %in% c("CR1","RD"))] <- NA

#' Spline fitted transition probabilities.
coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=data.frame(prdData, time0=as.numeric(clinicalData$Recurrence_date-clinicalData$CR_date)[prdData$index])) 
coxphOs <- coxph(Surv(time1,time2, status)~ pspline(time0, df=10), data=data.frame(osData, time0=pmin(500,cr[osData$index,1]))) 

#' Censoring distributions
censInd <- Surv(clinicalData$OS, 1-clinicalData$Status)[is.na(clinicalData$CR_date)]
censCr <- Surv(as.numeric(clinicalData$Date_LF - clinicalData$CR_date), 1-clinicalData$Status)[!is.na(clinicalData$CR_date) & is.na(clinicalData$Recurrence_date)]
censRel <- Surv(as.numeric(clinicalData$Date_LF - clinicalData$Recurrence_date), 1-clinicalData$Status)[!is.na(clinicalData$CR_date) & !is.na(clinicalData$Recurrence_date)]

#' Simulate outcomes
#+ simSurv5, cache=TRUE
set.seed(42)
simSurv5 <- SimSurv5(coxRFXNcdTD, coxRFXCrTD, coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, coxphOs, coxphPrs, censInd, censCr, censRel)

plot(survfit(Surv(apply(simSurv5[,1:5],1,max,na.rm=TRUE), simSurv5$status) ~ 1), xlim=c(0,5000))
lines(survfit(Surv(clinicalData$OS, clinicalData$Status) ~ 1), col='red')

#' #### Estimation based on simulated data
#' Now reestimate models in the scenario of a 10,000 patient cohort
set.seed(42)
simDataFrame$transplantCR1 <- rbinom(nrow(simDataFrame), 1, mean(data$transplantCR1))
simDataFrame$transplantRel <- rbinom(nrow(simDataFrame), 1, mean(data$transplantRel))
simDataSurv5 <- SimSurv5(coxRFXNcdTD, coxRFXCrTD, coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, simDataFrame, coxphOs, coxphPrs, censInd, censCr, censRel)

#' Estimate RFX transition rates
#+ simRfx, cache=TRUE
simCr <- Surv(ifelse(!is.na(simDataSurv5$crTime), simDataSurv5$crTime, simDataSurv5$edTime), !is.na(simDataSurv5$crTime))
simRfxCr <- CoxRFX(simDataFrame[names(crGroups)], simCr, groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))

simNcd <- Surv(ifelse(!is.na(simDataSurv5$edTime), simDataSurv5$edTime, simDataSurv5$crTime), simDataSurv5$status & !is.na(simDataSurv5$edTime))
simRfxEs <- CoxRFX(simDataFrame[names(crGroups)], simNcd, groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))

simRel <- Surv(ifelse(!is.na(simDataSurv5$relTime), simDataSurv5$relTime, simDataSurv5$nrdTime) - simDataSurv5$crTime, !is.na(simDataSurv5$relTime))
simRfxRel <- CoxRFX(simDataFrame[names(crGroups)], simRel, groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))

simNrd <- Surv(ifelse(!is.na(simDataSurv5$relTime), simDataSurv5$relTime, simDataSurv5$nrdTime) - simDataSurv5$crTime, simDataSurv5$status & !is.na(simDataSurv5$nrdTime))
simRfxNrs <- CoxRFX(simDataFrame[names(crGroups)], simNrd, groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))

simPrd <- Surv(simDataSurv5$prdTime - simDataSurv5$relTime, simDataSurv5$status)
simRfxPrs <- CoxRFX(simDataFrame[names(crGroups)], simPrd, groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))

plot(coef(coxRFXCrTD),coef(simRfxCr))
cor(coef(coxRFXCrTD),coef(simRfxCr))

plot(coef(coxRFXNcdTD),coef(simRfxEs))
cor(coef(coxRFXNcdTD),coef(simRfxEs))

plot(coef(coxRFXRelTD),coef(simRfxRel))
cor(coef(coxRFXRelTD),coef(simRfxRel))

plot(coef(coxRFXNrdTD),coef(simRfxNrs))
cor(coef(coxRFXNrdTD),coef(simRfxNrs))

plot(coef(coxRFXPrdTD),coef(simRfxPrs))
cor(coef(coxRFXPrdTD),coef(simRfxPrs))

#' Now compute the multistage RFX model
#+ simMultiRfx5, cache=TRUE 
xmax <- 2000
xx <- 0:ceiling(xmax)
simPrs <- coxph(Surv(prdTime-relTime, status)~ pspline(relTime-crTime, df=10), data=simDataSurv5) 
simPrsBaseline <- exp(predict(simPrs, newdata=data.frame(relTime=xx[-1], crTime=0))) ## Hazard (function of CR length)	

simOs <- coxph(Surv(pmax(nrdTime, prdTime, na.rm=TRUE)-crTime, status)~ pspline(crTime, df=5), data=simDataSurv5) 
simOsBaseline <- exp(predict(simOs, newdata=data.frame(crTime=xx[-1]))) ## Hazard (function of CR length)	

simMultiRfx5 <- MultiRFX5(simRfxEs, simRfxCr, simRfxNrs, simRfxRel, simRfxPrs, data, tdPrmBaseline = simPrsBaseline, tdOsBaseline = simOsBaseline, x=xmax)

plot(colSums(fiveStagePredicted[3*365,1:3,]), colSums(simMultiRfx5[3*365,1:3,]))
abline(0,1)

#' Also compute the predicted benefit of allografts with confidence intervals
#+ simMultiRFX3TplCi, cache=TRUE
d <- osData[1:nrow(dataFrame),]
d$transplantCR1 <- 0
d$transplantRel <- 0
simMultiRFX3TplCi <- MultiRFX3TplCi(simRfxNrs, simRfxRel, simRfxPrs, data=d[,colnames(coxRFXNrdTD$Z)], x=3*365, nSim=200, prsData=prsData) ## others with 200
plot(multiRFX3TplCi["dCr1Rel","hat","os",] , simMultiRFX3TplCi["dCr1Rel","hat","os",], xlab="Benefit 1,540 patients", ylab="Benefit 10,000 patients")
plot(multiRFX3TplCi["dCr1Rel","upper","os",] - multiRFX3TplCi["dCr1Rel","lower","os",], simMultiRFX3TplCi["dCr1Rel","upper","os",]-simMultiRFX3TplCi["dCr1Rel","lower","os",], xlab="CI width 1,540 patients", ylab="CI width 10,000 patients")
abline(0,0.5)

#' #### HSCTs
#' Here we reassess the effect of HSCTs, also cosidering the magnitude of prediction errors in the current data set and based on extrapolated errors.
#' 
#' #### Figure 4d
#' Benefit v number of allografts in CR1
#+ survNallo10000
par(bty="L")
fAlloRelapse <- sum(prdData$transplantRel & clinicalData$AOD[ !is.na(clinicalData$Recurrence_date)][prdData$index] < 60)/sum(relData$status & clinicalData$AOD[relData$index] < 60 ) # fraction of patients that have received a salvage transplant
benefitAllo <- multiRFX3LOO[,"CR1"] - (fAlloRelapse*multiRFX3LOO[,"Relapse"] +(1-fAlloRelapse)*multiRFX3LOO[,"None"])
o <- order(-benefitAllo + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>60,NA,0), na.last=NA)
pRelapse <- 1+multiRFX3TplCi[1:2,1,"aar",] - multiRFX3TplCi[1:2,1,"rs",] ## Relapse probabilities
fRelapse <- sapply(seq_along(o), function(i) mean(c(pRelapse[2,o[1:i]], pRelapse[1,o[-(1:i)]]), na.rm=TRUE)) # Personalised

s <- sapply(seq_along(o), function(i) mean(c(multiRFX3LOO[o[1:i],"CR1"], (1-fAlloRelapse)*multiRFX3LOO[o[-(1:i)],"None"] + fAlloRelapse*multiRFX3LOO[o[-(1:i)],"Relapse"]), na.rm=TRUE))
x <- seq_along(s)/length(s)
plot(x + (1-x)*fRelapse*fAlloRelapse,s, type='l', xlab="Total fraction of allografts", ylab="Survival of eligible patients 3yrs after CR", col=set1[1], xaxs="i", yaxs="i", lty=3)

#n <- c(5,10,15,20,25,30)
#u <- par("usr")
#mtext( at=(u[4] - s[1])*n, text=n,side=3, line=0)
#mtext( at=s[1]+u[2]/n, text=n,side=4, line=0)
#for(i in seq_along(n)) abline(s[1], 1/n[i], col='grey', lty=3)
#lines(seq_along(s)/length(s), s, type='l',col=set1[1], lty=3)

ci <- multiRFX3TplCi["dCr1Rel","upper","os",]-multiRFX3TplCi["dCr1Rel","lower","os",] # 1540 patients
sCi1540 <- rowMeans(sapply(1:10, function(foo){ set.seed(foo)
					o <- order(-benefitAllo + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>60,NA,0) + rnorm(1540,sd=ci/4), na.last=NA)
					s <- sapply(seq_along(o), function(i) mean(c(multiRFX3LOO[o[1:i],"CR1"], (1-fAlloRelapse)*multiRFX3LOO[o[-(1:i)],"None"] + fAlloRelapse*multiRFX3LOO[o[-(1:i)],"Relapse"]), na.rm=TRUE))
				}))
lines(x + (1-x)*fRelapse*fAlloRelapse, sCi1540, type='l',col=set1[1], lty=1)
#w <- max(which(abs(sCi1540-sCi1540[1] - 1/10 * seq_along(sCi1540)/length(sCi1540))<1e-5))
#points(w/length(sCi1540), sCi1540[w], pch=19, col=set1[1])

simCi <- simMultiRFX3TplCi["dCr1Rel","upper","os",]-simMultiRFX3TplCi["dCr1Rel","lower","os",]

sCi10000 <- rowMeans(sapply(1:10, function(foo){ set.seed(foo)
					o <- order(-benefitAllo + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>60,NA,0) + rnorm(1540,sd=simCi/4), na.last=NA)
					s <- sapply(seq_along(o), function(i) mean(c(multiRFX3LOO[o[1:i],"CR1"], (1-fAlloRelapse)*multiRFX3LOO[o[-(1:i)],"None"] + fAlloRelapse*multiRFX3LOO[o[-(1:i)],"Relapse"]), na.rm=TRUE))
				}))
lines(x + (1-x)*fRelapse*fAlloRelapse, sCi10000, type='l',col=set1[1], lty=2)
p <- order(na.zero(c(1,4,2,3)[clinicalData$M_Risk])  + dataFrame$AOD_10/20 + ifelse(is.na(clinicalData$CR_date),NA,0) + ifelse(clinicalData$AOD>60,NA,0), na.last=NA)
fRelapseEln <- sapply(seq_along(p), function(i) mean(c(pRelapse[2,p[1:i]], pRelapse[1,p[-(1:i)]]), na.rm=TRUE)) # ELN
sEln <- sapply(seq_along(p), function(i) mean(c(multiRFX3LOO[p[1:i],"CR1"], (1-fAlloRelapse)*multiRFX3LOO[p[-(1:i)],"None"] + fAlloRelapse*multiRFX3LOO[p[-(1:i)],"Relapse"]), na.rm=TRUE))
x <- seq_along(sEln)/length(sEln)

lines(x + (1-x)*fRelapseEln*fAlloRelapse,sEln, sEln, type='l', col=set1[2])
legend("bottomright", c("Personalised risk", "Idealised","10,000 patients","This cohort", "Standard risk","ELN and age"),  col=set1[c(NA,1,1,1,NA,2)],lty=c(NA,3,2,1,NA,1), bty="n", text.font=c(2,1,1,1,2,1))



#' # R session
#' This document was written entirely in R with markdown annotation. It was compiled with `knitr::spin()` [@Xie2015] and `pandoc` using the `rmarkdown` package [@Allaire2015]:
#+ compile, eval=FALSE
rmarkdown::render("SupplementaryMethodsCode.R")
#' The total runtime is approximately 24h using 10 cores. This excludes the extrapolations, which were run on a a computing grid.
#' 
#' The packages and specifics of the R session are:
#+ sessionInfo, eval=TRUE
library(devtools)
devtools::session_info()
sessionInfo()
#' 
#' # References
