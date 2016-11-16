# Code accompanying Gerstung *et al.*, Precision ocology for acute myeloid leukemia using a knowledge bank approach., *Nature Genetics*, (in press).

This repository contains all code used for running the analysis in the above manuscript. It was used and tested in `R-3.1.2`. 

The main code can be found in `doc/SupplementaryMethodsCode.R`. It was run in `R` using the command 

   > rmarkdown::render("SupplementaryMethodsCode.R")
   
Additional code chunks can be found in `code`, mainly for parallelising some heavy computations in an LSF environment, and also code for the `multistage` prediction webtool. 
Data is located in the `data` subfolder.