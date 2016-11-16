# Code for the AML-multistage calculator

The predictive tool is for **research use only**.

To run this code you'll need a recent R installation and the following packages:

* CoxHD (from http://github.com/mg14/CoxHD)
* shiny
  
In the parent directory of `aml-multistage`, run `R` and type:

   > shiny::runApp("aml-multistage")

This will launch a browser with the webportal. A stable web-based installation can be found at http://cancer.sanger.ac.uk/aml-multistage.

The multistage calculator is also available as a docker container `gerstunglab/aml-multistage`. See https://hub.docker.com/r/gerstunglab/aml-multistage for details. 