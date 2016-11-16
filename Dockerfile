FROM rocker/shiny
MAINTAINER Moritz Gerstung <moritz.gerstung@ebi.ac.uk>
RUN apt-get install -y libssl-dev
RUN R -e 'install.packages("devtools"); devtools::install_github("mg14/mg14"); devtools::install_github("mg14/CoxHD/CoxHD")' 
RUN rm -rf /srv/shiny-server/*
COPY ./aml-multistage /srv/shiny-server/aml-multistage
