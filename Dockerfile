# The suggested name for this image is: bioconductor/release_base.

FROM kevinwang09/scmerge

MAINTAINER kevin.wang@sydney.edu.au

ADD install.R /home/
ADD pbmc.R /home/rstudio/
ADD pbmc_mat_mult.R /home/rstudio/

# Running install
RUN R -f /home/install.R
RUN sudo apt-get update
RUN sudo apt-get install htop