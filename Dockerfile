# The suggested name for this image is: bioconductor/release_base.

FROM bioconductor/release_core2

MAINTAINER kevin.wang@sydney.edu.au

ADD install.R /home/
# ADD create_liver.R /home/rstudio/
ADD pbmc.R /home/rstudio/

# Running install
RUN R -f /home/install.R
RUN sudo apt-get update
RUN sudo apt-get install htop