# DO NOT EDIT FILES CALLED 'Dockerfile'; they are automatically
# generated. Edit 'Dockerfile.in' and generate the 'Dockerfile'
# with the 'rake' command.

# The suggested name for this image is: bioconductor/release_base.

FROM gcr.io/scmerge/scmerge_mem_docker:master

MAINTAINER kevin.wang@sydney.edu.au

ADD install.R /home/
ADD create_liver.R /home/rstudio/

# Running install
RUN R -f /home/install.R
RUN sudo apt-get update
RUN sudo apt-get install htop