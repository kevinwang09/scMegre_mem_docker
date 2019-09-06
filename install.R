# DO NOT EDIT 'install.R'; instead, edit 'install.R.in' and
# use 'rake' to generate 'install.R'.

install.packages("devtools", repos="https://cran.rstudio.com")
devtools::install_github("SydneyBioX/scMerge", ref = "master")