BiocManager::install(c("SingleCellExperiment", "DelayedArray", 
                       "DelayedMatrixStats", "BiocSingular",
                       "pryr", "profvis",
                       "tidyverse", "furrr", 
                       "Matrix", "TENxPBMCData"), version = "3.10", update = TRUE, ask = FALSE)
install.packages("devtools", repos="https://cran.rstudio.com")
devtools::install_github("SydneyBioX/scMerge", ref = "biocsing")