BiocManager::install(version = "3.11", update = TRUE, ask = FALSE)
BiocManager::install(c("SingleCellExperiment", "DelayedArray", 
                       "DelayedMatrixStats", "BiocSingular",
                       "pryr", "profvis",
                       "tidyverse", "furrr", 
                       "Matrix", "TENxPBMCData",
                       "scater", "M3Drop", 
                       "ruv", "devtools",
                       "hadley/lineprof", "shiny"), version = "3.11", update = TRUE, ask = FALSE)

BiocManager::install("SydneyBioX/scMerge", ref = "master")