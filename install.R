BiocManager::install(c("SingleCellExperiment", "DelayedArray", 
                       "DelayedMatrixStats", "BiocSingular",
                       "pryr", "profvis",
                       "tidyverse", "furrr", 
                       "Matrix", "TENxPBMCData",
                       "scater", "M3Drop", 
                       "ruv", "devtools",
                       "hadley/lineprof"), version = "3.10", update = TRUE, ask = FALSE)

BiocManager::install("SydneyBioX/scMerge", ref = "master")