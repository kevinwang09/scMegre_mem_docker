library(TENxPBMCData)
library(scMerge)
library(DelayedArray)
# NOTE: This turns off parallelisation, which we do for didactic purposes.
setAutoBPPARAM(SerialParam())
DelayedArray:::set_verbose_block_processing(TRUE)

pbmc3k <- TENxPBMCData('pbmc3k')
setRealizationBackend("HDF5Array")
getRealizationBackend()

d = counts(pbmc3k[1:11000,])
pryr::mem_change({m <- d %*% t(d)})
pryr::object_size(d)
rm(m)

s = as(d, "dgCMatrix")
pryr::mem_change({m <- s %*% t(s)})
pryr::object_size(s)
rm(m)
rm(s)

r = as(d, "matrix")
pryr::mem_change({m <- r %*% t(r)})
pryr::object_size(r)
rm(m)
rm(r)