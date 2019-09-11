library(scMerge)
library(SummarizedExperiment)
library(DelayedArray)

set.seed(1)
L = ruvSimulate(m = 200, n = 1000, nc = 200, 
                nCelltypes = 3, nBatch = 2, lambda = 0.1, k = 10, sce = TRUE)
mean(counts(L) == 0)
print(L)

# assay(L, "counts") = DelayedArray(assay(L, "counts"))
# assay(L, "logcounts") = DelayedArray(assay(L, "logcounts"))

da = L
assay(da, "counts") = DelayedArray(as(assay(da, "counts"), "dgCMatrix"))
assay(da, "logcounts") = DelayedArray(as(assay(da, "logcounts"), "dgCMatrix"))

pryr::object_size(L)
pryr::object_size(da)


profvis::profvis({
  matrix_output <- scMerge(sce_combine = L,
                           ctl = paste0('gene', 1:500),
                           cell_type = L$cellTypes,
                           ruvK = 10,
                           assay_name = 'scMerge')
})

profvis::profvis({
  da_output <- scMerge(sce_combine = da,
                       ctl = paste0('gene', 1:500),
                       cell_type = da$cellTypes,
                       ruvK = 10,
                       assay_name = 'scMerge')
})



## whenever there is a matrix operation, the speed associated with DA is very slow. Why is this?
## in comparison, the speed decrease associated with SVD is acceptable