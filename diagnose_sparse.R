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
L = ruvSimulate(m = 1000, n = 20000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = FALSE)
Y = L$Y
M = L$M
Y_sp = as(Y, "dgCMatrix")
M_sp = as(M, "dgCMatrix")
mean(Y == 0)
mean(M == 0)

profvis::profvis({
  t(Y) %*% M
  t(Y_sp) %*% M_sp
})


Y = matrix(rnorm(10000*1000) %>% as.numeric(), nrow = 1000)
M = matrix(rpois(10000*1000, lambda = 0.1) %>% as.numeric(), ncol = 1000)
dim(Y)
dim(M)
Y_sp = as(Y, "dgCMatrix")
M_sp = as(M, "dgCMatrix")
pryr::object_size(Y)
pryr::object_size(Y_sp)

profvis::profvis({
  Y %*% M
  # eigenResidop(Y, M)
  Y_sp %*% M_sp
  # eigenSpResidop(Y_sp, M_sp)
})

profvis::profvis({
  M %*% Y
  M_sp %*% Y_sp
})



## whenever there is a matrix operation, the speed associated with DA is very slow. Why is this?
## in comparison, the speed decrease associated with SVD is acceptable