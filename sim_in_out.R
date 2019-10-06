library(lineprof)
library(BiocSingular)
library(DelayedArray)
library(HDF5Array)
library(scMerge)
set.seed(12345)


standardize2 <- function(Y, batch) {
  num_cell <- ncol(Y)
  num_batch <- length(unique(batch))
  batch <- as.factor(batch)
  stand_mean <- DelayedArray::rowMeans(Y)
  design <- stats::model.matrix(~-1 + batch)
  B.hat = solve_axb(a = t(design) %*% design,
                    b = t(Y %*% design))
  
  stand_var <- DelayedArray::rowSums(((Y - t(B.hat) %*% t(design))^2))/(num_cell - num_batch)
  stand_Y <- (Y-stand_mean)/sqrt(stand_var)
  return(res = list(stand_Y = stand_Y, 
                    stand_mean = stand_mean, 
                    stand_var = stand_var))
}

solve_axb = function(a, b){
  x = solve(t(a) %*% a) %*% t(a) %*% b
  return(x)
}
##############################################
m = 100;
n = 1000;
L = ruvSimulate(m = m, n = n, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
Y = L$Y; M = L$M; ctl = L$ctl; batch = L$batch;
Y_hdf = as(Y, "HDF5Array")
Y_hdf_da = DelayedArray::DelayedArray(Y_hdf)

res_mat = lineprof({BiocSingular::runRandomSVD(t(standardize2(t(Y), batch)$stand_Y), k = 20)})
res_hdf = lineprof({BiocSingular::runRandomSVD(t(standardize2(t(Y_hdf), batch)$stand_Y), k = 20)})
# res_hdf_da = lineprof({BiocSingular::runRandomSVD(Y_hdf_da, k = 20)})

L = ruvSimulate(m = m, n = n, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = TRUE)
SummarizedExperiment::assay(L, "counts") = as(SummarizedExperiment::assay(L, "counts"), "HDF5Array")
SummarizedExperiment::assay(L, "logcounts") = as(SummarizedExperiment::assay(L, "logcounts"), "HDF5Array")

res_hdf_sce = lineprof({
  L_result1 <- scMerge(
    sce_combine = L,
    ctl = paste0("gene",1:100),
    cell_type = L$cellTypes,
    assay_name = 'scMerge')
})
View(print(res_mat, 3))
View(print(res_hdf, 3))
View(print(res_hdf_sce, 3))
