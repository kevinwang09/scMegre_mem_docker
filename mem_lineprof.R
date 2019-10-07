library(lineprof)
library(BiocSingular)
library(DelayedArray)
library(HDF5Array)
library(scMerge)
set.seed(12345)

# x_mat = matrix(rnorm(8734 * 6791), nrow = 8734)
x_mat = matrix(rnorm(100 * 1000), nrow = 100)
x_hdf = as(x_mat, "HDF5Array")
# x_mat_da = DelayedArray::DelayedArray(x_mat)
# x_hdf_da = DelayedArray::DelayedArray(x_hdf)

res_mat = lineprof({BiocSingular::runRandomSVD(x_mat, k = 20)})
# res_mat_da = lineprof({BiocSingular::runRandomSVD(x_mat_da, k = 20)})
res_hdf = lineprof({BiocSingular::runRandomSVD(x_hdf, k = 20)})
# res_hdf_da = lineprof({BiocSingular::runRandomSVD(x_hdf_da, k = 20)})

View(print(res_mat))
View(print(res_hdf))

lineprof(solve(x_mat %*% t(x_mat)))
lineprof(solve(x_hdf %*% t(x_hdf)))
pryr::object_size(x_mat)
pryr::object_size(x_hdf)
##############################################
# L = ruvSimulate(m = 20, n = 6000, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
# Y = L$Y; M = L$M; ctl = L$ctl
# 
# res_mat = lineprof(scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, BSPARAM = BiocSingular::RandomParam(), svd_k = 20))
# Y_hdf = as(L$Y, "HDF5Array"); M_hdf = as(L$M, "HDF5Array");
# res_hdf = lineprof(scMerge::fastRUVIII(Y = Y_hdf, M = M_hdf, ctl = ctl, k = 20, BSPARAM = BiocSingular::RandomParam(), svd_k = 20))

# profvis::profvis(ruv::residop(Y, M))
# profvis::profvis(ruv::residop(Y_hdf, M_hdf))
# 
# my_residop = function (A, B) {return(A - B %*% Matrix::solve(DelayedArray::t(B) %*% B) %*% DelayedArray::t(B) %*% A)}
# profvis::profvis(my_residop(Y, M))
# profvis::profvis(my_residop(Y_hdf, M_hdf))

# Y_hdf_da = DelayedArray::DelayedArray(Y_hdf); M_hdf_da = DelayedArray::DelayedArray(M_hdf)
# res_hdf_da = lineprof(scMerge::fastRUVIII(Y = Y_hdf_da, M = M, ctl = ctl, k = 20, BSPARAM = BiocSingular::RandomParam(), svd_k = 20))