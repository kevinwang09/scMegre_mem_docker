library(scMerge)
library(Matrix)
library(DelayedArray)
library(profvis)
library(microbenchmark)
library(dplyr)

set.seed(1)
# L = ruvSimulate(m = 100, n = 20000, nc = 400, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
# Y_mat = L$Y; M_mat = L$M; ctl = L$ctl
n = 100
p = 10000
Y_mat = matrix(rpois(n*p, lambda = 0.1) %>% as.numeric(), nrow = p)
M_mat = matrix(rpois(n*p, lambda = 0.1) %>% as.numeric(), nrow = n)
Y_sp = as(Y_mat, "dgCMatrix")
M_sp = as(M_mat, "dgCMatrix")
Y_mat_da = DelayedArray(Y_mat)
M_mat_da = DelayedArray(M_mat)
Y_sp_da = DelayedArray(Y_sp)
M_sp_da = DelayedArray(M_sp)

## Whenever DA is used, there is a slow down. 
# system.time({mat_output <- scMerge::fastRUVIII(Y = Y_mat, M = M_mat, ctl = ctl, k = 20, BSPARAM = BiocSingular::ExactParam())})
# system.time({mat_output <- scMerge::fastRUVIII(Y = Y_sp, M = M_sp, ctl = ctl, k = 20, BSPARAM = BiocSingular::ExactParam())})
# system.time({mat_output <- scMerge::fastRUVIII(Y = Y_mat_da, M = M_mat_da, ctl = ctl, k = 20, BSPARAM = BiocSingular::ExactParam())})
# system.time({mat_output <- scMerge::fastRUVIII(Y = Y_sp_da, M = M_sp_da, ctl = ctl, k = 20, BSPARAM = BiocSingular::ExactParam())})


## SVD takes the same time. So the slowdown is purely in the DA matrix multiplication 
# profvis::profvis({
#   scMerge::fastRUVIII(Y = Y_mat, M = M_mat, ctl = ctl, k = 20, BSPARAM = BiocSingular::ExactParam())
#   scMerge::fastRUVIII(Y = Y_mat_da, M = M_mat_da, ctl = ctl, k = 20, BSPARAM = BiocSingular::ExactParam())
# })

# microbenchmark(
#   eigenMatMult(t(M_mat), Y_mat),
#   DelayedArray::t(M_mat) %*% Y_mat,
#   DelayedArray::t(M_sp) %*% Y_sp,
#   DelayedArray::t(M_mat_da) %*% Y_mat_da, times = 5)

p = BiocParallel::bpstart(MulticoreParam(workers = 10))
# register(p)
# bpparam()


microbenchmark(
  eigenMatMult(M_mat, Y_mat),
  M_mat %*% Y_mat,
  M_sp %*% Y_sp,
  M_mat_da %*% Y_mat_da, times = 1)


setAutoBPPARAM(p)
getAutoBPPARAM()

microbenchmark(
  eigenMatMult(M_mat, Y_mat),
  M_mat %*% Y_mat,
  M_sp %*% Y_sp,
  M_mat_da %*% Y_mat_da, 
  M_sp_da %*% Y_sp_da, times = 5)


microbenchmark(
  eigenMatMult(Y_mat, M_mat),
  Y_mat %*% M_mat,
  Y_sp %*% M_sp,
  Y_mat_da %*% M_mat_da,
  Y_sp_da %*% M_sp_da,times = 5)

# DelayedArray::setDefaultBPPARAM(p)
  