## I want to test if DeferredArray with a SparseMatrix 
## Backend will speed up the calculations (since standardise2 is not needed and sparse residop could be made fast)

library(Matrix)
library(profvis)
set.seed(1234)
m = 1000
n = 20000
L = ruvSimulate(m = m, n = n, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
# Y = t(log2(L$Y + 1L)); M = L$M; ctl = L$ctl; batch = L$batch;
Y = as.matrix(t(rsparsematrix(m, n, density = 0.01))); 
M = L$M; 
ctl = L$ctl; 
batch = L$batch;
my_residop = function (A, B) 
{return(A - B %*% solve(DelayedArray::t(B) %*% B) %*% DelayedArray::t(B) %*% A)}

## Operation with only Y as sparse matrix
profvis::profvis({
  scale_res <- standardize2(Y, batch)
  normY <- scale_res$s.data
  geneSdMat <- sqrt(scale_res$stand.var)
  geneMeanVec <- scale_res$stand.mean
  res1 = my_residop(t(normY), M)
})



## Operation with both Y, M as sparse matrix
profvis::profvis({
  scale_res <- standardize2(Y, batch)
  normY <- scale_res$s.data
  geneSdMat <- sqrt(scale_res$stand.var)
  geneMeanVec <- scale_res$stand.mean
  sp_M = Matrix(M, sparse = TRUE)
  res2 = my_residop(t(normY), sp_M)
})


## Operation with Deferred Y and normal M
profvis::profvis({
  scale_res <- standardize2(Y, batch)
  normY <- scale_res$s.data
  geneSdMat <- sqrt(scale_res$stand.var)
  geneMeanVec <- scale_res$stand.mean
  
  dm_tY = BiocSingular::DeferredMatrix(
    Matrix(t(Y), sparse = TRUE), 
    center = geneMeanVec, scale = geneSdMat)
  
  sp_M = Matrix(M, sparse = TRUE)
  res3 = my_residop(dm_tY, sp_M)
  my_residop(dm_tY, M)
  eigenResidop(as.matrix(dm_tY), M)
})