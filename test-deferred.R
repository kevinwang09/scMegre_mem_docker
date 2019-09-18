# library(Matrix)
# set.seed(1234)
# m = 200
# n = 20000
# L = ruvSimulate(m = m, n = n, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
# # Y = t(log2(L$Y + 1L)); M = L$M; ctl = L$ctl; batch = L$batch;
# Y = as.matrix(t(rsparsematrix(m, n, density = 0.01))); M = L$M; ctl = L$ctl; batch = L$batch;
# 
# scale_res <- standardize2(Y, batch)
# normY <- as.matrix(scale_res$s.data)
# geneSdMat <- sqrt(scale_res$stand.var)
# geneMeanVec <- scale_res$stand.mean
# 
# all.equal((Y - geneMeanVec)/geneSdMat, 
#           normY)
# 
# dm_Y = t(BiocSingular::DeferredMatrix(Matrix(t(Y), sparse = TRUE), center = geneMeanVec, scale = geneSdMat))
# dm_M = BiocSingular::DeferredMatrix(Matrix(M, sparse = TRUE))
# t_dm_M = t(dm_M)
# 
# all.equal(as.matrix(dm_Y[1:5,1:5]), 
#           normY[1:5,1:5])
# 
# mean(Y == 0)
# mean(M == 0)
# pryr::object_size(Y)
# pryr::object_size(M)
# pryr::object_size(dm_Y)
# pryr::object_size(dm_M)
# 
# A = t(Y)
# B = M
# dm_A = t(dm_Y)
# dm_B = dm_M
# I = diag(m)
# dm_I = DelayedArray::DelayedArray(Matrix(I, sparse = TRUE))
# 
# profvis::profvis({ ## Why so slow??
#   B2 = (B %*% solve(t(B) %*% B) %*% t(B))
#   dm_B2 = (dm_B %*% solve(t(dm_B) %*% dm_B) %*% t(dm_B))
# })
# 
# dim(B2)
# dim(I)
# 
# dim(dm_B2)
# dim(dm_I)
# 
# profvis::profvis({
#   (I - B2) %*% A
#   BiocSingular::DeferredMatrix(Matrix::Matrix(dm_I - dm_B2, sparse = TRUE, m, m)) %*% dm_A
#   # eigenResidop(A, B)
# })
# 