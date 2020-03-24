library(Matrix)
library(BiocSingular)
set.seed(1234)
m = 100
n = 20000
L = scMerge::ruvSimulate(m = m, n = n, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
sp_Y = t(rsparsematrix(m, n, density = 0.5))
Y = as.matrix(sp_Y); M = L$M;
sp_M = Matrix::Matrix(M, sparse = TRUE)

mean(Y == 0)
mean(M == 0)
pryr::object_size(Y)
pryr::object_size(M)
pryr::object_size(sp_Y)
pryr::object_size(sp_M)

A = t(Y)
B = M
sp_A = t(sp_Y)
sp_B = sp_M
I = diag(m)
sp_I = Matrix(I, sparse = TRUE)

profvis::profvis({
  IB2 = I - (B %*% solve(t(B) %*% B) %*% t(B))
  IB2A = IB2 %*% A
  exact_svd = BiocSingular::runExactSVD(IB2A, k = 20)
  BiocSingular::runIrlbaSVD(IB2A, k = 20)
  BiocSingular::runRandomSVD(IB2A, k = 20)
})

profvis::profvis({
  sp_IB2 = sp_I - (sp_B %*% solve(t(sp_B) %*% sp_B) %*% t(sp_B))
  sp_IB2A = sp_IB2 %*% sp_A
  svd_sp_IB2 = sparsesvd::sparsesvd(sp_IB2, k = 20)
  svd_sp_A = sparsesvd::sparsesvd(sp_A, k = 20)
  svd_sp_IB2A = sparsesvd::sparsesvd(sp_IB2A, k = 20)
})

recon_mat = function(res){
  res$u %*% diag(res$d) %*% t(res$v)
}

appx_svd = recon_mat(svd_sp_IB2) %*% recon_mat(svd_sp_A)
exact_svd_mat = recon_mat(exact_svd)
colCor = function(j, x, y){
  cor(x[,j], y[,j])
}
rowCor = function(i, x, y){
  cor(x[i,], y[i,])
}


dim(exact_svd_mat)
dim(appx_svd)

rowCorResult = sapply(1:nrow(appx_svd), 
                      function(i){rowCor(i = i, x = exact_svd_mat, y = appx_svd)})
colCorResult = sapply(1:ncol(appx_svd), 
                      function(j){colCor(j = j, x = exact_svd_mat, y = appx_svd)})
boxplot(rowCorResult)
boxplot(colCorResult)
