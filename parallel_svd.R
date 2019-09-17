library(Matrix)
library(BiocSingular)
library(BiocParallel)
n = 1000
p = 20000
A = rsparsematrix(n*p, n, p, density = 0.01)
system.time(BiocSingular::runIrlbaSVD(A, k = 20))
# system.time(BiocSingular::runIrlbaSVD(A, k = 20, BPPARAM = BiocParallel::MulticoreParam(workers = 14)))
# install.packages("sparsesvd")
system.time(sparsesvd::sparsesvd(A, rank = 20))