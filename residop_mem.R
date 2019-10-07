library(lineprof)
library(DelayedArray)
library(HDF5Array)
set.seed(12345)
my_residop = function (A, B) {
  tBB = DelayedArray::t(B) %*% B
  tBB_inv = Matrix::solve(tBB)
  BtBB_inv = B %*% tBB_inv
  tBA = DelayedArray::t(B) %*% A
  BtBB_inv_tBA = BtBB_inv %*% tBA
  A - BtBB_inv_tBA
}

n = 100
p = 10000
Y = matrix(rnorm(n*p), nrow = n)
M = matrix(rnorm(n*10), nrow = n)
Y_hdf = as(Y, "HDF5Array")
M_hdf = as(M, "HDF5Array")
ctl = 1:100
profvis::profvis({
  A = Y
  B = M
  # tBB = DelayedArray::t(B) %*% B
  # tBB_inv = Matrix::solve(tBB)
  # tBA = DelayedArray::t(B) %*% A
  # BtBB_inv = B %*% tBB_inv
  # BtBB_inv_tBA = BtBB_inv %*% tBA
  # C = A - BtBB_inv_tBA
  C = my_residop(A, B)
  svdObj = BiocSingular::runRandomSVD(C, k = 20)
  fullalpha = t(svdObj$u[, seq_len(20), drop = FALSE]) %*% Y
  alpha <- fullalpha[seq_len(20), , drop = FALSE]
  ac <- alpha[, ctl, drop = FALSE]
  W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
  newY <- Y - W %*% alpha
  
  
  A = Y_hdf
  B = M
  # tBB = DelayedArray::t(B) %*% B
  # tBB_inv = Matrix::solve(tBB)
  # tBA = DelayedArray::t(B) %*% A
  # BtBB_inv = B %*% tBB_inv
  # BtBB_inv_tBA = BtBB_inv %*% tBA
  # C = A - BtBB_inv_tBA
  C = my_residop(A, B)
  svdObj = BiocSingular::runRandomSVD(C, k = 20)
  fullalpha = t(svdObj$u[, seq_len(20), drop = FALSE]) %*% Y
  alpha <- fullalpha[seq_len(20), , drop = FALSE]
  ac <- alpha[, ctl, drop = FALSE]
  W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
  newY <- Y - W %*% alpha
})