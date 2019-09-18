# Apologies, I am new to DA and this question may not be well-formated. 
# 
# The following example shows a small function `my_residop` that takes in two arguments of different matrix types. As expected, the return is a DelayedMatrix with two seeds. My understanding is that I need to run `realize` by manually setting a `BACKEND` to remove the conflicting seeds. Is there an automatic way to do this based on the seed of `hdf_A`? 
#   
#   
# ```
# library(DelayedArray)
# library(HDF5Array)
# 
# set.seed(1)
# A = matrix(rnorm(4), 2, 2)
# hdf_A = DelayedArray(as(A, "HDF5Array"))
# B = matrix(rnorm(4), 2, 2)
# my_residop = function(A,B){return(A - B %*% solve(t(B) %*% B) %*% t(B) %*% A)}
# C = my_residop(hdf_A, B) ## Works, return a DelayedMatrix
# seed(C) ## Error due to multiple seeds
# showtree(C)
# 
# 
# realize(C, BACKEND = class(seed(hdf_A))) ## This fails, hence my question
# ```