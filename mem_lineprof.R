library(scMerge)
library(lineprof)
set.seed(12345)
L = ruvSimulate(m = 6000, n = 10000, nc = 100, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

tmp = lineprof({
  L_result1 <- scMerge(
    sce_combine = L,
    ctl = paste0("gene",1:100),
    kmeansK = c(3, 3),
    assay_name = 'scMerge',
    svd_k = 50, 
    verbose = TRUE,
    BSPARAM = BiocSingular::RandomParam())
}, interval = 0.5)

View(print(tmp))
