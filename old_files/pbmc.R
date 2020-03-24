library(TENxPBMCData)
library(SummarizedExperiment)
library(scMerge)
library(scater)
pbmc3k <- TENxPBMCData('pbmc3k')
pbmc4k <- TENxPBMCData('pbmc4k')

unfiltered <- pbmc3k
is.mito <- grep("MT", rowData(pbmc3k)$Symbol_TENx)
stats <- perCellQCMetrics(pbmc3k, subsets=list(Mito=is.mito))
high.mito <- isOutlier(stats$subsets_Mito_percent, nmads=3, type="higher")
pbmc3k <- pbmc3k[,!high.mito]
pbmc3k <- logNormCounts(pbmc3k)

unfiltered <- pbmc4k
is.mito <- grep("MT", rowData(pbmc4k)$Symbol_TENx)
stats <- perCellQCMetrics(pbmc4k, subsets=list(Mito=is.mito))
high.mito <- isOutlier(stats$subsets_Mito_percent, nmads=3, type="higher")
pbmc4k <- pbmc4k[,!high.mito]
pbmc4k <- logNormCounts(pbmc4k)

pbmc_combine = scMerge::sce_cbind(list(pbmc3k = pbmc3k, pbmc4k = pbmc4k),
                                  method = "intersect")
dim(pbmc_combine)
data("segList_ensemblGeneID")
colnames(pbmc_combine) = paste0("gene", seq_len(ncol(pbmc_combine)))
##################################
p1 = lineprof::lineprof(code = {
# system.time(
  obj1 <- scMerge::scMerge(
    sce_combine = pbmc_combine,
    ctl = segList_ensemblGeneID$human$human_scSEG,
    kmeansK = c(10, 10),
    assay_name = "scMerge_unsupervised2", 
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    svd_k = 20,
    BSPARAM = BiocSingular::RandomParam(fold = Inf),
    BACKEND = "HDF5Array")
# )
}, interval = 0.5)


SummarizedExperiment::assay(pbmc_combine, "counts") = as.matrix(SummarizedExperiment::assay(pbmc_combine, "counts"))
SummarizedExperiment::assay(pbmc_combine, "logcounts") = as.matrix(SummarizedExperiment::assay(pbmc_combine, "logcounts"))

p2 = lineprof::lineprof(code = {
# system.time(
  obj2 <- scMerge::scMerge(
    sce_combine = pbmc_combine,
    ctl = segList_ensemblGeneID$human$human_scSEG,
    kmeansK = c(10, 10),
    assay_name = "scMerge_unsupervised2", 
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    svd_k = 20,
    BSPARAM = BiocSingular::RandomParam(fold = Inf))
# )
}, interval = 0.5)


obj1@metadata$timeRuv + obj1@metadata$timeReplicates
sum(print(p1)$time)
sum(print(p1)$alloc)


obj2@metadata$timeRuv + obj2@metadata$timeReplicates
sum(print(p2)$time)
sum(print(p2)$alloc)


scater::plotPCA(obj1, colour_by = 'batch',
                run_args = list(exprs_values = 'logcounts'))

scater::plotPCA(obj1, colour_by = 'batch',
                run_args = list(exprs_values = 'scMerge_unsupervised1'))

############################

scater::plotPCA(obj2, colour_by = 'batch',
                run_args = list(exprs_values = 'logcounts'))

scater::plotPCA(obj2, colour_by = 'batch',
                run_args = list(exprs_values = 'scMerge_unsupervised2'))


# save(p1, p2, file = "standard-4_pbmc_profile_k20_bc1e1bc.RData")