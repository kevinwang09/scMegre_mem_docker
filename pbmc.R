library(TENxPBMCData)
library(SummarizedExperiment)
library(scMerge)
library(scater)
# BiocManager::install("M3Drop", version = "3.10")
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
scMerge::scMerge(
  sce_combine = pbmc_combine,
  ctl = segList_ensemblGeneID$human$human_scSEG,
  kmeansK = c(10, 10),
  assay_name = "scMerge_unsupervised", 
  verbose = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(workers = 6),
  BSPARAM = BiocSingular::IrlbaParam(fold = 5))
