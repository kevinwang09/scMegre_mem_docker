library(SummarizedExperiment)
library(scMerge)
library(Matrix)
###################### Identical to scMerge::scMerge vignette #####################
data('example_sce', package = 'scMerge')
data('segList_ensemblGeneID', package = 'scMerge')
profvis::profvis({
  sce_mESC <- scMerge(
    sce_combine = example_sce,
    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
    kmeansK = c(3, 3),
    assay_name = 'scMerge')
})
scater::plotPCA(sce_mESC, colour_by = 'cellTypes', shape = 'batch',
                run_args = list(exprs_values = 'logcounts'))
scater::plotPCA(sce_mESC, colour_by = 'cellTypes', shape = 'batch',
                run_args = list(exprs_values = 'scMerge'))


#######################################
example_sp = example_sce
## Assuming that the counts and logcounts slots are dgCMatrix
assay(example_sp, "counts") = as(assay(example_sce, "counts"), "dgCMatrix")
assay(example_sp, "logcounts") = as(assay(example_sce, "logcounts"), "dgCMatrix")

profvis::profvis({
sce_mESC_sp <- scMerge(
  sce_combine = example_sp,
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = 'scMerge')
})

example_hd = example_sce
## Assuming that the counts and logcounts slots are HDF5Array
assay(example_hd, "counts") = as(assay(example_sce, "counts"), "HDF5Array")
assay(example_hd, "logcounts") = as(assay(example_sce, "logcounts"), "HDF5Array")

profvis::profvis({
  sce_mESC_hd <- scMerge(
    sce_combine = example_hd,
    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
    kmeansK = c(3, 3),
    assay_name = 'scMerge')
})