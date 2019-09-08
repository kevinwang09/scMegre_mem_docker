library(scMerge)
library(SummarizedExperiment)
data('example_sce', package = 'scMerge')
data('segList_ensemblGeneID', package = 'scMerge')
exprs_mat = counts(example_sce)
batch = example_sce$batch
this_batch_list = "batch2"
batch == this_batch_list 
marker = sample(rownames(example_sce))
t(exprs_mat[marker, batch == this_batch_list])
t(exprs_sp[marker, batch == this_batch_list])

exprs_sp = as(exprs_mat, "dgCMatrix")