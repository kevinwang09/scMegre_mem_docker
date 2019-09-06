# source("~/setup.R")
library(rhdf5)
library(HDF5Array)
download.file(url = "https://storage.googleapis.com/scp_data/se.rds",
              destfile = "su_yang_data/se.rds")
download.file(url = "https://storage.googleapis.com/scp_data/assays.h5",
              destfile = "su_yang_data/assays.h5")

## top -p 920 -b -n40 > top.txt

# su = readRDS(paste0(datapath, "sce_GSE87795.rds"))
# yang = readRDS(paste0(datapath,"sce_GSE90047.rds"))
# 
# su <- readRDS("~/Dropbox (Sydney Uni)/Single Cell Reserach/SingleCellPlus/data/sce_GSE87795.rds")
# yang <- readRDS("~/Dropbox (Sydney Uni)/Single Cell Reserach/SingleCellPlus/data/sce_GSE90047.rds")
# 
# sce_list = list(
#   su = su,
#   yang = yang
# )
# 
# sce_list
# 
# sce_combine = scMerge::sce_cbind(sce_list = sce_list,
#                                  method = "union",
#                                  colData_names = c("cellTypes", "stage"),
#                                  batch_names = c("Su", "Yang"))
# 
# sce_combine
# 
# sce_combine = sce_combine[rowSums(SingleCellExperiment::counts(sce_combine)) != 0,
#                           colSums(SingleCellExperiment::counts(sce_combine)) != 0]
# 
# 
# saveHDF5SummarizedExperiment(sce_combine, "./su_yang_data")

sce_combine_da = loadHDF5SummarizedExperiment("./su_yang_data")
sce_combine_WASTE = realize(sce_combine_da)


library(scMerge)
data("segList_ensemblGeneID", package = "scMerge")
(t1 = Sys.time())
scMerge_supervised = scMerge(
  sce_combine = sce_combine_da,
  ctl = which(rownames(sce_combine_da) %in% segList_ensemblGeneID$mouse$mouse_scSEG),
  cell_type = sce_combine_da$cellTypes,
  replicate_prop = 1,
  assay_name = "scMerge_supervised",
  verbose = TRUE)
t2 = Sys.time()
t2 - t1

sce_combine_real = realize(sce_combine_da)
(t1 = Sys.time())
scMerge_supervised = scMerge(
  sce_combine = sce_combine_real,
  ctl = which(rownames(sce_combine_real) %in% segList_ensemblGeneID$mouse$mouse_scSEG),
  cell_type = sce_combine_real$cellTypes,
  replicate_prop = 1,
  assay_name = "scMerge_supervised",
  verbose = TRUE)
t2 = Sys.time()
t2 - t1