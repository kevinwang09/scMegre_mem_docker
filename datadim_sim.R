library(scMerge)
library(purrr)
library(furrr)
library(SummarizedExperiment)
library(pryr)
library(HDF5Array)
library(Matrix)
library(tidyverse)
# library(DelayedArray)
# library(DelayedMatrixStats)


set.seed(1234)
num_cells = c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
plan(multisession, workers = length(num_cells))

listSim = purrr::map(.x = num_cells,
                     .f = ~ scMerge::ruvSimulate(m = .x, n = 20000, lambda = 1, sce = TRUE))

purrr::map_dbl(listSim, ~mean(assay(.x, "counts") == 0))
#####################################################
one_run = function(obj){
  mem = pryr::mem_change(
    scmerge <-  scMerge::scMerge(sce_combine = obj,
                                 ctl = paste0('gene', 1:500),
                                 cell_type = obj$cellTypes,
                                 ruvK = 10, assay_name = "scMerge_sim")
  )
  result = list(scmerge = scmerge,
                mem = mem)
  return(result)
}
##########################################
convert_sparse = function(obj){
  result = obj
  assay(result, "counts") = DelayedArray(as(assay(result, "counts"), "dgCMatrix"))
  assay(result, "logcounts") = DelayedArray(as(assay(result, "logcounts"), "dgCMatrix"))
  return(result)
}
##########################################
convert_hdf = function(obj){
  result = obj
  assay(result, "counts") = DelayedArray(as(assay(result, "counts"), "HDF5Array"))
  assay(result, "logcounts") = DelayedArray(as(assay(result, "logcounts"), "HDF5Array"))
  return(result)
}
##########################################
make_output = function(output){
  time = purrr::map(.x = output,
                    .f = ~ .x$scmerge@metadata$timeRuv) %>% 
    purrr::map_dbl(.f = ~ as.numeric(.x, units = "secs"))
  
  mem = purrr::map_dbl(.x = output, .f = ~ .x$mem)/1e9
  
  result = tibble::tibble(
    num_cells = num_cells,
    time = time,
    mem = mem
  )
  return(result)
}
##########################################
list_matrix_output = furrr::future_map(
  .x = listSim,
  .f = ~ one_run(.x),
  .progress = TRUE)

matrix_df = make_output(list_matrix_output) %>% 
  dplyr::mutate(type = "matrix")

##########################################
listSim_sparse = listSim %>% purrr::map(convert_sparse)
list_sparse_output = furrr::future_map(
  .x = listSim_sparse,
  .f = ~ one_run(.x),
  .progress = TRUE)

sparse_df = make_output(list_sparse_output) %>%
  dplyr::mutate(type = "sparse")
##########################################
listSim_hdf = listSim %>% purrr::map(convert_hdf)
list_hdf_output = furrr::future_map(
  .x = listSim_hdf,
  .f = ~ one_run(.x),
  .progress = TRUE)

hdf_df = make_output(list_hdf_output) %>%
  dplyr::mutate(type = "hdf5")
##########################################
combine_df = dplyr::bind_rows(
  matrix_df, 
  sparse_df,
  hdf_df
)

combine_df %>%
  ggplot(aes(x = as.factor(num_cells),
             y = time,
             colour = type,
             group = type)) +
  geom_point(size = 2) +
  geom_line(size = 1.2) +
  # scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Number of cells",
       y = "Time (s)",
       title = "Computational time of scMerge",
       subtitle = "Number of genes fixed at 20,000") +
  theme_bw(14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())


combine_df %>%
  ggplot(aes(x = as.factor(num_cells),
             y = mem,
             colour = type,
             group = type)) +
  geom_point(size = 2) +
  geom_line(size = 1.2) +
  # scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Number of cells",
       y = "Memory (GB)",
       title = "Computational memory usage of scMerge",
       subtitle = "Number of genes fixed at 20,000") +
  theme_bw(14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())