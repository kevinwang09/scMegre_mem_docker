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
# num_cells = c(100, 200, 500, 1000, 2000, 5000, 10000)
num_cells = c(50, 100, 200, 500)
plan(multisession, workers = length(num_cells))
# options(future.globals.maxSize = 50000*1024^2)

listSim = furrr::future_map(
  .x = num_cells,
  .f = ~ scMerge::ruvSimulate(m = .x, n = 20000, lambda = 0.1, sce = TRUE),
  .progress = TRUE)

purrr::map_dbl(listSim, ~mean(assay(.x, "counts") == 0))
#####################################################
one_run = function(obj){
  mem_comp = pryr::mem_change(
    scmerge <-  scMerge::scMerge(sce_combine = obj,
                                 ctl = paste0('gene', 1:500),
                                 cell_type = obj$cellTypes,
                                 ruvK = 10, assay_name = "scMerge_sim")
  )
  mem_obj = pryr::object_size(scmerge)
  mem_total = mem_comp + mem_obj
  result = list(scmerge = scmerge,
                mem_comp = mem_comp,
                mem_obj = mem_obj,
                mem_total = mem_total)
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
  
  mem_total = purrr::map_dbl(.x = output, .f = ~ .x$mem_total)/1e9
  mem_comp = purrr::map_dbl(.x = output, .f = ~ .x$mem_comp)/1e9
  mem_obj = purrr::map_dbl(.x = output, .f = ~ .x$mem_obj)/1e9
  
  result = tibble::tibble(
    num_cells,
    time,
    mem_comp,
    mem_obj,
    mem_total
  )
  return(result)
}
##########################################
# list_matrix_output = furrr::future_map(
#   .x = listSim[1:2],
#   .f = ~ one_run(.x),
#   .progress = TRUE)

list_matrix_output = mclapply(listSim, one_run, mc.cores = length(num_cells))

matrix_df = make_output(list_matrix_output) %>% 
  dplyr::mutate(type = "matrix")

# mem_change(rm(list_matrix_output))
##########################################
listSim_sparse = listSim %>% purrr::map(convert_sparse)
list_sparse_output = furrr::future_map(
  .x = listSim_sparse,
  .f = ~ one_run(.x),
  .progress = TRUE)

# list_sparse_output = mclapply(listSim_sparse, one_run, mc.cores = length(num_cells))
# one_run(listSim_sparse[[3]])

sparse_df = make_output(list_sparse_output) %>%
  dplyr::mutate(type = "sparse")

# mem_change(rm(listSim_sparse))
##########################################
listSim_hdf = listSim %>% purrr::map(convert_hdf)
list_hdf_output = furrr::future_map(
  .x = listSim_hdf,
  .f = ~ one_run(.x),
  .progress = TRUE)

# list_hdf_output = mclapply(listSim_hdf, one_run, mc.cores = length(num_cells))

hdf_df = make_output(list_hdf_output) %>%
  dplyr::mutate(type = "hdf5")

# mem_change(rm(list_hdf_output))
##########################################
combine_df = dplyr::bind_rows(
  matrix_df, 
  sparse_df,
  hdf_df) %>% 
  tidyr::gather(key = mem_type, 
                value = mem_value, 
                dplyr::contains("mem_"))

combine_df %>%
  ggplot(aes(x = mem_value,
             y = time,
             colour = type,
             group = type)) +
  geom_point(size = 2) +
  geom_line(size = 1.2) +
  geom_text(aes(label = num_cells), nudge_y = 1) +
  facet_wrap(~ mem_type, scales = "free_x") +
  # scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Memory (GB)",
       y = "Time (s)",
       title = "Computational time and memory of scMerge",
       caption = "Total memory = object memory + computation memory",
       subtitle = "Number of genes fixed at 20,000") +
  theme_bw(14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

# combine_df %>%
#   ggplot(aes(x = as.factor(num_cells),
#              y = time,
#              colour = type,
#              group = type)) +
#   geom_point(size = 2) +
#   geom_line(size = 1.2) +
#   # scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(x = "Number of cells",
#        y = "Time (s)",
#        title = "Computational time of scMerge",
#        subtitle = "Number of genes fixed at 20,000") +
#   theme_bw(14) +
#   theme(legend.position = "bottom",
#         panel.grid.minor = element_blank())
# 
# 
# combine_df %>%
#   ggplot(aes(x = as.factor(num_cells),
#              y = mem,
#              colour = type,
#              group = type)) +
#   geom_point(size = 2) +
#   geom_line(size = 1.2) +
#   # scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(x = "Number of cells",
#        y = "Memory (GB)",
#        title = "Computational memory usage of scMerge",
#        subtitle = "Number of genes fixed at 20,000") +
#   theme_bw(14) +
#   theme(legend.position = "bottom",
#         panel.grid.minor = element_blank())