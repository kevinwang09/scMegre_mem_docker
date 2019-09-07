library(scMerge)
library(purrr)
library(furrr)
library(SummarizedExperiment)
library(pryr)
library(HDF5Array)
library(Matrix)

plan(multisession, workers = 2)
set.seed(1234)
nCells = c(100, 200)

listSim = purrr::map(.x = nCells,
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
  assay(result, "counts") = as(assay(result, "counts"), "dgCMatrix")
  assay(result, "logcounts") = as(assay(result, "logcounts"), "dgCMatrix")
  return(result)
}
##########################################
convert_hd = function(obj){
  result = obj
  assay(result, "counts") = as(assay(result, "counts"), "HDF5Array")
  assay(result, "logcounts") = as(assay(result, "logcounts"), "HDF5Array")
  return(result)
}
##########################################
make_output = function(output){
  time = purrr::map(.x = output,
                    .f = ~ .x$scmerge@metadata$timeRuv) %>% 
    purrr::map_dbl(.f = ~ as.numeric(.x, units = "secs"))
  
  mem = purrr::map_dbl(.x = output, .f = ~ .x$mem)/1e9
  
  result = tibble::tibble(
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

# 
# 
# list_scMerge_fast = furrr::future_map(
#   .x = listSim,
#   .f = ~ scMerge::scMerge(sce_combine = .x,
#                           ctl = paste0('gene', 1:500),
#                           cell_type = .x$cellTypes,
#                           ruvK = 10, assay_name = "scMerge_sim", fast_svd = TRUE),
#   .progress = TRUE)
# 
# list_scMerge_fast_time = purrr::map(.x = list_scMerge_fast, 
#                                .f = ~ .x@metadata$timeRuv)
# save(list_scMerge_time, list_scMerge_fast_time, file = "~/Desktop/scMerge_time.RData")

# load("~/Desktop/scMerge_time.RData")
# library(tidyverse)
# scMerge_slow_time = tibble(
#   num_cells = c(100, 200, 500, 1000, 2000, 5000),
#   time = purrr::map_dbl(.x = list_scMerge_time, .f = ~ as.numeric(.x, units = "secs")), 
#   type = "default"
# )
# 
# scMerge_fast_time = tibble(
#   num_cells = c(100, 200, 500, 1000, 2000, 5000),
#   time = purrr::map_dbl(.x = list_scMerge_fast_time, .f = ~ as.numeric(.x, units = "secs")),
#   type = "fast"
# )
# 
# 
# timeplotdf = bind_rows(scMerge_slow_time,
#                        scMerge_fast_time)
# 
# 
# timeplotdf %>% 
#   dplyr::filter(num_cells >= 200) %>% 
#   ggplot(aes(x = as.factor(num_cells), 
#              y = time, 
#              colour = type, 
#              group = type)) +
#   geom_point(size = 2) +
#   geom_line(size = 1.2) +
#   scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(x = "Number of cells", 
#        y = "Time (s)",
#        title = "Computational time of scMerge", 
#        subtitle = "Number of genes fixed at 20,000") +
#   theme_bw(14) +
#   theme(legend.position = "bottom",
#         panel.grid.minor = element_blank())
# ggsave(filename = "scMerge_time.png", height = 6, width = 6)
