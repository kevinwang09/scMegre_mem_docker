library(DelayedArray)
library(DelayedMatrixStats)
library(parallel)
library(dplyr)
library(purrr)
library(ggplot2)
library(tibble)
# library(scMerge)
set.seed(1234)
ncores = 30

p = c(20000)
n = c(100, 200, 500, 1000, 2000, 5000, 10000)

# p = c(20000)
# n = c(500, 1000, 2000)

(grid = expand.grid(n = n, p = p) %>% 
    as_tibble() %>% 
    dplyr::mutate(
      A = purrr::map2(
        .x = n,
        .y = p,
        .f = ~ matrix(rpois(.x*.y, lambda = 0.1) %>% as.numeric(), nrow = .x)),
      B = purrr::map2(
        .x = n,
        .y = p,
        .f = ~ matrix(rpois(.x*.y, lambda = 0.1) %>% as.numeric(), nrow = .y)),
      A_da = purrr::map(A, DelayedArray),
      B_da = purrr::map(B, DelayedArray),
      A_sp = purrr::map(A, ~ as(.x, "dgCMatrix")),
      B_sp = purrr::map(B, ~ as(.x, "dgCMatrix")),
      A_da_sp = purrr::map(A_sp, DelayedArray),
      B_da_sp = purrr::map(B_sp, DelayedArray),
    )
)


A_times_B = function(A, B){
  t1 = Sys.time()
  mem_comp = pryr::mem_change(A %*% B)
  t2 = Sys.time()
  
  result = list(time = as.numeric(t2 - t1, units = "secs"),
                mem_comp = mem_comp)
}

A_times_B_cpp = function(A, B){
  t1 = Sys.time()
  mem_comp = pryr::mem_change(scMerge::eigenMatMult(A, B))
  t2 = Sys.time()

  result = list(time = as.numeric(t2 - t1, units = "secs"),
                mem_comp = mem_comp)
}

# mat_result = furrr::future_map2(.x = grid$A, .y = grid$B, .f = A_times_B)
mat_result = parallel::mcmapply(FUN = A_times_B, grid$A, grid$B, mc.cores = ncores, SIMPLIFY = FALSE)
cpp_result = parallel::mcmapply(FUN = A_times_B_cpp, grid$A, grid$B, mc.cores = ncores, SIMPLIFY = FALSE)
da_result = mapply(FUN = A_times_B, grid$A_da, grid$B_da, SIMPLIFY = FALSE)
sp_result = parallel::mcmapply(FUN = A_times_B, grid$A_sp, grid$B_sp, mc.cores = ncores, SIMPLIFY = FALSE)
da_sp_result = mapply(FUN = A_times_B, grid$A_da_sp, grid$B_da_sp, SIMPLIFY = FALSE)

comp_result = grid %>% 
  dplyr::mutate(
    mat_result,
    cpp_result,
    da_result,
    sp_result,
    da_sp_result) %>% 
  dplyr::select(n, p, dplyr::contains("result"))

comp_time_df = comp_result %>% 
  dplyr::mutate(mat_time = purrr::map_dbl(mat_result, "time"),
                cpp_time = purrr::map_dbl(cpp_result, "time"),
                da_time = purrr::map_dbl(da_result, "time"),
                sp_time = purrr::map_dbl(sp_result, "time"),
                da_sp_time = purrr::map_dbl(da_sp_result, "time")) %>% 
  tidyr::gather(key = time_key, 
                value = time_value,
                dplyr::contains("time"))


comp_time_df %>% 
  ggplot(aes(x = factor(n), y = time_value/60, 
             colour = time_key,
             group = time_key)) +
  geom_point() +
  geom_line() +
  geom_text(aes(label = signif(time_value/60, 2)), nudge_y = 0.05) +
  facet_wrap(~p) +
  # scale_x_log10() +
  scale_y_log10(breaks = c(1, 2, 3, 4, 5, 10, 15, 20 )) +
  labs(x = "n", 
       y = "time (s)") +
  theme_bw(18) +
  theme(legend.position = "bottom")

ggsave(filename = "mat_multiply_time.pdf", plot = last_plot(), height = 8, width = 8)

####################################


# mem_time_df = comp_result %>% 
#   dplyr::select(n, p, dplyr::contains("result")) %>% 
#   dplyr::mutate(mat_mem = purrr::map_dbl(mat_result, "mem_comp"),
#                 da_mem = purrr::map_dbl(da_result, "mem_comp"),
#                 sp_mem = purrr::map_dbl(sp_result, "mem_comp"),
#                 da_sp_mem = purrr::map_dbl(da_sp_result, "mem_comp")) %>% 
#   tidyr::gather(key = mem_key, 
#                 value = mem_value,
#                 dplyr::contains("mem"))
# 
# 
# mem_time_df %>% 
#   ggplot(aes(x = n, y = mem_value, 
#              colour = mem_key)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~p) +
#   scale_x_log10() +
#   labs(x = "n", 
#        y = "Memory (GB)")



# profvis::profvis({
#   n = 10000
#   p = 500
#   A = matrix(rpois(n*p, lambda = 0.1), nrow = n)
#   B = matrix(rpois(n*p, lambda = 0.1), nrow = n)
#   
#   A_da = DelayedArray(A)
#   B_da = DelayedArray(B)
#   
#   A_sp = as(A, "dgCMatrix")
#   B_sp = as(B, "dgCMatrix")
#   
#   A %*% t(B)
#   A_da %*% t(B_da)
#   A_sp %*% t(B_sp)
# })
