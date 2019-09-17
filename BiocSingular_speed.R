library(DelayedArray)
library(tibble)
library(dplyr)
library(BiocSingular)
library(ggplot2)
library(tidyr)

set.seed(1234)

p = 20000
n = c(500, 1000, 2000, 5000)
k = c(10)

BSPARAM = c(ExactParam(fold = Inf),
            IrlbaParam(fold = Inf), 
            RandomParam(fold = Inf))

oneset_grid = expand.grid(n = n, p = p, k = k) %>% 
  as_tibble() %>% 
  dplyr::mutate(
    A_mat = purrr::map2(
      .x = n,
      .y = p,
      .f = ~ matrix(rnorm(.x*.y), nrow = .x)),
    A_sp = purrr::map(A_mat, ~ as(.x, "dgeMatrix")),
    A_mat_da = purrr::map(A_mat, DelayedArray),
    A_sp_da = purrr::map(A_sp, DelayedArray))

(grid = rep(list(oneset_grid), length(BSPARAM)) %>% 
    purrr::map2_dfr(
      .x = ., .y = BSPARAM, 
      .f = ~ dplyr::mutate(.x, bsparam = list(.y))) %>% 
    tidyr::gather(key = mat_key, 
                  value = mat_value, 
                  dplyr::contains("A_"))
)

comp_result = mcmapply(
  FUN = function(A, k, BSPARAM){
    t1 = Sys.time()
    svd_obj = BiocSingular::runSVD(x = A, k = k, BSPARAM = BSPARAM)
    t2 = Sys.time()
    result = list(
      time = as.numeric(t2 - t1, units = "secs"),
      svd_obj = svd_obj
    )
    return(result)
  },
  A = grid$mat_value,
  k = grid$k, 
  BSPARAM = grid$bsparam,
  SIMPLIFY = FALSE,
  mc.cores = 10)

# comp_result

comp_tibble = grid %>% 
  dplyr::mutate(time = purrr::map_dbl(comp_result, "time"),
                bsparam_class = purrr::map_chr(bsparam, class),
                k = factor(k))

comp_tibble %>% 
  ggplot(aes(x = factor(n), y = time,
             colour = bsparam_class, group = interaction(k, bsparam_class))) + 
  geom_point() +
  geom_line() +
  scale_y_log10(breaks = c(10, 100, 200, 500, 700, 800, 1000, 1200, 1500, 2000, 2200, 3000)) +
  facet_grid(~mat_key)

# ggsave(filename = "BiocSingular_speed_fold2.pdf")
  ## What about accuracy of the decomp?