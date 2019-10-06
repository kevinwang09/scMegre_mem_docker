library(tidyverse)
load("standard-4_pbmc_profile_k20_96b9319.RData")
df1 = print(p1, depth = 3) %>% 
  dplyr::mutate(mat_type = "DelayedArray+HDF5Array",
                machine_type = "15GB mem")
df2 = print(p2, depth = 3) %>% 
  dplyr::mutate(mat_type = "matrix",
                machine_type = "15GB mem")

mem15_df = bind_rows(df1, df2) 


load("standard-8_pbmc_profile_k20_96b9319.RData")
df1 = print(p1, depth = 3) %>% 
  dplyr::mutate(mat_type = "DelayedArray+HDF5Array",
                machine_type = "30GB mem")
df2 = print(p2, depth = 3) %>% 
  dplyr::mutate(mat_type = "matrix",
                machine_type = "30GB mem")

mem30_df = bind_rows(df1, df2) 

rm(df1, df2)

df = bind_rows(mem15_df, mem30_df)

df = df %>% dplyr::mutate(
    src = trimws(src),
    function_label = purrr::map_chr(
      .x = src,
      .f = ~ ifelse(str_detect(.x, "Delayed"), str_split(.x, "/")[[1]], .x)))

mem_plot = df %>% 
  ggplot(aes(x = fct_reorder(function_label, alloc), y = alloc,
             fill = mat_type)) +
  geom_col(position = "dodge") +
  labs(x = "Function call (depth = 3)",
       y = "Memory allocated (Mb)",
       caption = "The pbmc data is 8734 rows and 6791 cols, 3.91MB as Delayed+HDF5Array and 713MB as matrix") +
  facet_wrap(~machine_type) +
  coord_flip() +
  theme_bw(18) +
  theme(legend.position = "bottom")


time_plot = df %>% 
  ggplot(aes(x = fct_reorder(function_label, time), y = time,
             fill = mat_type)) +
  geom_col(position = "dodge") +
  labs(x = "Function call (depth = 3)",
       y = "Time (s)") +
  facet_wrap(~machine_type) +
  coord_flip() +
  theme_bw(18) +
  theme(legend.position = "bottom")



ggsave(filename = "./figures/mem_plot.pdf", plot = mem_plot,
       width = 12, height = 6)
ggsave(filename = "./figures/time_plot.pdf", plot = time_plot,
       width = 12, height = 6)

