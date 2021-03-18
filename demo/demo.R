rm(list=ls())
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

p1 <- 100; p2 <- 20; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
head(df$df_x); head(df$df_y)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
dim(mat_g)
mat_g[1:5,1:5]
image(.rotate(mat_g))
