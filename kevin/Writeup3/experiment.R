set.seed(10)
p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
timepoints <- 20
mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                 list(mat_traj), verbose = F)
dat <- generate_data(obj_next, number_runs = 10, sample_perc = 0.9, verbose = F)

vec_start <- which(dat$df_info$time <= 0.1)
list_end <- list(which(dat$df_info$time >= 0.9))

set.seed(10)
res <- chromatin_potential_prepare(dat$obs_x, dat$obs_y, df$df_x, df$df_y,
                                   vec_start, list_end)
#############

mat_x = dat$obs_x
mat_y = dat$obs_y
df_x = df$df_x
df_y = df$df_y
dim_method = "pca"
nn_method = "annoy"
form_method = "literal"
est_method = "glmnet"
cand_method = "nn_any"
rec_method = "distant_cor"
options = list()

stopifnot(nrow(mat_x) == nrow(mat_y), ncol(mat_x) == nrow(df_x), 
          ncol(mat_y) == nrow(df_y), is.list(options))
stopifnot(all(mat_x >= 0), all(mat_y >= 0))
n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y); cell_name <- rownames(mat_x)

# check all the options
full_options <- .chrom_options(dim_method, nn_method,
                               form_method, est_method, 
                               cand_method, rec_method, 
                               options)
dim_options <- full_options$dim_options; nn_options <- full_options$nn_options

# compute the dimension reduction
dim_reduc_obj <- vector("list", 0)
tmp <- dimension_reduction(mat_x, mode = "x", dim_options)
