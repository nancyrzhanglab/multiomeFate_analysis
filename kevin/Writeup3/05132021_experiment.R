rm(list=ls())
source("../multiomeFate_analysis/kevin/Writeup3/nancy_simulation_05042021.R")

set.seed(10)
dat <- generate_nancy_data()
names(dat)
mat_x <- dat$X; mat_y <- dat$Y
mat_y <- mat_y - min(mat_y)

mat_x <- t(mat_x); mat_y <- t(mat_y)
colnames(mat_x) <- paste0("p", 1:ncol(mat_x))
colnames(mat_y) <- paste0("g", 1:ncol(mat_y))
rownames(mat_x) <- paste0("n", 1:nrow(mat_x))
rownames(mat_y) <- paste0("n", 1:nrow(mat_y))

################

n <- nrow(mat_x); p1 <- ncol(mat_x); p2 <- ncol(mat_y)
gene_loc <- 100*(1:p2)
peak_loc <- unlist(lapply(gene_loc, function(x){x+(1:dat$nRegionsPerGene)}))
df_x <- data.frame(name = paste0("p",1:p1), location = peak_loc)
df_y <- data.frame(name = paste0("g",1:p2), location = gene_loc, baseline = 0)
stopifnot(nrow(df_x) == ncol(mat_x), nrow(df_y) == ncol(mat_y))

vec_start <- order(dat$trueTime, decreasing = F)[1:50]
list_end <- lapply(1:2, function(x){
  idx <- which(dat$trueBranch == x)
  idx[order(dat$trueTime[idx], decreasing = T)[1:50]]
})

# zz <- svd(mat_x); plot(zz$d)
rank_x <- 5
# zz <- svd(mat_y); plot(zz$d)
rank_y <- 3
########

set.seed(10)
prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        options = list(nn_nn = 20, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 10))
prep_obj$options
table(prep_obj$df_res$init_state)
table(prep_obj$df_res$order_rec)

set.seed(10)
res <- chromatin_potential(prep_obj, verbose = T)






