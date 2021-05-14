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

######################################
plot(res$df_res$order_rec, dat$trueTime)

set.seed(10)
plot_umap.chromatin_potential(res, multiple_to = "umap_avg", col_vec = dat$trueBranch+1, k =3)

set.seed(10)
X.umap=umap(mat_x, n_neighbors=30)
plot(X.umap$layout[,1], X.umap$layout[,2], asp = T, col = dat$trueBranch+1, pch = 16)
for(i in 1:nrow(mat_x)){
  flip <- rbinom(1, 1, 0.3)
  if(flip == 1){
    vec_from <- X.umap$layout[i,]
    idx_to <- res$ht_neighbor[[as.character(i)]]
    vec_to <- colMeans(X.umap$layout[idx_to,])
    
    graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                     x1 = vec_to[1], y1 = vec_to[2], length = 0.05)
  }
}

time_start <- dat$trueTime
time_end <- sapply(1:nrow(mat_x), function(i){
  idx_to <- res$ht_neighbor[[as.character(i)]]
  mean(dat$trueTime[idx_to])
})
time_start <- time_start[order(res$df_res$order_rec)]
time_end <- time_end[order(res$df_res$order_rec)]

plot(NA, xlim = c(0,nrow(mat_x)), ylim = range(dat$trueTime))
for(i in 1:length(time_start)){
  if(time_start[i] <= time_end[i]) col = "green" else col = "red"
  graphics::arrows(x0 = i, y0 = time_start[i],
                   x1 = i, y1 = time_end[i], length = 0.05, col = col)
}


