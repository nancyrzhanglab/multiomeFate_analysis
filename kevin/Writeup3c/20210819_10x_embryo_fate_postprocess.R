rm(list=ls())

load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed.RData")
load("../../../../out/kevin/Writeup3c/20210819_10x_embryo_result.RData")

# for now, a small percentage of cells haven't been matched yet
# let's remove those cells
zz <- which(apply(res$df_res, 1, function(x){is.na(x["init_state"]) & is.na(x["order_rec"])}))
stopifnot(length(zz) == 0)

#####################################

# let's do the naive fate calculations first
list_diagnos <- res$list_diagnos
n <- nrow(mat_x)
celltype <- as.numeric(as.character(celltype))
end_states <- list("16", "6", c("1","2","4"), "9")

# select endstate
prob_mat <- sapply(end_states, function(end_state){
  prob_vec <- rep(NA, n)
  prob_vec[which(celltype %in% end_state)] <- 1
  
  # set all other endstates to prob = 0
  prob_vec[which(celltype %in% end_states[which(!end_states %in% end_state)])] <- 0
  
  # iterate backward
  total_iter <- 1
  while(any(is.na(prob_vec)) && total_iter <= 10){
    print(sum(is.na(prob_vec)))
    
    for(iter in 1:length(list_diagnos)){
      for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
        idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
        idx_to <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
        
        ## ignore all "another nodes" that don't have a value assigned yet
        valid_to <- idx_to[which(!is.na(prob_vec[idx_to]))]
        valid_from <- idx_from[which(is.na(prob_vec[idx_from]))]
        
        if(length(valid_from) == 0 | length(valid_to) == 0) next()
        ## do multiplication: prob of a node going into another node * prob of that node going to end state
        prob_vec[valid_from] <- sum(prob_vec[valid_to])/length(idx_to)
        prob_vec[setdiff(idx_to, valid_to)] <- sum(prob_vec[valid_to])/length(idx_to)
      }
    }
    
    total_iter <- total_iter + 1
  }
  
  prob_vec
})

sum(is.na(prob_mat))
prob_mat[is.na(prob_mat)] <- 0
prob_mat <- t(apply(prob_mat, 1, function(x){
  if(all(x == 0)) return(x)
  x/sum(x)
}))

#####################################
tmp <- mbrain2[["RNA"]]@data[,cell_idx]
mbrain3 <- Seurat::CreateSeuratObject(counts = tmp)
mbrain3[["celltype"]] <- mbrain2@meta.data$celltype[cell_idx]
mbrain3[["seurat_cluster"]] <- mbrain2@meta.data$seurat_cluster[cell_idx]
mbrain3[["umap"]] <- Seurat::CreateDimReducObject(embedding = mbrain2[["umap"]]@cell.embeddings[cell_idx,], key = "UMAP_")
mbrain3[["umap.atac"]] <- Seurat::CreateDimReducObject(embedding = mbrain2[["umap.atac"]]@cell.embeddings[cell_idx,], key = "atacUMAP_")
mbrain3[["wnn.umap"]] <- Seurat::CreateDimReducObject(embedding = mbrain2[["wnn.umap"]]@cell.embeddings[cell_idx,], key = "wnnUMAP_")

rownames(prob_mat) <- rownames(mbrain3@meta.data)
colnames(prob_mat) <- paste0("prob_", 1:ncol(prob_mat))
mbrain3[["prob"]] <- Seurat::CreateDimReducObject(embedding = prob_mat, key = "prob_")

plot1 <- Seurat::FeaturePlot(mbrain3, features = "prob_1", reduction = "wnn.umap")
plot1 <- plot1 + ggplot2::ggtitle("To Oligodendrocyte")
plot2 <- Seurat::FeaturePlot(mbrain3, features = "prob_2", reduction = "wnn.umap")
plot2 <- plot2 + ggplot2::ggtitle("To Forebrain gabaergic")
plot3 <- Seurat::FeaturePlot(mbrain3, features = "prob_3", reduction = "wnn.umap")
plot3 <- plot3 + ggplot2::ggtitle("To one of the Corticals")
plot4 <- Seurat::FeaturePlot(mbrain3, features = "prob_4", reduction = "wnn.umap")
plot4 <- plot4 + ggplot2::ggtitle("To the other Cortical")

tmp2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_fateprob_naive.png",
                   tmp2, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

######################################

# let's deploy my fancier markov method
list_diagnos <- res$list_diagnos
max_to <- 10*max(sapply(list_diagnos, function(lis){
  max(sapply(lis$recruit$postprocess, function(x){
    length(x$to)
  }))
}))

n <- nrow(mat_x)
count_mat <- matrix(0, n, n)

# handle all the estimated transitions
for(iter in 1:length(list_diagnos)){
  for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
    from_idx <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
    to_idx <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
    
    count_mat[from_idx, to_idx] <- count_mat[from_idx, to_idx] + round(max_to/length(to_idx))
  }
}

# now handle all the terminal states
df_res <- res$df_res
tmp <- df_res$init_state[!is.na(df_res$init_state)]
val <- median(count_mat[count_mat != 0])
uniq_terminal <- unique(tmp[tmp > 0])
for(i in uniq_terminal){
  idx <- which(df_res$init_state == i)
  count_mat[idx, idx] <- val
}

# [[NOTE TO SELF:]] some columns are all-0's.
set.seed(10)
col_sum <- matrixStats::colSums2(count_mat)
idx <- which(col_sum == 0)
for(i in idx){
  count_mat[sample(1:n,1),i] <- 1
}


######################################
# run it piece-by-piece
celltype <- as.numeric(as.character(celltype))
end_states <- list("16", "6", c("1","2","4"), "9")

K <- 9
fixed_clustering = lapply(end_states, function(x){which(celltype %in% x)})
m = round(min(c(10*K, nrow(count_mat)/5)))
K0 = min(c(2*K, m-1))
delta = 0.9
num_restart = 10
max_tries = 500

# [[NOTE TO SELF: Let .mult_mat_vec take in a sparse matrix]]
col_sum <- sparseMatrixStats::colSums2(count_mat)
stopifnot(all(col_sum > 0))
tmp <- multiomeFate:::.mult_mat_vec(count_mat, 1/col_sum^(0.5))
svd_res <- multiomeFate:::.svd_truncated(tmp, K = K, K_full_rank = F, 
                          vec_mean = NULL, vec_sd = NULL)
d_mat <- multiomeFate:::.mult_vec_mat(1/svd_res$v[,1], svd_res$v[,-1])

set.seed(10)
# vertices <-  multiomeFate:::.vertex_hunting(d_mat, fixed_clustering = fixed_clustering,
#                             m = m, K0 = K0, 
#                             num_restart = num_restart, 
#                             max_tries = max_tries)
verbose <- T
mat <- d_mat
all_fixed_cell <- unlist(fixed_clustering)
has_fixed <- length(all_fixed_cell) > 0
if(verbose) print("Initializing cluster centers")
K <- ncol(mat)+1
if(has_fixed){
  mat2 <- mat[-all_fixed_cell,,drop = F]
}
res <- stats::kmeans(mat2, m, iter.max = 100, nstart = num_restart)
theta <- res$centers
if(has_fixed){
  extra_means <- t(sapply(fixed_clustering, function(x){
    colMeans(mat[x,,drop = F])
  }))
  theta <- rbind(theta, extra_means)
}

# [[note to self: Why are some entries of theta so large?
# is d_mat supposed to have such large entries? Why are some negative?]]
scaling_factor <- 100/max(abs(theta))
if(max(abs(theta)) > 100) theta <- theta*scaling_factor
dist_mat <- as.matrix(stats::dist(theta, method = "euclidean"))

if(verbose) print("Greedy selecting")
if(has_fixed) {
  len_fixed <- length(fixed_clustering)
  fixed_idx <- (nrow(theta)-len_fixed+1):nrow(theta)
  idx <- fixed_idx
} else {
  len_fixed <- 0
  idx <- unique(as.numeric(which(dist_mat == max(dist_mat), arr.ind = T)))
}
#length(idx)-len_fixed is the number of "free vertices"
if(K0 > length(idx)-len_fixed){
  while(length(idx)-len_fixed <= K0) {
    remaining_idx <- c(1:nrow(theta))[-idx]
    tmp <- remaining_idx[which.max(apply(dist_mat[remaining_idx,idx,drop = F], 1, mean))]
    idx <- c(idx, tmp)
  }
}

if(verbose) print("Finding best vertices")
combn_mat <- utils::combn(1:K0, K-len_fixed)
if(has_fixed) {
  for(i in fixed_idx) combn_mat <- rbind(combn_mat, i)
}
if(max_tries < ncol(combn_mat)){
  combn_mat <- combn_mat[,sample(1:ncol(combn_mat), max_tries)]
}
max_values <- rep(0, ncol(combn_mat))
for (i in 1:ncol(combn_mat)){
  if(verbose) print(i)
  for (j in setdiff(1:m, combn_mat[,i])){
    max_values[i] <- max(multiomeFate:::.simplex_dist(theta[j,], theta[combn_mat[,i],])$value, max_values[i])
  }
}

min_idx <- which.min(max_values)
vertices <-theta[combn_mat[,min_idx],]

vertices <- vertices/scaling_factor
w_mat <- t(sapply(1:nrow(d_mat), function(i){
  multiomeFate:::.weighting_optimization(d_mat[i,], vertices)
}))
w_mat[w_mat < 0] <- 0
w_mat <- t(apply(w_mat, 1, function(x){x/sum(x)}))
v_est <- multiomeFate:::.mult_vec_mat(svd_res$v[,1]*col_sum^(1/2), w_mat)
v_est <- apply(v_est, 2, function(x){x/sum(x)})

row_sum <- matrixStats::rowSums2(count_mat)
p_mat <- multiomeFate:::.mult_vec_mat(1/row_sum, count_mat)
u_est <- p_mat %*% v_est %*% solve(crossprod(v_est))
# [[note to self: why are some entries of U negative???]]
# [[not to self: Is it because I need to collapse some points to just one point 
#  to make it an "anchor point"?]]
u_est <- t(apply(u_est, 1, function(x){x[x< 0] <- 0; x/sum(x)}))

t(v_est) %*% u_est

# [[not to self: Visualize count_mat to see if I'm crazy]]

