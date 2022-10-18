rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

all_data[["spliced"]] <- NULL
all_data[["unspliced"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6/Writeup6_timeAll_simplified.RData")

###############

load("../../../../out/kevin/Writeup6/Writeup6_timeAll_simplified.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

lin_mat <- all_data[["Lineage"]]@counts
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]

###################

lin_mat_t <- Matrix::t(lin_mat)
p <- nrow(lin_mat)
num_cells <- sapply(1:p, function(j){
  length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
})
lin_mat_t <- lin_mat_t[,which(num_cells >= 100)]
cor_mat <- stats::cor(as.matrix(lin_mat_t))

diag(cor_mat) <- 0
off_diag_vec <- cor_mat[lower.tri(cor_mat, diag = F)]
tmp <- hist(off_diag_vec, plot = F)
tmp$counts <- log10(tmp$counts+1)
png("../../../../out/figures/Writeup6/Writeup6_lineage-doubleinfection-correlation_histogram.png",
    height = 1500, width = 2000, units = "px", res = 300)
plot(tmp, xlab = "Correlation between two lineage counts (across cells)",
     main = paste0("Only lineages with 100 or more cells\n(Among all ", ncol(lin_mat_t), " lineages)"),
     ylab = "Log10 Frequency", col = "gray")
graphics.off()

###################

res_clusters <- barcode_clustering(lin_mat)
table(sapply(res_clusters$lineage_clusters, length))

cluster_size <- sapply(res_clusters$lineage_clusters, length)
for(i in 1:length(cluster_size)){
  if(cluster_size[i] >= 3){
    tmp <- lin_mat[which(rownames(lin_mat) %in% res_clusters$lineage_clusters[[i]]), ]
    tmp <- as.matrix(Matrix::t(tmp))
    print(paste0("On cluster ", i))
    print(round(cor(tmp), 2))
    print("====")
  }
}

##################

# combining lineages
lin_mat <- as.matrix(lin_mat)
lin_mat2 <- lin_mat
lin_mat <- barcode_combine(lin_mat = lin_mat,
                           lineage_clusters = res_clusters$lineage_clusters,
                           verbose = 1)

# check:
k <- which.max(sapply(res_clusters$lineage_clusters, length))
idx_prev <- which(rownames(lin_mat2) %in% res_clusters$lineage_clusters[[k]])
idx_new <- which(rownames(lin_mat) %in% res_clusters$lineage_clusters[[k]])
stopifno(length(idx_new) == 1)
zz <- colSums(lin_mat2[idx_prev,])
yy <- lin_mat[idx_new,]
stopifnot(sum(abs(zz-yy)) == 0)

lin_mat_t <- Matrix::t(lin_mat)
p <- nrow(lin_mat)
num_cells <- sapply(1:p, function(j){
  length(which(lin_mat[j,] != 0))
})
lin_mat_t <- lin_mat_t[,which(num_cells >= 100)]
if(ncol(lin_mat_t) <= 1000) cor_mat <- stats::cor(as.matrix(lin_mat_t))
off_diag_vec <- cor_mat[lower.tri(cor_mat, diag = F)]
round(quantile(off_diag_vec), 2)

#########################

barcoding_res <- barcoding_assignment(lin_mat = lin_mat,
                                      verbose = 1)
barcoding_res$Bhat[1:5,1:5]
table(!is.na(barcoding_res$beta1), barcoding_res$lineage_num_winner > 0)

maxBhat <- apply(barcoding_res$Bhat,2,max) #chosen lineage 
argmaxBhat <- apply(barcoding_res$Bhat, 2, which.max)
argmaxX <- apply(lin_mat, 2, which.max) #largest count
maxX <- apply(lin_mat,2,max)
n <- ncol(lin_mat)
is_unique_vec_x <- sapply(1:n, function(i){length(which(lin_mat[,i] == maxX[i])) == 1})
is_unique_vec_b <- sapply(1:n, function(i){length(which(abs(barcoding_res$Bhat[,i] - maxBhat[i]) <= 1e-3)) == 1})
table(is_unique_vec_b)
table(is_unique_vec_x)
disagrees <- which(argmaxX!=argmaxBhat & maxX>5& maxBhat>0.9 & is_unique_vec_x)
length(disagrees)
length(disagrees)/length(which(maxX>5& maxBhat>0.9 & is_unique_vec_x))
length(which(maxX>5& maxBhat>0.9 & is_unique_vec_x))

disagrees2 <- which(argmaxX!=argmaxBhat & maxX>5 & is_unique_vec_x)
length(disagrees2)

disagrees3 <- which(argmaxX!=argmaxBhat & maxX>0 & is_unique_vec_x)
length(disagrees3)
length(which(maxX>0 & is_unique_vec_x))

############################################

# print some stats
quantile(maxX)
length(which(maxX == 0))/length(maxX)
naive_idx <- which(all_data$dataset == "day0")
quantile(maxX[naive_idx])
length(which(maxX[naive_idx] == 0))/length(naive_idx)
length(which(maxX[naive_idx] == 1))/length(naive_idx)

difference_Xvec <- apply(lin_mat, 2, function(x){
  tmp <- sort(x, decreasing = T)[1:2]
  -diff(tmp)
})
quantile(difference_Xvec)
length(which(difference_Xvec > 1/3-1e-3))/length(difference_Xvec)
quantile(difference_Xvec[naive_idx])
length(which(difference_Xvec[naive_idx] > 1))/length(naive_idx)
length(which(difference_Xvec[naive_idx] > 2))/length(naive_idx)

difference_Bvec <- apply(barcoding_res$Bhat, 2, function(x){
  -diff(sort(x, decreasing = T)[1:2])
})
quantile(difference_Bvec)
length(which(difference_Bvec >= 0.1))/length(difference_Bvec)
length(which(difference_Bvec >= 0.1 & is_unique_vec_b))/length(difference_Bvec)
quantile(difference_Bvec[naive_idx])
length(which(difference_Bvec[naive_idx] >= 0.05))/length(naive_idx)
length(which(difference_Bvec[naive_idx] >= 0.2))/length(naive_idx)
length(which(difference_Bvec[naive_idx] >= 0.1 & is_unique_vec_b[naive_idx]))/length(naive_idx)
quantile(difference_Bvec[naive_idx[which(maxX[naive_idx] == 0)]])

table(is_unique_vec_b[naive_idx])

table(table(argmaxBhat))
table(table(argmaxBhat[is_unique_vec_b]))
table(table(argmaxBhat[naive_idx]))
table(table(argmaxBhat[intersect(naive_idx,which(is_unique_vec_b))]))
which(table(argmaxBhat[naive_idx]) == 3527)
length(which(argmaxBhat == 1))
head(which(argmaxBhat == 1))
quantile(maxX[which(argmaxBhat == 1)]) # oohh.....
quantile(barcoding_res$Bhat[,2])

table(is_unique_vec_b[naive_idx])
length(which(is_unique_vec_b[naive_idx]))/length(naive_idx)

zz <- which(!is_unique_vec_b)
sapply(zz[1:5], function(i){
  quantile(barcoding_res$Bhat[,i])
})
yy <- sapply(zz, function(i){
  sum(lin_mat[,i])
})
head(which(yy != 0))
zz2 <- zz[which(yy != 0)[1:5]]
sapply(zz2, function(i){
  print(quantile(barcoding_res$Bhat[,i]))
  print(quantile(lin_mat[,i]))
  print(length(which(lin_mat[,i] > 0)))
  print("===")
})
round(barcoding_res$beta0[zz2],8)
round(barcoding_res$beta1[zz2],8)

i <- naive_idx[which(!is_unique_vec_b[naive_idx])[1]]
quantile(lin_mat[,i])
all_data$dataset[i]
barcoding_res$beta0[i]; barcoding_res$beta1[i]
any_count <- colSums(lin_mat)>0
table(is_unique_vec_b, any_count)
table(is_unique_vec_b[naive_idx], any_count[naive_idx])
idx2 <- naive_idx[intersect(which(!is_unique_vec_b[naive_idx]), which(any_count[naive_idx]))]
i <- idx2[1]
length(which(lin_mat[,i] > 0))
barcoding_res$beta0[which(lin_mat[,i] > 0)]
barcoding_res$beta1[which(lin_mat[,i] > 0)]
barcoding_res$Bhat[which(lin_mat[,i] > 0),i]
quantile(barcoding_res$Bhat[,i])
beta0 <- barcoding_res$beta0
tol <- 1e-8
beta0_thresh <- pmax(beta0, max(quantile(beta0[beta0 > tol], 0.02)))
quantile(beta0[beta0 > tol])
beta0_thresh[which(lin_mat[,i] > 0)]

##

assignment_vec <- apply(barcoding_res$Bhat, 2, function(vec){
  ord_vec <- sort(vec, decreasing = T)[1:2]
  if(abs(diff(ord_vec)) >= 1e-3) return(which.max(vec)) else return(NA)
})
table(table(assignment_vec))
length(unique(assignment_vec))
length(setdiff(1:nrow(lin_mat), unique(assignment_vec)))
length(which(is.na(assignment_vec[naive_idx])))
length(unique(assignment_vec[naive_idx]))
length(setdiff(1:nrow(lin_mat), unique(assignment_vec[naive_idx])))
quantile(table(assignment_vec[naive_idx]))
quantile(table(assignment_vec))

assignment_vec_simple <- apply(lin_mat, 2, function(x){
  tmp <- sort(x, decreasing = T)[1:2]
  if(abs(diff(tmp)) >= 2) return(which.max(x)) else return(NA)
})
length(unique(assignment_vec_simple))
length(setdiff(1:nrow(lin_mat), unique(assignment_vec_simple)))
length(which(is.na(assignment_vec_simple[naive_idx])))
length(unique(assignment_vec_simple[naive_idx]))
length(setdiff(1:nrow(lin_mat), unique(assignment_vec_simple[naive_idx])))
quantile(table(assignment_vec_simple[naive_idx]))

################################

# for comparison

barcoding_res2 <- barcoding_assignment(lin_mat = lin_mat2,
                                       verbose = 1)
argmaxX2 <- apply(lin_mat2, 2, which.max) #largest count
maxX2 <- apply(lin_mat2,2,max)
maxBhat2 <- apply(barcoding_res2$Bhat,2,max) #chosen lineage 
argmaxBhat2 <- apply(barcoding_res2$Bhat, 2, which.max)
is_unique_vec_x2 <- sapply(1:n, function(i){length(which(lin_mat2[,i] == maxX2[i])) == 1})
is_unique_vec_b2 <- sapply(1:n, function(i){length(which(abs(barcoding_res2$Bhat[,i] - maxBhat2[i]) <= 1e-3)) == 1})

disagrees <- which(argmaxX2!=argmaxBhat2 & maxX2>5 & maxBhat2>0.9 & is_unique_vec_x2)
length(disagrees)
length(disagrees)/length(which(maxX2>5& maxBhat2>0.9 & is_unique_vec_x2))

# let's take a look at the first one
cell_idx <- disagrees[1]
zz <- sort(lin_mat2[,cell_idx],decreasing = T)[1:10]
round(barcoding_res2$Bhat[names(zz),cell_idx],2)

disagrees3 <- which(argmaxX2!=argmaxBhat2 & maxX2>0 & is_unique_vec_x2)
length(disagrees3)
length(which(maxX2>0 & is_unique_vec_x2))

difference_Bvec2 <- apply(barcoding_res2$Bhat, 2, function(x){
  -diff(sort(x, decreasing = T)[1:2])
})
length(which(difference_Bvec2[naive_idx] >= 0.05))/length(naive_idx)
length(which(difference_Bvec2[naive_idx] >= 0.1))/length(naive_idx)
table(is_unique_vec_b2[naive_idx])


assignment_vec2 <- apply(barcoding_res2$Bhat, 2, function(vec){
  ord_vec <- sort(vec, decreasing = T)[1:2]
  if(abs(diff(ord_vec)) >= 1e-3) return(which.max(vec)) else return(NA)
})
table(table(assignment_vec2))
length(unique(assignment_vec2))
length(setdiff(1:nrow(lin_mat2), unique(assignment_vec2)))
length(unique(assignment_vec2[naive_idx]))
length(setdiff(1:nrow(lin_mat2), unique(assignment_vec2[naive_idx])))
quantile(table(assignment_vec2[naive_idx]))

k <- which.max(sapply(res_clusters$lineage_clusters, length))
idx_prev <- which(rownames(lin_mat2) %in% res_clusters$lineage_clusters[[k]])
idx_new <- which(rownames(lin_mat) %in% res_clusters$lineage_clusters[[k]])
table(assignment_vec2[assignment_vec2 %in% idx_prev])
table(assignment_vec[assignment_vec %in% idx_new])

#################################

dataset_vec <- sort(unique(all_data$dataset))
lin_total_count <- sapply(1:nrow(lin_mat), function(j){
  if(j %% floor(nrow(lin_mat)/10) == 0) cat('*')
  idx <- which(assignment_vec == j)
  sapply(dataset_vec, function(dataset){
    length(which(all_data$dataset[idx] == dataset))
  })
})
lin_total_count <- t(lin_total_count)
lin_sum <- rowSums(lin_total_count)
lin_total_count <- lin_total_count[order(lin_sum, decreasing = T),]
lin_total_count[1:10,]

save(assignment_vec, lin_mat, barcoding_res,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6/Writeup6_timeAll_lineage-assignments.RData")
