rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep4_lineage.RData"))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

lin_mat <- SeuratObject::LayerData(all_data,
                                   assay = "Lineage",
                                   layer = "counts")
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]

res_clusters <- multiomeFate:::barcode_clustering(lin_mat)
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

# combining lineages
lin_mat <- as.matrix(lin_mat)
lin_mat2 <- lin_mat
lin_mat <- multiomeFate:::barcode_combine(lin_mat = lin_mat,
                                         lineage_clusters = res_clusters$lineage_clusters,
                                         verbose = 1)

# check:
k <- which.max(sapply(res_clusters$lineage_clusters, length))
idx_prev <- which(rownames(lin_mat2) %in% res_clusters$lineage_clusters[[k]])
idx_new <- which(rownames(lin_mat) %in% res_clusters$lineage_clusters[[k]])
stopifnot(length(idx_new) == 1)
zz <- colSums(lin_mat2[idx_prev,])
yy <- lin_mat[idx_new,]
stopifnot(sum(abs(zz-yy)) == 0)

#########################

posterior_res <- multiomeFate:::barcoding_posterior(lin_mat = lin_mat,
                                                   verbose = 1)

maxBhat <- apply(posterior_res$posterior_mat,2,max) #chosen lineage 
argmaxBhat <- apply(posterior_res$posterior_mat, 2, which.max)
n <- ncol(lin_mat)
is_unique_vec_b <- sapply(1:n, function(i){length(which(abs(posterior_res$posterior_mat[,i] - maxBhat[i]) <= 1e-3)) == 1})
naive_idx <- which(all_data$dataset == "day0")
lineage_sum <- colSums(lin_mat)
table(is_unique_vec_b, lineage_sum > 0)
table(is_unique_vec_b[naive_idx], lineage_sum[naive_idx] > 0)

assignment_vec <- multiomeFate:::barcoding_assignment(posterior_mat = posterior_res$posterior_mat,
                                                     difference_val = 0.2,
                                                     verbose = 1)

SeuratObject::LayerData(all_data,
                        assay = "Lineage",
                        layer = "data") <- posterior_res$posterior_mat
all_data$assigned_lineage <- assignment_vec
all_data$assigned_posterior <- maxBhat

print("Before subsetting")
print(all_data)
print(stats::quantile(all_data$assigned_posterior))

keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
keep_vec[which(all_data$assigned_posterior > 0.5)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

print("After subsetting")
print(all_data)

save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep5_barcode-assignment.RData"))

###################################

# do a simple check
posterior_res$posterior_mat[assignment_vec[1],1]
sort(posterior_res$posterior_mat[,2], decreasing = T)[1:5]

tab_mat <- table(assignment_vec, all_data$dataset)
length(which(tab_mat[,1]>=1))
quantile(tab_mat[tab_mat[,1]>=1,1])

tab_mat[order(rowSums(tab_mat), decreasing = T)[1:10],]
tab_mat[order(rowSums(tab_mat), decreasing = F)[1:10],]

tab_mat[order(rowSums(tab_mat[,2:4]), decreasing = T)[1:10],]
tab_mat[order(rowSums(tab_mat[,2:4]), decreasing = F)[1:10],]

table(tab_mat[,1], rowSums(tab_mat[,5:7]) > 10)
table(tab_mat[,1], cut(rowSums(tab_mat[,5:7]), breaks = c(-0.5,0.5,5.5,10.5,20.5,max(tab_mat))))
zz <- table(tab_mat[,1], rowSums(tab_mat[,5:7]) > 10)
1-zz["0","TRUE"]/sum(zz[,2])
sum(zz[,2])/sum(zz)

table(tab_mat[,1], rowSums(tab_mat[,2:4]) > 10)
table(tab_mat[,1], cut(rowSums(tab_mat[,2:4]), breaks = c(-0.5,0.5,5.5,10.5,20.5,max(tab_mat))))
zz <- table(tab_mat[,1], rowSums(tab_mat[,2:4]) > 10)
1-zz["0","TRUE"]/sum(zz[,2])
sum(zz[,2])/sum(zz)

table(cut(rowSums(tab_mat[,2:4]), breaks = c(-0.5,0.5,5.5,10.5,20.5,max(tab_mat))), 
      cut(rowSums(tab_mat[,5:7]), breaks = c(-0.5,0.5,5.5,10.5,20.5,max(tab_mat))))
zz <- table(cut(rowSums(tab_mat[,2:4]), breaks = c(-0.5,0.5,5.5,10.5,20.5,max(tab_mat))), 
            cut(rowSums(tab_mat[,5:7]), breaks = c(-0.5,0.5,5.5,10.5,20.5,max(tab_mat))))
sum(zz[1,-1])/sum(zz)
sum(zz[-1,1])/sum(zz)


tab_mat2 <- tab_mat[rowSums(tab_mat[,5:7])>20,]
tab_mat2[order(apply(tab_mat2[,5:7],1,function(x){
  tmp <- sort(x, decreasing = T)[1:2]
  (tmp[1]-tmp[2])
}), decreasing = F)[1:10],]
tab_mat2[order(apply(tab_mat2[,5:7],1,function(x){
  tmp <- sort(x, decreasing = T)[1:2]
  (tmp[1]-tmp[2])
}), decreasing = T)[1:10],]


tab_mat2 <- tab_mat[rowSums(tab_mat[,2:4])>20,]
tab_mat2[order(apply(tab_mat2[,2:4],1,function(x){
  tmp <- sort(x, decreasing = T)[1:2]
  (tmp[1]-tmp[2])/tmp[2]
}), decreasing = F)[1:10],]
tab_mat2[order(apply(tab_mat2[,2:4],1,function(x){
  tmp <- sort(x, decreasing = T)[1:2]
  (tmp[1]-tmp[2])
}), decreasing = T)[1:10],]

