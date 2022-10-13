rm(list=ls())
library(Seurat)
library(Signac)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

all_data[["spliced"]] <- NULL
all_data[["unspliced"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6/Writeup6_timeAll_simplified.RData")

###############

lin_mat <- all_data[["Lineage"]]@counts
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]

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

# combining lineages
lin_mat <- as.matrix(lin_mat)
lin_mat2 <- lin_mat
lin_mat <- barcode_combine(lin_mat = lin_mat,
                           lineage_clusters = res_clusters$lineage_clusters)

barcoding_res <- barcoding_assignment(lin_mat = lin_mat,
                                      verbose = 1)

maxBhat <- apply(barcoding_res$Bhat,2,max) #chosen lineage 
argmaxBhat <- apply(barcoding_res$Bhat, 2, which.max)
argmaxX <- apply(lin_mat, 2, which.max) #largest count
maxX <- apply(lin_mat,2,max)
is_unique_vec_x <- sapply(1:n, function(i){length(which(lin_mat[,i] == maxX[i])) == 1})
is_unique_vec_b <- sapply(1:n, function(i){length(which(abs(barcoding_res$Bhat[,i] - maxBhat[i]) <= 1e-3)) == 1})
table(is_unique_vec_b)
disagrees <- which(argmaxX!=argmaxBhat & maxX>5& maxBhat>0.9 & is_unique_vec_x)
length(disagrees)
length(disagrees)/length(argmaxX)

disagrees2 <- which(argmaxX!=argmaxBhat & maxX>5 & is_unique_vec_x)
length(disagrees2)

disagrees3 <- which(argmaxX!=argmaxBhat & maxX>0 & is_unique_vec_x)
length(disagrees3)

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
  -diff(tmp)/sum(tmp+1)
})
quantile(difference_Xvec)
length(which(difference_Xvec > 1/3-1e-3))/length(difference_Xvec)
quantile(difference_Xvec[naive_idx])
length(which(difference_Xvec[naive_idx] >= 1/3+1e-3))/length(naive_idx)
length(which(difference_Xvec[naive_idx] >= 1/3-1e-3))/length(naive_idx)
difference_Bvec <- apply(barcoding_res$Bhat, 2, function(x){
  -diff(sort(x, decreasing = T)[1:2])
})
quantile(difference_Bvec)
length(which(difference_Bvec >= 0.1))/length(difference_Bvec)
length(which(difference_Bvec >= 0.1 & is_unique_vec_b))/length(difference_Bvec)
quantile(difference_Bvec[naive_idx])
length(which(difference_Bvec[naive_idx] >= 0.1 & is_unique_vec_b[naive_idx]))/length(naive_idx)
quantile(difference_Bvec[naive_idx[which(maxX[naive_idx] == 0)]])

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

