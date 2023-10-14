rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6n/Writeup6n_topics.RData")
tab_mat <- t(tab_mat)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cocl2_d10 <- cocl2_fasttopics@cell.embeddings
cell_idx <- grep("day10_COCL2", names(lineage_assignments))
cocl2_d10 <- cocl2_d10[cell_idx,]

# pick the top 3 topics
sum_vec <- colSums(cocl2_d10)
col_idx <- order(sum_vec, decreasing = T)[1:3]
sd_vec <- apply(cocl2_d10, 2, sd)
round(sort(sd_vec*100, decreasing = T))

cocl2_d10 <- cocl2_d10[,col_idx]
sum_vec <- rowSums(cocl2_d10)
min_val <- 0.1
cocl2_d10 <- cocl2_d10[sum_vec >= min_val,]
for(i in 1:nrow(cocl2_d10)){
  cocl2_d10[i,] <- cocl2_d10[i,]/sum(cocl2_d10[i,])
}

lineage_vec <- lineage_assignments[rownames(cocl2_d10)]
tab_vec <- table(lineage_vec)
lineage_names <- names(tab_vec)[which(tab_vec >= 20)]

save(cocl2_d10, lineage_names, lineage_vec, tab_mat,
     file = "../../../../out/kevin/Writeup6n/Writeup6n_COCL2_topic-simplex.RData")


