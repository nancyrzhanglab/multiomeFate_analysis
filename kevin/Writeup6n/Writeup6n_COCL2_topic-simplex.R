rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6n/Writeup6n_topics.RData")
tab_mat <- t(tab_mat)
load("../../../../out/kevin/Writeup6n/Writeup6n_COCL2_day10_lineage-imputation_stepdown_concise-postprocessed.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cocl2_d10 <- cocl2_fasttopics@cell.embeddings
cell_idx <- grep("day10_COCL2", names(lineage_assignments))
cocl2_d10 <- cocl2_d10[cell_idx,]
growth_potential <- cell_imputed_count[rownames(cocl2_d10)]

# pick the top 3 topics
# cor_vec <- sapply(1:ncol(cocl2_d10), function(j){
#   abs(stats::cor(cocl2_d10[,j], growth_potential))
# })
# col_idx <- order(cor_vec, decreasing = T)[1:3]
# 
# cocl2_d10 <- cocl2_d10[,col_idx]
# sum_vec <- rowSums(cocl2_d10)
# min_val <- 0.1
# cocl2_d10 <- cocl2_d10[sum_vec >= min_val,]
# for(i in 1:nrow(cocl2_d10)){
#   cocl2_d10[i,] <- cocl2_d10[i,]/sum(cocl2_d10[i,])
# }

lineage_vec <- lineage_assignments[rownames(cocl2_d10)]
tab_vec <- table(lineage_vec)
lineage_names <- names(tab_vec)[which(tab_vec >= 20)]

save(cocl2_d10, lineage_names, lineage_vec, tab_mat, growth_potential,
     file = "../../../../out/kevin/Writeup6n/Writeup6n_COCL2_topic-simplex.RData")


