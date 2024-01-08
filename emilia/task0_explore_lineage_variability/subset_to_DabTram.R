library(Seurat)

load("/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup6b/Writeup6b_all-data.RData")
day0_idx <- intersect(which(all_data$dataset == "day10_COCL2"),
                       which(!is.na(all_data$assigned_lineage)))
day0_idx <- intersect(day0_idx,which(all_data$assigned_posterior >= 0.5))
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[day0_idx] <- TRUE
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep == TRUE)
saveRDS(all_data2, '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/data/day10_COCL2.rds')