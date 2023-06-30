rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment <- "COCL2"

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

keep_vec <- rep(NA, ncol(all_data))
idx1 <- which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))
idx2 <- which(all_data$assigned_posterior >= 0.5)
idx3 <- which(!is.na(all_data$assigned_lineage))
keep_vec[intersect(intersect(idx1, idx2), idx3)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_COCL2",
                            dims = 1:30)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("day10_", treatment)] >= 3),
                                              which(tab_mat[,paste0("day10_", treatment)] <= 19))]
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] <= 2)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

tier3_idx <- intersect(
  which(all_data$assigned_lineage %in% tier3_lineages),
  which(all_data$dataset == paste0("day0"))
)
tier2_idx <- intersect(
  which(all_data$assigned_lineage %in% tier2_lineages),
  which(all_data$dataset == paste0("day0"))
)
tier1_idx <- intersect(
  which(all_data$assigned_lineage %in% tier1_lineages),
  which(all_data$dataset == paste0("day0"))
)
keep_vec <- rep(NA, ncol(all_data))
keep_vec[tier1_idx] <- paste0("1loser_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("3high_winner_", treatment)
table(keep_vec)
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep %in% c(paste0("3high_winner_", treatment),
                                          paste0("2mid_winner_", treatment),
                                          paste0("1loser_", treatment)))

save(all_data2, all_data,
     date_of_run, session_info, treatment,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_COCL2-day0_extracted.RData")


