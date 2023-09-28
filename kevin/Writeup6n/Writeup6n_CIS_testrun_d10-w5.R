rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)
treatment <- "CIS"

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

# assign tiers for the upcoming feature-selection for ATAC
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 20)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 19))]
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]

tier3_idx <- intersect(
  which(all_data$assigned_lineage %in% tier3_lineages),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier2_idx <- intersect(
  which(all_data$assigned_lineage %in% tier2_lineages),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier1_idx <- intersect(
  which(all_data$assigned_lineage %in% tier1_lineages),
  which(all_data$dataset == paste0("day10_", treatment))
)
keep_vec <- rep(NA, ncol(all_data))
keep_vec[tier1_idx] <- paste0("1loser_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("3high_winner_", treatment)
all_data$tier <- keep_vec
# keep only the relevant cells for this particular analysis
all_data2 <- subset(all_data, tier %in% c(paste0("3high_winner_", treatment),
                                          paste0("2mid_winner_", treatment),
                                          paste0("1loser_", treatment)))

# construct cell_features matrix
topic_mat <- all_data2[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_mat <- all_data2[["lsi"]]@cell.embeddings

log10pval_vec <- sapply(1:ncol(atac_mat), function(j){
  x <- atac_mat[,j]
  y <- as.factor(all_data2$tier)
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(atac_mat)
print(sort(log10pval_vec))



