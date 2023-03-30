rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

################

treatment <- "DABTRAM"
Seurat::DefaultAssay(all_data) <- "ATAC"

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 10)]
cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
cell_names2 <- colnames(all_data)[which(all_data$dataset == "day0")]
cell_names_winning <- intersect(cell_names1, cell_names2)
cell_names_losing <- setdiff(cell_names2, cell_names1)
ident_vec <- rep(NA, ncol(all_data))
names(ident_vec) <- colnames(all_data)
ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
ident_vec[cell_names_losing] <- paste0("day0_lose_", treatment)
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"
table(Seurat::Idents(all_data))

# from https://stuartlab.org/signac/articles/pbmc_vignette.html#find-differentially-accessible-peaks-between-clusters
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = paste0("day0_win_", treatment),
  ident.2 = paste0("day0_lose_", treatment),
  test.use = 'LR',
  latent.vars = 'peak_region_fragments',
  verbose = T
)

save(date_of_run, session_info, de_res,
     file = "../../../../out/kevin/Writeup6d/Writeup6d_DABTRAM_differential-peak_day0-week5.RData")
