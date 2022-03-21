rm(list=ls())
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_saver.RData")
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_exploration_alldata_onlyGEX.RData")
library(Seurat)
library(SAVER)
source("saver_helpers.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
gene_idx <- which(rownames(saver_res$estimate) %in% jackpot_genes)

cell_names <- rownames(all_data@meta.data)[which(all_data$original_dataset == "time0")]
cell_idx <- which(colnames(saver_res$estimate) %in% cell_names)
cor_time0 <- compute_cor_genes(saver_res, cell_idx = cell_idx,
                               gene_idx = gene_idx)
diag(cor_time0) <- 0

cell_names <- rownames(all_data@meta.data)[which(all_data$original_dataset == "time10_cis")]
cell_idx <- which(colnames(saver_res$estimate) %in% cell_names)
cor_time10_cis <- compute_cor_genes(saver_res, cell_idx = cell_idx,
                               gene_idx = gene_idx)
diag(cor_time10_cis) <- 0

cell_names <- rownames(all_data@meta.data)[which(all_data$original_dataset == "time10_cocl2")]
cell_idx <- which(colnames(saver_res$estimate) %in% cell_names)
cor_time10_cocl2 <- compute_cor_genes(saver_res, cell_idx = cell_idx,
                                    gene_idx = gene_idx)
diag(cor_time10_cocl2) <- 0

cell_names <- rownames(all_data@meta.data)[which(all_data$original_dataset == "time10_dabtram")]
cell_idx <- which(colnames(saver_res$estimate) %in% cell_names)
cor_time10_dabtram <- compute_cor_genes(saver_res, cell_idx = cell_idx,
                                      gene_idx = gene_idx)
diag(cor_time10_dabtram) <- 0

cor_list <- list(time0 = cor_time0,
                 time10_cis = cor_time10_cis,
                 time10_cocl2 = cor_time10_cocl2,
                 time10_dabtram = cor_time10_dabtram)

max_val <- max(sapply(cor_list, max))
min_val <- min(sapply(cor_list, min))
abs_val <- max(abs(c(max_val, min_val)))
break_vec <- seq(-abs_val, abs_val, length.out = 101)
  
# see https://stackoverflow.com/questions/60538007/heatmap2-graph-doesnt-display-complete-dataset
for(i in 1:length(cor_list)){
  png(paste0("../../../../out/figures/Writeup4b/Writeup4b_time0time10_saver_correlation_", names(cor_list)[i], ".png"),
      height = 1500, width = 1500, units = "px", res = 300)
  gplots::heatmap.2(cor_list[[i]], scale = "none", 
                    symm = T,
                    na.rm	= F,
                    col = gplots::bluered(100), 
                    breaks = break_vec,
                    symbreaks = T,
                    cexRow = 0.5,
                    cexCol = 0.5,
                    trace = "none", density.info = "none")
  graphics.off()
}

####################

zz <- t(saver_res$estimate)
zz <- scale(zz)
svd_res <- irlba::irlba(zz, nv = 53)
dimred <- svd_res$u[,1:50] %*% diag(svd_res$d[1:30])
quantile(dimred, probs = seq(0,1,length.out=11))
rownames(dimred) <- colnames(saver_res$estimate)
keep_vec <- rep(0, ncol(all_data))
keep_vec[which(colnames(all_data) %in% colnames(saver_res$estimate))] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)

set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)

all_data[["saver"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings)

plot1 <-Seurat::DimPlot(all_data, reduction = "saver",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (SAVER),\n", length(all_data[["RNA"]]@var.features), " genes, using 50 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_saver_umap_original_dataset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


