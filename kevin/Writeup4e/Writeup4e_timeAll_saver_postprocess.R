rm(list=ls())
library(Seurat)
library(Signac)
library(SAVER)
load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_saver.RData")

# modified from https://github.com/mohuangx/SAVER/blob/master/R/cor_adjust.R
compute_cor_genes <- function(x, cell_idx = NULL,
                              gene_idx = NULL) {
  
  if(all(is.null(cell_idx))){
    cell_idx <- 1:ncol(x$estimate)
  }
  if(all(is.null(gene_idx))){
    gene_idx <- 1:nrow(x$estimate)
  }
  ngenes <- length(gene_idx)
  adj_vec <- rep(0, ngenes)
  
  cor_mat <- cor(t(x$estimate[gene_idx,cell_idx]))
  for (i in 1:ngenes) {
    adj_vec[i] <- sqrt(var(x$estimate[gene_idx[i],cell_idx], na.rm = TRUE)/
                         (var(x$estimate[gene_idx[i],cell_idx], na.rm = TRUE) +
                            mean(x$se[gene_idx[i],cell_idx]^2, na.rm = TRUE)))
  }
  adj_mat <- outer(adj_vec, adj_vec)
  cor_adj <- adj_mat*cor_mat
  
  cor_adj
}

jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")
gene_idx <- which(rownames(saver_res$estimate) %in% jackpot_genes)

dataset_vec <- unique(all_data$dataset)
cor_list <- lapply(dataset_vec, function(dataset){
  print(dataset)
  
  cell_names <- colnames(all_data)[which(all_data$dataset == dataset)]
  cell_idx <- which(colnames(saver_res$estimate) %in% cell_names)
  cor_mat <- compute_cor_genes(saver_res, cell_idx = cell_idx,
                               gene_idx = gene_idx)
  diag(cor_mat) <- 0
  cor_mat
})
cor_list[[1]][1:5,1:5]

max_val <- max(sapply(cor_list, max))
min_val <- min(sapply(cor_list, min))
abs_val <- max(abs(c(max_val, min_val)))
break_vec <- seq(-abs_val, abs_val, length.out = 101)

for(i in 1:length(cor_list)){
  print(dataset_vec[i])
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_saver_correlation_", dataset_vec[i], ".png"),
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
  
  pdf(paste0("../../../../out/figures/Writeup4e/Writeup4e_saver_correlation_", dataset_vec[i], ".pdf"),
      height = 5, width = 5)
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


