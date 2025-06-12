rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"

load(paste0(out_folder, "Writeup10a_ppStep4_lineage.RData"))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

########

# plot histogram of correlation on log-10 scale
cell_lower_limit <- 100
lin_mat <- SeuratObject::LayerData(all_data,
                                   assay = "Lineage",
                                   layer = "counts")
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]
lin_mat_t <- Matrix::t(lin_mat)
nlineages <- nrow(lin_mat)
num_cells <- sapply(1:nlineages, function(j){
  length(multiomeFate:::.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
})
lin_mat_t <- lin_mat_t[,which(num_cells >= cell_lower_limit)]
# compute a correlation matrix (# rows/columsn = number of lineages)
cor_mat <- multiomeFate:::.custom_correlation(lin_mat_t)
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA

cor_vec <- as.numeric(cor_mat)
cor_vec <- cor_vec[!is.na(cor_vec)]
cor_vec <- pmax(cor_vec, 0)

hist_res <- graphics::hist(cor_vec, 
                           breaks = seq(0, 1, length.out = 21),
                           plot = FALSE)
hist_res$counts <- log10(hist_res$counts+1)

png(paste0(plot_folder, "Writeup10a_ppStep5b_lineage-correlation.png"),
    height = 1200, width = 1800, res = 300, units = "px")
plot(hist_res, 
     xlab = "Correlation between two lineages",
     ylab = "Frequency (Log10+1)",
     main = "Lineage barcoding duplicates")
graphics.off()

# plot histogram of (maximum - second_maximum)/(maximum), with a bar reserved for specifically zero
ratio_vec <- sapply(1:ncol(lin_mat), function(j){
  barcode_vals <- multiomeFate:::.nonzero_col(lin_mat, 
                                              col_idx = j, 
                                              bool_value = TRUE)
  if(length(barcode_vals) < 2) return(NA)
  max_vals <- sort(barcode_vals, decreasing = TRUE)[1:2]
  (max_vals[1] - max_vals[2])/max_vals[2]
})
ratio_vec <- ratio_vec[!is.na(ratio_vec)]
ratio_vec <- pmin(ratio_vec, 10)
min_nonzero <- min(ratio_vec[ratio_vec > 0])
ratio_vec <- ratio_vec + 0.5

custom_breaks <- c(min_nonzero/2-1, seq(min_nonzero/2, 11, by = 1)) + 0.5

png(paste0(plot_folder, "Writeup10a_ppStep5b_lineage-ratio_maximum.png"),
    height = 1200, width = 1800, res = 300, units = "px")
hist_res <- graphics::hist(ratio_vec, 
                           freq = FALSE,
                           breaks = custom_breaks,
                           xlab = "Ratio between the maximum and second maximum lineage count",
                           ylab = "Frequency",
                           main = "Lineage barcoding shared maximum")
graphics.off()

