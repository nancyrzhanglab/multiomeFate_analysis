rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

lin_mat <- all_data[["Lineage"]]@counts
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]

n <- ncol(lin_mat)
dominant_assignment_pairs <- lapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  idx <- .nonzero_col(lin_mat, col_idx = i, bool_value = F)
  val <- .nonzero_col(lin_mat, col_idx = i, bool_value = T)
  
  if(length(idx) == 0) return(numeric(0))
  if(length(idx) == 1) return(c(idx[1],i))
  if(length(which(val == max(val))) == 1) return(c(idx[which.max(val)],i))
  return(numeric(0))
})
dominant_assignment_pairs <- do.call(rbind, dominant_assignment_pairs)
dominant_mat <- Matrix::sparseMatrix(i = dominant_assignment_pairs[,1],
                                     j = dominant_assignment_pairs[,2],
                                     x = rep(1, nrow(dominant_assignment_pairs)),
                                     dims = dim(lin_mat))
dominant_mat_t <- Matrix::t(dominant_mat)

p <- nrow(lin_mat)
library_vec <- all_data$nCount_RNA
cor_vec <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  cell_idx <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
  
  if(length(cell_idx) < 5) return(NA)
  library_tmp <- library_vec[cell_idx]
  lin_tmp <- lin_mat[j,cell_idx]
  
  df <- data.frame(lineage_count = lin_tmp, library_count = library_tmp)
  lm_res <- stats::lm(lineage_count ~ . - 1, data = df)
  summary(lm_res)$r.squared
})
length(which(is.na(cor_vec)))
quantile(cor_vec, na.rm = T)

png("../../../../out/figures/Writeup4f/Writeup4f_lineage-library-correlation_histogram.png",
    height = 1500, width = 2000, units = "px", res = 300)
hist(cor_vec, xlab = "Correlation between lineage and GEX library count",
     main = "Only lineages with more than 5 winning cells")
graphics.off()

num_dominant <- sapply(1:p, function(j){
  length(.nonzero_col(dominant_mat_t, col_idx = j, bool_value = F))
})
idx <- which(!is.na(cor_vec))
stats::cor(cor_vec[idx], log(num_dominant[idx]))
  
non_na_idx <- which(!is.na(cor_vec))
idx <- non_na_idx[which.min(abs(cor_vec[non_na_idx] - median(cor_vec[non_na_idx])))]

p <- nrow(lin_mat)
library_vec <- all_data$nCount_RNA
lin_mat_t <- Matrix::t(lin_mat)
cor_vec2 <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  cell_idx <- .nonzero_col(lin_mat_t, col_idx = j, bool_value = F)
  cell_idx2 <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
  cell_idx <- setdiff(cell_idx, cell_idx2)
  
  if(length(cell_idx) < 5) return(NA)
  library_tmp <- library_vec[cell_idx]
  lin_tmp <- lin_mat[j,cell_idx]
  
  df <- data.frame(lineage_count = lin_tmp, library_count = library_tmp)
  lm_res <- stats::lm(lineage_count ~ . - 1, data = df)
  summary(lm_res)$r.squared
})
length(which(is.na(cor_vec2)))
quantile(cor_vec2, na.rm = T)

png("../../../../out/figures/Writeup4f/Writeup4f_lineage-library-correlation_histogram_nonwinning.png",
    height = 1500, width = 2000, units = "px", res = 300)
hist(cor_vec2, xlab = "Correlation between lineage and GEX library count",
     main = "Only lineages with more than 5 non-winning cells")
graphics.off()

####

p <- nrow(lin_mat)
dominant_count <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  idx <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
  if(length(idx) == 0) return(0)
  sum(lin_mat[j,idx])
})
background_count <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  idx <- .nonzero_col(lin_mat_t, col_idx = j, bool_value = F)
  idx2 <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
  idx <- setdiff(idx, idx2)
  
  if(length(idx) == 0) return(0)
  sum(lin_mat[j,idx])
})

dominant_count <- log10(dominant_count+1)
background_count <- log10(background_count+1)
cor_val <- cor(dominant_count,background_count)

png("../../../../out/figures/Writeup4f/Writeup4f_lineage_correlation_winners-vs-all.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot(x = dominant_count, y = background_count,
     xlab = "Log10 winner count", ylab = "Log10 background count",
     main = paste0("Correlation: ", round(cor_val, 2)),
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1), asp = T)
lines(c(-5,5), c(-5,5), col = 2, lwd = 2, lty = 2)
graphics.off()


dominant_num <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  idx <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
  length(idx)
})
background_num <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  idx <- .nonzero_col(lin_mat_t, col_idx = j, bool_value = F)
  idx2 <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
  idx <- setdiff(idx, idx2)
  length(idx)
})
dominant_num <- log10(dominant_num+1)
background_num <- log10(background_num+1)
cor_val <- cor(dominant_num,background_num)

png("../../../../out/figures/Writeup4f/Writeup4f_lineage_correlation_winners-vs-all_number.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot(x = dominant_num, y = background_num,
     xlab = "Log10 winner number of cells", ylab = "Log10 background number of cells",
     main = paste0("Correlation: ", round(cor_val, 2)),
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1), asp = T)
lines(c(-5,5), c(-5,5), col = 2, lwd = 2, lty = 2)
graphics.off()



###################

library_vec <- log10(all_data$nCount_RNA+1)
lineage_count <- log10(all_data$nCount_Lineage+1)

png(paste0("../../../../out/figures/Writeup4f/Writeup4f_lineage-library_scatterplot.png"),
    height = 2250, width = 2000, units = "px", res = 300)
par(mfrow = c(3,3), mar = c(4,4,4,0.5))
dataset_vec <- c("day0", "day10_CIS", "week5_CIS",
                 NA, "day10_COCL2", "week5_COCL2",
                 NA, "day10_DABTRAM", "week5_DABTRAM")
for(dataset in dataset_vec){
  
  if(is.na(dataset)) {
    plot(NA, xaxt = "n", yaxt = "n", xlim = c(0,1), ylim = c(0,1), 
         bty = "n", main = "", xlab = "", ylab = "")
    
  } else {
    cell_idx <- which(all_data$dataset == dataset)
    plot(library_vec[cell_idx], jitter(lineage_count[cell_idx]),
         xlab = "Log10 GEX library count", ylab = " Log10 lineage count",
         main = paste0(dataset, "\nCorrelation: ", round(cor(library_vec[cell_idx], lineage_count[cell_idx]), 2)),
         pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1))
  }
}
graphics.off()

png(paste0("../../../../out/figures/Writeup4f/Writeup4f_maximum-ratio_histogram.png"),
    height = 2250, width = 2000, units = "px", res = 300)
par(mfrow = c(3,3), mar = c(4,4,4,0.5))
dataset_vec <- c("day0", "day10_CIS", "week5_CIS",
                 NA, "day10_COCL2", "week5_COCL2",
                 NA, "day10_DABTRAM", "week5_DABTRAM")
for(dataset in dataset_vec){
  print(dataset)
  
  if(is.na(dataset)) {
    plot(NA, xaxt = "n", yaxt = "n", xlim = c(0,1), ylim = c(0,1), 
         bty = "n", main = "", xlab = "", ylab = "")
    
  } else {
    cell_idx <- which(all_data$dataset == dataset)
    ratio <- sapply(cell_idx, function(i){
      val <- .nonzero_col(lin_mat, col_idx = i, bool_value = T)
      if(length(val) == 0) return(NA)
      if(length(val) == 1) return(0)
      val <- sort(val, decreasing = T)
      val[2]/val[1]
    })
    ratio <- ratio[!is.na(ratio)]
    hist(ratio, breaks = 50,
         xlab = "Ratio of counts: 2nd/1st",
         main = paste0(dataset, "\n(", length(ratio), " of ", length(cell_idx), " cells)"))
  }
}
graphics.off()

#################################

lin_mat_t <- Matrix::t(lin_mat)
p <- nrow(lin_mat)
num_cells <- sapply(1:p, function(j){
  length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
})
lin_mat_t <- lin_mat_t[,which(num_cells >= 10)]
cor_mat <- stats::cor(as.matrix(lin_mat_t))

cor_mat2 <- cor_mat
diag(cor_mat2) <- 0
off_diag_vec <- cor_mat2[lower.tri(cor_mat2, diag = F)]
tmp <- hist(off_diag_vec, plot = F)
tmp$counts <- log10(tmp$counts+1)
png("../../../../out/figures/Writeup4f/Writeup4f_lineage-doubleinfection-correlation_histogram.png",
    height = 1500, width = 2000, units = "px", res = 300)
plot(tmp, xlab = "Correlation between two lineage counts (across cells)",
     main = paste0("Only lineages with 10 or more cells\n(Among all ", ncol(lin_mat_t), " lineages)"),
     ylab = "Log10 Frequency", col = "gray")
graphics.off()


lin_mat_t <- Matrix::t(lin_mat)
p <- nrow(lin_mat)
num_cells <- sapply(1:p, function(j){
  length(.nonzero_col(lin_mat_t, col_idx = j, bool_value = F))
})
lin_mat_t <- lin_mat_t[,which(num_cells >= 100)]
cor_mat <- stats::cor(as.matrix(lin_mat_t))

cor_mat2 <- cor_mat
diag(cor_mat2) <- 0
off_diag_vec <- cor_mat2[lower.tri(cor_mat2, diag = F)]
tmp <- hist(off_diag_vec, plot = F)
tmp$counts <- log10(tmp$counts+1)
png("../../../../out/figures/Writeup4f/Writeup4f_lineage-doubleinfection-correlation_histogram2.png",
    height = 1500, width = 2000, units = "px", res = 300)
plot(tmp, xlab = "Correlation between two lineage counts (across cells)",
     main = paste0("Only lineages with 100 or more cells\n(Among all ", ncol(lin_mat_t), " lineages)"),
     ylab = "Log10 Frequency", col = "gray")
graphics.off()

ordering <- order(rowSums(cor_mat2), decreasing = T)
cor_mat2 <- cor_mat2[ordering,ordering]

.rotate <- function(mat){t(mat)[,nrow(mat):1]}
png("../../../../out/figures/Writeup4f/Writeup4f_lineage-doubleinfection-correlation_heatmap.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(.rotate(abs(cor_mat2)), xaxt = "n", yaxt = "n", main = "", bty = "n", xlab = "", ylab = "",
      col = c("white", hcl.colors(12, "YlOrRd", rev = TRUE)[-c(1:4)]),
      breaks = c(-.5, seq(min(off_diag_vec), max(off_diag_vec)+1, length = 10)[-1]))
graphics.off()
