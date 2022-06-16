rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_idx <- which(all_data$dataset == "day0")
mat <- all_data[["Lineage"]]@counts[,cell_idx]
lineage_total_count <- Matrix::rowSums(mat)
mat <- mat[which(lineage_total_count > 0),]
cell_total_count <- Matrix::colSums(mat)
mat <- mat[,which(cell_total_count > 0)]
quantile(mat@x)
mat_log <- mat
mat_log@x <- log1p(mat_log@x)
quantile(mat_log@x)
length(mat_log@x)/prod(dim(mat_log))*100
mat_log2 <- mat_log[order(Matrix::rowSums(mat_log)), order(Matrix::colSums(mat_log))]
mat_log2 <- as.matrix(mat_log2)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png("../../../../out/figures/Writeup4f/Writeup4f_day0_barcode_heatmap.png",
    height = 2500, width = 3500, units = "px", res = 300)
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(.rotate(mat_log2), xaxt = "n", yaxt = "n", main = "", bty = "n", xlab = "", ylab = "",
      col = c("white", hcl.colors(12, "YlOrRd", rev = TRUE)[-c(1:4)]),
      breaks = c(-.5, .1, seq(min(mat_log@x), max(mat_log@x), length = 9)[-1]))
graphics.off()

####################

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

num_lineages_per_cell <- sapply(1:ncol(mat), function(j){
  length(.nonzero_col(mat, col_idx = j, bool_value = F))
})

quantile(num_lineages_per_cell)

########################################
# now the same but for day10_COCL2

cell_idx <- which(all_data$dataset == "day10_COCL2")
mat <- all_data[["Lineage"]]@counts[,cell_idx]
lineage_total_count <- Matrix::rowSums(mat)
mat <- mat[which(lineage_total_count > 0),]
cell_total_count <- Matrix::colSums(mat)
mat <- mat[,which(cell_total_count > 0)]
mat_log <- mat
mat_log@x <- log1p(mat_log@x)
quantile(mat_log@x)
length(mat_log@x)/prod(dim(mat_log))*100
mat_log2 <- mat_log[order(Matrix::rowSums(mat_log)), order(Matrix::colSums(mat_log))]
mat_log2 <- as.matrix(mat_log2)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png("../../../../out/figures/Writeup4f/Writeup4f_day10-COCL2_barcode_heatmap.png",
    height = 2500, width = 3500, units = "px", res = 300)
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(.rotate(mat_log2), xaxt = "n", yaxt = "n", main = "", bty = "n", xlab = "", ylab = "",
      col = c("white", hcl.colors(12, "YlOrRd", rev = TRUE)[-c(1:4)]),
      breaks = c(-.5, .1, seq(min(mat_log@x), max(mat_log@x), length = 9)[-1]))
graphics.off()


num_lineages_per_cell <- sapply(1:ncol(mat), function(j){
  length(.nonzero_col(mat, col_idx = j, bool_value = F))
})

quantile(num_lineages_per_cell)

