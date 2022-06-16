rm(list=ls())

library(Seurat)
library(Signac)
source("barcodePlot_functions.R")

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

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

lineage_mat <- all_data[["Lineage"]]@counts
lineage_frequency <- Matrix::rowSums(all_data[["Lineage"]]@counts)
lineage_mat <- lineage_mat[which(lineage_frequency >= 1),]

dataset_unique <- unique(all_data$dataset)
lineage_frequency_list <- lapply(dataset_unique, function(dataset){
  cell_idx <- which(all_data$dataset == dataset)
  tmp <- Matrix::t(lineage_mat[,cell_idx])
  sapply(1:ncol(tmp), function(j){
    length(.nonzero_col(tmp, col_idx = j, bool_value = F))
  })
})
sapply(lineage_frequency_list, function(x){quantile(x)})
sapply(lineage_frequency_list, function(x){length(which(x>0))})

lineage_mat2 <- Matrix::t(lineage_mat)

idx <- which(lineage_frequency_list[[1]] == 30)
colnames(lineage_mat2)[idx]
cell_idx <- .nonzero_col(lineage_mat2, col_idx = idx, bool_value = F)
table(all_data$dataset[cell_idx])

##########

p <- ncol(lineage_mat2)
num_cells_per_lineage <- sapply(1:p, function(j){
  length(.nonzero_col(lineage_mat2, col_idx = j, bool_value = F))
})
n <- ncol(lineage_mat)
cell_anybarcodes <- which(sapply(1:n, function(i){
  length(.nonzero_col(lineage_mat, col_idx = i, bool_value = F)) > 0
}))
lineage_pure_ingredients <- lapply(1:n, function(j){
  if(j %% floor(n/10) == 0) cat('*')
  
  lineage_idx <- .nonzero_col(lineage_mat, col_idx = j, bool_value = F)
  lineage_val <- .nonzero_col(lineage_mat, col_idx = j, bool_value = T)
  
  if(length(lineage_idx) == 0) return(numeric(0))
  
  sorted_vals <- sort(lineage_val, decreasing = T)
  if(length(lineage_idx) == 1 || sorted_vals[1] > sorted_vals[2]){
    c(lineage_idx[which.max(lineage_val)], j, 1)
    
  } else if(sorted_vals[1] == sorted_vals[2]){
    lineage_idx2 <- lineage_idx[which(lineage_val == max(lineage_val))]
    lineage_vals2 <- num_cells_per_lineage[lineage_idx2]
    
    sorted_vals2 <- sort(lineage_vals2, decreasing = T)
    if(sorted_vals2[1] > sorted_vals2[2]){
      c(lineage_idx2[which.max(lineage_vals2)], j, 1)
    } else {
      print(j)
      max_val <- max(lineage_vals2)
      tmp_idx <- which(lineage_vals2 == max_val)
      c(lineage_idx2[sample(tmp_idx, 1)], j, 1)
    }
    
  } else {
    numeric(0)
  }
})
tmp <- do.call(rbind, lineage_pure_ingredients)
lineage_pure_mat <- Matrix::sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3],
                                         dims = dim(lineage_mat))
rownames(lineage_pure_mat) <- rownames(lineage_mat)
colnames(lineage_pure_mat) <- colnames(lineage_mat)
idx <- which(Matrix::rowSums(lineage_pure_mat) > 0)
lineage_pure_mat <- lineage_pure_mat[idx,]
quantile(Matrix::rowSums(lineage_pure_mat))
quantile(Matrix::colSums(lineage_pure_mat))

all_data[["Lineage"]]@data <- lineage_pure_mat
col_palette <- scales::hue_pal()(nrow(all_data[["Lineage"]]@data))
col_vec <- sapply(1:ncol(all_data), function(i){
  if(i %% floor(ncol(all_data)/10) == 0) cat('*')
  
  idx <- .nonzero_col(all_data[["Lineage"]]@data, col_idx = i, bool_value = F)
  if(length(idx) > 0){
    col_palette[idx[1]]
  } else {
    NA
  }
})
length(which(is.na(col_vec)))/length(col_vec)

png("../../../../out/figures/Writeup4e/Writeup4e_rna_umap-assigned_lineage.png",
    height = 3000, width = 3000, units = "px", res = 300)
plot(x = all_data[["umap"]]@cell.embeddings[,1],
     y = all_data[["umap"]]@cell.embeddings[,2],
     xlab = "umap_1", ylab = "umap_2",
     pch = 16,
     col = col_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     main = "Cells by assigned lineage\n(Color for barcodes randomly assigned)")
axis(1)
axis(2)
graphics.off()

#######################
#######################
#######################

mat <- all_data[["Lineage"]]@counts
mat <- mat[which(Matrix::rowSums(mat) != 0),]
mat@x <- rep(1, length(mat@x))
tmp <- mat[,which(all_data$dataset %in% c("week5_CIS", "week5_COCL2", "week5_DABTRAM"))]
lineage_order <- order(Matrix::rowSums(tmp), decreasing = T)
dataset_vec <- c("day0", "day10_CIS", "week5_CIS", 
                 NA, "day10_COCL2", "week5_COCL2", 
                 NA, "day10_DABTRAM", "week5_DABTRAM")
png("../../../../out/figures/Writeup4e/Writeup4e_barcode_size.png",
    height = 3500, width = 5000, units = "px", res = 300)
par(mfrow = c(3,3), mar = c(4,4,4,0.5))
for(dataset in dataset_vec){
  if(is.na(dataset)){
    plot(NA, xaxt = "n", yaxt = "n", xlim = c(0,1), ylim = c(0,1), 
         bty = "n", main = "", xlab = "", ylab = "")
  } else {
    mat2 <- mat[lineage_order,which(all_data$dataset == dataset)]
    vec <- log10(Matrix::rowSums(mat2)+1)
    names(vec) <- NULL
    barplot(vec, xlab = "Lineage barcode", main = dataset,
            ylab = "Number of cells (log10+1)", ylim = c(0,3.8))
  }
  
}
graphics.off()

for(selection_dataset in c("week5_CIS", "week5_COCL2", "week5_DABTRAM")){
  mat <- all_data[["Lineage"]]@counts
  mat <- mat[which(Matrix::rowSums(mat) != 0),]
  mat@x <- rep(1, length(mat@x))
  tmp <- mat[,which(all_data$dataset == selection_dataset)]
  lineage_order <- order(Matrix::rowSums(tmp), decreasing = T)
  
  dataset_vec <- c("day0", "day10_CIS", "week5_CIS", 
                   NA, "day10_COCL2", "week5_COCL2", 
                   NA, "day10_DABTRAM", "week5_DABTRAM")
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_size-", selection_dataset, ".png"),
      height = 3500, width = 5000, units = "px", res = 300)
  par(mfrow = c(3,3), mar = c(4,4,4,0.5))
  for(dataset in dataset_vec){
    if(is.na(dataset)){
      plot(NA, xaxt = "n", yaxt = "n", xlim = c(0,1), ylim = c(0,1), 
           bty = "n", main = "", xlab = "", ylab = "")
    } else {
      mat2 <- mat[lineage_order,which(all_data$dataset == dataset)]
      vec <- log10(Matrix::rowSums(mat2)+1)
      names(vec) <- NULL
      barplot(vec, xlab = paste0("Lineage barcode (Ordered by ", selection_dataset), 
              main = dataset,
              ylab = "Number of cells (log10+1)", ylim = c(0,3.8))
    }
    
  }
  graphics.off()
}

dataset_vec <- unique(all_data$dataset)
num_counts <- sapply(dataset_vec, function(dataset){
  mat <- all_data[["Lineage"]]@counts
  mat <- mat[which(Matrix::rowSums(mat) != 0),]
  mat@x <- rep(1, length(mat@x))
  tmp <- mat[,which(all_data$dataset == dataset)]
  Matrix::rowSums(tmp)
})
max_vec <- rowSums(num_counts)
num_counts2 <- num_counts[order(max_vec, decreasing = T),]
num_counts2[1:15,]

num_counts2 <- num_counts[order(num_counts[,"week5_CIS"], decreasing = T),]
num_counts2[1:15,]

num_counts2 <- num_counts[order(num_counts[,"week5_COCL2"], decreasing = T),]
num_counts2[1:15,]

num_counts2 <- num_counts[order(num_counts[,"week5_DABTRAM"], decreasing = T),]
num_counts2[1:15,]

#######################
#######################
#######################

tmp <- all_data[["Lineage"]]@counts[,which(all_data$dataset %in% c("week5_CIS", "week5_COCL2", "week5_DABTRAM"))]
tmp@x <- rep(1, length(tmp@x))
lineage_idx <- order(Matrix::rowSums(tmp), decreasing = T)[1:10]
lineage_name <- rownames(all_data[["Lineage"]]@counts)[lineage_idx]
lineage_mat <- all_data[["Lineage"]]@counts

col_palette_enrichment <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(100)
col_transparency <- rgb(1,1,1,0.2)
last_two_letters <- substr(col_transparency, 
                           start = nchar(col_transparency)-1,
                           stop = nchar(col_transparency))
col_palette_enrichment <- sapply(col_palette_enrichment, function(col){
  paste0(col, last_two_letters)
})

for(lineage in lineage_name){
  print(lineage)
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, ".png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["umap"]]@cell.embeddings,
                   mode = "indicator")
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, "-enrichment.png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["umap"]]@cell.embeddings,
                   col_palette_enrichment = col_palette_enrichment,
                   mode = "enrichment")
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, "_atac.png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["atac.umap"]]@cell.embeddings,
                   mode = "indicator")
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, "-enrichment_atac.png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["atac.umap"]]@cell.embeddings,
                   col_palette_enrichment = col_palette_enrichment,
                   mode = "enrichment")
  graphics.off()
}

################

# now let's select lineages with large number of cells at the terminal states, but also
# highly enriched day0 cells
tmp <- all_data[["Lineage"]]@counts[,which(all_data$dataset %in% c("week5_CIS", "week5_COCL2", "week5_DABTRAM"))]
lineage_idx <- order(Matrix::rowSums(tmp), decreasing = T)[1:200]
lineage_name <- rownames(all_data[["Lineage"]]@counts)[lineage_idx]
lineage_mat <- all_data[["Lineage"]]@counts
lineage_mat2 <- all_data[["Lineage"]]@counts[lineage_name,which(all_data$dataset == "day0")]
lineage_mat2 <- Matrix::t(lineage_mat2)
day0_cellcount <- sapply(1:ncol(lineage_mat2), function(j){
  length(.nonzero_col(lineage_mat2, col_idx = j, bool_value = F))
})
lineage_name <- lineage_name[which(day0_cellcount >= 4)]
enrichment_score <- sapply(lineage_name, function(lineage){
  j <- which(colnames(lineage_mat2) == lineage)
  j2 <- which(rownames(lineage_mat) == lineage)
  cell_idx <- .nonzero_col(lineage_mat2, col_idx = j, bool_value = F)
  cell_names <- rownames(lineage_mat2)[cell_idx]
  cell_enrichment <- sapply(cell_names, function(cell_name){
    i <- which(colnames(lineage_mat) == cell_name)
    tmp_idx <- .nonzero_col(lineage_mat, col_idx = i, bool_value = F)
    tmp_vec <- .nonzero_col(lineage_mat, col_idx = i, bool_value = T)
    stopifnot(j2 %in% tmp_idx)
    if(length(tmp_idx) == 1){
      denominator <- tmp_vec[1]
    } else {
      denominator <- sum(sort(tmp_vec, decreasing = T)[1:2])
    }
    
    tmp_val <- tmp_vec[which(tmp_idx == j2)]/denominator
  })
  
  mean(cell_enrichment) - sd(cell_enrichment)
})
lineage_name <- names(enrichment_score)[order(enrichment_score, decreasing = T)[1:10]]
for(lineage in lineage_name){
  print(lineage)
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, ".png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["umap"]]@cell.embeddings,
                   mode = "indicator")
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, "-enrichment.png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["umap"]]@cell.embeddings,
                   col_palette_enrichment = col_palette_enrichment,
                   mode = "enrichment")
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, "_atac.png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["atac.umap"]]@cell.embeddings,
                   mode = "indicator")
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage, "-enrichment_atac.png"),
      height = 3000, width = 3000, units = "px", res = 300)
  barplot_function(lineage_name = lineage,
                   lineage_mat = all_data[["Lineage"]]@counts,
                   membership_vec = all_data$dataset,
                   umap_mat = all_data[["atac.umap"]]@cell.embeddings,
                   col_palette_enrichment = col_palette_enrichment,
                   mode = "enrichment")
  graphics.off()
}


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

# make PDF versions of the plots
library(png)
lineage_vec <- c("Lin1550", "Lin81746", "Lin56370", "Lin89477")
lineage_mat <- all_data[["Lineage"]]@counts
lineage_mat2 <- Matrix::t(lineage_mat)
membership_vec <- all_data$dataset
col_palette_enrichment <- paste0(grDevices::colorRampPalette(c('lightgrey', 'blue'))(100), "33")
mode <- "enrichment"

for(lineage_name in lineage_vec){
  print(lineage_name)
  
  for(modality in c(1,2)){
    if(modality == 1){
      umap_mat <- all_data[["atac.umap"]]@cell.embeddings
      modality_name <- "_atac"
    } else {
      umap_mat <- all_data[["umap"]]@cell.embeddings
      modality_name <- ""
    }
    
    j <- which(colnames(lineage_mat2) == lineage_name)
    cell_idx <- .nonzero_col(lineage_mat2, col_idx = j, bool_value = F)
    
    # https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
    pdf(paste0("../../../../out/figures/Writeup4e/Writeup4e_barcode_", lineage_name, "-enrichment", modality_name, ".pdf"),
        height = 10, width = 10)
    dataset_vec <- c("day0", "day10_CIS", "week5_CIS",
                     NA, "day10_COCL2", "week5_COCL2",
                     NA, "day10_DABTRAM", "week5_DABTRAM")
    
    par(mfrow = c(3,3), mar = c(4,4,4,0.5))
    for(dataset in dataset_vec){
      print(dataset)
      cell_idx_all <- which(membership_vec == dataset)
      
      if(is.na(dataset)) {
        plot(NA, xaxt = "n", yaxt = "n", xlim = c(0,1), ylim = c(0,1), 
             bty = "n", main = "", xlab = "", ylab = "")
        
      } else {
        cell_idx2 <- intersect(cell_idx, cell_idx_all)
        
        plot(x = umap_mat[,1], y = umap_mat[,2],
             xlab = colnames(umap_mat)[1], ylab = colnames(umap_mat)[2],
             xaxt = "n", yaxt = "n", bty = "n",
             main = paste0(lineage_name, " for ", dataset, "\n(Present in ", 
                           length(cell_idx2), " of ", length(cell_idx_all), " cells)"),
             type = "n")
        axis(1); axis(2)
        
        # Extract plot area in both user and physical coordinates
        coords <- par("usr")
        gx <- grconvertX(coords[1:2], "user", "inches")
        gy <- grconvertY(coords[3:4], "user", "inches")
        width <- max(gx) - min(gx)
        height <- max(gy) - min(gy)
        
        # Get a temporary file name for our rasterized plot area
        tmp <- tempfile()
        
        # Can increase resolution from 300 if higher quality is desired.
        png(tmp, width = width, height = height, units = "in", res = 300, bg = "transparent")
        par(mar = c(0,0,0,0))
        plot(NA, xaxt = "n", yaxt = "n", xlim = coords[1:2], ylim = coords[3:4], 
             bty = "n", main = "", xlab = "", ylab = "")
        cell_enrichment <- sapply(cell_idx2, function(i){
          tmp_idx <- .nonzero_col(lineage_mat, col_idx = i, bool_value = F)
          tmp_vec <- .nonzero_col(lineage_mat, col_idx = i, bool_value = T)
          stopifnot(j %in% tmp_idx)
          if(length(tmp_idx) == 1){
            denominator <- tmp_vec[1]
          } else {
            denominator <- sum(sort(tmp_vec, decreasing = T)[1:2])
          }
          
          tmp_val <- tmp_vec[which(tmp_idx == j)]/denominator
        })
        cell_col <- sapply(cell_enrichment, function(val){
          col_palette_enrichment[ceiling(val*length(col_palette_enrichment))]
        })
        
        points(x = umap_mat[,1], y = umap_mat[,2],
               pch = 16, col = "gray")
        points(x = umap_mat[cell_idx2,1], y = umap_mat[cell_idx2,2],
               pch = 16, col = "white", cex = 2)
        points(x = umap_mat[cell_idx2,1], y = umap_mat[cell_idx2,2],
               pch = 16, col = cell_col, cex = 1.5)
        dev.off()
        
        # Windows users may have trouble with transparent plot backgrounds; if this is the case,
        # set bg = "white" above and move the legend plot command below the raster plot command.
        panel <- readPNG(tmp)
        rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
      }
    }
    dev.off()
  }
}




