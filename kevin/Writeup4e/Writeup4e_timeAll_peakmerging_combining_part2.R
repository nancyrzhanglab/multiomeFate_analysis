rm(list=ls())

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete_tmp.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
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

# integrate datasets one at time
file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_05_19_GEXandLin_time0", 
                  "2022_05_19_GEXandLin_day10_CIS",
                  "2022_05_19_GEXandLin_day10_COCL2", "2022_05_19_GEXandLin_day10_DABTRAM", 
                  "2022_05_19_GEXandLin_week5_CIS",
                  "2022_05_19_GEXandLin_week5_COCL2", "2022_05_19_GEXandLin_week5_DABTRAM")
dataset_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
                 "week5_CIS", "week5_COCL2", "week5_DABTRAM")
stopifnot(length(file_folders) == length(dataset_vec))

barcode_matrix_list <- lapply(1:length(file_folders), function(i){
  print(paste0("Working on ", dataset_vec[i]))
  file_folder <- file_folders[i]
  barcode_input <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  
  original_barcodes <- colnames(all_data)[which(all_data$dataset == dataset_vec[i])]
  ncell_original <- length(original_barcodes)
  ncell_barcode <- ncol(barcode_input[["Custom"]])
  colnames(barcode_input[["Custom"]]) <- paste0(dataset_vec[i], "_", colnames(barcode_input[["Custom"]]))
  intersect_barcodes <- intersect(colnames(barcode_input[["Custom"]]), colnames(all_data))
  ncell_intersect <- length(intersect_barcodes)
  ncell_all <- length(unique(c(colnames(barcode_input[["Custom"]]), original_barcodes)))
  print(paste0("Number of original cells: ", ncell_original))
  print(paste0("Number of barcoded cells: ", ncell_barcode))
  print(paste0("Percentage overlap: ", round(ncell_intersect/ncell_all, 2)))
  
  barcode_list <- lapply(1:length(original_barcodes), function(j){
    if(j %% floor(length(original_barcodes)/10) == 0) cat('*')
    
    barcode <- original_barcodes[j]
    barcode_idx <- which(colnames(barcode_input[["Custom"]]) == barcode)
    
    if(!barcode %in% intersect_barcodes) return()
    lineage_idx <- .nonzero_col(mat = barcode_input[["Custom"]],
                                col_idx = barcode_idx,
                                bool_value = F)
    lineage_val <- .nonzero_col(mat = barcode_input[["Custom"]],
                                col_idx = barcode_idx,
                                bool_value = T)
    stopifnot(length(lineage_idx) == length(lineage_val))
    
    cbind(lineage_idx,
          rep(j, length(lineage_idx)),
          lineage_val)
  })
  
  tmp <- do.call(rbind, barcode_list)
  barcode_mat <- Matrix::sparseMatrix(i = tmp[,1],
                                      j = tmp[,2],
                                      x = tmp[,3],
                                      dims = c(nrow(barcode_input[["Custom"]]), length(original_barcodes)))
  colnames(barcode_mat) <- original_barcodes
  rownames(barcode_mat) <- rownames(barcode_input[["Custom"]])
  
  barcode_mat
})

sapply(barcode_matrix_list, dim)
barcode_matrix_all <- do.call(cbind, barcode_matrix_list)
all_data[["Lineage"]] <- Seurat::CreateAssayObject(counts = barcode_matrix_all)

########

# also port in the fastTopics
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  load(paste0("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_fasttopics_",
              treatment, ".RData"))
  
  topic_mat <- matrix(NA, nrow = ncol(all_data), ncol = ncol(topic_res$L))
  rownames(topic_mat) <- colnames(all_data)
  topic_res$L <- topic_res$L[rownames(topic_res$L) %in% colnames(all_data),]
  
  for(i in 1:nrow(topic_res$L)){
    if(i %% floor(nrow(topic_res$L)/10) == 0) cat('*')
    
    topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
  }
  colnames(topic_mat) <- paste0("fastTopic", treatment, "_", 1:ncol(topic_mat))
  colnames(topic_res$F) <- paste0("fastTopic", treatment, "_", 1:ncol(topic_mat))
  
  all_data[[paste0("fasttopic_", treatment)]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                                              loadings =  topic_res$F,
                                                                              assay = "RNA",
                                                                              key =  paste0("fastTopic", treatment, "_"))
}

########################

# also port in the spliced and unspliced
treatment_vec <- c("CIS", "COCL2", "DABTRAM")
splicedunspliced_list <- lapply(treatment_vec, function(treatment){
  print(paste0("Working on ", treatment))
  
  load(paste0("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_", treatment, ".RData"))
  list(spliced_count = all_data_subset[["spliced"]]@counts,
       unspliced_count = all_data_subset[["unspliced"]]@counts)
})
sapply(splicedunspliced_list, function(x){
  dim(x$spliced_count)
})
sapply(splicedunspliced_list, function(x){
  head(colnames(x$spliced_count))
})
sum(abs(as.numeric(splicedunspliced_list[[1]]$spliced_count[,1]) - as.numeric(splicedunspliced_list[[2]]$spliced_count[,1])))
sum(abs(as.numeric(splicedunspliced_list[[1]]$spliced_count[,1]) - as.numeric(splicedunspliced_list[[3]]$spliced_count[,1])))

splicedunspliced_list[[2]]$spliced_count <- splicedunspliced_list[[2]]$spliced_count[,-grep("day0", colnames(splicedunspliced_list[[2]]$spliced_count))]
splicedunspliced_list[[3]]$spliced_count <- splicedunspliced_list[[3]]$spliced_count[,-grep("day0", colnames(splicedunspliced_list[[3]]$spliced_count))]
spliced_original <- cbind(splicedunspliced_list[[1]]$spliced_count,
                          splicedunspliced_list[[2]]$spliced_count,
                          splicedunspliced_list[[3]]$spliced_count)
spliced_list <- lapply(1:ncol(all_data), function(j){
  if(j %% floor(ncol(all_data)/10) == 0) cat('*')
  
  barcode <- colnames(all_data)[j]
  barcode_idx <- which(colnames(spliced_original) == barcode)
  if(length(barcode_idx) != 1) return()
  
  lineage_idx <- .nonzero_col(mat = spliced_original,
                              col_idx = barcode_idx,
                              bool_value = F)
  lineage_val <- .nonzero_col(mat = spliced_original,
                              col_idx = barcode_idx,
                              bool_value = T)
  stopifnot(length(lineage_idx) == length(lineage_val))
  
  cbind(lineage_idx,
        rep(j, length(lineage_idx)),
        lineage_val)
})

tmp <- do.call(rbind, spliced_list)
spliced_mat <- Matrix::sparseMatrix(i = tmp[,1],
                                    j = tmp[,2],
                                    x = tmp[,3],
                                    dims = c(nrow(spliced_original), ncol(all_data)))
colnames(spliced_mat) <- colnames(all_data)
rownames(spliced_mat) <- rownames(spliced_original)
all_data[["spliced"]] <- Seurat::CreateAssayObject(counts = spliced_mat)

splicedunspliced_list[[2]]$unspliced_count <- splicedunspliced_list[[2]]$unspliced_count[,-grep("day0", colnames(splicedunspliced_list[[2]]$unspliced_count))]
splicedunspliced_list[[3]]$unspliced_count <- splicedunspliced_list[[3]]$unspliced_count[,-grep("day0", colnames(splicedunspliced_list[[3]]$unspliced_count))]
unspliced_original <- cbind(splicedunspliced_list[[1]]$unspliced_count,
                            splicedunspliced_list[[2]]$unspliced_count,
                            splicedunspliced_list[[3]]$unspliced_count)
unspliced_list <- lapply(1:ncol(all_data), function(j){
  if(j %% floor(ncol(all_data)/10) == 0) cat('*')
  
  barcode <- colnames(all_data)[j]
  barcode_idx <- which(colnames(unspliced_original) == barcode)
  if(length(barcode_idx) != 1) return()
  
  lineage_idx <- .nonzero_col(mat = unspliced_original,
                              col_idx = barcode_idx,
                              bool_value = F)
  lineage_val <- .nonzero_col(mat = unspliced_original,
                              col_idx = barcode_idx,
                              bool_value = T)
  stopifnot(length(lineage_idx) == length(lineage_val))
  
  cbind(lineage_idx,
        rep(j, length(lineage_idx)),
        lineage_val)
})

tmp <- do.call(rbind, unspliced_list)
unspliced_mat <- Matrix::sparseMatrix(i = tmp[,1],
                                      j = tmp[,2],
                                      x = tmp[,3],
                                      dims = c(nrow(unspliced_original), ncol(all_data)))
colnames(unspliced_mat) <- colnames(all_data)
rownames(unspliced_mat) <- rownames(unspliced_original)
all_data[["unspliced"]] <- Seurat::CreateAssayObject(counts = unspliced_mat)

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

###############################

Seurat::DefaultAssay(all_data) <- "Saver"
plot1 <-Seurat::DimPlot(all_data, reduction = "saverumap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (Saver),\n", length(all_data[["Saver"]]@var.features), " genes, using 50 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_saver_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


