rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep3_combine.RData"))

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
file_prefix <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2024_08/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2024_08_21_GEXandLin_time0", "2024_08_21_GEXandLin_day10_CIS", 
                  "2024_08_21_GEXandLin_day10_COCL2", "2024_08_21_GEXandLin_day10_DABTRAM",
                  "2024_08_21_GEXandLin_week5_CIS", "2024_08_21_GEXandLin_week5_COCL2",
                  "2024_08_21_GEXandLin_week5_DABTRAM")
name_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
              "week5_CIS", "week5_COCL2", "week5_DABTRAM")
stopifnot(length(file_folders) == length(name_vec))

barcode_matrix_list <- lapply(1:length(file_folders), function(i){
  print(paste0("Working on ", name_vec[i]))
  file_folder <- file_folders[i]
  barcode_input <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  
  barcode_input <- barcode_input[["Custom"]]
  print(barcode_input[1:5,1:5])
  
  original_barcodes <- Seurat::Cells(all_data)[which(all_data$dataset == name_vec[i])]
  colnames(barcode_input) <- paste0(name_vec[i], "_", colnames(barcode_input))
  intersect_barcodes <- intersect(colnames(barcode_input), original_barcodes)
  
  barcode_list <- lapply(1:length(original_barcodes), function(j){
    if(j %% floor(length(original_barcodes)/10) == 0) cat('*')
    
    barcode <- original_barcodes[j]
    barcode_idx <- which(colnames(barcode_input) == barcode)
    
    if(!barcode %in% intersect_barcodes) return()
    lineage_idx <- .nonzero_col(mat = barcode_input,
                                col_idx = barcode_idx,
                                bool_value = FALSE)
    lineage_val <- .nonzero_col(mat = barcode_input,
                                col_idx = barcode_idx,
                                bool_value = TRUE)
    stopifnot(length(lineage_idx) == length(lineage_val))
    
    cbind(lineage_idx,
          rep(j, length(lineage_idx)),
          lineage_val)
  })
  
  tmp <- do.call(rbind, barcode_list)
  barcode_mat <- Matrix::sparseMatrix(i = tmp[,1],
                                      j = tmp[,2],
                                      x = tmp[,3],
                                      dims = c(nrow(barcode_input), length(original_barcodes)))
  colnames(barcode_mat) <- original_barcodes
  rownames(barcode_mat) <- rownames(barcode_input)
  
  barcode_mat
})

barcode_matrix_all <- do.call(cbind, barcode_matrix_list)
barcode_matrix_all <- barcode_matrix_all[,Seurat::Cells(all_data)]

print(barcode_matrix_all[1:5,1:5])

all_data[["Lineage"]] <- Seurat::CreateAssayObject(counts = barcode_matrix_all)
print(all_data)

print("Saving")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep4_lineage.RData"))

