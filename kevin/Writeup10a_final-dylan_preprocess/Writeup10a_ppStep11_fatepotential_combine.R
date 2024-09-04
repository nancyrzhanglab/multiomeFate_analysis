rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"

all_data <- multiomeFate:::data_loader(which_files = c("lineage"))

file_vec <- c(CIS_d0_d10 = "Writeup10a_CIS-from-d0_fatepotential.RData",
              CIS_d10_w5 = "Writeup10a_CIS-from-d10_fatepotential.RData",
              COCL2_d0_d10 = "Writeup10a_COCL2-from-d0_fatepotential.RData",
              COCL2_d10_w5 = "Writeup10a_COCL2-from-d10_fatepotential.RData",
              DABTRAM_d0_d10 = "Writeup10a_DABTRAM-from-d0_fatepotential.RData",
              DABTRAM_d10_w5 = "Writeup10a_DABTRAM-from-d10_fatepotential.RData")

for(kk in 1:length(file_vec)){
  filepath <- file_vec[kk]
  filename <- names(file_vec)[kk]
  
  print(paste0("Working on ", filename))
  # load the fate potential
  load(paste0(out_folder, filepath))
  
  # make the training-testing plots
  train_mat <- sapply(fit_res, function(x){
    x$train_loglik
  })
  
  test_mat <- sapply(fit_res, function(x){
    x$test_loglik
  })
  
  train_quantile <- apply(train_mat, 1, function(vec){stats::quantile(vec, probs = c(0.1, 0.5, 0.9))})
  test_quantile <- apply(test_mat, 1, function(vec){stats::quantile(vec, probs = c(0.1, 0.5, 0.9))})
  
  lambda_sequence <- fit_res[[1]]$train_fit$lambda_sequence
  lambda <- lambda_sequence[which.min(test_quantile[2,])]
  
  lambda_sequence2 <- lambda_sequence+1
  max_base10 <- max(ceiling(log10(lambda_sequence+1)))
  xaxt_values <- unlist(lapply(1:max_base10, function(x){
    seq(10^(x-1), 10^x, length.out = 10)
  }))
  xaxt_values <- xaxt_values[!duplicated(xaxt_values)]
  
  png(paste0(plot_folder, "Writeup10a_", filename, "_train-test-curve.png"),
      height = 1500, width = 3000, units = "px", res = 300)
  par(mfrow = c(1,2), mar = c(4,4,4,4))
  plot(lambda_sequence2, train_quantile[2,], 
       log = "x", type = "n",
       main = paste0(filename, " growth potential (Training)"),
       xlab = "Lambda+1 (Log-scale tickmarks)", 
       ylab = "Negative loglikelihood (Training)",
       ylim = range(c(train_quantile[1,], rev(train_quantile[3,]))))
  polygon(x = c(lambda_sequence2, rev(lambda_sequence2)),
          y = c(train_quantile[1,], rev(train_quantile[3,])),
          col = "gray", 
          border = "black")
  points(lambda_sequence2, train_quantile[2,], pch = 16)
  lines(lambda_sequence2, train_quantile[2,], lwd = 2)
  
  test_upper <- test_quantile[3,]
  test_upper <- pmin(test_upper, stats::quantile(test_upper, probs = 0.95))
  test_lower <- test_quantile[1,]
  test_upper <- pmax(test_upper, stats::quantile(test_upper, probs = 0.05))
  
  plot(lambda_sequence2, test_quantile[2,], 
       log = "x", type = "n",
       main = paste0(filename, " growth potential (Testing)\nLambda = ", 
                     round(lambda,3)),
       xlab = "Lambda+1 (Log-scale tickmarks)", 
       ylab = "Negative loglikelihood (Testing)",
       ylim = range(c(test_upper, rev(test_lower))))
  polygon(x = c(lambda_sequence2, rev(lambda_sequence2)),
          y = c(test_upper, rev(test_lower)),
          col = "gray", 
          border = "black")
  points(lambda_sequence2, test_quantile[2,], pch = 16)
  lines(lambda_sequence2, test_quantile[2,], lwd = 2)
  lines(rep(lambda+1, 2), c(-1e6, 1e6), col = "red", lty = 2)
  
  graphics.off()
  
  # preparing the objects to be in seurat format
  fate_vec <- final_fit$cell_imputed_score
  cellnames <- Seurat::Cells(all_data)
  full_vec <- rep(NA, length(cellnames))
  names(full_vec) <- cellnames
  full_vec[names(fate_vec)] <- fate_vec
  all_data@meta.data[,paste0("fatepotential_", filename)] <- full_vec
  
  # put it into misc
  all_data@misc <- c(all_data@misc,
                     list(final_fit))
  names(all_data@misc)[length(all_data@misc)] <- paste0("fatepotential_", filename)
}

# save all the fate potential miscs and metadata
date_of_run <- Sys.time()
session_info <- devtools::session_info()

all_data_fatepotential <- all_data@misc
all_data_fatepotential$dataset_colors <- NULL

save(all_data_fatepotential, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_fatepotential.RData"))

# Updating the Writeup10a_data_empty.RData
all_data_tmp <- all_data
load(paste0(out_folder, "Writeup10a_data_empty.RData"))
metadata_full <- all_data_tmp@meta.data
metadata_empty <- all_data@meta.data
cellnames_full <- Seurat::Cells(all_data_tmp)
cellnames_empty <- Seurat::Cells(all_data)

var_setdiff <- setdiff(colnames(metadata_full), colnames(metadata_empty))
for(variable in var_setdiff){
  vec <- rep(NA, length(cellnames_empty))
  names(vec) <- cellnames_empty
  vec[rownames(metadata_full)] <- metadata_full[,variable]
  metadata_empty[,variable] <- vec
}

all_data@meta.data <- metadata_empty

save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_empty.RData"))




