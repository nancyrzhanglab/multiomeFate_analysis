rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup10a/"

#######

all_data <- multiomeFate:::data_loader(which_files = c("rna", "saver", "peakvi", "fasttopics", "wnn"))

dimred_vec <- names(all_data@reductions)
dimred_vec <- dimred_vec[grep("umap", dimred_vec)]

pdf(paste0(plot_folder, "Writeup10a_all-umaps.pdf"),
    onefile = TRUE, width = 8, height = 5)

for(kk in 1:length(dimred_vec)){
  plot1 <- scCustomize::DimPlot_scCustom(all_data, 
                                         reduction = dimred_vec[kk],
                                         group.by = "dataset",
                                         colors_use = all_data@misc$dataset_colors)
  plot1 <- plot1 + ggplot2::ggtitle(dimred_vec[kk])
  print(plot1)
}

dev.off()

##########

all_data$keep <- as.numeric(is.na(all_data$assigned_lineage))

table(all_data$dataset, is.na(all_data$assigned_lineage))

dimred_vec <- names(all_data@reductions)
dimred_vec <- dimred_vec[grep("umap", dimred_vec)]

pdf(paste0(plot_folder, "Writeup10a_all-umaps_diagnostic-missing-lineage.pdf"),
    onefile = TRUE, width = 8, height = 5)

for(kk in 1:length(dimred_vec)){
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                             reduction = dimred_vec[kk],
                                             features = "keep")
  plot1 <- plot1 + ggplot2::ggtitle(dimred_vec[kk])
  print(plot1)
}

dev.off()

#####################

treatment_week5_vec <- c("week5_CIS", "week5_COCL2", "week5_DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
min_val <- sapply(1:nrow(tab_mat), function(i){
  min(tab_mat[i,treatment_week5_vec])
})
lineage_name <- rownames(tab_mat)[which.max(min_val)]
umap_res <- all_data[["wnn.umap"]]@cell.embeddings

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  png(paste0(plot_folder, "Writeup10a_ppStep9_lineage-week5_", treatment,".png"),
      height = 1500, width = 1500, res = 300, units = "px")
  
  plot(umap_res[,1], umap_res[,2],
       pch = 16, 
       cex = 1,
       col = "lightgray",
       main = paste0(lineage_name, ": ", treatment),
       xlab = "UMAP 1",
       ylab = "UMAP 2",
       xaxt = "n",
       yaxt = "n",
       bty = "n")
  axis(1)
  axis(2)
  
  idx <- intersect(
    which(all_data$assigned_lineage == lineage_name),
    which(all_data$dataset == "day0")
  )
  
  points(umap_res[idx,1], umap_res[idx,2],
         pch = 16,
         cex = 2.5,
         col = "white")
  
  points(umap_res[idx,1], umap_res[idx,2],
         pch = 16,
         cex = 2,
         col = all_data@misc$dataset_colors["day0"])
  
  
  idx <- intersect(
    which(all_data$assigned_lineage == lineage_name),
    which(all_data$dataset == paste0("day10_", treatment))
  )
  
  points(umap_res[idx,1], umap_res[idx,2],
         pch = 16,
         cex = 2.5,
         col = "white")
  
  points(umap_res[idx,1], umap_res[idx,2],
         pch = 16,
         cex = 2,
         col = all_data@misc$dataset_colors[paste0("day10_", treatment)])
  
  
  idx <- intersect(
    which(all_data$assigned_lineage == lineage_name),
    which(all_data$dataset == paste0("week5_", treatment))
  )
  
  points(umap_res[idx,1], umap_res[idx,2],
         pch = 16,
         cex = 2.5,
         col = "white")
  
  points(umap_res[idx,1], umap_res[idx,2],
         pch = 16,
         cex = 2,
         col = all_data@misc$dataset_colors[paste0("week5_", treatment)])
  
  dev.off()
}

###############

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- log10(tab_mat + 1)

# Create data
data <- data.frame(
  pair = c("day10_CIS to day10_COCL2",
           "day10_CIS to day10_DABTRAM",
           "day10_COCL2 to day10_DABTRAM"), 
  correlation_log10 = c(stats::cor(tab_mat[,"day10_CIS"],
                                   tab_mat[,"day10_COCL2"]),
                        stats::cor(tab_mat[,"day10_CIS"],
                                   tab_mat[,"day10_DABTRAM"]),
                        stats::cor(tab_mat[,"day10_COCL2"],
                                   tab_mat[,"day10_DABTRAM"]))
)

# Barplot
plot1 <- ggplot2::ggplot(data, ggplot2::aes(x=pair, y=correlation_log10)) + 
  ggplot2::geom_bar(stat = "identity") + ggplot2::ylim(c(0,0.75))
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_ppStep9_barplot-correlations_day10.png"),
                plot1,
                height = 4, 
                width = 7,
                units = "in")


library(tidyr)
library(dplyr)
pivot_longer_heatmap <- function(cor_mat){
  diag(cor_mat) <- NA
  # Convert the matrix to a tibble (a type of data frame) for easier manipulation
  mat_df <- tidyr::as_tibble(cor_mat)
  
  # Add row numbers as a new column, since pivot_longer() will melt all existing columns
  mat_df <- as.data.frame(cor_mat) %>%
    dplyr::mutate(Row = rownames(cor_mat))
  
  # Use pivot_longer() to convert from wide to long format
  mat_long <- mat_df %>%
    tidyr::pivot_longer(cols = !Row, 
                        names_to = "Column", 
                        values_to = "Value") %>% 
    dplyr::filter(!is.na(Value)) %>%
    dplyr::mutate(Row = factor(Row, levels = rownames(df_cor)),
                  Column = factor(Column, levels = rownames(df_cor)))
  
  # Filter for the lower triangle (including the diagonal)
  mat_lower_triangle <- mat_long %>%
    dplyr::filter(as.numeric(Row) >= as.numeric(Column)) 
  
  
  return(mat_lower_triangle)
}

# https://stackoverflow.com/questions/51697101/how-to-do-a-triangle-heatmap-in-r-using-ggplot2-reshape2-and-hmisc
df_cor <- stats::cor(as.matrix(tab_mat[,c("day10_CIS", "day10_COCL2", "day10_DABTRAM")]))
df_long <- pivot_longer_heatmap(df_cor)
plot1 <- ggplot2::ggplot(df_long, ggplot2::aes(x = Row, y = Column, fill = Value)) +
  ggplot2::geom_raster() +
  ggplot2::geom_text(ggplot2::aes(label = round(Value, 2)), angle = 225, size = 15) + # Round values for better display
  ggplot2::scale_fill_gradient2(low = "white", high = "gray20", limits = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_ppStep9_heatmap-correlations_day10.png"),
                plot1,
                height = 5, 
                width = 6,
                units = "in")

###

# Create data
data <- data.frame(
  pair = c("week5_CIS to week5_COCL2",
           "week5_CIS to week5_DABTRAM",
           "week5_COCL2 to week5_DABTRAM"), 
  correlation_log10 = c(stats::cor(tab_mat[,"week5_CIS"],
                                   tab_mat[,"week5_COCL2"]),
                        stats::cor(tab_mat[,"week5_CIS"],
                                   tab_mat[,"week5_DABTRAM"]),
                        stats::cor(tab_mat[,"week5_COCL2"],
                                   tab_mat[,"week5_DABTRAM"]))
)

# Barplot
plot1 <- ggplot2::ggplot(data, ggplot2::aes(x=pair, y=correlation_log10)) + 
  ggplot2::geom_bar(stat = "identity") + ggplot2::ylim(c(0,0.75))
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_ppStep9_barplot-correlations_week5.png"),
                plot1,
                height = 4, 
                width = 7,
                units = "in")

# https://stackoverflow.com/questions/51697101/how-to-do-a-triangle-heatmap-in-r-using-ggplot2-reshape2-and-hmisc
df_cor <- stats::cor(as.matrix(tab_mat[,c("week5_CIS", "week5_COCL2", "week5_DABTRAM")]))
df_long <- pivot_longer_heatmap(df_cor)
plot1 <- ggplot2::ggplot(df_long, ggplot2::aes(x = Row, y = Column, fill = Value)) +
  ggplot2::geom_raster() +
  ggplot2::geom_text(ggplot2::aes(label = round(Value, 2)), angle = 225, size = 15) + # Round values for better display
  ggplot2::scale_fill_gradient2(low = "white", high = "gray20", limits = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_ppStep9_heatmap-correlations_week5.png"),
                plot1,
                height = 5, 
                width = 6,
                units = "in")

###

# Create data
data <- data.frame(
  pair = c("day10_CIS to week5_CIS",
           "day10_COCL2 to week5_COCL2",
           "day10_DABTRAM to week5_DABTRAM"), 
  correlation_log10 = c(stats::cor(tab_mat[,"day10_CIS"],
                                   tab_mat[,"week5_CIS"]),
                        stats::cor(tab_mat[,"day10_COCL2"],
                                   tab_mat[,"week5_COCL2"]),
                        stats::cor(tab_mat[,"day10_DABTRAM"],
                                   tab_mat[,"week5_DABTRAM"]))
)

# Barplot
plot1 <- ggplot2::ggplot(data, ggplot2::aes(x=pair, y=correlation_log10)) + 
  ggplot2::geom_bar(stat = "identity") + ggplot2::ylim(c(0,0.75))
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_ppStep9_barplot-correlations_treatment.png"),
                plot1,
                height = 4, 
                width = 7,
                units = "in")

###############################

# Gini plots

unique_lineages <- sort(unique(all_data$assigned_lineage))

lineage_idx_list <- lapply(unique_lineages, function(lineage){
  which(all_data$assigned_lineage == lineage)
})

dataset_colors <- all_data@misc$dataset_colors

treatment_list <- list(CIS = c("day0", "day10_CIS", "week5_CIS"),
                       COCL2 = c("day0", "day10_COCL2", "week5_COCL2"),
                       DABTRAM = c("day0", "day10_DABTRAM", "week5_DABTRAM"))

pdf(paste0(plot_folder, "Writeup10a_gini-index.pdf"),
    onefile = T, width = 5, height = 5)

for(kk in 1:length(treatment_list)){
  print(paste0("Working on: ", names(treatment_list)[kk]))
  treatment_vec <- treatment_list[[kk]]
  
  # count
  cell_treatment_idx_list <- lapply(treatment_vec, function(treatment){
    which(all_data$dataset == treatment)
  })
  
  vec_list <- lapply(cell_treatment_idx_list, function(cell_treatment_idx){
    sapply(lineage_idx_list, function(lineage_idx){
      length(intersect(cell_treatment_idx,
                       lineage_idx))
    })
  })
  
  gini_vec <- sapply(vec_list, function(vec){
    dineq::gini.wtd(vec)
  })
  
  plot(NA,
       xlim = c(0,1),
       ylim = c(0,1),
       asp = TRUE,
       xlab = "Cumulative share of lineages (smallest to largest)",
       ylab = "Cumulative share of cells",
       main = paste0("Gini index for ", names(treatment_list)[kk], ": ",
                     paste0(round(gini_vec, 2), collapse = ", ")))
  for(i in 1:length(treatment_vec)){
    lines(x = seq(0, 1, length.out = length(vec_list[[i]])), 
          y = cumsum(sort(vec_list[[i]], decreasing = FALSE))/sum(vec_list[[i]]),
          col = dataset_colors[treatment_vec[i]],
          lwd = 3)
  }
  
  lines(c(0,1), 
        c(0,1), 
        col = "red",
        lwd = 2, 
        lty = 2)
}

dev.off()

####################

all_data <- multiomeFate:::data_loader(which_files = c("lineage"))

lin_mat <- SeuratObject::LayerData(all_data,
                                   layer = "counts",
                                   assay = "Lineage")
nonzero_val <- sapply(1:ncol(lin_mat), function(i){
  length(multiomeFate:::.nonzero_col(lin_mat, col_idx = i, bool_value = FALSE))
})
nonzero_val <- pmin(nonzero_val, stats::quantile(nonzero_val, prob = 0.99))

png(paste0(plot_folder, "Writeup10a_ppStep5b_lineage-nonzero.png"),
    height = 1000, width = 2400, res = 300, units = "px")
graphics::hist(nonzero_val, 
               breaks = seq(0.5, max(nonzero_val)+1, by = 1),
               xlab = "Number of detected lineage barcodes",
               ylab = "Frequency",
               main = "Multiple lineage barcodes are detected in a cell")
graphics.off()

#######################################
#######################################
#######################################

all_data <- multiomeFate:::data_loader(which_files = c("lineage"))

set.seed(10)
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- log10(tab_mat+1)
tab_mat_jittered <- tab_mat + matrix(
  stats::runif(prod(dim(tab_mat)), min = -0.05, max = 0.05),
  nrow = nrow(tab_mat),
  ncol = ncol(tab_mat)
)
tab_mat_jittered <- pmax(tab_mat_jittered, 0)
max_val <- max(tab_mat) + 0.05

k <- ncol(tab_mat)
pdf(paste0(plot_folder, "Writeup10a_lineage-size_separate.pdf"),
    onefile = T, width = 5, height = 5)
for(i in 1:(k-1)){
  for(j in (i+1):k){
    x_vec <- tab_mat_jittered[,i]
    y_vec <- tab_mat_jittered[,j]
    
    plot(x_vec,
         y_vec,
         xlab = paste0("Size of ", colnames(tab_mat)[i], " (Log10+1, jittered)"),
         ylab = paste0("Size of ", colnames(tab_mat)[j], " (Log10+1, jittered)"),
         main = paste0("Between ", colnames(tab_mat)[i], " and ",
                       colnames(tab_mat)[j], "\n(Cor: ",
                       round(stats::cor(tab_mat[,i], tab_mat[,j]), 2), ")"),
         xlim = c(0, max_val),
         ylim = c(0, max_val),
         asp = TRUE,
         pch = 16,
         col = rgb(0.5, 0.5, 0.5, 0.5))
  }
}
graphics.off()

