rm(list=ls())

library(Seurat)
library(Signac)
library(fastTopics)
library(ggplot2)
library(dplyr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")

gene_list <- list(housekeeping = read.csv("~/project/eSVD/data/housekeeping/housekeeping.txt", header = F)[,1],
                  cellcycle =  c(cc.genes$s.genes, cc.genes$g2m.genes),
                  jackpot = c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                              "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                              "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                              "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                              "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"))

for(treatment in treatment_vec){
  print(treatment)
  load(paste0("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_fasttopics_", treatment, ".RData"))
  load(paste0("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver_onlyRNA_", treatment, "-subset.RData"))
  
  print("Plotting")
  # UMAPs
  Seurat::DefaultAssay(all_data_subset) <- "Saver"
  plot1 <-Seurat::DimPlot(all_data_subset, reduction = "saverumap",
                          group.by = "dataset", label = TRUE,
                          repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (SAVER, ", treatment, ")"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  plot1 <- Seurat::FeaturePlot(all_data_subset, 
                               reduction = "saverumap",
                               slot = "scale.data",
                               features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E", 
                                                 "CD44", "LOXL2", "ID3")),
                               ncol = 3)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap_jackpot1.png"),
                  plot1, device = "png", width = 12, height = 12, units = "in")
  
  # Density
  dataset_vec <- tolower(all_data_subset$dataset)
  slice_vec <- sapply(all_data_subset$dataset, function(x){
    if(x %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment))) {
      "barcoded"
    } else {"test"}})
  df <- data.frame(all_data_subset[["saverumap"]]@cell.embeddings, dataset_vec, dataset_vec, slice_vec)
  colnames(df) = c("umap1", "umap2", "dataset", "dataset2", "slice")
  df <- tibble::as_tibble(df)
  df_median <- df %>% dplyr::group_by(dataset2) %>%
    dplyr::summarize(umap1 = stats::median(umap1), umap2 = stats::median(umap2))
  plot1 <- ggplot2::ggplot(subset(df, slice == "barcoded"), ggplot2::aes(x = umap1, y = umap2)) 
  plot1 <- plot1 + ggplot2::stat_density_2d(aes(fill = ..ndensity..), 
                                            geom = "raster", 
                                            contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option="inferno") 
  plot1 <- plot1 + ggplot2::geom_text(df_median, mapping = ggplot2::aes(label = dataset2), size = 2, color = "white",
                                      nudge_y = 0, nudge_x = -1)
  plot1 <- plot1 + ggplot2::geom_density_2d(size = 0.3, alpha = 0.8) 
  plot1 <- plot1 + ggplot2::facet_grid(slice ~ dataset)
  plot1 <- plot1 + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap_density-barcoded.png"),
                  plot1, device = "png", width = 15, height = 5, units = "in")
  
  plot1 <- ggplot2::ggplot(subset(df, slice == "test"), ggplot2::aes(x = umap1, y = umap2)) 
  plot1 <- plot1 + ggplot2::stat_density_2d(aes(fill = ..ndensity..), 
                                            geom = "raster", 
                                            contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option="inferno") 
  plot1 <- plot1 + ggplot2::geom_text(df_median, mapping = ggplot2::aes(label = dataset2), size = 2, color = "white",
                                      nudge_y = 0, nudge_x = -1)
  plot1 <- plot1 + ggplot2::geom_density_2d(size = 0.3, alpha = 0.8) 
  plot1 <- plot1 + ggplot2::facet_grid(slice ~ dataset)
  plot1 <- plot1 + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap_density-test.png"),
                  plot1, device = "png", width = 10, height = 5, units = "in")
  
  # Topics
  all(rownames(topic_res$L) == colnames(all_data_subset))
  all_data_subset[["fasttopic"]] <- Seurat::CreateDimReducObject(topic_res$L, 
                                                                 assay = "RNA",
                                                                 key = "fastTopic_")
  plot1 <- Seurat::FeaturePlot(all_data_subset, 
                               reduction = "saverumap",
                               features = paste0("fastTopic_", 1:ncol(topic_res$L)),
                               ncol = 6)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap_fasttopics.png"),
                  plot1, device = "png", width = 20, height = 15, units = "in")
  
  # Violin plots
  Seurat::Idents(all_data_subset) <- "dataset"
  plot1 <- Seurat::VlnPlot(all_data_subset, 
                           features = paste0("fastTopic_", 1:ncol(topic_res$L)),
                           ncol = 6,
                           pt.size = 0)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_fasttopics-violin.png"),
                  plot1, device = "png", width = 20, height = 15, units = "in")
  
  # Custom heatmap
  col_vec <- viridis::viridis(25)
  break_vec <- seq(0, 1, length.out = 26)
  for(i in 1:3){
    mat <- topic_res$F
    idx <- which(rownames(mat) %in% gene_list[[i]])
    mat <- mat[idx,]
    rowname_vec <- rownames(mat)
    l1_vec <- apply(mat, 1, sum)
    mat <- diag(1/l1_vec) %*% mat
    rownames(mat) <- rowname_vec
    
    hclust_res <- stats::hclust(stats::dist(mat))
    mat <- mat[hclust_res$order,]
    
    ratio <- min(2000*(nrow(mat)-1)/(ncol(mat)-1), 6500)/(2000)
    png(paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_fasttopics-scores_", names(gene_list)[i], ".png"),
        height = 2000*ratio, width = 2000, res = 500, units = "px")
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    image(tiltedCCA:::.rotate(mat), 
          asp = ratio, ylim = c(0,1),
          main = "", xlab = "", ylab = "",
          xaxt = "n", yaxt = "n", bty = "n",
          breaks = break_vec, col = col_vec)
    
    # draw lines
    x_indent_val <- 1/(2*(ncol(mat)-1))
    x_vec <- seq(-x_indent_val, 1+x_indent_val, by = 2*x_indent_val)
    for(x in x_vec[seq(31,6,by=-5)[-1]]){
      graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 2, col = "white")
      graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 1.5, lty = 2)
    }
    
  
    # label genes
    if(names(gene_list)[i] != "housekeeping"){
      y_indent_val <- 1/(2*(nrow(mat)-1))
      y_vec <- seq(-y_indent_val, 1+y_indent_val, by = 2*y_indent_val)
      for(y in y_vec[seq(6,length(y_vec), by = 5)]){
        graphics::lines(x = c(0,1), y = rep(y, 2), lwd = 2, col = "white")
      }
      
      x_vec_label <- rep(x_vec[c(3,5)], times = ceiling(nrow(mat)/2))[1:nrow(mat)]
      y_vec_label <- seq(1, 0, by = -2*y_indent_val)
      graphics::text(x = x_vec_label, y = y_vec_label, labels = rownames(mat),
                     col = "white", cex = 0.5)
    }
    
    graphics.off()
  }
  
  # Genes in each topic
  mat <- topic_res$F
  break_vec <- seq(-5, 0, length.out = 50)
  gene_list2 <- lapply(gene_list, function(x){x[x %in% rownames(mat)]})
  col_vec <- c(4,3,2)
  idx_list <-  lapply(gene_list2, function(x){which(rownames(mat) %in% x)})
  png(paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_fasttopics-geneweights.png"),
      height = 4000, width = 5000, res = 500, units = "px")
  par(mfrow = c(5,6), mar = c(4,2,1,0.5))
  for(i in 1:ncol(mat)){
    hist_obj <- hist(pmax(-5, log10(mat[,i])), 
                     breaks = break_vec, 
                     plot = F)
    hist_obj$counts <- log10(hist_obj$counts+1)
    plot(hist_obj, 
         xlab = "", 
         ylab = "", 
         main = paste0("Topic ",i), col = "gray")
    for(k in 1:3){
      rug(x = pmax(-5, log10(mat[idx_list[[k]],i])), col = col_vec[k], lwd = 2)
    }
  }
  graphics.off()
  
  # Cell-cycling phases
  plot1 <-Seurat::DimPlot(all_data_subset, reduction = "saverumap",
                          group.by = "Phase", label = TRUE,
                          repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (SAVER, ", treatment, "),\nCell-cycle phase"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap-cellcycle.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  phase_vec <- tolower(all_data_subset$Phase)
  slice_vec <- rep("main", length(phase_vec))
  df <- data.frame(all_data_subset[["saverumap"]]@cell.embeddings, phase_vec, phase_vec, slice_vec)
  colnames(df) = c("umap1", "umap2", "phase", "phase2", "slice")
  df <- tibble::as_tibble(df)
  df_median <- df %>% dplyr::group_by(phase2) %>%
    dplyr::summarize(umap1 = stats::median(umap1), umap2 = stats::median(umap2))
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = umap1, y = umap2)) 
  plot1 <- plot1 + ggplot2::stat_density_2d(aes(fill = ..ndensity..), 
                                            geom = "raster", 
                                            contour = FALSE) +
    ggplot2::scale_fill_viridis_c(option="inferno") 
  plot1 <- plot1 + ggplot2::geom_text(df_median, mapping = ggplot2::aes(label = phase2), size = 2, color = "white",
                                      nudge_y = 0, nudge_x = -1)
  plot1 <- plot1 + ggplot2::geom_density_2d(size = 0.3, alpha = 0.8) 
  plot1 <- plot1 + ggplot2::facet_grid(slice ~ phase)
  plot1 <- plot1 + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_saver-", treatment, "_umap_cellcycle-density.png"),
                  plot1, device = "png", width = 15, height = 5, units = "in")
}
