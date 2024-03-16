rm(list=ls())
library(Seurat)
library(multiomeFate)

load("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step3_fasttopics.RData")

treatment_vec <- sort(unique(all_data$OG_condition))
treatment_vec <- treatment_vec[treatment_vec != "naive"]
day_early <- "naive"

for(treatment in treatment_vec){
  day_treatment <- paste0(day_early, "_", treatment)
  print(day_treatment)
  
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step6_growth-potential_", treatment, "_postprocess.RData"))
  # has the cell_imputed_score
  
  growth_vector <- rep(NA, ncol(all_data))
  names(growth_vector) <- SeuratObject::Cells(all_data)
  growth_vector[names(cell_imputed_score)] <- cell_imputed_score
  
  all_data@meta.data[,paste0(treatment, "_gp")] <- growth_vector
}

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(all_data, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step8_updating-seurat.RData")

gp_df <- all_data@meta.data
gp_df <- gp_df[,grep("^.*gp", colnames(gp_df))]
gp_df <- gp_df[!is.na(gp_df[,1]),]
gp_mat <- as.matrix(gp_df)

gp_umap <- Seurat::RunUMAP(gp_mat)
gp_umap_mat <- gp_umap@cell.embeddings
gp_umap_mat_full <- matrix(NA, nrow = ncol(all_data), ncol = 2)
rownames(gp_umap_mat_full) <- SeuratObject::Cells(all_data)
colnames(gp_umap_mat_full) <- colnames(gp_umap_mat)
gp_umap_mat_full[rownames(gp_umap_mat),] <- gp_umap_mat
all_data[["gp_umap"]] <- Seurat::CreateDimReducObject(gp_umap_mat_full)

plot_list <- lapply(treatment_vec, function(treatment){
  vec <- all_data@meta.data[,paste0(treatment, "_gp")]
  vec <- pmin(vec, quantile(vec, probs = 0.999, na.rm = T))
  all_data@meta.data[,paste0(treatment, "_gp-thres")] <- vec
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                             colors_use = list("red", "lightgray", "blue"),
                                             na_cutoff = quantile(vec, probs = 0.05, na.rm = T),
                                             na_color = "bisque",
                                             reduction = "gp_umap", 
                                             features = paste0(treatment, "_gp-thres"))
  plot1 <- plot1 + ggplot2::ggtitle(paste0(treatment, " GP"))
  plot1
})

plot_all <- cowplot::plot_grid(
  plotlist = plot_list, ncol = 3
)
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_dylan_step8_gp-umap.png",
                plot_all, device = "png", width = 12, height = 6, units = "in")


