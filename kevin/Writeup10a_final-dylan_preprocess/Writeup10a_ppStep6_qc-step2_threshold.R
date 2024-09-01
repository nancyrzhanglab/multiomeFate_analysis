rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc.RData"))

Seurat::DefaultAssay(all_data) <- "RNA"
all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, 
                           verbose = FALSE)

set.seed(10)
all_data <- Seurat::FindNeighbors(all_data, dims = 1:30)
all_data <- Seurat::FindClusters(all_data, resolution = 0.5)
set.seed(10)
all_data <- Seurat::FindClusters(all_data, resolution = 1)

set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            dims = 1:30)

categorical_vars <- c("dataset", "Phase", "high.tss",
                      "RNA_snn_res.0.5",
                      "RNA_snn_res.1")
numerical_vars <- c("nCount_ATAC", "nFeature_ATAC", "nCount_RNA", "nFeature_RNA",
                    "nCount_Lineage", "nFeature_Lineage", "assigned_posterior",
                    "S.Score", "G2M.Score", "percent.mt", "percent.rb",
                    "nucleosome_signal", "nucleosome_percentile",
                    "TSS.enrichment", "TSS.percentile",
                    "blacklist_fraction")

plot_list_categorial <- lapply(categorical_vars, function(variable){
  plot1 <- Seurat::DimPlot(all_data, 
                           reduction = "umap",
                           group.by = variable, 
                           label = TRUE,
                           repel = TRUE, 
                           label.size = 2.5)
  plot1
})

plot_list_numerical <- lapply(numerical_vars, function(variable){
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                             features = variable,
                                             reduction = "umap")
  plot1
})

plot_list <- c(plot_list_categorial, plot_list_numerical)

pdf(paste0(plot_folder, "Writeup10a_exploration_umaps.pdf"),
    width = 6, height = 5, onefile = TRUE)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
graphics.off()

###############################

# compute the average percent.rb and nFeature_RNA for each cluster
cluster_names <- sort(unique(all_data$RNA_snn_res.1))
tmp <- sapply(cluster_names, function(cluster){
  idx <- which(all_data$RNA_snn_res.1 == cluster)
  c(mean(all_data$percent.rb[idx]), 
    mean(all_data$nFeature_RNA[idx]))
})
colnames(tmp) <- cluster_names
rownames(tmp) <- c("percent.rb", "nFeature_RNA")
round(tmp,1)

# remove clusters
keep_vec <- rep(TRUE, length(Seurat::Cells(all_data)))
keep_vec[which(all_data$RNA_snn_res.1 %in% c("7", "10", "14", "18", "19", "22"))] <- FALSE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

# recompute the visualization
set.seed(10)
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, 
                           verbose = FALSE)
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            dims = 1:30)

categorical_vars <- c("dataset", "Phase", "high.tss",
                      "RNA_snn_res.0.5",
                      "RNA_snn_res.1")
numerical_vars <- c("nCount_ATAC", "nFeature_ATAC", "nCount_RNA", "nFeature_RNA",
                    "nCount_Lineage", "nFeature_Lineage", "assigned_posterior",
                    "S.Score", "G2M.Score", "percent.mt", "percent.rb",
                    "nucleosome_signal", "nucleosome_percentile",
                    "TSS.enrichment", "TSS.percentile",
                    "blacklist_fraction")

plot_list_categorial <- lapply(categorical_vars, function(variable){
  plot1 <- Seurat::DimPlot(all_data, 
                           reduction = "umap",
                           group.by = variable, 
                           label = TRUE,
                           repel = TRUE, 
                           label.size = 2.5)
  plot1
})

plot_list_numerical <- lapply(numerical_vars, function(variable){
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                             features = variable,
                                             reduction = "umap")
  plot1
})

plot_list <- c(plot_list_categorial, plot_list_numerical)

pdf(paste0(plot_folder, "Writeup10a_exploration_umaps_after-qc.pdf"),
    width = 6, height = 5, onefile = TRUE)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
graphics.off()

all_data
table(all_data$dataset)
length(unique(all_data$assigned_lineage))


print("Finished preprocessing data")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep6_qc-step2.RData"))

all_data_rna <- all_data[["RNA"]]
all_data_pca <- all_data[["pca"]]
all_data_umap <- all_data[["umap"]]
all_data_atac <- all_data[["ATAC"]]
all_data_lineage <- all_data[["Lineage"]]
all_data_metadata <- all_data@meta.data

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(all_data_rna, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_rna.RData"))
save(all_data_pca, all_data_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_rna_dimred.RData"))
save(all_data_atac, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_atac.RData"))
save(all_data_lineage, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_lineage.RData"))
save(all_data_metadata, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_metadata.RData"))

print("Done! :)")
