rm(list=ls())
library(Seurat)
source("../Writeup3d/funcs.R")
source("select_cells.R")

load("../../../../out/kevin/Writeup3e/10192021_sydney_preprocess.RData")

tabulate_mat <- .tabulate_lineages(all_data)
max_val <- sapply(1:nrow(tabulate_mat), function(i){max(tabulate_mat[i,-5])})
tabulate_mat2 <- tabulate_mat[order(max_val, decreasing = T),]
head(tabulate_mat2)

##############

naive_mat <- all_data[["RNA"]]@counts[,which(all_data@meta.data[,"Original_condition"] == "naive")]
naive_data <- Seurat::CreateSeuratObject(counts = naive_mat)
Seurat::DefaultAssay(naive_data) <- "RNA"
naive_data[["Original_condition"]] <- rep("naive", ncol(naive_data))
lin_mat <-all_data[["lineage"]]@counts
lin_mat <- lin_mat[,which(colnames(lin_mat) %in% colnames(naive_data))]
naive_data[["lineage"]] <- Seurat::CreateAssayObject(counts = lin_mat)

# naive_data2 <- .collapse_by_lineage(naive_data, assay = "RNA", slot = "counts")
# naive_data2[["Original_condition"]] <- rep("naive", ncol(naive_data2))
# Seurat::DefaultAssay(naive_data2) <- "RNA"
# naive_data2 <- Seurat::SCTransform(naive_data2, clip.range = c(-30,30))
# quantile(naive_data2[["SCT"]]@scale.data)
# naive_data2 <- Seurat::RunPCA(naive_data2, verbose = F)
# set.seed(10)
# naive_data2 <- Seurat::RunUMAP(naive_data2, dims = 1:50)
# Seurat::DefaultAssay(naive_data2) <- "RNA"
# naive_data2 <- Seurat::NormalizeData(naive_data2)
# naive_data2 <- Seurat::ScaleData(naive_data2, 
#                                  features = naive_data2[["SCT"]]@var.features,
#                                  scale.max = 30)

##################

# gene_vec <- c("EGFR", "NGFR", "NT5E")
# which(gene_vec %in% naive_data2[["SCT"]]@var.features)
# 
# treatment_vec <- c("dab", "tram")
# 
# # when looking at lineages
# for(treatment in treatment_vec){
#   res <- .select_expansion_naives(naive_data2, 
#                                   tabulate_mat, 
#                                   treatment = treatment,
#                                   threshold = 3, 
#                                   type = 2)
#   naive_deathtoall <- .select_naive_deathtoall(naive_data2,
#                                                tabulate_mat,
#                                                threshold = 1, 
#                                                type = 2)
#   newgroup <- rep("none", ncol(naive_data2))
#   names(newgroup) <- rownames(naive_data2@meta.data)
#   newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
#   newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
#   print(table(newgroup))
#   naive_data2[[paste0(treatment, "Status")]] <- newgroup
#   Seurat::Idents(naive_data2) <- paste0(treatment, "Status")
#   
#   marker_genes <- c("EGFR", "NGFR", "NT5E")
#   plot_suffix <- c("sctransform-scaledata", "log1p", "lognormalized")
#   assay_vec <- c("SCT", "SCT", "RNA")
#   slot_vec <- c("scale.data", "data", "scale.data")
#   
#   for(i in 1:length(marker_genes)){
#     idx1 <- which(Seurat::Idents(naive_data2) == paste0("naive_no", treatment))
#     idx2 <- which(Seurat::Idents(naive_data2) == paste0("naive_survive", treatment))
#     
#     for(j in 1:length(plot_suffix)){
#       plot1 <- Seurat::VlnPlot(naive_data2, features = marker_genes[i],
#                                idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
#                                assay = assay_vec[j],
#                                slot = slot_vec[j],
#                                pt.size = 1.5) + ggplot2::theme_classic()
#       plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
#         paste0("naive - will not survive\n(", length(idx1), " lineages)"),
#         paste0("naive - will survive\n(", length(idx2), " lineages)")))
#       plot1 <- plot1 + Seurat::NoLegend()
#       ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10272021_sydney_naive_within", treatment, "_vln_", marker_genes[i], "_avglineage_", plot_suffix[j], ".png"),
#                       plot1, device = "png", width = 6, height = 5, units = "in")
#       
#     }
#   }
# }

########################################

Seurat::DefaultAssay(naive_data) <- "RNA"
naive_data <- Seurat::SCTransform(naive_data, clip.range = c(-30,30))
Seurat::DefaultAssay(naive_data) <- "RNA"
naive_data <- Seurat::NormalizeData(naive_data)
naive_data <- Seurat::ScaleData(naive_data, 
                                 features = naive_data[["SCT"]]@var.features,
                                 scale.max = 30)
Seurat::DefaultAssay(naive_data) <- "SCT"
naive_data <- Seurat::RunPCA(naive_data, verbose = F)
set.seed(10)
naive_data <- Seurat::RunUMAP(naive_data, dims = 1:50)

gene_vec <- c("EGFR", "NGFR", "NT5E", "FN1", "AXL")
which(gene_vec %in% naive_data[["SCT"]]@var.features)

treatment_vec <- c("dab", "tram")

# when looking at lineages
for(treatment in treatment_vec){
  res <- .select_expansion_naives(naive_data, 
                                  tabulate_mat, 
                                  treatment = treatment,
                                  threshold = 3, 
                                  type = 1)
  naive_deathtoall <- .select_naive_deathtoall(naive_data,
                                               tabulate_mat,
                                               threshold = 1, 
                                               type = 1)
  newgroup <- rep("none", ncol(naive_data))
  names(newgroup) <- rownames(naive_data@meta.data)
  newgroup[res$naive_terminal] <- paste0("naive_survive", treatment)
  newgroup[naive_deathtoall] <- paste0("naive_no", treatment)
  print(table(newgroup))
  naive_data[[paste0(treatment, "Status")]] <- newgroup
  Seurat::Idents(naive_data) <- paste0(treatment, "Status")

  plot_suffix <- c("sctransform-scaledata", "log1p", "lognormalized")
  assay_vec <- c("SCT", "SCT", "RNA")
  slot_vec <- c("scale.data", "data", "scale.data")
  
  for(i in 1:length(gene_vec)){
    idx1 <- which(Seurat::Idents(naive_data) == paste0("naive_no", treatment))
    idx2 <- which(Seurat::Idents(naive_data) == paste0("naive_survive", treatment))
    
    plot1 <- Seurat::FeaturePlot(naive_data,
                                 features = gene_vec[i],
                                 reduction = "umap")
    plot1 <- plot1 + ggplot2::ggtitle(paste0(treatment, ": ", gene_vec[i]))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10272021_sydney_naive_within", treatment, "_umap_", gene_vec[i], "_singlecell.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
    
    col_vec <- c(rev(scales::hue_pal()(2)), "gray")
    names(col_vec) <- c(paste0("naive_no", treatment), paste0("naive_survive", treatment), "none")
    plot1 <- Seurat::DimPlot(naive_data, group.by = paste0(treatment, "Status"),
                             reduction = "umap", 
                             cols = col_vec
    )
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    plot1 <- plot1 + ggplot2::ggtitle(paste0(treatment, ": Cell fate"))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10272021_sydney_naive_within", treatment, "_umap_fate.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
    
    
    for(j in 1:length(plot_suffix)){
      plot1 <- Seurat::VlnPlot(naive_data, features = gene_vec[i],
                               idents = c(paste0("naive_no", treatment), paste0("naive_survive", treatment), treatment),
                               assay = assay_vec[j],
                               slot = slot_vec[j],
                               pt.size = 1.5) + ggplot2::theme_classic()
      plot1 <- plot1 + ggplot2::scale_x_discrete(labels = c(
        paste0("naive - will not survive\n(", length(idx1), " cells)"),
        paste0("naive - will survive\n(", length(idx2), " cells)")))
      plot1 <- plot1 + Seurat::NoLegend()
      ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10272021_sydney_naive_within", treatment, "_vln_", gene_vec[i], "_singlecell_", plot_suffix[j], ".png"),
                      plot1, device = "png", width = 4, height = 5, units = "in")
    }
  }
}


