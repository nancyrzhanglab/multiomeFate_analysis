rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(ggplot2) 

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

n <- ncol(all_data)
keep_vec <- rep(FALSE, n)
keep_vec[which(all_data$dataset %in% c("day0", "week5_COCL2"))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

motif_focus <- "FOS"
motif_name_vec <- names(data.use)
motif_idx <- grep(motif_focus, motif_name_vec)

chromvar_mat <- all_data[["chromvar"]]@data
rownames(chromvar_mat) <- motif_name_vec
chromvar_mat <- chromvar_mat[motif_idx,]

status_vec <- all_data$dataset

# let's try a dimension-reduction
seurat_obj <- Seurat::CreateSeuratObject(counts = chromvar_mat, 
                                         meta.data = data.frame(status_vec))
seurat_obj[["RNA"]]@var.features <- rownames(chromvar_mat)
seurat_obj <- Seurat::ScaleData(seurat_obj)
set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj, verbose = F, dims = 1:5)

plot1 <-Seurat::DimPlot(seurat_obj, reduction = "umap",
                        group.by = "status_vec", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Day0 and COCL week5 cells:\nChromvar (original) scores for ", motif_focus, " motifs"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6m/Writeup6m_day0-vs-COCL2week5_chromvar_", motif_focus, "_umap.png"), 
                plot1, device = "png", width = 7, height = 5, units = "in")






