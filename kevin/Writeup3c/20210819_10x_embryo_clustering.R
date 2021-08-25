rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

mbrain <- readRDS("../../../../data/Embrian_mouse/obj_seurat.rds")

# port over labels
mbrain_jane <- readRDS("../../../../data/Embrian_mouse/data_tenx_labels_Jane.rds")
cell_id <- rownames(mbrain@meta.data)
idx <- sapply(1:length(cell_id), function(x){
  which(rownames(mbrain_jane@meta.data) == cell_id[x])[1]
})
mbrain@meta.data$celltype <- mbrain_jane$savercatLable[idx]

dim(mbrain[["SCT"]]@data[Seurat::VariableFeatures(mbrain, assay = "SCT"),])
dim(mbrain[["ATAC"]]@data)

###############################

# select the desired cells
idx <- which(mbrain@meta.data$celltype %in% c("Oligodendrocyte", "Radial glia", 
                                              "Forebrain GABAergic", "Neuroblast", 
                                              "Glioblast", "Cortical or hippocampal glutamatergic"))

mat_x <- mbrain[["ATAC"]]@data[,idx]
mat_y <- mbrain[["RNA"]]@data[,idx]

#################################

# start a new Seurat object and preprocess as usual

mbrain2 <- Seurat::CreateSeuratObject(counts = mat_y)
Seurat::DefaultAssay(mbrain2) <- "RNA"
mbrain2[["celltype"]] <- mbrain@meta.data$celltype[idx]
set.seed(10)
mbrain2 <- Seurat::SCTransform(mbrain2, verbose = T)
mbrain2 <- Seurat::FindVariableFeatures(mbrain2)
mbrain2 <- Seurat::RunPCA(mbrain2, verbose = FALSE)
set.seed(10)
mbrain2 <- Seurat::RunUMAP(mbrain2, dims = 1:50)

####################

mbrain2[["ATAC"]] <- Seurat::CreateAssayObject(counts = mat_x)
Seurat::DefaultAssay(mbrain2) <- "ATAC"
mbrain2 <- Signac::RunTFIDF(mbrain2)
mbrain2 <-  Signac::FindTopFeatures(mbrain2, min.cutoff="q10")
mbrain2 <-  Signac::RunSVD(mbrain2)  # question: is this svd run with only the top variable features?
set.seed(10)
mbrain2 <- Seurat::RunUMAP(mbrain2, reduction="lsi", dims=2:50, reduction.name="umap.atac", 
                             reduction.key="atacUMAP_")

######################

set.seed(10)
mbrain2 <- Seurat::FindMultiModalNeighbors(
  mbrain2, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:50), modality.weight.name = "RNA.weight"
)
set.seed(10)
mbrain2 <- Seurat::RunUMAP(mbrain2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
set.seed(10)
mbrain2 <- Seurat::FindClusters(mbrain2, graph.name = "wsnn", algorithm = 3, verbose = T)

save(mbrain2, file = "../../../../out/kevin/Writeup3c/10x_mbrain_subset.RData")

###################

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nRNA")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3b_10x_embryo_umap_seurat_rna.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nRNA (Seurat clusters)")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3b_10x_embryo_umap_seurat_rna_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap.atac", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nATAC")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3b_10x_embryo_umap_seurat_atac.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nATAC (Seurat clusters)")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3b_10x_embryo_umap_seurat_atac_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain2, reduction = "wnn.umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nWNN")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3b_10x_embryo_umap_seurat_wnn.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nWNN (Seurat clusters)")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3b_10x_embryo_umap_seurat_wnn_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")



################################

# now we want to remove certain cells in "Forebrain GABAergic" or "Cortical or hippocampal glutamatergic"
# NOTE: be sure to check the specific clusters prior to removing these clusters
#   I've noticed that sometimes setting the seed isn't enough to guarantee the same exact cluster-labelings

cell_idx <- which(!mbrain2@meta.data$seurat_clusters %in% c(5, 7, 10, 12, 13, 14))

#################################

# now find the relevant genes/peaks
# NOTE: replace this with DE analyses

getRegion <- function(vicinity.size, TSS){
  upstream.size <- vicinity.size[1]
  downstream.size <- vicinity.size[2]
  if (upstream.size == 0) {
    upstream <- TSS
  } else if (upstream.size > 0) {
    upstream <- GenomicRanges::flank(TSS, width=upstream.size, start=T) # upstream excluding the TSS
  }
  if (downstream.size == 0) {
    downstream <- TSS
  } else if (downstream.size > 0) {
    downstream <- GenomicRanges::flank(TSS, width=downstream.size, start=F) # downstream excluding the TSS
  }
  vicinity <- GenomicRanges::punion(upstream, downstream, fill.gap=T)
  BiocGenerics::start(vicinity) <- pmax(0, BiocGenerics::start(vicinity))
  vicinity$gene_name <- TSS$gene_name
  return(vicinity)
}

###

gene_names <- rownames(mbrain2[["SCT"]]@data[Seurat::VariableFeatures(mbrain2, assay = "SCT"),])
p2 <- length(gene_names)
# NOTE: we're using mbrain here since we need to use a ChromatinAssay object
peaks_gr <- GenomicRanges::granges(mbrain@assays$ATAC)
annotations <- mbrain@assays$ATAC@annotation

relevant_peaks <- lapply(1:p2, function(j){
  if(j %% floor(p2/10) == 0) cat('*')
  
  temp <- Signac::LookupGeneCoords(mbrain, gene_names[j])
  if(is.null(temp)) return(numeric(0))
  GenomeInfoDb::seqlevelsStyle(temp)="UCSC"  # just so that the seqlevelsStyle matches up with peaks_gr
  temp2 <- annotations[which(annotations$gene_name==gene_names[j] & annotations$type=="cds")] #cds stands for "coding sequence"
  if(length(temp2)==0) return(numeric(0))
  gene_strand <- BiocGenerics::strand(temp2)[1]
  TSS_position <- ifelse(gene_strand == "+", BiocGenerics::start(temp), BiocGenerics::end(temp))
  TSS <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(temp),
                                ranges = IRanges::IRanges(start = TSS_position, width = 1),
                                strand = BiocGenerics::strand(temp),
                                gene_name = gene_names[j])
  upstream_gr <- getRegion(c(50000, 0), TSS)
  downstream_gr <- getRegion(c(0, 50000), TSS)
  
  # populate the peaks counts array.
  peaks_upstream <- GenomicRanges::countOverlaps(peaks_gr, upstream_gr)>0
  peaks_downstream <- GenomicRanges::countOverlaps(peaks_gr, downstream_gr)>0
  
  c(sort(unique(c(which(peaks_upstream), which(peaks_downstream)))))
})

# reparameterize the data
gene_reparam_mat <- data.frame(org_idx = which(sapply(relevant_peaks, length) > 0))
gene_reparam_mat$new_idx <- 1:nrow(gene_reparam_mat)

peak_reparam_mat <- data.frame(org_idx = sort(unique(unlist(relevant_peaks))))
peak_reparam_mat$new_idx <- 1:nrow(peak_reparam_mat)

relevant_peaks2 <- relevant_peaks[gene_reparam_mat$org_idx]
for(i in 1:length(relevant_peaks2)){
  relevant_peaks2[[i]] <- peak_reparam_mat$new_idx[which(peak_reparam_mat$org_idx %in% relevant_peaks2[[i]])]
}


save(mbrain2, gene_reparam_mat, peak_reparam_mat, relevant_peaks2, 
     file = "../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed.RData")

############################3

# to extract the relevant objects:
mat_y <- Matrix::t(mbrain2[["SCT"]]@data[Seurat::VariableFeatures(mbrain2, assay = "SCT"),])[cell_idx,gene_reparam_mat$org_idx]
mat_x <- Matrix::t(mbrain2[["ATAC"]]@data)[cell_idx,peak_reparam_mat$org_idx]
dim(mat_x); dim(mat_y)
p1 <- ncol(mat_x); p2 <- ncol(mat_y)
ht_map <- hash::hash()
for(i in 1:p2){
  ht_map[[as.character(i)]] <- relevant_peaks2[[i]]
}
