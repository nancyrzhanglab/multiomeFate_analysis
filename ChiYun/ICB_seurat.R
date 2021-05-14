library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(Matrix)


dir.create("./plots")

## load data

rna_counts=readRDS("~/project/Multiome_fate/data/ICB_mouse/gmat_RNA_counts_ALL_dim35.rds")
atac_counts=readRDS("~/project/Multiome_fate/data/ICB_mouse/ALL_processed_FM.rds")
cell_meta=readRDS("~/project/Multiome_fate/data/ICB_mouse/metadat.rds")

rna_counts=rna_counts[,match(rownames(cell_meta), colnames(rna_counts))]
atac_counts=atac_counts[,match(rownames(cell_meta), colnames(atac_counts))]

ind_2000_each=c()
for(ii in  unique(cell_meta$Sample)){
  ind_2000_each=c(ind_2000_each,
                  sample(which(cell_meta$Sample==ii),2000, replace=F))
}

ind_2000_each=sort(ind_2000_each)
rna_counts=rna_counts[, ind_2000_each]
atac_counts=atac_counts[, ind_2000_each]


# Create Seurat object
myeloid <- CreateSeuratObject(counts = rna_counts)
myeloid[["percent.mt"]] <- PercentageFeatureSet(myeloid, pattern = "^mt-")


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  #fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
myeloid[["ATAC"]] <- chrom_assay


pdf("plots/vlnplot.pdf")
VlnPlot(myeloid, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
dev.off()

# RNA analysis
DefaultAssay(myeloid) <- "RNA"
myeloid <- SCTransform(myeloid, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(myeloid) <- "ATAC"
myeloid <- RunTFIDF(myeloid)
myeloid <- FindTopFeatures(myeloid, min.cutoff = 'q0')
myeloid <- RunSVD(myeloid)
myeloid <- RunUMAP(myeloid, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


myeloid <- FindMultiModalNeighbors(myeloid, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50),k.nn = 10)
myeloid <- RunUMAP(myeloid, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
myeloid <- FindClusters(myeloid, graph.name = "wsnn", algorithm = 1, verbose = FALSE)


#  plotting
pdf("./plots/Seurat_umap.pdf", width=10, height=6)
p1 <- DimPlot(myeloid, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(myeloid, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(myeloid, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

myeloid$celltype= myeloid$orig.ident

pdf("./plots/Seurat_umap_celltype.pdf", width=10, height=6)
p1 <- DimPlot(myeloid, reduction = "umap.rna",group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(myeloid, reduction = "umap.atac",group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(myeloid, reduction = "wnn.umap",group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 & theme(plot.title = element_text(hjust = 0.5))
p2 & theme(plot.title = element_text(hjust = 0.5))
p3 & theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf("./plots/Seurat_genes.pdf", width=10, height=6)
gene_plot1 <- FeaturePlot(myeloid, features = c("sct_Stat1"),reduction = 'umap.rna', ncol=1, pt.size = 0.1)
gene_plot2 <- FeaturePlot(myeloid, features = c("sct_Stat1"),reduction = 'umap.atac', ncol=1, pt.size = 0.1)
gene_plot3 <- FeaturePlot(myeloid, features = c("sct_Stat1"),reduction = 'wnn.umap', ncol=1, pt.size = 0.1)

print(gene_plot1)
print(gene_plot2)
print(gene_plot3)
#gene_plot1 + gene_plot2 + gene_plot3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

