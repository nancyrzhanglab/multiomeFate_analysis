
# Initialize ----
rm(list = ls())
gc()
setwd('/Volumes/GoogleDrive/My Drive/Schaff_Shared/Cloud/Experiment_IDs/DLS005/Segmented_Analysis_Pipeline/preprocess_GEX')
library(dplyr)
library(Seurat)
library(xlsx)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)

# Load in the GEX data ----
cocl2.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_CoCl2_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
acid.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Acid_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
dab.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Dab_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
tram.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Tram_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
cis.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Cis_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
dox.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Dox_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
naive1.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Naive1_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
naive2.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Naive2_BC_count_fromUnited/outs/filtered_feature_bc_matrix')
naive3.data <- Read10X(data.dir = '/Volumes/GoogleDrive/My Drive/Schaff_Shared/Super_cells/2022_01_14_RunCellRanger/DLS005_Naive3_BC_count_fromUnited/outs/filtered_feature_bc_matrix')


# Convert GEX into a seurat object ----
cocl2 <- CreateSeuratObject(counts = cocl2.data$`Gene Expression`)
acid <- CreateSeuratObject(counts = acid.data$`Gene Expression`)
dab <- CreateSeuratObject(counts = dab.data$`Gene Expression`)
tram <- CreateSeuratObject(counts = tram.data$`Gene Expression`)
cis <- CreateSeuratObject(counts = cis.data$`Gene Expression`)
dox <- CreateSeuratObject(counts = dox.data$`Gene Expression`)
naive1 <- CreateSeuratObject(counts = naive1.data$`Gene Expression`)
naive2 <- CreateSeuratObject(counts = naive2.data$`Gene Expression`)
naive3 <- CreateSeuratObject(counts = naive3.data$`Gene Expression`)

# Find the number of cells in each condition prior to filtering
starting_num_cells <- c(cocl2 = ncol(cocl2),
                        acid = ncol(acid),
                        dab = ncol(dab),
                        tram = ncol(tram),
                        cis = ncol(cis),
                        dox = ncol(dox),
                        naive1 = ncol(naive1),
                        naive2 = ncol(naive2),
                        naive3 = ncol(naive3))


# Add lineage data to the seurat object ----
cocl2[['lineage']] <- CreateAssayObject(counts = cocl2.data$Custom, min.cells = 1)
acid[['lineage']] <- CreateAssayObject(counts = acid.data$Custom, min.cells = 1)
dab[['lineage']] <- CreateAssayObject(counts = dab.data$Custom, min.cells = 1)
tram[['lineage']] <- CreateAssayObject(counts = tram.data$Custom, min.cells = 1)
cis[['lineage']] <- CreateAssayObject(counts = cis.data$Custom, min.cells = 1)
dox[['lineage']] <- CreateAssayObject(counts = dox.data$Custom, min.cells = 1)
naive1[['lineage']] <- CreateAssayObject(counts = naive1.data$Custom, min.cells = 1)
naive2[['lineage']] <- CreateAssayObject(counts = naive2.data$Custom, min.cells = 1)
naive3[['lineage']] <- CreateAssayObject(counts = naive3.data$Custom, min.cells = 1)

# Determine the amount of mitochondrial and ribosomal genes in each condition --------
cocl2[["percent.mt"]] <- PercentageFeatureSet(object = cocl2, pattern = "^MT-")
cocl2[["percent.rb"]] <- PercentageFeatureSet(object = cocl2, pattern = "^RPS")

acid[["percent.mt"]] <- PercentageFeatureSet(object = acid, pattern = "^MT-")
acid[["percent.rb"]] <- PercentageFeatureSet(object = acid, pattern = "^RPS")

dab[["percent.mt"]] <- PercentageFeatureSet(object = dab, pattern = "^MT-")
dab[["percent.rb"]] <- PercentageFeatureSet(object = dab, pattern = "^RPS")

tram[["percent.mt"]] <- PercentageFeatureSet(object = tram, pattern = "^MT-")
tram[["percent.rb"]] <- PercentageFeatureSet(object = tram, pattern = "^RPS")

cis[["percent.mt"]] <- PercentageFeatureSet(object = cis, pattern = "^MT-")
cis[["percent.rb"]] <- PercentageFeatureSet(object = cis, pattern = "^RPS")

dox[["percent.mt"]] <- PercentageFeatureSet(object = dox, pattern = "^MT-")
dox[["percent.rb"]] <- PercentageFeatureSet(object = dox, pattern = "^RPS")

naive1[["percent.mt"]] <- PercentageFeatureSet(object = naive1, pattern = "^MT-")
naive1[["percent.rb"]] <- PercentageFeatureSet(object = naive1, pattern = "^RPS")

naive2[["percent.mt"]] <- PercentageFeatureSet(object = naive2, pattern = "^MT-")
naive2[["percent.rb"]] <- PercentageFeatureSet(object = naive2, pattern = "^RPS")

naive3[["percent.mt"]] <- PercentageFeatureSet(object = naive3, pattern = "^MT-")
naive3[["percent.rb"]] <- PercentageFeatureSet(object = naive3, pattern = "^RPS")

#Visualize QC metrics as a violin plot ----

VlnPlot(object = cocl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = cocl2, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = cocl2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = acid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = acid, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = acid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = dab, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = dab, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = dab, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = tram, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = tram, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = tram, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = cis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = cis, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = cis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = dox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = dox, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = dox, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = naive1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = naive1, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = naive1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = naive2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = naive2, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = naive2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(object = naive3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(object = naive3, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(object = naive3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# From QC plots above set filter parameters, and replot after filtration ----

#subset based on filter parameters
cocl2_onc <- cocl2@assays$lineage@counts@Dim[2]
cocl2 <- subset(x = cocl2, subset = nFeature_RNA > 2500 & nCount_RNA < 50000 & percent.mt < 15) 
VlnPlot(object = cocl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
cocl2_f1 <- cocl2@assays$lineage@counts@Dim[2]/cocl2_onc
print(cocl2_f1*100)

#subset based on filter parameters
acid_onc <- acid@assays$lineage@counts@Dim[2]
acid <- subset(x = acid, subset = nFeature_RNA > 1500 & nCount_RNA < 30000 & percent.mt < 15) #  
VlnPlot(object = acid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
acid_f1 <- acid@assays$lineage@counts@Dim[2]/acid_onc
print(acid_f1*100)

#subset based on filter parameters
dab_onc <- dab@assays$lineage@counts@Dim[2]
dab <- subset(x = dab, subset = nFeature_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 15) #  
VlnPlot(object = dab, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
dab_f1 <- dab@assays$lineage@counts@Dim[2]/dab_onc
print(dab_f1*100)

#subset based on filter parameters
tram_onc <- tram@assays$lineage@counts@Dim[2]
tram <- subset(x = tram, subset = nFeature_RNA > 1500 & nCount_RNA < 20000 & percent.mt < 15) #  
VlnPlot(object = tram, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
tram_f1 <- tram@assays$lineage@counts@Dim[2]/tram_onc
print(tram_f1*100)

#subset based on filter parameters
cis_onc <- cis@assays$lineage@counts@Dim[2]
cis <- subset(x = cis, subset = nFeature_RNA > 1500 & nCount_RNA < 20000 & percent.mt < 20) 
VlnPlot(object = cis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
cis_f1 <- cis@assays$lineage@counts@Dim[2]/cis_onc
print(cis_f1*100)

#subset based on filter parameters
dox_onc <- dox@assays$lineage@counts@Dim[2]
dox <- subset(x = dox, subset = nFeature_RNA > 1500 & nCount_RNA < 20000 & percent.mt < 20) 
VlnPlot(object = dox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
dox_f1 <- dox@assays$lineage@counts@Dim[2]/dox_onc
print(dox_f1*100)

#subset based on filter parameters
naive1_onc <- naive1@assays$lineage@counts@Dim[2]
naive1 <- subset(x = naive1, subset = nFeature_RNA > 3000 & nCount_RNA < 75000 & percent.mt < 20) # 
VlnPlot(object = naive1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
naive1_f1 <- naive1@assays$lineage@counts@Dim[2]/naive1_onc
print(naive1_f1*100)

#subset based on filter parameters
naive2_onc <- naive2@assays$lineage@counts@Dim[2]
naive2 <- subset(x = naive2, subset = nFeature_RNA > 3000 & nCount_RNA < 60000 & percent.mt < 20) # 
VlnPlot(object = naive2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
naive2_f1 <- naive2@assays$lineage@counts@Dim[2]/naive2_onc
print(naive2_f1*100)

#subset based on filter parameters
naive3_onc <- naive3@assays$lineage@counts@Dim[2]
naive3 <- subset(x = naive3, subset = nFeature_RNA > 2500 & nCount_RNA < 50000 & percent.mt < 15) # 
VlnPlot(object = naive3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
#percent of cells kept after this filter
naive3_f1 <- naive3@assays$lineage@counts@Dim[2]/naive3_onc
print(naive3_f1*100)

# Find the number of cells per condition after filtering
postfilt_num_cells <- c(cocl2 = ncol(cocl2),
                        acid = ncol(acid),
                        dab = ncol(dab),
                        tram = ncol(tram),
                        cis = ncol(cis),
                        dox = ncol(dox),
                        naive1 = ncol(naive1),
                        naive2 = ncol(naive2),
                        naive3 = ncol(naive3))

# Write xlsx file with pre and post filtered cell numbers
df <- data.frame(starting = starting_num_cells, post_filt = postfilt_num_cells)
write.xlsx(df, file = 'cellsPreandPostFilt.xlsx')
rm(df)


# Remove unnecessary files ----
rm(list = c('acid.data',
            'cis.data',
            'cocl2.data',
            'dab.data',
            'tram.data',
            'dox.data',
            'naive1.data',
            'naive2.data',
            'naive3.data', 
            'acid_f1', 
            'acid_onc', 
            'cis_f1',
            'cis_onc',
            'cocl2_f1',
            'cocl2_onc', 
            'dab_f1',
            'dab_onc',
            'dox_f1',
            'dox_onc',
            'naive1_f1',
            'naive1_onc',
            'naive2_f1',
            'naive2_onc', 
            'naive3_f1',
            'naive3_onc',
            'tram_f1',
            'tram_onc'))



# allow R to use more memory ----
options(future.globals.maxSize= 35 * 1024^2)

# Normalize and identify variable features ----

cocl2 <- cocl2 %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
acid <- acid%>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
dab <- dab %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
tram <- tram %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
cis <- cis %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
dox <- dox %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
naive1 <- naive1 %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
naive2 <- naive2 %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
naive3 <- naive3 %>% Seurat::NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 20000) %>%  ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)

# Add cell cycle scores to cells ----

cocl2 <- CellCycleScoring(object = cocl2, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
acid <- CellCycleScoring(object = acid, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
dab <- CellCycleScoring(object = dab, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
tram <- CellCycleScoring(object = tram, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
cis <- CellCycleScoring(object = cis, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
dox <- CellCycleScoring(object = dox, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
naive1 <- CellCycleScoring(object = naive1, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
naive2 <- CellCycleScoring(object = naive2, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
naive3 <- CellCycleScoring(object = naive3, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# Run UMAP ----
cocl2 <- FindNeighbors(object = cocl2, dims = c(1:20))
cocl2 <- FindClusters(object = cocl2, resolution = .5)
cocl2 <- RunUMAP(object = cocl2, dims = c(1:20))

acid <- FindNeighbors(object = acid, dims = c(1:20))
acid <- FindClusters(object = acid, resolution = .5)
acid <- RunUMAP(object = acid, dims = c(1:20))

dab <- FindNeighbors(object = dab, dims = c(1:20))
dab <- FindClusters(object = dab, resolution = .5)
dab <- RunUMAP(object = dab, dims = c(1:20))

tram <- FindNeighbors(object = tram, dims = c(1:20))
tram <- FindClusters(object = tram, resolution = .5)
tram <- RunUMAP(object = tram, dims = c(1:20))

cis <- FindNeighbors(object = cis, dims = c(1:20))
cis <- FindClusters(object = cis, resolution = .5)
cis <- RunUMAP(object = cis, dims = c(1:20))

dox <- FindNeighbors(object = dox, dims = c(1:20))
dox <- FindClusters(object = dox, resolution = .5)
dox <- RunUMAP(object = dox, dims = c(1:20))

naive1 <- FindNeighbors(object = naive1, dims = c(1:20))
naive1 <- FindClusters(object = naive1, resolution = .5)
naive1 <- RunUMAP(object = naive1, dims = c(1:20))

naive2 <- FindNeighbors(object = naive2, dims = c(1:20))
naive2 <- FindClusters(object = naive2, resolution = .5)
naive2 <- RunUMAP(object = naive2, dims = c(1:20))

naive3 <- FindNeighbors(object = naive3, dims = c(1:20))
naive3 <- FindClusters(object = naive3, resolution = .5)
naive3 <- RunUMAP(object = naive3, dims = c(1:20))

# Plot UMAP clusters ----
pdf('CoCl2_UMAP_clusters.pdf')
DimPlot(cocl2, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters CoCl2')
dev.off()

('Acid_UMAP_clusters.pdf')
DimPlot(acid, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Acid')
dev.off()

('Dab_UMAP_clusters.pdf')
DimPlot(dab, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Dabrafenib')
dev.off()

pdf('Tram_UMAP_clusters.pdf')
DimPlot(tram, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Trametinib')
dev.off()

pdf('Cis_UMAP_clusters.pdf')
DimPlot(cis, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Cisplatin')
dev.off()

pdf('Dox_UMAP_clusters.pdf')
DimPlot(dox, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Doxorubicin')
dev.off()

pdf('Naive1_UMAP_clusters.pdf')
DimPlot(naive1, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Naive1')
dev.off()

pdf('Naive2_UMAP_clusters.pdf')
DimPlot(naive2, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Naive2')
dev.off()

pdf('Naive3_UMAP_clusters.pdf')
DimPlot(naive3, reduction = "umap",dims = c(1,2), group.by = 'seurat_clusters', pt.size = 1) + 
  DarkTheme() + ggtitle('Seurat Clusters Naive3')
dev.off()

# Make sure assay is RNA ----
DefaultAssay(object = cocl2) <- "RNA"
DefaultAssay(object = acid) <- "RNA"
DefaultAssay(object = dab) <- "RNA"
DefaultAssay(object = tram) <- "RNA"
DefaultAssay(object = cis) <- "RNA"
DefaultAssay(object = dox) <- "RNA"
DefaultAssay(object = naive1) <- "RNA"
DefaultAssay(object = naive2) <- "RNA"
DefaultAssay(object = naive3) <- "RNA"

# Plot UMAP markers ----
pdf('CoCl2_UMAP_markers.pdf')
FeaturePlot(object = cocl2, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

('Acid_UMAP_markers.pdf')
FeaturePlot(object = acid, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & 
  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

('Dab_UMAP_markers.pdf')
FeaturePlot(object = dab, reduction = 'umap',features = c("NGFR","EGFR","AXL","VGF","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & 
  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

pdf('Tram_UMAP_markers.pdf')
FeaturePlot(object = tram, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & 
  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

pdf('Cis_UMAP_markers.pdf')
FeaturePlot(object = cis, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & 
  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

pdf('Dox_UMAP_markers.pdf')
FeaturePlot(object = dox, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & 
  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

pdf('Naive1_UMAP_markers.pdf')
FeaturePlot(object = naive1, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) &  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

pdf('Naive2_UMAP_markers.pdf')
FeaturePlot(object = naive2, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

pdf('Naive3_UMAP_markers.pdf')
FeaturePlot(object = naive3, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) &  DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu')))
dev.off()

# Save objects prior to merging
save(acid,cis,cocl2,dab,dox,naive1,naive2,naive3,tram, file = 'Objects_premerged.RData')

# Merge the datasets ----
all_naive <- merge(naive1, y = c(naive2, naive3), add.cell.ids = c("naive1", "naive2", "naive3"), project = "Naive", merge.data = T)
all_resistant <- merge(cocl2, y = c(acid, dab, tram, cis, dox), add.cell.ids = c("cocl2", "acid", "dab", "tram", "cis", "dox"), project = "Resistant", merge.data = T)
all_data <- merge(naive1, y = c(naive2, naive3, cocl2, acid, dab, tram, cis, dox), add.cell.ids = c("naive1", "naive2", "naive3", "cocl2", "acid", "dab", "tram", "cis", "dox"), project = "All_Data", merge.data = T)

# Save the merged objects and remove the excess objects ----
rm(acid,cis,cocl2,dab,tram,naive1,naive2,naive3,dox,postfilt_num_cells,starting_num_cells)   
save(all_naive,all_resistant,all_data, file = 'postmerged_objects.RData')

# scale merged data and generate UMAPs

all_naive <- FindVariableFeatures(all_naive) %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunUMAP(dims = 1:20, verbose = FALSE) %>% FindNeighbors(dims = 1:20, verbose = FALSE) %>%  FindClusters(resolution = 0.5)
DimPlot(all_naive)
save(all_naive, file = 'all_naive_merged.RData')
rm(all_naive)

all_resistant <- FindVariableFeatures(all_resistant) %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunUMAP(dims = 1:20, verbose = FALSE) %>% FindNeighbors(dims = 1:20, verbose = FALSE) %>%  FindClusters(resolution = 0.5)
DimPlot(all_resistant)
save(all_resistant, file = 'all_resistant_merged.RData')
rm(all_resistant)

all_data <- FindVariableFeatures(all_data) %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunUMAP(dims = 1:20, verbose = FALSE) %>% FindNeighbors(dims = 1:20, verbose = FALSE) %>%  FindClusters(resolution = 0.5)
DimPlot(all_data)
save(all_data, file = 'all_data_merged.RData')
rm(all_data)

# load in data and assign the original condition to metadata ----

# all_naive
load('all_naive_merged.RData')
condition <- c()
for (i in 1:length(all_naive@assays$RNA@counts@Dimnames[2][[1]])){
  condition <- rbind(condition, strsplit(all_naive@assays$RNA@counts@Dimnames[2][[1]][[i]], '_')[[1]][1])
}
all_naive$OG_condition <- condition
save(all_naive, file = 'all_naive_merged.RData')
pdf('all_naive_merged.pdf')
print(DimPlot(all_naive, label = T, group.by = 'OG_condition') + DarkTheme())
print(FeaturePlot(object = all_naive, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu'))))
dev.off()
rm(all_naive)

# all_resistant
load('all_resistant_merged.RData')
condition <- c()
for (i in 1:length(all_resistant@assays$RNA@counts@Dimnames[2][[1]])){
  condition <- rbind(condition, strsplit(all_resistant@assays$RNA@counts@Dimnames[2][[1]][[i]], '_')[[1]][1])
}
all_resistant$OG_condition <- condition
save(all_resistant, file = 'all_resistant_merged.RData')
pdf('all_resistant_merged.pdf')
print(DimPlot(all_resistant, label = T, group.by = 'OG_condition') + DarkTheme())
print(FeaturePlot(object = all_resistant, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu'))))
dev.off()
rm(all_resistant)

# all_data
load('all_data_merged.RData')
condition <- c()
for (i in 1:length(all_data@assays$RNA@counts@Dimnames[2][[1]])){
  condition <- rbind(condition, strsplit(all_data@assays$RNA@counts@Dimnames[2][[1]][[i]], '_')[[1]][1])
}
all_data$OG_condition <- condition
all_data$OG_condition[all_data$OG_condition %in% c('naive1', 'naive2', 'naive3')] <- 'naive'
save(all_data, file = 'all_data_merged.RData')
pdf('all_data_merged.pdf')
print(DimPlot(all_data, label = T, group.by = 'OG_condition') + DarkTheme())
print(FeaturePlot(object = all_data, reduction = 'umap',features = c("NGFR","EGFR","AXL","NT5E","FN1","SERPINE2","MITF","nFeature_RNA","nCount_RNA"), pt.size = 1, combine = T, order = TRUE) & DarkTheme() & scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = 'RdBu'))))
dev.off()
rm(all_data, condition)

# Plot the all_data umap with each condition in its color separately - (ground truth not barcodes) ----

load('all_data_merged.RData')

# Plotting overlay
pdf('CoCl2_res.pdf')
cocl2_res <- names(all_data$orig.ident[all_data$OG_condition == 'cocl2']) # Gets cocl2 res cells
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(cocl2_res),
        cols.highlight = 'forestgreen') + DarkTheme() + ggtitle('CoCl2')
dev.off()

pdf('Acid_res.pdf')
acid_res <- names(all_data$orig.ident[all_data$OG_condition == 'acid']) # Gets resistant cell identifiers
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(acid_res),
        cols.highlight = 'darksalmon') + DarkTheme() + ggtitle('Acid')
dev.off()

pdf('Dab_res.pdf')
dab_res <- names(all_data$orig.ident[all_data$OG_condition == 'dab']) # Gets resistant cell identifiers
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(dab_res),
        cols.highlight = 'cyan') + DarkTheme() + ggtitle('Dabrafenib')
dev.off()

pdf('Tram_res.pdf')
tram_res <- names(all_data$orig.ident[all_data$OG_condition == 'tram']) # Gets resistant cell identifiers
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(tram_res),
        cols.highlight = 'hotpink') + DarkTheme() + ggtitle('Trametinib')
dev.off()

pdf('Cis_res.pdf')
cis_res <- names(all_data$orig.ident[all_data$OG_condition == 'cis']) # Gets resistant cell identifiers
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(cis_res),
        cols.highlight = 'darkgoldenrod') + DarkTheme() + ggtitle('Cisplatin')
dev.off()

pdf('Dox_res.pdf')
dox_res <- names(all_data$orig.ident[all_data$OG_condition == 'dox']) # Gets resistant cell identifiers
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(dox_res),
        cols.highlight = 'darkmagenta') + DarkTheme() + ggtitle('Doxorubicin')
dev.off()

pdf('Naive_res.pdf')
naive_res <- names(all_data$orig.ident[all_data$OG_condition == 'naive']) # Gets resistant cell identifiers
DimPlot(all_data, reduction = "umap",dims = c(1,2),
        group.by = 'OG_condition', pt.size = .1, cells.highlight = list(naive_res),
        cols.highlight = '#00A9FF') + DarkTheme() + ggtitle('Naive')
dev.off()