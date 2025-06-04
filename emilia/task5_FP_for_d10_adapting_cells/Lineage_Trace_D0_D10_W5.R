rm(list = rm())
library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

remove_unassigned_cells <- TRUE
# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

df.bias <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
fp.d0_d10 <- as.data.frame(all_data_fatepotential[[paste0("fatepotential_DABTRAM_d0_d10")]][["cell_imputed_score"]])
# ==============================================================================
# Calcualte module scores 
# ==============================================================================
custom.cis <- c('NUCKS1','TUBB','HMGB2', 'ICMT', 'CBX5', 'TUBA1B', 'ANP32B','TYMS',
                'GMNN', 'USP1', 'NASP', 'TMPO', 'NCAPH', 'TK1', 'TUBG1', 'PRC1',
                'PBK', 'SMC3', 'RRM2', 'RAD51AP1')
custom.dabtram <- c('ACTB', 'TMEM43', 'TPM4', 'CALM2', 'FN1', 'PALLD', 'LMO7',
                    'ACTN1', 'HSPG2', 'MYOF', 'TNFRSF12A', 'TUBB', 'RCN1', 'CRIM1',
                    'COL5A2', 'SAMD5', 'TPM1', 'OXSR1', 'CBX5')
custom.cocl2 <- c('GXYLT2', 'ANTXR1', 'CADM1', 'ITGB3', 'BICC1', 'SLC1A4',
                  'CADPS', 'HMGA2', 'TIMP3', 'PTPRG', 'SERPINE2', 'IMMP2L', 'LRMDA',
                  'MFSD12', 'SOX5', 'EPHA3', 'PRKG2', 'IL1RAP', 'SLC44A1', 'KCNQ5')

metadat <- all_data@meta.data
scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=c('custom.cis' = list(custom.cis), 
                                           'custom.dabtram' = list(custom.dabtram), 
                                           'custom.cocl2' = list(custom.cocl2)))
all_data <- AddMetaData(all_data, metadata = scores)

gavish.mp.df <- read.csv(paste0(result_dir, 'GavishMP_UCell_scores.csv'), row.names = 1)

# ==============================================================================
# Subset
# ==============================================================================
dabtram_d10_w5 <- subset(all_data, subset = (dataset %in% c('day10_DABTRAM', 'week5_DABTRAM')))

metadat_dabtram_d10_w5 <- dabtram_d10_w5@meta.data
lin.size_dabtram_d10_w5 <- table(metadat_dabtram_d10_w5$assigned_lineage, metadat_dabtram_d10_w5$dataset)
lin.size_dabtram_d10_w5 <- as.data.frame(lin.size_dabtram_d10_w5)
lin.size_dabtram_d10_w5 <- pivot_wider(lin.size_dabtram_d10_w5, names_from = Var2, values_from = Freq)
lin.present.dabtram_d10_w5 <- lin.size_dabtram_d10_w5[lin.size_dabtram_d10_w5$day10_DABTRAM > 0 & 
                                                      lin.size_dabtram_d10_w5$week5_DABTRAM > 0, ]   
cell.high.fp <- metadat_dabtram_d10_w5[metadat_dabtram_d10_w5$fatepotential_DABTRAM_d10_w5 > 0, ]
cell.w5 <- metadat_dabtram_d10_w5[metadat_dabtram_d10_w5$dataset == 'week5_DABTRAM', ]

dabtram_d10_w5_both <- subset(dabtram_d10_w5, subset = assigned_lineage %in% lin.present.dabtram_d10_w5$Var1)
dabtram_d10_w5_both <- dabtram_d10_w5_both[, c(rownames(cell.high.fp), rownames(cell.w5))]
metadat.dabtram_d10_w5_both <- dabtram_d10_w5_both@meta.data

# recluster
dabtram_d10_w5_both@active.assay <- 'Saver'
dabtram_d10_w5_both <- RunPCA(dabtram_d10_w5_both, features = VariableFeatures(object = dabtram_d10_w5_both), verbose = FALSE)
dabtram_d10_w5_both  <- FindNeighbors(dabtram_d10_w5_both, dims = 1:20)
dabtram_d10_w5_both  <- FindClusters(dabtram_d10_w5_both, resolution = 0.5)
dabtram_d10_w5_both  <- RunUMAP(dabtram_d10_w5_both, dims = 1:20)
DimPlot(dabtram_d10_w5_both, reduction = "umap", group.by = "dataset", label = TRUE) + NoLegend() 
FeaturePlot(dabtram_d10_w5_both, features = c('fatepotential_DABTRAM_d10_w5'), 
            cols = c("blue", "red"), pt.size = 0.1, order = T) + NoLegend()
  
ggplot(metadat.dabtram_d10_w5_both, aes(x = dataset, y = signature_2_UCell)) +
  geom_violin(aes(fill = dataset), scale = 'width') +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', label.y = 0.95) +
  scale_fill_manual(values = dataset_colors) +
  ylab('DABTRAM_adaptation_UCell') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gavish.mp.df.dabtram <- gavish.mp.df[rownames(metadat.dabtram_d10_w5_both), ]
gavish.mp.df.dabtram$dataset <- ifelse(grepl('day10', rownames(gavish.mp.df.dabtram)), 'day10_DABTRAM', 'week5_DABTRAM')
ggplot(gavish.mp.df.dabtram, aes(x = dataset, y = Interferon.MHC.II..I._UCell)) +
  geom_violin(aes(fill = dataset), scale = 'width') +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  stat_compare_means(method = 'wilcox.test', label.y = 0.95) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gavish.mp.df.dabtram, aes(x = dataset, y = Stress_UCell)) +
  geom_violin(aes(fill = dataset), scale = 'width') +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  stat_compare_means(method = 'wilcox.test', label.y = 0.95) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ===============================================================================
# Linke d0 d10 w5
# ===============================================================================
colnames(fp.d0_d10) <- c("fatepotential_d0_d10")
fp.d0_d10$cell_id <- rownames(fp.d0_d10)

df <- merge(df.bias, fp.d0_d10, by = 'cell_id')
df$imputed_cell_count <- 10**df$fatepotential_d0_d10
df$adapting_cell_count <- 10**df$adaptingFP
df$frac_adapting <- df$adapting_cell_count / df$imputed_cell_count * 100
df$frac_adapting <- ifelse(df$frac_adapting > 100, 100, df$frac_adapting)

df$fp_dabtram_hi_lo <- ifelse(df$fatepotential_d0_d10 > quantile(df$fatepotential_d0_d10)[4], 'high', 'Other')
df$fp_dabtram_hi_lo <- ifelse(df$fatepotential_d0_d10 < quantile(df$fatepotential_d0_d10)[2], 'low', df$fp_dabtram_hi_lo)
df$fb_dabtram_hi_lo <- ifelse(df$frac_adapting > quantile(df$frac_adapting)[4], 'high', 'Other')
df$fb_dabtram_hi_lo <- ifelse(df$frac_adapting < quantile(df$frac_adapting)[2], 'low', df$fb_dabtram_hi_lo)

day0.cells.high.FP.FB <- df[df$fp_dabtram_hi_lo == 'high' & df$fb_dabtram_hi_lo == 'high', ]

mat.saver <- all_data@assays[["Saver"]]@data
mat.saver.use <- mat.saver[, c(day0.cells.high.FP.FB$cell_id,
                               colnames(dabtram_d10_w5_both))]
mat.saver.use <- as.data.frame(t(mat.saver.use))
mat.saver.use$cell_id <- rownames(mat.saver.use)
mat.saver.use$dataset <- ifelse(grepl('day0', mat.saver.use$cell_id), 'day0',
                                ifelse(grepl('day10', mat.saver.use$cell_id), 'day10_DABTRAM',
                                       'week5_DABTRAM'))
mat.saver.summary <- mat.saver.use %>% 
  select(-cell_id) %>% 
  group_by(dataset) %>% 
  summarise_all( mean)
mat.saver.summary <- as.data.frame(mat.saver.summary)
rownames(mat.saver.summary) <- mat.saver.summary$dataset
mat.saver.summary <- mat.saver.summary[,-1]

row.var <- matrixStats::colVars(as.matrix(mat.saver.summary))
gene.highvar <- names(row.var[row.var > 1])

pheatmap::pheatmap(mat.saver.summary[, gene.highvar],
         cluster_rows = F,
         cluster_cols = T,
         fontsize_col = 6,
         width = 20,
         filename = '~/Downloads/heatmap.png')


# ===============================================================================
# Get d0 and d10 cells
# ===============================================================================
metadat.dabtram_d10 <- metadat.dabtram_d10_w5_both[metadat.dabtram_d10_w5_both$dataset == 'day10_DABTRAM', ]
metadat.dabtram_d10 <- metadat.dabtram_d10[order(metadat.dabtram_d10$fatepotential_DABTRAM_d10_w5), ]

day0.cells.high.FP.FB <- day0.cells.high.FP.FB[order(day0.cells.high.FP.FB$frac_adapting), ]

mat.saver.use <- mat.saver.use[c(day0.cells.high.FP.FB$cell_id, rownames(metadat.dabtram_d10)), ]
mat.saver.use$order <- seq(1, nrow(mat.saver.use), 1)

genes <- c('TPM4', 'CALM2', 'FN1', 'LMO7',
           'AXL', 'NGFR', 'JUN', 'WNT5A',
           'IFI27', 'IFIT3', 'HLA-A', 'MITF')

mat.saver.use.genes <- mat.saver.use[, c(genes, 'order')]
mat.saver.use.genes.melt <- melt(mat.saver.use.genes, id.vars = c('order'))

ggplot(mat.saver.use, aes(x = order, y = log2(AXL))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'loess', se = T) +
  geom_vline(xintercept = nrow(day0.cells.high.FP.FB), linetype = 'dashed', color = 'red') +
  theme_bw()

ggplot(mat.saver.use.genes.melt, aes(x = order, y = log2(value))) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = variable),lwd=2, method = 'loess', se = T) +
  geom_vline(xintercept = nrow(day0.cells.high.FP.FB), linetype = 'dashed', color = 'red') +
  facet_wrap(~ variable, ncol = 4, scale = 'free_y') +
  theme_bw()


row.var <- matrixStats::colVars(as.matrix(mat.saver.use[, seq(1, ncol(mat.saver.use)-3)]))
gene.highvar <- names(row.var[row.var > 0.1])

row.mean <- colMeans(as.matrix(mat.saver.use[, seq(1, ncol(mat.saver.use)-3)]))
gene.highmean <- names(row.mean[row.mean > 0.1])

genes.to.plot <- intersect(gene.highvar, gene.highmean)

pheatmap::pheatmap(mat.saver.use[, seq(1, ncol(mat.saver.use)-3)],
         cluster_rows = F,
         cluster_cols = T,
         fontsize_col = 6,
         width = 20,
         filename = '~/Downloads/heatmap_ordered_by_adaptation_d0_d10.png')

mat.saver.use.toplot <- t(mat.saver.use[, genes.to.plot])
pheatmap::pheatmap(mat.saver.use.toplot,
                   cluster_rows = T,
                   cluster_cols = F,
                   fontsize_col = 6,
                   width = 25,
                   height = 10)

set.seed(123)
g <- pheatmap::pheatmap(mat.saver.use.toplot,
                   kmeans_k = 3,
                   cluster_rows = T,
                   cluster_cols = F,
                   fontsize_col = 6,
                   width = 25,
                   height = 10)
clusters <- g$kmeans$cluster
clusters <- as.data.frame(clusters)
clusters$gene <- rownames(clusters)

genes.to.plot <- clusters[clusters$clusters %in% c(1, 3), ]$gene
pheatmap::pheatmap(mat.saver.use.toplot[genes.to.plot, ],
                   cluster_rows = T,
                   cluster_cols = F,
                   fontsize_col = 6,
                   width = 25,
                   height = 10)
