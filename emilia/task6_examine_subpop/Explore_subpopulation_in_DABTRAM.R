rm(list = ls())

library(tidyverse)
library(ggplot2)
library(GGally)

remove_unassigned_cells <- TRUE
# =============================================================================
# reading data
# =============================================================================
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_week5_DABTRAM.RData'))

all_data[[paste0("chromVar.day10_DABTRAM")]] <- all_data_chromVar_day10_DABTRAM
all_data[[paste0("chromVar.week5_DABTRAM")]] <- all_data_chromVar_week5_DABTRAM

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# =============================================================================
# wrangle
# =============================================================================

# day10
all_data.day10_DABTRAM <- subset(all_data, dataset == "day10_DABTRAM")


cv.day10_DABTRAM <- as.data.frame(t(all_data.day10_DABTRAM@assays[["chromVar.day10_DABTRAM"]]@data))

tfs <- colnames(cv.day10_DABTRAM)
tfs_of_interest <- tfs[grepl("JUN", tfs)]
tfs_of_interest <- c(tfs_of_interest, tfs[grepl("STAT", tfs)])

tfs_of_interest <- c('JUN', 'FOS::JUN',  'STAT1')

cv.day10_DABTRAM.toplot <- cv.day10_DABTRAM[, tfs_of_interest]

ggpairs(cv.day10_DABTRAM.toplot)

# week5
all_data.week5_DABTRAM <- subset(all_data, dataset == "week5_DABTRAM")

cv.week5_DABTRAM <- as.data.frame(t(all_data.week5_DABTRAM@assays[["chromVar.week5_DABTRAM"]]@data))


tfs_of_interest <- c('JUN', 'FOS::JUN', 'STAT1')

cv.week5_DABTRAM.toplot <- cv.week5_DABTRAM[, tfs_of_interest]

ggpairs(cv.week5_DABTRAM.toplot)

# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

rm(list = ls())

library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggdensity)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '~/Downloads/'

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))


all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[['FastTopic_DABTRAM']] <- all_data_fasttopic_DABTRAM

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

gavish.mp <- read_csv(paste0(ref_dir, 'Gavish_Malignant_Meta_Programs.csv'))
isg_rs <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/ISG.RS.txt', col_names = F)
isg_mem <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Persistent IFN ISG Groups/Memory ISGs Human.csv', col_names = F)

colnames(isg_rs) <- 'Gene'
colnames(isg_mem) <- 'Gene'

# ==============================================================================
# Wrangle data
# ==============================================================================
gavish.mp.list <- list()
for (mp in colnames(gavish.mp)) {
  gavish.mp.list[[mp]] <- c(gavish.mp[[mp]] %>% as.character())
}

gavish.mp.list <- lapply(1:length(gavish.mp.list), function(x) {
  gs <- gavish.mp.list[[x]]
  gs <- unique(na.omit(gs[gs %in% all_data_saver@data@Dimnames[[1]]]))
  return(gs)
})
names(gavish.mp.list) <- colnames(gavish.mp)

isg_rs <- unique(na.omit(isg_rs$Gene[isg_rs$Gene %in% all_data_saver@data@Dimnames[[1]]]))
isg_mem <- unique(na.omit(isg_mem$Gene[isg_mem$Gene %in% all_data_saver@data@Dimnames[[1]]]))

# ==============================================================================
# Calculate module scores 
# ==============================================================================

metadat <- all_data@meta.data
metadat.week5_dabtram <- metadat[metadat$dataset == 'week5_DABTRAM', ]

scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=c(gavish.mp.list, 
                                           list(isg_rs), 
                                           list(isg_mem)))

colnames(scores) <- c(names(gavish.mp.list), 'ISG.RS', 'ISG.Mem')
scores.df <- as.data.frame(scores) 

scores.df.week5_dabtram <- scores.df[rownames(metadat.week5_dabtram), ]

mps <- c('Stress', 'Interferon/MHC-II (I)', 'Interferon/MHC-II (II)', 'ISG.RS', 'ISG.Mem')

scores.df.week5_dabtram.toplot <- scores.df.week5_dabtram[, mps]


ggpairs(scores.df.week5_dabtram.toplot)


ft_DABTRAM_umap <- as.data.frame(all_data_ft_DABTRAM_umap@cell.embeddings)
ft_DABTRAM_umap$cell_id <- rownames(ft_DABTRAM_umap)

scores.df.week5_dabtram$cell_id <- rownames(scores.df.week5_dabtram)

ft_DABTRAM_umap <- merge(ft_DABTRAM_umap, scores.df.week5_dabtram, by = 'cell_id', all = T)

ft_DABTRAM_umap <- ft_DABTRAM_umap[order(ft_DABTRAM_umap$Stress, decreasing = F), ]
ft_DABTRAM_umap <- ft_DABTRAM_umap[order(ft_DABTRAM_umap$`Interferon/MHC-II (I)`, decreasing = F), ]

ft_DABTRAM_umap$`Interferon/MHC-II (I) SCALED` <- scale(ft_DABTRAM_umap$`Interferon/MHC-II (I)`)
ft_DABTRAM_umap$`Stress SCALED` <- scale(ft_DABTRAM_umap$`Stress`)
ggplot(ft_DABTRAM_umap, 
       aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = `Stress SCALED`)) +
  geom_point() +
  xlim(-1, 10) +
  ylim(-12, 0) +
  scale_color_gradient2(low = 'blue', high = 'red')



# ==============================================================================================
# ==============================================================================================
# ==============================================================================================

rm(list = ls())

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))

load(paste0(data_dir, 'Writeup10a_data_chromVar_week5_DABTRAM.RData'))


all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[['FastTopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[["ft.DABTRAM.umap"]] <- all_data_ft_DABTRAM_umap


# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# ==============================================================================
# Wrangle data
# ==============================================================================
all_data@active.assay <- 'Saver'

all_data.dabtram <- subset(all_data, dataset != 'week5_COCL2')
all_data.dabtram <- subset(all_data.dabtram, dataset != 'week5_CIS')
all_data.dabtram <- subset(all_data.dabtram, dataset != 'day10_COCL2')
all_data.dabtram <- subset(all_data.dabtram, dataset != 'day10_CIS')

all_data.dabtram <- RunPCA(all_data.dabtram, features = VariableFeatures(object = all_data.dabtram))
all_data.dabtram <- FindNeighbors(all_data.dabtram, dims = 1:10)
all_data.dabtram <- FindClusters(all_data.dabtram, resolution = 0.2)
all_data.dabtram <- RunUMAP(all_data.dabtram, dims = 1:10)

DimPlot(all_data.dabtram, reduction = 'ft.DABTRAM.umap')

metadat.all_data_dabtram <- all_data.dabtram@meta.data
metadat.all_data_dabtram$cell_id <- rownames(metadat.all_data_dabtram)

# write.csv(metadat.all_data_dabtram, paste0(data_dir, 'metadat.all_data_dabtram_recluster.csv'),row.names = T)

# ==============================================================================
# Plot
# ==============================================================================

scores.df.week5_dabtram <- merge(scores.df.week5_dabtram, metadat.all_data_dabtram[, c('cell_id', 'seurat_clusters')], by = 'cell_id')
scores.df.week5_dabtram$Stress_Scaled <- scale(scores.df.week5_dabtram$Stress)
scores.df.week5_dabtram$`Interferon/MHC-II (I) Scaled` <- scale(scores.df.week5_dabtram$`Interferon/MHC-II (I)`)
ggplot(scores.df.week5_dabtram, 
       aes(x = seurat_clusters, y = `ISG.RS`)) +
  geom_violin(scale = 'width') +
  geom_jitter(alpha = 0.1, size = 1, width = 0.1) +
  geom_boxplot(width = 0.1) +
  # coord_cartesian(ylim = c(-3, 3)) +
  theme_bw()

# ==============================================================================
# Add in TF
# ==============================================================================
cv.week5_dabtram <- all_data_chromVar_week5_DABTRAM@data
cv.week5_dabtram <- cv.week5_dabtram[, scores.df.week5_dabtram$cell_id]
cv.week5_dabtram <- as.data.frame(t(cv.week5_dabtram))
cv.week5_dabtram$cell_id <- rownames(cv.week5_dabtram)

scores.df.week5_dabtram <- merge(scores.df.week5_dabtram[, c('cell_id', 'seurat_clusters')], cv.week5_dabtram[, c('cell_id',  'JUN')], by = 'cell_id')

tfs <- colnames(cv.week5_dabtram)
tfs_of_interest <- tfs[grepl("FOS", tfs)]

ggplot(scores.df.week5_dabtram, 
       aes(x = seurat_clusters, y = `JUN`)) +
  geom_violin(scale = 'width') +
  geom_jitter(alpha = 0.1, size = 1, width = 0.1) +
  geom_boxplot(width = 0.1) +
  # coord_cartesian(ylim = c(-3, 3)) +
  theme_bw()
