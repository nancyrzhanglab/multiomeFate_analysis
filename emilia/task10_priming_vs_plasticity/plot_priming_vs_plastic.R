rm(list=ls())
library(Seurat)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))


# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

metadat.DABTRAM <- read.csv(paste0(data_dir, 'metadat.all_data_dabtram_recluster.csv'), row.names = 1)
metadat.week5_DABTRAM <- metadat.DABTRAM[metadat.DABTRAM$dataset == 'week5_DABTRAM', ]
metadat.week5_DABTRAM.clust0 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 0)
metadat.week5_DABTRAM.clust3 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 3)

metadat.week5_DABTRAM.clust0$category <- ifelse(metadat.week5_DABTRAM.clust0$assigned_lineage %in% lin.priming,
                                                'Priming', 'Other')
metadat.week5_DABTRAM.clust0$category <- ifelse(metadat.week5_DABTRAM.clust0$assigned_lineage %in% lin.plasticity,
                                                'Plastic', metadat.week5_DABTRAM.clust0$category)
metadat.week5_DABTRAM.clust3$category <- ifelse(metadat.week5_DABTRAM.clust3$assigned_lineage %in% lin.priming,
                                                'Priming', 'Other')
metadat.week5_DABTRAM.clust3$category <- ifelse(metadat.week5_DABTRAM.clust3$assigned_lineage %in% lin.plasticity,
                                                'Plastic', metadat.week5_DABTRAM.clust3$category)


df.c0 <- metadat.week5_DABTRAM.clust0 %>%
  group_by(category) %>%
  summarise(n = n()) %>% 
  mutate(cluster = 'Cluster 0')
df.c0$freq <- df.c0$n / nrow(metadat.week5_DABTRAM)

df.c3 <- metadat.week5_DABTRAM.clust3 %>%
  group_by(category) %>%
  summarise(n = n()) %>% 
  mutate(cluster = 'Cluster 3')
df.c3$freq <- df.c3$n / nrow(metadat.week5_DABTRAM)

df <- rbind(df.c0, df.c3)
df$category <- factor(df$category, levels = c('Priming', 'Plastic', 'Other'))

p.w5 <- ggplot(df, aes(x = cluster, y = freq, fill = category)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  labs(x = '', y = 'Fraction of cells (out of w5 cells)') +
  scale_fill_manual(values = c('Priming' = '#92C7CF', 'Plastic' = '#FCDC94', 'Other' = '#FBF9F1')) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


umap <- as.data.frame(all_data_ft_DABTRAM_umap@cell.embeddings)
umap$cell_id <- rownames(umap)
umap <- merge(umap, metadat.DABTRAM[, c('cell_id', 'seurat_clusters', 'dataset')], by = 'cell_id')
umap$dataset <- factor(umap$dataset, levels = c('day0', 'day10_DABTRAM', 'week5_DABTRAM'))
umap <- umap[order(umap$dataset, decreasing = F), ]

# captain = c("0" = '#FF7A5C', "1"="grey","2"="#ccebc5","3"="#4858A7","4"="#d4c3a7")
captain = c("0" = '#FF7A5C', "1"="#E8E8E8","2"="#E8E8E8","3"="#4858A7","4"="#E8E8E8")
umap$seurat_clusters <- as.character(umap$seurat_clusters)
p.umap <- ggplot(umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = as.factor(seurat_clusters))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = captain) +
  labs(x = '', y = '', color = 'Cluster') +
  theme_Publication() +
  theme(legend.position = 'right')

ggarrange(p.umap, p.w5, ncol = 2, nrow = 1, widths = c(1, 0.7), align = 'hv')
ggsave(paste0(out_dir, 'DABTRAM_w5_cluster0_3.pdf'), width = 6, height = 2.8)

out_dir <- '~/Downloads/'







gavish.mp <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/saver.GavishMP.UCellScores.csv')
rownames(gavish.mp) <- gavish.mp$cell_id
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))

cv.mat <- all_data_chromVar_day10_DABTRAM@data
cv.mat <- t(cv.mat)
cv.mat <- cv.mat[rownames(gavish.mp), ]
cv.mat <- as.data.frame(cv.mat)
cv.mat <- cv.mat[rownames(metadat.day10_DABTRAM), ]

cv.mat$cell_id <- rownames(cv.mat)

ap1 <- grep('JUN|FOS', colnames(cv.mat), value = T)
cv.mat.ap1 <- cv.mat[, ap1]
cv.mat.ap1$cell_id <- rownames(cv.mat.ap1)

cv.mat.ap1 <- merge(cv.mat.ap1, gavish.mp, by = 'cell_id')

ggplot(cv.mat.ap1, aes(x = JUN, y = Stress)) +
  geom_point() +
  stat_cor() +
  geom_smooth(method = 'lm')

