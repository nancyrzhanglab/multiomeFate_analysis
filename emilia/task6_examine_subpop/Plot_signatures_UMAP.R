rm(list = rm())
library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggdensity)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig5/'

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            # axis.ticks = element_blank(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}
# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
score.df <- read.csv(paste0(data_dir, 'saver.GavishMP.UCellScores.csv'))

all_data@misc <- all_data_fatepotential
all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
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
# wrangle
# ==============================================================================

# ft umap
ft_umap <- all_data@reductions[["ft.DABTRAM.umap"]]@cell.embeddings
ft_umap <- as.data.frame(ft_umap)
ft_umap$cell_id <- rownames(ft_umap)

ft_umap <- merge(ft_umap, score.df, by = 'cell_id')

# metadata
metadata <- all_data@meta.data
metadata$cell_id <- rownames(metadata)
metadata.day10_DABTRAM <- metadata[metadata$dataset == 'day10_DABTRAM', ]
metadata.rest <- metadata[metadata$dataset != 'day10_DABTRAM', ]

# ==============================================================================
# plot
# ==============================================================================

ft_umap <- ft_umap[order(ft_umap$fatepotential_DABTRAM_d10_w5),]
ft_umap$fatepotential_DABTRAM_d10_w5 <- ifelse(ft_umap$fatepotential_DABTRAM_d10_w5 < -1, -1, ft_umap$fatepotential_DABTRAM_d10_w5)
ft_umap$fatepotential_DABTRAM_d10_w5 <- ifelse(ft_umap$fatepotential_DABTRAM_d10_w5 > 1, 1, ft_umap$fatepotential_DABTRAM_d10_w5)
rownames(ft_umap) <- ft_umap$cell_id
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.day10_DABTRAM$cell_id), ]
ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = fatepotential_DABTRAM_d10_w5)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "red", mid = "bisque", high = "blue", midpoint = -0.2, na.value = 'gray') +
  theme_Publication()

ft_umap$ISG.RS.scale <- scale(ft_umap$ISG.RS)
ft_umap$ISG.RS.scale <- ifelse(ft_umap$cell_id %in% metadata.day10_DABTRAM$cell_id, ft_umap$ISG.RS.scale, NA)
ft_umap$ISG.RS.scale <- ifelse(ft_umap$ISG.RS.scale < -1, -1, ft_umap$ISG.RS.scale)
ft_umap$ISG.RS.scale <- ifelse(ft_umap$ISG.RS.scale > 1, 1, ft_umap$ISG.RS.scale)

ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = ISG.RS.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('DABTRAM(day10)') +
  theme_Publication()

ft_umap <- ft_umap[order(ft_umap$ISG.RS.scale),]
ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = ISG.RS.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('DABTRAM (all)') +
  theme_Publication()



ft_umap$Interferon.MHC.II..I..scale <- scale(ft_umap$Interferon.MHC.II..I.)
ft_umap$Interferon.MHC.II..I..scale <- ifelse(ft_umap$Interferon.MHC.II..I..scale < -1, -1, ft_umap$Interferon.MHC.II..I..scale)
ft_umap$Interferon.MHC.II..I..scale <- ifelse(ft_umap$Interferon.MHC.II..I..scale  > 1, 1, ft_umap$Interferon.MHC.II..I..scale )
ft_umap$Interferon.MHC.II..I..scale <- ifelse(ft_umap$cell_id %in% metadata.day10_DABTRAM$cell_id, ft_umap$Interferon.MHC.II..I..scale, NA)
ft_umap <- ft_umap[order(ft_umap$Interferon.MHC.II..I..scale),]
ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = Interferon.MHC.II..I..scale )) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('DABTRAM(day10)') +
  theme_Publication()

ggplot(ft_umap, aes(x = fatepotential_DABTRAM_d10_w5, y = ISG.RS)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', alpha = 0.2, se = T) +
  stat_cor() +
  ggtitle('DABTRAM (day10)') +
  theme_Publication()

metadata.day0 <- metadata[metadata$dataset == 'day0', ]
ft_umap.filter <- ft_umap[ft_umap$cell_id %in% metadata.day0$cell_id, ]
ggplot(ft_umap.filter, aes(x = fatepotential_DABTRAM_d0_d10, y = ISG.RS)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', alpha = 0.2) +
  stat_cor(method = 'spearman') +
  ggtitle('DABTRAM (day0)') +
  # ylim(0, 0.45) +
  theme_Publication()

ft_umap$Stress.scale <- scale(ft_umap$Stress)
ft_umap$Stress.scale <- ifelse(ft_umap$Stress.scale < -2, -2, ft_umap$Stress.scale)
ft_umap$Stress.scale <- ifelse(ft_umap$Stress.scale  > 2, 2, ft_umap$Stress.scale)
# ft_umap$Stress.scale <- ifelse(ft_umap$cell_id %in% metadata.day10_DABTRAM$cell_id, ft_umap$Interferon.MHC.II..I..scale, NA)
ft_umap <- ft_umap[order(ft_umap$Stress.scale),]
ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = Stress.scale )) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  theme_Publication()

ggsave(paste0(figure_dir, 'Fig5E.UMAP.Stress.pdf'), width = 4.4, height = 3)


p1 <- ggplot(ft_umap, aes(x = fatepotential_DABTRAM_d0_d10, y = Stress)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', alpha = 0.2) +
  stat_cor(method = 'spearman') +
  ggtitle('DABTRAM (day0)') +
  theme_Publication()

p2 <- ggplot(ft_umap, aes(x = fatepotential_DABTRAM_d10_w5, y = Stress)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', alpha = 0.2) +
  stat_cor(method = 'spearman') +
  ggtitle('DABTRAM (day10)') +
  theme_Publication()

p3 <- grid.arrange(p1, p2, ncol = 1)
ggsave(paste0(figure_dir, 'Fig5E.Stress.fatepotential.pdf'), p3, width = 3.5, height = 6.5)




# ==========================================================================================

remove_unassigned_cells <- TRUE
# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
score.df <- read.csv(paste0(data_dir, 'saver.GavishMP.UCellScores.csv'))

all_data[['fasttopic_COCL2']] <- all_data_fasttopic_COCL2
all_data[["ft.COCL2.umap"]] <- all_data_ft_COCL2_umap

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
# wrangle
# ==============================================================================

# ft umap
ft_umap <- all_data@reductions[["ft.COCL2.umap"]]@cell.embeddings
ft_umap <- as.data.frame(ft_umap)
ft_umap$cell_id <- rownames(ft_umap)

ft_umap <- merge(ft_umap, score.df, by = 'cell_id')

# metadata
metadata <- all_data@meta.data
metadata$cell_id <- rownames(metadata)
metadata.day10_COCL2 <- metadata[metadata$dataset == 'day10_COCL2', ]
metadata.week5_COCL2 <- metadata[metadata$dataset == 'week5_COCL2', ]
metadata.rest <- metadata[metadata$dataset != 'day10_COCL2', ]

# ==============================================================================
# plot
# ==============================================================================

ft_umap <- ft_umap[order(ft_umap$fatepotential_COCL2_d10_w5),]
ft_umap$fatepotential_COCL2_d10_w5 <- ifelse(ft_umap$fatepotential_COCL2_d10_w5 < -1, -1, ft_umap$fatepotential_COCL2_d10_w5)
ft_umap$fatepotential_COCL2_d10_w5 <- ifelse(ft_umap$fatepotential_COCL2_d10_w5 > 1, 1, ft_umap$fatepotential_COCL2_d10_w5)
rownames(ft_umap) <- ft_umap$cell_id
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.day10_COCL2$cell_id), ]
ggplot(ft_umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = fatepotential_COCL2_d10_w5)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "red", mid = "bisque", high = "blue", midpoint = -0.2, na.value = 'gray') +
  theme_Publication()

ft_umap$EMT.IV.scale <- scale(ft_umap$EMT.IV)
ft_umap$EMT.IV.scale <- ifelse(ft_umap$EMT.IV.scale < -1, -1, ft_umap$EMT.IV.scale)
ft_umap$EMT.IV.scale <- ifelse(ft_umap$EMT.IV.scale > 1, 1, ft_umap$EMT.IV.scale)

ft_umap <- ft_umap[order(ft_umap$EMT.IV.scale),]
ggplot(ft_umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = EMT.IV.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('COCL2 (all)') +
  theme_Publication()

ft_umap$Hypoxia.scale <- scale(ft_umap$Hypoxia)
ft_umap$Hypoxia.scale <- ifelse(ft_umap$Hypoxia.scale < -1, -1, ft_umap$Hypoxia.scale)
ft_umap$Hypoxia.scale <- ifelse(ft_umap$Hypoxia.scale > 1, 1, ft_umap$Hypoxia.scale)
ft_umap$Hypoxia.w5.scale <- ifelse(ft_umap$cell_id %in% metadata.week5_COCL2$cell_id, ft_umap$Hypoxia.scale, NA)
ft_umap$Hypoxia.d10.scale <- ifelse(ft_umap$cell_id %in% metadata.day10_COCL2$cell_id, ft_umap$Hypoxia.scale, NA)

ft_umap <- ft_umap[order(ft_umap$Hypoxia.scale),]
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.week5_COCL2$cell_id), ]
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.day10_COCL2$cell_id), ]
ggplot(ft_umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = Hypoxia.d10.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('COCL2 (all)') +
  theme_Publication()


ft_umap$Interferon.MHC.II..II..scale <- scale(ft_umap$Interferon.MHC.II..II.)
ft_umap$Interferon.MHC.II..II..scale <- ifelse(ft_umap$Interferon.MHC.II..II..scale < -1, -1, ft_umap$Interferon.MHC.II..II..scale)
ft_umap$Interferon.MHC.II..II..scale <- ifelse(ft_umap$Interferon.MHC.II..II..scale > 1, 1, ft_umap$Interferon.MHC.II..II..scale)
ft_umap$Interferon.MHC.II..II..scale <- ifelse(ft_umap$cell_id %in% metadata.day10_COCL2$cell_id, ft_umap$Interferon.MHC.II..II..scale, NA)

ft_umap <- ft_umap[order(ft_umap$Interferon.MHC.II..II..scale),]
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.day10_COCL2$cell_id), ]
ggplot(ft_umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = Interferon.MHC.II..II..scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('COCL2 (day10)') +
  theme_Publication()


ft_umap$Secreted.II.scale <- scale(ft_umap$Secreted.II)
ft_umap$Secreted.II.scale <- ifelse(ft_umap$Secreted.II.scale < -1, -1, ft_umap$Secreted.II.scale)
ft_umap$Secreted.II.scale <- ifelse(ft_umap$Secreted.II.scale > 1, 1, ft_umap$Secreted.II.scale)
ft_umap$Secreted.II.scale <- ifelse(ft_umap$cell_id %in% metadata.day10_COCL2$cell_id, ft_umap$Secreted.II.scale, NA)

ft_umap <- ft_umap[order(ft_umap$Secreted.II.scale),]
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.day10_COCL2$cell_id), ]
ggplot(ft_umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = Secreted.II.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('COCL2 (day10)') +
  theme_Publication()

ft_umap$Stress.scale <- ifelse(ft_umap$cell_id %in% metadata.day10_COCL2$cell_id, ft_umap$Stress, NA)
ft_umap$Stress.scale <- scale(ft_umap$Stress.scale)
ft_umap$Stress.scale <- ifelse(ft_umap$Stress.scale < -1, -1, ft_umap$Stress.scale)
ft_umap$Stress.scale <- ifelse(ft_umap$Stress.scale > 1, 1, ft_umap$Stress.scale)

ft_umap <- ft_umap[order(ft_umap$Stress.scale),]
ft_umap <- ft_umap[c(metadata.rest$cell_id, metadata.day10_COCL2$cell_id), ]
ggplot(ft_umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = Stress.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  ggtitle('COCL2 (day10)') +
  theme_Publication()


# =====================================================================================

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
score.df <- read.csv(paste0(data_dir, 'saver.GavishMP.UCellScores.csv'))

remove_unassigned_cells <- TRUE

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

mps <- c('Cell.Cycle...G2.M', 'Cell.Cycle...G1.S', 'Cell.Cycle.HMG.rich', 'Chromatin',
         'Stress', 'Hypoxia', 'Stress..in.vitro.', 'Proteasomal.degradation', 'Unfolded.protein.response',
         'Protein.maturation', 'Translation.initiation', 'EMT.I', 'EMT.II', 'EMT.III', 'EMT.IV', 
         'Interferon.MHC.II..I.', 'Interferon.MHC.II..II.', 'Epithelial.Senescence', 'MYC', 'Respiration',
         'Secreted.I', 'Secreted.II', 'Skin.pigmentation')

# ==============================================================================
# Correlate
# ==============================================================================
metadata <- all_data@meta.data
metadata$cell_id <- rownames(metadata)

metadat.day0 <- metadata[metadata$dataset == 'day0',]
metadat.day10_DABTRAM <- metadata[metadata$dataset == 'day10_DABTRAM',]
metadat.day10_COCL2 <- metadata[metadata$dataset == 'day10_COCL2',]
metadat.day10_CIS <- metadata[metadata$dataset == 'day10_CIS',]
metadat.week5_DABTRAM <- metadata[metadata$dataset == 'week5_DABTRAM',]
metadat.week5_COCL2 <- metadata[metadata$dataset == 'week5_COCL2',]
metadat.week5_CIS <- metadata[metadata$dataset == 'week5_CIS',]

rownames(score.df) <- score.df$cell_id
score.df <- score.df[, -1]
score.df <- score.df[, mps]
score.df.day0 <- score.df[metadat.day0$cell_id, ]
score.df.day10_DABTRAM <- score.df[metadat.day10_DABTRAM$cell_id, ]
score.df.day10_COCL2 <- score.df[metadat.day10_COCL2$cell_id, ]
score.df.day10_CIS <- score.df[metadat.day10_CIS$cell_id, ]
score.df.week5_DABTRAM <- score.df[metadat.week5_DABTRAM$cell_id, ]
score.df.week5_COCL2 <- score.df[metadat.week5_COCL2$cell_id, ]
score.df.week5_CIS <- score.df[metadat.week5_CIS$cell_id, ]

# perform pairwise correlation
correlation_day0 <- cor(score.df.day0)
pheatmap(correlation_day0, cluster_rows = T, cluster_cols = T,  main = 'Day0')


correlation_day10_COCL2 <- cor(score.df.day10_COCL2)
pheatmap(correlation_day10_COCL2, cluster_rows = T, cluster_cols = T,  main = 'Day10 COCL2')

correlation_day10_CIS <- cor(score.df.day10_CIS)
pheatmap(correlation_day10_CIS, cluster_rows = T, cluster_cols = T,  main = 'Day10 CIS')

correlation_week5_COCL2 <- cor(score.df.week5_COCL2)
pheatmap(correlation_week5_COCL2, cluster_rows = T, cluster_cols = T,  main = 'Week5 COCL2')

correlation_week5_CIS <- cor(score.df.week5_CIS)
pheatmap(correlation_week5_CIS, cluster_rows = T, cluster_cols = T,  main = 'Week5 CIS')


correlation_day10_DABTRAM <- cor(score.df.day10_DABTRAM)
pdf(paste0(figure_dir, 'Supp_correlation_day10_DABTRAM.pdf'), width = 5.5, height = 5)
pheatmap(correlation_day10_DABTRAM, cluster_rows = T, cluster_cols = T,  main = 'Day10 DABTRAM',
         treeheight_row = 0, treeheight_col = 0)
dev.off()

correlation_week5_DABTRAM <- cor(score.df.week5_DABTRAM)
pdf(paste0(figure_dir, 'Supp_correlation_week5_DABTRAM.pdf'), width = 5.5, height = 5)
pheatmap(correlation_week5_DABTRAM, cluster_rows = T, cluster_cols = T,  main = 'Week5 DABTRAM',
         treeheight_row = 0, treeheight_col = 0)
dev.off()
