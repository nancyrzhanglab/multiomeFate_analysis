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
  
  