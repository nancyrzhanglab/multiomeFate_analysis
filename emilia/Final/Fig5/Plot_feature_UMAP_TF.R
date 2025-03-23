rm(list=ls())
library(Seurat)
library(clusterProfiler)
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
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_week5_DABTRAM.RData'))

all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[["ft.DABTRAM.umap"]] <- all_data_ft_DABTRAM_umap

all_data[['chromVar.day0']] <- all_data_chromVar_day0
all_data[['chromVar.day10_DABTRAM']] <- all_data_chromVar_day10_DABTRAM
all_data[['chromVar.week5_DABTRAM']] <- all_data_chromVar_week5_DABTRAM


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
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM', ]

umap <- as.data.frame(all_data[['ft.DABTRAM.umap']]@cell.embeddings)
umap$cell_id <- rownames(umap)

cv.day0 <- as.data.frame(all_data@assays[["chromVar.day0"]]@data)[c('FOSL1'), ]
cv.day10 <- as.data.frame(all_data@assays[["chromVar.day10_DABTRAM"]]@data)[c('FOSL1'), ]
cv.week5 <- as.data.frame(all_data@assays[["chromVar.week5_DABTRAM"]]@data)[c('FOSL1'), ]

cv.day0 <- as.data.frame(t(cv.day0)) %>% drop_na()
cv.day10 <- as.data.frame(t(cv.day10)) %>% drop_na()
cv.week5 <- as.data.frame(t(cv.week5)) %>% drop_na()

cv.day0$cell_id <- rownames(cv.day0)
cv.day10$cell_id <- rownames(cv.day10)
cv.week5$cell_id <- rownames(cv.week5)

cv.day0 <- cv.day0[metadat.day0$cell_id, ]
cv.day10 <- cv.day10[metadat.day10_DABTRAM$cell_id, ]
cv.week5 <- cv.week5[metadat.week5_DABTRAM$cell_id, ]

cv <- rbind(cv.day0, cv.day10, cv.week5)
cv$FOSL1.scale <- scale(cv$FOSL1)

cv <- merge(cv, umap, by = 'cell_id')

cv <- cv[order(cv$FOSL1.scale), ]
quantile(cv$FOSL1.scale)
# cv$JUNB.scale <- ifelse(cv$JUNB.scale < -2, 0, cv$JUNB.scale)
ggplot(cv, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = FOSL1.scale)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', midpoint = 4) +
  theme_minimal()
