rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

treatment <- 'DABTRAM'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_wnn.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_CIS.RData'))

all_data[["wnn.umap"]] <- all_data_wnn
all_data[[paste0("ft.DABTRAM.umap")]] <- eval(parse(text = paste0("all_data_ft_DABTRAM_umap")))
all_data[[paste0("ft.COCL2.umap")]] <- eval(parse(text = paste0("all_data_ft_COCL2_umap")))
all_data[[paste0("ft.CIS.umap")]] <- eval(parse(text = paste0("all_data_ft_CIS_umap")))

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

final_fit.dabtram <- readRDS('~/Downloads/final_fit_d0_w5_DABTRAM_scaled.rds')
fp.dabtram <- as.data.frame(final_fit.dabtram[["cell_imputed_score"]])
colnames(fp.dabtram) <- c('DABTRAM_FP')
fp.dabtram$cell_id <- rownames(fp.dabtram)

final_fit.cocl2 <- readRDS('~/Downloads/final_fit_d0_w5_COCL2_scaled.rds')
fp.cocl2 <- as.data.frame(final_fit.cocl2[["cell_imputed_score"]])
colnames(fp.cocl2) <- c('COCL2_FP')
fp.cocl2$cell_id <- rownames(fp.cocl2)

final_fit.cis <- readRDS('~/Downloads/final_fit_d0_w5_CIS_scaled.rds')
fp.cis <- as.data.frame(final_fit.cis[["cell_imputed_score"]])
colnames(fp.cis) <- c('CIS_FP')
fp.cis$cell_id <- rownames(fp.cis)

# =============================================================================
# Subset to day0
# =============================================================================
all_data.day0 <- subset(all_data, dataset == 'day0')

wnn.umap.day0 <- as.data.frame(all_data.day0[["wnn.umap"]]@cell.embeddings)
wnn.umap.day0$cell_id <- rownames(wnn.umap.day0)

wnn.umap.day0 <- merge(wnn.umap.day0, fp.dabtram, by = 'cell_id')
wnn.umap.day0 <- merge(wnn.umap.day0, fp.cocl2, by = 'cell_id')
wnn.umap.day0 <- merge(wnn.umap.day0, fp.cis, by = 'cell_id')

# =============================================================================
# Plot
# =============================================================================

wnn.umap.day0 <- wnn.umap.day0[order(wnn.umap.day0$DABTRAM_FP), ]
p1 <- ggplot(wnn.umap.day0, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = DABTRAM_FP)) +
  geom_point() +
  scale_color_gradient2(low = 'red', mid = 'gray', high = 'blue') +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank()) +
  labs(title = 'DABTRAM_FP')

wnn.umap.day0 <- wnn.umap.day0[order(wnn.umap.day0$COCL2_FP), ]
p2 <- ggplot(wnn.umap.day0, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = COCL2_FP)) +
  geom_point() +
  scale_color_gradient2(low = 'red', mid = 'gray', high = 'blue') +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank()) +
  labs(title = 'COCL2_FP')

wnn.umap.day0 <- wnn.umap.day0[order(wnn.umap.day0$CIS_FP), ]
p3 <- ggplot(wnn.umap.day0, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = CIS_FP)) +
  geom_point() +
  scale_color_gradient2(low = 'red', mid = 'gray', high = 'blue') +
  theme_bw() +
  theme(panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank()) +
  labs(title = 'CIS_FP')

grid.arrange(p1, p2, p3, ncol = 3)
