rm(list=ls())
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
metadat <- seurat_object@meta.data
metadat$cell_barcode <- rownames(metadat)

load(paste0(out_dir, 'PC9_time_course_d7_to_d14_high.RData'))
non_cycling_fate <- final_fit[['cell_imputed_score']]
non_cycling_fate <- 10**non_cycling_fate
non_cycling_fate <- as.data.frame(non_cycling_fate)
non_cycling_fate$cell_barcode <- rownames(non_cycling_fate)

load(paste0(out_dir, 'PC9_time_course_d7_to_d14_mid_low.RData'))
cycling_fate <- final_fit[['cell_imputed_score']]
cycling_fate <- 10**cycling_fate
cycling_fate <- as.data.frame(cycling_fate)
cycling_fate$cell_barcode <- rownames(cycling_fate)

fate_potentials <- merge(non_cycling_fate, cycling_fate, by = 'cell_barcode')
fate_potentials$non_cycling_bias <- fate_potentials$non_cycling_fate / (fate_potentials$non_cycling_fate + fate_potentials$cycling_fate)
fate_potentials$cycling_bias <- fate_potentials$cycling_fate / (fate_potentials$non_cycling_fate + fate_potentials$cycling_fate)

umap <- as.data.frame(seurat_object@reductions[["umap"]]@cell.embeddings)
umap$cell_barcode <- rownames(umap)

metadat <- merge(umap, metadat, by = 'cell_barcode')
metadat <- merge(metadat, fate_potentials, by = 'cell_barcode', all = T)
metadat <- metadat %>% drop_na(lineage_barcode)


# ==============================================================================
# plotting
# ==============================================================================

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(metadat, aes(x = umap_1, y = umap_2, color = cycling_bias)) +
  geom_point() +
  scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1)) +
  theme_bw()
  # scale_color_gradient2(low = 'blue', mid = 'lightyellow', high = 'red')

ggplot(metadat, aes(x = majority_fate, y = cycling_bias)) +
  geom_violin() +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  ggtitle('D7 -> D14') +
  theme_bw()


