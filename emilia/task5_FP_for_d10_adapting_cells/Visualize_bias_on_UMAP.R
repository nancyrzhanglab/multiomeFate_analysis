rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE


treatment <- 'CIS'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_', treatment, '.RData'))

all_data[[paste0("fasttopic.", treatment)]] <- eval(parse(text = paste0("all_data_fasttopic_", treatment)))
all_data[[paste0("ft.", treatment, ".umap")]] <- eval(parse(text = paste0("all_data_ft_", treatment, "_umap")))

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

df.bias <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_', treatment, '.csv'))

metadata <- all_data@meta.data
metadata$cell_id <- rownames(metadata)
metadata.day0 <- subset(metadata, dataset == 'day0')
metadata.nonday0 <- subset(metadata, dataset != 'day0')
# =============================================================================
# Wrangle data
# =============================================================================

ft.umap <- all_data@reductions[[paste0("ft.", treatment, ".umap")]]@cell.embeddings
ft.umap <- as.data.frame(ft.umap)
ft.umap$cell_id <- rownames(ft.umap)

ft.umap <- merge(ft.umap, df.bias, by = 'cell_id', all = T)

# ft.umap <- ft.umap[order(ft.umap$bias, decreasing = F),]
ft.umap.nonday0 <- subset(ft.umap, cell_id %in% metadata.nonday0$cell_id)
ft.umap.day0 <- subset(ft.umap, cell_id %in% metadata.day0$cell_id)
ft.umap.day0 <- ft.umap.day0[order(ft.umap.day0$bias, decreasing = F),]
ft.umap <- rbind(ft.umap.nonday0, ft.umap.day0)
ggplot(ft.umap, aes(x = ftCISumap_1, y = ftCISumap_2, color = bias)) +
  geom_point(size = 0.5) +
  scale_color_gradient(low = '#fed7b4', high = 'blue', na.value = "#DCDCDC") +
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave(paste0(out_dir, 'adapting_bias_umap_', treatment, '.png'), width = 4, height = 3.6, dpi = 300)
 