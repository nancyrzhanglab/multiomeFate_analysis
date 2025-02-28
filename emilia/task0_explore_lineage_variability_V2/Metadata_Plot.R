rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

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
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

remove_unassigned_cells <- TRUE

treatment <- 'DABTRAM'
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_', treatment, '.RData'))

all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[['ft.DABTRAM.UMAP']] <- all_data_ft_DABTRAM_umap

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
# Wrangle
# =============================================================================
umap <- all_data@reductions[["ft.DABTRAM.UMAP"]]@cell.embeddings
umap <- as.data.frame(umap)
umap$cell_id <- rownames(umap)
umap <- umap %>% drop_na()

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

umap <- merge(umap, metadat[, c('cell_id', 'dataset')], by = 'cell_id')

# =============================================================================
# Plot
# =============================================================================
umap$dataset <- str_split_fixed(umap$dataset, '_', 2)[,1]
ggplot(umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = dataset)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(day0 = "darkgray",
                                day10 = "#9D85BE",
                                week5 = "#623594")) +
  theme_Publication()

    