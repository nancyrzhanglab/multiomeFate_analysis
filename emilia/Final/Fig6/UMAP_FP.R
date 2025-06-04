rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'


remove_unassigned_cells <- TRUE


theme_Clean<- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            axis.text = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour=NA,fill=NA),
            panel.background = element_blank(),
            strip.text = element_text(face="bold")
    ))
}
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))

all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap
all_data@misc <- all_data_fatepotential

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
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.week5 <- metadat[grepl('week5', metadat$dataset), ]
lin.size.week5 <- as.data.frame(table(metadat.week5$assigned_lineage, metadat.week5$dataset))
colnames(lin.size.week5) <- c('assigned_lineage', 'dataset', 'size')
lin.size.week5 <- lin.size.week5[lin.size.week5$size > 0, ]
lin.size.week5 <- spread(lin.size.week5, key = 'dataset', value = 'size')

# DABTRAM
umap.dabtram <- all_data[[paste0("ft.DABTRAM.umap")]]@cell.embeddings
umap.dabtram <- as.data.frame(umap.dabtram) %>% drop_na()
umap.dabtram$cell_id <- rownames(umap.dabtram)

fp.day10_dabtram <- as.data.frame(all_data@misc[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(fp.day10_dabtram) <- 'fatepotential'
fp.day10_dabtram$cell_id <- rownames(fp.day10_dabtram)

max_val <- stats::quantile(fp.day10_dabtram$fatepotential, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(fp.day10_dabtram$fatepotential, 
                           probs = 0.01, 
                           na.rm = TRUE)

fp.day10_dabtram$fatepotential <- scales::rescale(
  pmax(pmin(fp.day10_dabtram$fatepotential, 
            max_val), 
       min_val),
  to = c(-1, 1)
)


fp.day0_dabtram <- as.data.frame(all_data@misc[["fatepotential_DABTRAM_d0_d10"]][["cell_imputed_score"]])
colnames(fp.day0_dabtram) <- 'fatepotential'
fp.day0_dabtram$cell_id <- rownames(fp.day0_dabtram)

max_val <- stats::quantile(fp.day0_dabtram$fatepotential, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(fp.day0_dabtram$fatepotential, 
                           probs = 0.01, 
                           na.rm = TRUE)

fp.day0_dabtram$fatepotential <- scales::rescale(
  pmax(pmin(fp.day0_dabtram$fatepotential, 
            max_val), 
       min_val),
  to = c(-1, 1)
)


fp.dabtram <- rbind(fp.day0_dabtram, fp.day10_dabtram)

umap.dabtram <- merge(umap.dabtram, fp.dabtram, by = 'cell_id', all = T)

umap.dabtram <- merge(umap.dabtram, metadat[, c('cell_id', 'assigned_lineage', 'dataset')], by = 'cell_id', all = T)
umap.dabtram$dataset <- factor(umap.dabtram$dataset, levels = c('day0', 'week5_DABTRAM', 'day10_DABTRAM'))
umap.dabtram <- umap.dabtram[order(umap.dabtram$dataset), ]


ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential), size = 0.5) +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Clean() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

ggsave(paste0(figure_dir, 'fatepotential.png'), width = 3.5, height = 3.5)
