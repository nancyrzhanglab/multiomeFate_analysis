library(stringr)
library(gridExtra)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
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

dataset_colors <- c(day0 = "darkgray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

umap.rna <- readRDS(paste0(data_dir, 'Writeup6b_umap_rna.rds'))@cell.embeddings
umap.atac <- readRDS(paste0(data_dir, 'Writeup6b_umap_atac.rds'))@cell.embeddings
metadat <- read.csv(paste0(data_dir, 'Writeup6b_all-data_metadata.csv'), row.names = 1)
metadat$cell_id <- rownames(metadat)

umap.rna <- as.data.frame(umap.rna)
umap.rna$cell_id <- rownames(umap.rna)

umap.atac <- as.data.frame(umap.atac)
umap.atac$cell_id <- rownames(umap.atac)

umap.rna <- merge(umap.rna, metadat, by = 'cell_id')

p1 <- ggplot() +
      geom_point(data = subset(umap.rna, dataset %in% c('day0', 'day10_CIS', 'week5_CIS', 'day10_COCL2', 'week5_COCL2')), 
                 aes(x = UMAP_1, y = UMAP_2), size = 0.1, color = '#FFEDFA') +
      geom_point(data = subset(umap.rna, dataset %in% c('day0', 'day10_DABTRAM', 'week5_DABTRAM')),
                 aes( x = UMAP_1, y = UMAP_2, color = dataset), size = 0.1) +
      scale_color_manual(values = dataset_colors) +
      ggtitle('DABTRAM') +
      theme_Publication()

p2 <- ggplot() +
  geom_point(data = subset(umap.rna, dataset %in% c('day0', 'day10_CIS', 'week5_CIS', 'day10_DABTRAM', 'week5_DABTRAM')), 
             aes(x = UMAP_1, y = UMAP_2), size = 0.1, color = '#FFEDFA') +
  geom_point(data = subset(umap.rna, dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2')),
             aes( x = UMAP_1, y = UMAP_2, color = dataset), size = 0.1) +
  scale_color_manual(values = dataset_colors) +
  ggtitle('COCL2') +
  theme_Publication()

p3 <- ggplot() +
  geom_point(data = subset(umap.rna, dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2', 'day10_DABTRAM', 'week5_DABTRAM')), 
             aes(x = UMAP_1, y = UMAP_2), size = 0.1, color = '#FFEDFA') +
  geom_point(data = subset(umap.rna, dataset %in% c('day0', 'day10_CIS', 'week5_CIS')),
             aes( x = UMAP_1, y = UMAP_2, color = dataset), size = 0.1) +
  scale_color_manual(values = dataset_colors) +
  ggtitle('CIS') +
  theme_Publication()

p4 <- grid.arrange(p1, p2, p3, ncol = 3)

ggsave('~/Downloads/UMAP_RNA_sm.png', p4, width = 15, height = 3.5, dpi = 300)


# ATAC
umap.atac <- merge(umap.atac, metadat, by = 'cell_id')

p5 <- ggplot() +
  geom_point(data = subset(umap.atac, dataset %in% c('day10_CIS', 'week5_CIS', 'day10_COCL2', 'week5_COCL2')), 
             aes(x = atacumap_1, y = atacumap_2), size = 0.1, color = '#FFEDFA') +
  geom_point(data = subset(umap.atac, dataset %in% c('day0', 'day10_DABTRAM', 'week5_DABTRAM')),
             aes(x = atacumap_1, y = atacumap_2, color = dataset), size = 0.1) +
  scale_color_manual(values = dataset_colors) +
  ggtitle('DABTRAM') +
  theme_Publication()

p6 <- ggplot() +
  geom_point(data = subset(umap.atac, dataset %in% c('day10_CIS', 'week5_CIS', 'day10_DABTRAM', 'week5_DABTRAM')), 
             aes(x = atacumap_1, y = atacumap_2), size = 0.1, color = '#FFEDFA') +
  geom_point(data = subset(umap.atac, dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2')),
             aes( x = atacumap_1, y = atacumap_2, color = dataset), size = 0.1) +
  scale_color_manual(values = dataset_colors) +
  ggtitle('COCL2') +
  theme_Publication()

p7 <- ggplot() +
  geom_point(data = subset(umap.atac, dataset %in% c('day10_COCL2', 'week5_COCL2', 'day10_DABTRAM', 'week5_DABTRAM')), 
             aes(x = atacumap_1, y = atacumap_2), size = 0.1, color = '#FFEDFA') +
  geom_point(data = subset(umap.atac, dataset %in% c('day0', 'day10_CIS', 'week5_CIS')),
             aes( x = atacumap_1, y = atacumap_2, color = dataset), size = 0.1) +
  scale_color_manual(values = dataset_colors) +
  ggtitle('CIS') +
  theme_Publication()

p8 <- grid.arrange(p5, p6, p7, ncol = 3)

ggsave('~/Downloads/UMAP_ATAC_sm.png', p8, width = 15, height = 3.5, dpi = 300)



