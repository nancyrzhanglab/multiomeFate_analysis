rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'


remove_unassigned_cells <- TRUE

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
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_CIS.RData'))

all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap
all_data[[paste0("fasttopic.COCL2")]] <- all_data_fasttopic_COCL2
all_data[[paste0("ft.COCL2.umap")]] <- all_data_ft_COCL2_umap
all_data[[paste0("fasttopic.CIS")]] <- all_data_fasttopic_CIS
all_data[[paste0("ft.CIS.umap")]] <- all_data_ft_CIS_umap
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
colnames(fp.day10_dabtram) <- 'fatepotential_DABTRAM_d10_w5'
fp.day10_dabtram$cell_id <- rownames(fp.day10_dabtram)

umap.dabtram <- merge(umap.dabtram, fp.day10_dabtram, by = 'cell_id', all = T)

umap.dabtram <- merge(umap.dabtram, metadat[, c('cell_id', 'assigned_lineage', 'dataset')], by = 'cell_id', all = T)
umap.dabtram$dataset <- factor(umap.dabtram$dataset, levels = c('day0', 'week5_DABTRAM', 'day10_DABTRAM'))
umap.dabtram <- umap.dabtram[order(umap.dabtram$dataset), ]

# COCL2
umap.cocl2 <- all_data[[paste0("ft.COCL2.umap")]]@cell.embeddings
umap.cocl2 <- as.data.frame(umap.cocl2) %>% drop_na()
umap.cocl2$cell_id <- rownames(umap.cocl2)

fp.day10_cocl2 <- as.data.frame(all_data@misc[["fatepotential_COCL2_d10_w5"]][["cell_imputed_score"]])
colnames(fp.day10_cocl2) <- 'fatepotential_COCL2_d10_w5'
fp.day10_cocl2$cell_id <- rownames(fp.day10_cocl2)

umap.cocl2 <- merge(umap.cocl2, fp.day10_cocl2, by = 'cell_id', all = T)

umap.cocl2 <- merge(umap.cocl2, metadat[, c('cell_id','assigned_lineage', 'dataset')], by = 'cell_id', all = T)
umap.cocl2$dataset <- factor(umap.cocl2$dataset, levels = c('day0', 'week5_COCL2', 'day10_COCL2'))
umap.cocl2 <- umap.cocl2[order(umap.cocl2$dataset), ]

# =============================================================================
# Plot UMAP
# =============================================================================
# DABTRAM
max_val <- stats::quantile(umap.dabtram$fatepotential_DABTRAM_d10_w5, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(umap.dabtram$fatepotential_DABTRAM_d10_w5, 
                           probs = 0.01, 
                           na.rm = TRUE)

umap.dabtram$fatepotential_DABTRAM_d10_w5_scaled <- scales::rescale(
  pmax(pmin(umap.dabtram$fatepotential_DABTRAM_d10_w5, 
            max_val), 
       min_val),
  to = c(-1, 1)
)

p1 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 0.5) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
p1

ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap.pdf'), p1, width = 4.5, height = 4.5)



# COCL2
max_val <- stats::quantile(umap.cocl2$fatepotential_COCL2_d10_w5, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(umap.cocl2$fatepotential_COCL2_d10_w5,
                           probs = 0.01, 
                           na.rm = TRUE)                           
umap.cocl2$fatepotential_COCL2_d10_w5_scaled <- scales::rescale(
  pmax(pmin(umap.cocl2$fatepotential_COCL2_d10_w5, 
            max_val), 
       min_val),
  to = c(-1, 1)
)

p2 <- ggplot(umap.cocl2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(aes(color = fatepotential_COCL2_d10_w5_scaled), size = 0.5) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
p2
ggsave(paste0(figure_dir, 'fatepotential_COCL2_d10_w5_umap.pdf'), p2, width = 4.5, height = 4.5)


# =============================================================================
# pick specific lineages to plot
# =============================================================================

# DABTRAM
fp.summary.dabtram <- umap.dabtram %>% 
  filter(dataset == 'day10_DABTRAM') %>% 
  group_by(assigned_lineage) %>% 
  summarize(fp.mean = mean(fatepotential_DABTRAM_d10_w5, na.rm = T),
            fp.var = var(fatepotential_DABTRAM_d10_w5, na.rm = T),
            fp.max = max(fatepotential_DABTRAM_d10_w5, na.rm = T),
            n = n())
fp.summary.dabtram <- merge(fp.summary.dabtram, lin.size.week5[, c('assigned_lineage', 'week5_DABTRAM')], by = 'assigned_lineage', all.x = T)


umap.dabtram.lin1 <- umap.dabtram %>% 
  filter(dataset == 'day10_DABTRAM') %>%
  filter(assigned_lineage == 'Lin77715')

p3 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.dabtram.lin1, aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap_lin_die.pdf'), p3, width = 4.5, height = 4.5)


umap.dabtram.lin2 <- umap.dabtram %>% 
  filter(dataset == 'day10_DABTRAM') %>%
  filter(assigned_lineage == 'Lin104681') #Lin130951

p4 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.dabtram.lin2, aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap_lin_variance.pdf'), p4, width = 4.5, height = 4.5)


umap.dabtram.lin3 <- umap.dabtram %>% 
  filter(dataset == 'day10_DABTRAM') %>%
  filter(assigned_lineage == 'Lin130951')

p5 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.dabtram.lin3, aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap_lin_priming.pdf'), p5, width = 4.5, height = 4.5)

umap.dabtram.lin4 <- umap.dabtram %>% 
  filter(dataset == 'day10_DABTRAM') %>%
  filter(assigned_lineage == 'Lin5351')

p6 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.dabtram.lin4, aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap_lin_size.pdf'), p6, width = 4.5, height = 4.5)


# COCL2
fp.summary.cocl2 <- umap.cocl2 %>% 
  filter(dataset == 'day10_COCL2') %>% 
  group_by(assigned_lineage) %>% 
  summarize(fp.mean = mean(fatepotential_COCL2_d10_w5_scaled, na.rm = T),
            fp.var = var(fatepotential_COCL2_d10_w5_scaled, na.rm = T),
            fp.max = max(fatepotential_COCL2_d10_w5_scaled, na.rm = T),
            n = n())
fp.summary.cocl2 <- merge(fp.summary.cocl2, lin.size.week5[, c('assigned_lineage', 'week5_COCL2')], by = 'assigned_lineage', all.x = T)


umap.cocl2.lin1 <- umap.cocl2 %>% 
  filter(dataset == 'day10_COCL2') %>%
  filter(assigned_lineage == 'Lin120586')

p7 <- ggplot(umap.cocl2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(aes(color = fatepotential_COCL2_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.cocl2.lin1, aes(color = fatepotential_COCL2_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_COCL2_d10_w5_umap_lin_die.pdf'), p7, width = 4.5, height = 4.5)




umap.cocl2.lin2 <- umap.cocl2 %>% 
  filter(dataset == 'day10_COCL2') %>%
  filter(assigned_lineage == 'Lin43999')

p8 <- ggplot(umap.cocl2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(aes(color = fatepotential_COCL2_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.cocl2.lin2, aes(color = fatepotential_COCL2_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_COCL2_d10_w5_umap_lin_variance.pdf'), p8, width = 4.5, height = 4.5)


umap.cocl2.lin3 <- umap.cocl2 %>% 
  filter(dataset == 'day10_COCL2') %>%
  filter(assigned_lineage == 'Lin104509')

p9 <- ggplot(umap.cocl2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(aes(color = fatepotential_COCL2_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.cocl2.lin3, aes(color = fatepotential_COCL2_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_COCL2_d10_w5_umap_lin_priming.pdf'), p9, width = 4.5, height = 4.5)


umap.cocl2.lin4 <- umap.cocl2 %>% 
  filter(dataset == 'day10_COCL2') %>%
  filter(assigned_lineage == 'Lin82418')

p10 <- ggplot(umap.cocl2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(aes(color = fatepotential_COCL2_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.cocl2.lin4, aes(color = fatepotential_COCL2_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'red', mid = 'bisque', high = 'blue', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
ggsave(paste0(figure_dir, 'fatepotential_COCL2_d10_w5_umap_lin_size.pdf'), p10, width = 4.5, height = 4.5)

