rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'


remove_unassigned_cells <- TRUE

theme_Publication<- function(base_size=14, base_family="sans") {
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

# CIS
umap.cis <- all_data[[paste0("ft.CIS.umap")]]@cell.embeddings
umap.cis <- as.data.frame(umap.cis) %>% drop_na()
umap.cis$cell_id <- rownames(umap.cis)
fp.day10_cis <- as.data.frame(all_data@misc[["fatepotential_CIS_d10_w5"]][["cell_imputed_score"]])
colnames(fp.day10_cis) <- 'fatepotential_CIS_d10_w5'
fp.day10_cis$cell_id <- rownames(fp.day10_cis)
umap.cis <- merge(umap.cis, fp.day10_cis, by = 'cell_id', all = T)
umap.cis <- merge(umap.cis, metadat[, c('cell_id','assigned_lineage', 'dataset')], by = 'cell_id', all = T)
umap.cis$dataset <- factor(umap.cis$dataset, levels = c('day0', 'week5_CIS', 'day10_CIS'))
umap.cis <- umap.cis[order(umap.cis$dataset), ]

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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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

ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap.pdf'), p1, height = 3, width = 3)

p1 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 0.5) +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'right',
        legend.direction = "vertical")
p1
ggsave(paste0(figure_dir, 'fatepotential_DABTRAM_d10_w5_umap_legend.pdf'), p1, height = 3, width = 3)

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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
ggsave(paste0(figure_dir, 'fatepotential_COCL2_d10_w5_umap.pdf'), p2, width = 3, height = 3)


# CIS
max_val <- stats::quantile(umap.cis$fatepotential_CIS_d10_w5, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(umap.cis$fatepotential_CIS_d10_w5,
                           probs = 0.01, 
                           na.rm = TRUE)                           
umap.cis$fatepotential_CIS_d10_w5_scaled <- scales::rescale(
  pmax(pmin(umap.cis$fatepotential_CIS_d10_w5, 
            max_val), 
       min_val),
  to = c(-1, 1)
)
p2 <- ggplot(umap.cis, aes(x = ftCISumap_1, y = ftCISumap_2)) +
  geom_point(aes(color = fatepotential_CIS_d10_w5_scaled), size = 0.5) +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
                        na.value = "#E0E0E0", midpoint = 0) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = 'none')
p2
ggsave(paste0(figure_dir, 'Supp_fatepotential_CIS_d10_w5_umap.pdf'), p2, width = 3, height = 3)

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
            fp.min = min(fatepotential_DABTRAM_d10_w5, na.rm = T),
            n = n())
fp.summary.dabtram <- merge(fp.summary.dabtram, lin.size.week5[, c('assigned_lineage', 'week5_DABTRAM')], by = 'assigned_lineage', all.x = T)


umap.dabtram.lin1 <- umap.dabtram %>% 
  filter(dataset == 'day10_DABTRAM') %>%
  filter(assigned_lineage == 'Lin77715')

p3 <- ggplot(umap.dabtram, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.dabtram.lin1, aes(color = fatepotential_DABTRAM_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
            fp.min = min(fatepotential_COCL2_d10_w5_scaled, na.rm = T),
            n = n())
fp.summary.cocl2 <- merge(fp.summary.cocl2, lin.size.week5[, c('assigned_lineage', 'week5_COCL2')], by = 'assigned_lineage', all.x = T)


umap.cocl2.lin1 <- umap.cocl2 %>% 
  filter(dataset == 'day10_COCL2') %>%
  filter(assigned_lineage == 'Lin120586')

p7 <- ggplot(umap.cocl2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(aes(color = fatepotential_COCL2_d10_w5_scaled), shape = NA) +
  geom_point(color = '#E0E0E0') +
  geom_point(data = umap.cocl2.lin1, aes(color = fatepotential_COCL2_d10_w5_scaled), size = 3) +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', 
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

# =============================================================================
# Plot violin
# =============================================================================

# DABTRAM

priming.lins.dabtram <- fp.summary.dabtram %>% 
  filter(n > 5) %>% 
  filter(week5_DABTRAM > 0) %>% 
  filter(fp.min > -0.5) %>% 
  arrange(desc(fp.mean)) %>%
  head(4)

variance.lins.dabtram <- fp.summary.dabtram %>% 
  filter(n > 5) %>% 
  filter(week5_DABTRAM > 0) %>% 
  filter(!assigned_lineage %in% priming.lins.dabtram$assigned_lineage) %>%
  arrange(desc(fp.var)) %>%
  head(4)

size.lins.dabtram <- fp.summary.dabtram %>% 
  filter(n > 5) %>% 
  filter(week5_DABTRAM > 0) %>% 
  arrange(desc(n)) %>%
  head(4)

death.lins.dabtram <- fp.summary.dabtram %>% 
  filter(n > 5) %>% 
  filter(is.na(week5_DABTRAM)) %>% 
  arrange(fp.mean) %>%
  head(4)

other.lins.dabtram <- fp.summary.dabtram %>% 
  filter(!assigned_lineage %in% c(priming.lins.dabtram$assigned_lineage, 
                                  variance.lins.dabtram$assigned_lineage, 
                                  size.lins.dabtram$assigned_lineage, 
                                  death.lins.dabtram$assigned_lineage))

fp.day10_dabtram$cell_id <- rownames(fp.day10_dabtram)
df <- merge(fp.day10_dabtram, metadat[, c('cell_id', 'assigned_lineage', 'dataset')], by = 'cell_id')
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% priming.lins.dabtram$assigned_lineage, 
                                   df$assigned_lineage, 'Other')
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% variance.lins.dabtram$assigned_lineage, 
                                   df$assigned_lineage, df$assigned_lineage_plot)
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% size.lins.dabtram$assigned_lineage,
                                   df$assigned_lineage, df$assigned_lineage_plot)
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% death.lins.dabtram$assigned_lineage,
                                   df$assigned_lineage, df$assigned_lineage_plot)


df$assigned_lineage_plot <- factor(df$assigned_lineage_plot, levels = c(priming.lins.dabtram$assigned_lineage, 
                                                                        variance.lins.dabtram$assigned_lineage, 
                                                                        size.lins.dabtram$assigned_lineage, 
                                                                        death.lins.dabtram$assigned_lineage, 
                                                                        'Other'))
col_vec <- c(rep("lightgray", length(unique(df$assigned_lineage_plot))-1), 'darkgray')
names(col_vec) <- levels(df$assigned_lineage_plot)
plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=assigned_lineage_plot, y=fatepotential_DABTRAM_d10_w5, fill = assigned_lineage_plot))
plot1 <- plot1 + ggplot2::geom_violin(trim=T, scale = "width")
plot1 <- plot1 + ggplot2::scale_fill_manual(values = col_vec) 
plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                      position=ggplot2::position_jitter(0.2), alpha = 0.5, size = 1)
plot1 <- plot1 + Seurat::NoLegend()
# plot1 <- plot1 + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))
plot1 <- plot1 + ggplot2::stat_summary(fun = median, geom = "crossbar", lwd = 3,
                                       width = 0.75, color = "#633895")
plot1 <- plot1 + geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5), linetype = 'dashed', color = 'gray', linewidth = 0.8)
plot1 <- plot1 + 
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = 'none',
                 axis.title = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.border = element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(color = "gray", linetype = 'dashed', linewidth = 0.8))

plot1

df2 <- rbind(priming.lins.dabtram, variance.lins.dabtram, size.lins.dabtram, death.lins.dabtram)
df2[nrow(df2)+1, ] <- c('Other', 0, 0, 0, 0, 0, NA)
df2$week5_DABTRAM <- as.numeric(df2$week5_DABTRAM)
df2$week5_DABTRAM[is.na(df2$week5_DABTRAM)] <- 0
df2$assigned_lineage <- factor(df2$assigned_lineage, levels = c(priming.lins.dabtram$assigned_lineage, 
                                                                 variance.lins.dabtram$assigned_lineage, 
                                                                 size.lins.dabtram$assigned_lineage, 
                                                                 death.lins.dabtram$assigned_lineage,
                                                                 'Other'))


inset_plot <- ggplot(df2, aes(x = assigned_lineage, y = 0, fill = log10(week5_DABTRAM))) +
  # geom_bar(stat = 'identity', fill = 'gray') +
  geom_tile(width = 0.8, height = 0.8) +
  coord_equal() +
  scale_fill_gradient(low = "#FFCDB2", high = "#FC2947") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.border = element_blank())
inset_plot

plot2 <- ggarrange(plot1, inset_plot, ncol = 1, nrow = 2, heights = c(4, 1), align = 'v')
plot2

ggplot2::ggsave(filename = paste0(figure_dir, "fatepotential-violinplot-DABTRAM.png"),
                plot2, width = 10, height = 4)
 

plot1 <- plot1 + 
  theme(axis.text.y = element_blank())
plot2 <- ggarrange(plot1, inset_plot, ncol = 1, nrow = 2, heights = c(4, 1), align = 'v')
plot2
ggplot2::ggsave(filename = paste0(figure_dir, "fatepotential-violinplot-DABTRAM_clean.pdf"),
                plot2, width = 11, height = 4)

inset_plot <- inset_plot +
  theme(legend.position = 'right')
ggplot2::ggsave(filename = paste0(figure_dir, "fatepotential-violinplot-DABTRAM_inset_legend.pdf"),
                inset_plot, width = 2, height = 4)
# COCL2

priming.lins.cocl2 <- fp.summary.cocl2 %>% 
  filter(n > 5) %>% 
  filter(week5_COCL2 > 0) %>% 
  filter(fp.min > -0.5) %>% 
  arrange(desc(fp.mean)) %>%
  head(4)

variance.lins.cocl2 <- fp.summary.cocl2 %>% 
  filter(n > 5) %>% 
  filter(week5_COCL2 > 0) %>% 
  filter(!assigned_lineage %in% priming.lins.cocl2$assigned_lineage) %>%
  arrange(desc(fp.var)) %>%
  head(4)

size.lins.cocl2 <- fp.summary.cocl2 %>% 
  filter(n > 5) %>% 
  filter(week5_COCL2 > 0) %>% 
  arrange(desc(n)) %>%
  head(4)

death.lins.cocl2 <- fp.summary.cocl2 %>% 
  filter(n > 5) %>% 
  filter(is.na(week5_COCL2)) %>% 
  arrange(fp.mean) %>%
  head(4)

other.lins.cocl2 <- fp.summary.cocl2 %>% 
  filter(!assigned_lineage %in% c(priming.lins.cocl2$assigned_lineage, 
                                  variance.lins.cocl2$assigned_lineage, 
                                  size.lins.cocl2$assigned_lineage, 
                                  death.lins.cocl2$assigned_lineage))

fp.day10_cocl2$cell_id <- rownames(fp.day10_cocl2)
df <- merge(fp.day10_cocl2, metadat[, c('cell_id', 'assigned_lineage', 'dataset')], by = 'cell_id')
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% priming.lins.cocl2$assigned_lineage, 
                                   df$assigned_lineage, 'Other')
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% variance.lins.cocl2$assigned_lineage, 
                                   df$assigned_lineage, df$assigned_lineage_plot)
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% size.lins.cocl2$assigned_lineage,
                                   df$assigned_lineage, df$assigned_lineage_plot)
df$assigned_lineage_plot <- ifelse(df$assigned_lineage %in% death.lins.cocl2$assigned_lineage,
                                   df$assigned_lineage, df$assigned_lineage_plot)


df$assigned_lineage_plot <- factor(df$assigned_lineage_plot, levels = c(priming.lins.cocl2$assigned_lineage, 
                                                                        variance.lins.cocl2$assigned_lineage, 
                                                                        size.lins.cocl2$assigned_lineage, 
                                                                        death.lins.cocl2$assigned_lineage, 
                                                                        'Other'))
col_vec <- c(rep("lightgray", length(unique(df$assigned_lineage_plot))-1), 'darkgray')
names(col_vec) <- levels(df$assigned_lineage_plot)
plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=assigned_lineage_plot, y=fatepotential_COCL2_d10_w5, fill = assigned_lineage_plot))
plot1 <- plot1 + ggplot2::geom_violin(trim=T, scale = "width")
plot1 <- plot1 + ggplot2::scale_fill_manual(values = col_vec) 
plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                      position=ggplot2::position_jitter(0.2), alpha = 0.1, size = 0.5)
plot1 <- plot1 + Seurat::NoLegend()
# plot1 <- plot1 + ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))
plot1 <- plot1 + ggplot2::stat_summary(fun = median, geom = "crossbar", lwd = 1,
                                       width = 0.75, color = "#6DC49C")
plot1 <- plot1 + geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5), linetype = 'dashed', color = 'gray', linewidth = 0.8)
plot1 <- plot1 + 
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = 'none',
                 axis.title = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.border = element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(color = "gray", linetype = 'dashed', linewidth = 0.8))

plot1

df2 <- rbind(priming.lins.cocl2, variance.lins.cocl2, size.lins.cocl2, death.lins.cocl2)
df2[nrow(df2)+1, ] <- c('Other', 0, 0, 0, 0, 0, NA)
df2$week5_COCL2 <- as.numeric(df2$week5_COCL2)
df2$week5_COCL2[is.na(df2$week5_COCL2)] <- 0
df2$assigned_lineage <- factor(df2$assigned_lineage, levels = c(priming.lins.cocl2$assigned_lineage, 
                                                                variance.lins.cocl2$assigned_lineage, 
                                                                size.lins.cocl2$assigned_lineage, 
                                                                death.lins.cocl2$assigned_lineage,
                                                                'Other'))


inset_plot <- ggplot(df2, aes(x = assigned_lineage, y = 0, fill = log10(week5_COCL2))) +
  # geom_bar(stat = 'identity', fill = 'gray') +
  geom_tile(width = 0.8, height = 0.8) +
  coord_equal() +
  scale_fill_gradient(low = "#FFCDB2", high = "#FC2947") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.border = element_blank())
inset_plot

plot2 <- ggarrange(plot1, inset_plot, ncol = 1, nrow = 2, heights = c(4, 1))

ggplot2::ggsave(filename = paste0(figure_dir, "fatepotential-violinplot-COCL2.png"),
                plot2, width = 8, height = 3)

plot1 <- plot1 + 
  theme(axis.text.y = element_blank())
plot2 <- ggarrange(plot1, inset_plot, ncol = 1, nrow = 2, heights = c(4, 1), align = 'v')
plot2
ggplot2::ggsave(filename = paste0(figure_dir, "fatepotential-violinplot-COCL2_clean.pdf"),
                plot2, width = 11, height = 4)

inset_plot <- inset_plot +
  theme(legend.position = 'right')
ggplot2::ggsave(filename = paste0(figure_dir, "fatepotential-violinplot-COCL2_inset_legend.pdf"),
                inset_plot, width = 2, height = 4)



# =============================================================================
# Lineage size comparison
# =============================================================================

imputed.lin.w5.size.CIS <- as.data.frame(all_data_fatepotential[["fatepotential_CIS_d10_w5"]][["lineage_imputed_count"]])
imputed.lin.w5.size.COCL2 <- as.data.frame(all_data_fatepotential[["fatepotential_COCL2_d10_w5"]][["lineage_imputed_count"]])
imputed.lin.w5.size.DABTRAM <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["lineage_imputed_count"]])

colnames(imputed.lin.w5.size.CIS) <- c('Imputed Week5 Lineage Size')
colnames(imputed.lin.w5.size.COCL2) <- c('Imputed Week5 Lineage Size')
colnames(imputed.lin.w5.size.DABTRAM) <- c('Imputed Week5 Lineage Size')

imputed.lin.w5.size.CIS$assigned_lineage <- rownames(imputed.lin.w5.size.CIS)
imputed.lin.w5.size.COCL2$assigned_lineage <- rownames(imputed.lin.w5.size.COCL2)
imputed.lin.w5.size.DABTRAM$assigned_lineage <- rownames(imputed.lin.w5.size.DABTRAM)

df.DABTRAM <- merge(imputed.lin.w5.size.DABTRAM, lin.size.week5[, c('assigned_lineage', 'week5_DABTRAM')], by = 'assigned_lineage', all.x = T)
df.COCL2 <- merge(imputed.lin.w5.size.COCL2, lin.size.week5[, c('assigned_lineage', 'week5_COCL2')], by = 'assigned_lineage', all.x = T)
df.CIS <- merge(imputed.lin.w5.size.CIS, lin.size.week5[, c('assigned_lineage', 'week5_CIS')], by = 'assigned_lineage', all.x = T)

ggplot(df.DABTRAM, aes(x = log10(week5_DABTRAM+1), y = log10(`Imputed Week5 Lineage Size`+1))) +
  geom_point(size = 2) +
  # stat_cor(size = 8, color = 'red', label.sep = '\n') +
  annotate("text", x = 0.8, y = 3, label = "R^2 == 0.95", color = 'red', size = 7, parse = TRUE) +
  xlab('log10(Observed Week5 Lineage Size)') +
  ylab('log10(Imputed Week5 Lineage Size)') +
  theme_Publication()


ggsave(filename = paste0(figure_dir, "fatepotential-lineagesize-DABTRAM.png"),
       width = 4, height = 4)

ggplot(df.DABTRAM, aes(x = log10(week5_DABTRAM+1), y = log10(`Imputed Week5 Lineage Size`+1))) +
  geom_point(size = 2) +
  # stat_cor(size = 8, color = 'red', label.sep = '\n') +
  # annotate("text", x = 0.8, y = 3, label = "R^2 == 0.95", color = 'red', size = 7, parse = TRUE) +
  xlab('') +
  ylab('') +
  theme_Clean()

ggsave(filename = paste0(figure_dir, "fatepotential-lineagesize-DABTRAM_clean.pdf"),
       width = 2, height = 2)

ggplot(df.COCL2, aes(x = log10(week5_COCL2+1), y = log10(`Imputed Week5 Lineage Size`+1))) +
  geom_point() +
  # stat_cor(size = 8, color = 'red', label.sep = '\n') +
  annotate("text", x = 0.8, y = 3, label = "R^2 == 0.99", color = 'red', size = 7, parse = TRUE) +
  xlab('log10(Observed Week5 Lineage Size)') +
  ylab('log10(Imputed Week5 Lineage Size)') +
  theme_Publication()  
ggsave(filename = paste0(figure_dir, "fatepotential-lineagesize-COCL2.png"),
       width = 4, height = 4)

ggplot(df.COCL2, aes(x = log10(week5_COCL2+1), y = log10(`Imputed Week5 Lineage Size`+1))) +
  geom_point() +
  # stat_cor(size = 8, color = 'red', label.sep = '\n') +
  # annotate("text", x = 0.8, y = 3, label = "R^2 == 0.99", color = 'red', size = 7, parse = TRUE) +
  xlab('') +
  ylab('') +
  theme_Clean()  
ggsave(filename = paste0(figure_dir, "fatepotential-lineagesize-COCL2_clean.pdf"),
       width = 2, height = 2)

ggplot(df.CIS, aes(x = log10(week5_CIS+1), y = log10(`Imputed Week5 Lineage Size`+1))) +
  geom_point(size = 2) +
  # stat_cor(size = 8, color = 'red', label.sep = '\n') +
  # annotate("text", x = 0.5, y = 3, label = "R^2 == 0.39", color = 'red', size = 7, parse = TRUE) +
  xlab('log10(Observed Week5 Lineage Size)') +
  ylab('log10(Imputed Week5 Lineage Size)') +
  theme_Publication()
ggsave(filename = paste0(figure_dir, "Supp_fatepotential-lineagesize-CIS.pdf"),
       width = 2, height = 2)
