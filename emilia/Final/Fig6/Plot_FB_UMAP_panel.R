rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'


remove_unassigned_cells <- TRUE

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

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


# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_CIS.RData'))

all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap
all_data[[paste0("fasttopic.COCL2")]] <- all_data_fasttopic_COCL2
all_data[[paste0("ft.COCL2.umap")]] <- all_data_ft_COCL2_umap
all_data[[paste0("fasttopic.CIS")]] <- all_data_fasttopic_CIS
all_data[[paste0("ft.CIS.umap")]] <- all_data_ft_CIS_umap

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

df.bias.DABTRAM <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
df.bias.COCL2 <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_COCL2.csv'))
df.bias.CIS <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_CIS.csv'))

metadata <- all_data@meta.data
metadata$cell_id <- rownames(metadata)
metadata.day0 <- subset(metadata, dataset == 'day0')
metadata.nonday0 <- subset(metadata, dataset != 'day0')
# =============================================================================
# Wrangle data
# =============================================================================

ft.umap.DABTRAM <- all_data@reductions[[paste0("ft.DABTRAM.umap")]]@cell.embeddings
ft.umap.DABTRAM <- as.data.frame(ft.umap.DABTRAM)
ft.umap.DABTRAM$cell_id <- rownames(ft.umap.DABTRAM)

ft.umap.DABTRAM <- merge(ft.umap.DABTRAM, df.bias.DABTRAM, by = 'cell_id', all = T)
ft.umap.DABTRAM <- merge(ft.umap.DABTRAM, metadata[, c('cell_id', 'dataset')], by = 'cell_id')

# ft.umap <- ft.umap[order(ft.umap$bias, decreasing = F),]
ft.umap.DABTRAM.nonday0 <- subset(ft.umap.DABTRAM, cell_id %in% metadata.nonday0$cell_id)
ft.umap.DABTRAM.day0 <- subset(ft.umap.DABTRAM, cell_id %in% metadata.day0$cell_id)
ft.umap.DABTRAM.day0 <- ft.umap.DABTRAM.day0[order(ft.umap.DABTRAM.day0$bias, decreasing = F),]
ft.umap.DABTRAM <- rbind(ft.umap.DABTRAM.nonday0, ft.umap.DABTRAM.day0)

midpoint.DABTRAM <- (max(ft.umap.DABTRAM.day0$bias) + min(ft.umap.DABTRAM.day0$bias))/2

p1 <- ggplot(ft.umap.DABTRAM.nonday0, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 1, color = "#DCDCDC") +
  geom_point(data = ft.umap.DABTRAM.day0, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = bias), size = 1) +
  scale_color_gradient2(low = "#3EDBF0",
                        mid = "#FFB433",
                        high = "#A31D1D",
                        midpoint = midpoint.DABTRAM,
                        na.value = "#DCDCDC") +
  scale_fill_manual(values = dataset_colors) +
  theme_Publication()

# COCL2
ft.umap.COCL2 <- all_data@reductions[[paste0("ft.COCL2.umap")]]@cell.embeddings
ft.umap.COCL2 <- as.data.frame(ft.umap.COCL2)
ft.umap.COCL2$cell_id <- rownames(ft.umap.COCL2)

ft.umap.COCL2 <- merge(ft.umap.COCL2, df.bias.COCL2, by = 'cell_id', all = T)
ft.umap.COCL2 <- merge(ft.umap.COCL2, metadata[, c('cell_id', 'dataset')], by = 'cell_id')

ft.umap.COCL2.nonday0 <- subset(ft.umap.COCL2, cell_id %in% metadata.nonday0$cell_id)
ft.umap.COCL2.day0 <- subset(ft.umap.COCL2, cell_id %in% metadata.day0$cell_id)
ft.umap.COCL2.day0 <- ft.umap.COCL2.day0[order(ft.umap.COCL2.day0$bias, decreasing = F),]
ft.umap.COCL2 <- rbind(ft.umap.COCL2.nonday0, ft.umap.COCL2.day0)

midpoint.COCL2 <- (max(ft.umap.COCL2.day0$bias) + min(ft.umap.COCL2.day0$bias))/2

p2 <- ggplot(ft.umap.COCL2.nonday0, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 1, color = "#DCDCDC") +
  geom_point(data = ft.umap.COCL2.day0, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = bias), size = 1) +
  scale_color_gradient2(low = "#3EDBF0",
                        mid = "#FFB433",
                        high = "#A31D1D",
                        midpoint = midpoint.COCL2,
                        na.value = "#DCDCDC") +
  theme_Publication()

# CIS
ft.umap.CIS <- all_data@reductions[[paste0("ft.CIS.umap")]]@cell.embeddings
ft.umap.CIS <- as.data.frame(ft.umap.CIS)
ft.umap.CIS$cell_id <- rownames(ft.umap.CIS)

ft.umap.CIS <- merge(ft.umap.CIS, df.bias.CIS, by = 'cell_id', all = T)
ft.umap.CIS <- merge(ft.umap.CIS, metadata[, c('cell_id', 'dataset')], by = 'cell_id')

ft.umap.CIS.nonday0 <- subset(ft.umap.CIS, cell_id %in% metadata.nonday0$cell_id)
ft.umap.CIS.day0 <- subset(ft.umap.CIS, cell_id %in% metadata.day0$cell_id)
ft.umap.CIS.day0 <- ft.umap.CIS.day0[order(ft.umap.CIS.day0$bias, decreasing = F),]
ft.umap.CIS <- rbind(ft.umap.CIS.nonday0, ft.umap.CIS.day0)

midpoint.CIS <- (max(ft.umap.CIS.day0$bias) + min(ft.umap.CIS.day0$bias))/2

p3 <- ggplot(ft.umap.CIS.nonday0, aes(x = ftCISumap_1, y = ftCISumap_2)) +
  geom_point(size = 1, color = "#DCDCDC") +
  geom_point(data = ft.umap.CIS.day0, aes(x = ftCISumap_1, y = ftCISumap_2, color = bias), size = 1) +
  scale_color_gradient2(low = "#3EDBF0",
                        mid = "#FFB433",
                        high = "#A31D1D",
                        midpoint = midpoint.CIS,
                        na.value = "#DCDCDC") +
  theme_Publication()

p4 <- ggarrange(p1, p2, p3, ncol = 3)
p4

ggsave(paste0(figure_dir, 'adapting_bias_umap_panel.pdf'), p4, width = 9.5, height = 2.3)
 
# ==================================================================================
# =============================================================================
# Wrangle
# =============================================================================
metadata.day0 <- metadata[metadata$dataset == 'day0',][, c('cell_id', 'assigned_lineage')]
metadata.week5.DABTRAM <- metadata[metadata$dataset == 'week5_DABTRAM',][, c('cell_id', 'assigned_lineage')]
metadata.week5.COCL2 <- metadata[metadata$dataset == 'week5_COCL2',][, c('cell_id', 'assigned_lineage')]
metadata.week5.CIS <- metadata[metadata$dataset == 'week5_CIS',][, c('cell_id', 'assigned_lineage')]

metadata.day0$lineage_in_week5_DABTRAM <- metadata.day0$assigned_lineage %in% metadata.week5.DABTRAM$assigned_lineage
metadata.day0$lineage_in_week5_COCL2 <- metadata.day0$assigned_lineage %in% metadata.week5.COCL2$assigned_lineage
metadata.day0$lineage_in_week5_CIS <- metadata.day0$assigned_lineage %in% metadata.week5.CIS$assigned_lineage


fb.DABTRAM.day0 <- ft.umap.DABTRAM.day0[, c('cell_id', 'bias')]
colnames(fb.DABTRAM.day0) <- c('cell_id', 'bias.DABTRAM')
fb.COCL2.day0 <- ft.umap.COCL2.day0[, c('cell_id', 'bias')]
colnames(fb.COCL2.day0) <- c('cell_id', 'bias.COCL2')
fb.CIS.day0 <- ft.umap.CIS.day0[, c('cell_id', 'bias')]
colnames(fb.CIS.day0) <- c('cell_id', 'bias.CIS')

metadata.day0 <- merge(metadata.day0, fb.DABTRAM.day0, by = 'cell_id')
metadata.day0 <- merge(metadata.day0, fb.COCL2.day0, by = 'cell_id')
metadata.day0 <- merge(metadata.day0, fb.CIS.day0, by = 'cell_id')

lin.count.day0.for.week5.DABTRAM <- metadata.day0 %>% 
  select(assigned_lineage, lineage_in_week5_DABTRAM) %>%
  distinct() %>%
  group_by(lineage_in_week5_DABTRAM) %>% 
  summarize(count = n())

lin.count.day0.for.week5.DABTRAM.high.bias <- metadata.day0 %>% 
  filter(bias.DABTRAM > 0.5) %>%
  select(assigned_lineage, lineage_in_week5_DABTRAM) %>%
  distinct() %>%
  group_by(lineage_in_week5_DABTRAM) %>% 
  summarize(count = n())


# Compute percentages
lin.count.day0.for.week5.DABTRAM$fraction = lin.count.day0.for.week5.DABTRAM$count / sum(lin.count.day0.for.week5.DABTRAM$count)

# Compute the cumulative percentages (top of each rectangle)
lin.count.day0.for.week5.DABTRAM$ymax = cumsum(lin.count.day0.for.week5.DABTRAM$fraction)

# Compute the bottom of each rectangle
lin.count.day0.for.week5.DABTRAM$ymin = c(0, head(lin.count.day0.for.week5.DABTRAM$ymax, n=-1))

# Compute label position
lin.count.day0.for.week5.DABTRAM$labelPosition <- (lin.count.day0.for.week5.DABTRAM$ymax + lin.count.day0.for.week5.DABTRAM$ymin) / 2

lin.count.day0.for.week5.DABTRAM$label <- ifelse(lin.count.day0.for.week5.DABTRAM$lineage_in_week5_DABTRAM, 'Present', 'Absent')
# Make the plot
p5 <- ggplot(lin.count.day0.for.week5.DABTRAM, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=c("#3EDBF0", "#FFB433")) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()

# Compute percentages
lin.count.day0.for.week5.DABTRAM.high.bias$fraction = lin.count.day0.for.week5.DABTRAM.high.bias$count / sum(lin.count.day0.for.week5.DABTRAM.high.bias$count)

# Compute the cumulative percentages (top of each rectangle)
lin.count.day0.for.week5.DABTRAM.high.bias$ymax = cumsum(lin.count.day0.for.week5.DABTRAM.high.bias$fraction)

# Compute the bottom of each rectangle
lin.count.day0.for.week5.DABTRAM.high.bias$ymin = c(0, head(lin.count.day0.for.week5.DABTRAM.high.bias$ymax, n=-1))

# Compute label position
lin.count.day0.for.week5.DABTRAM.high.bias$labelPosition <- (lin.count.day0.for.week5.DABTRAM.high.bias$ymax + lin.count.day0.for.week5.DABTRAM.high.bias$ymin) / 2

lin.count.day0.for.week5.DABTRAM.high.bias$label <- ifelse(lin.count.day0.for.week5.DABTRAM.high.bias$lineage_in_week5_DABTRAM, 'Present', 'Absent')
# Make the plot
p6 <- ggplot(lin.count.day0.for.week5.DABTRAM.high.bias, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=c("#3EDBF0", "#FFB433")) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()

p7 <- ggarrange(p5, p6, ncol = 1, common.legend = T, legend.position = "none")

ggsave(paste0(figure_dir, 'day0_for_week5_DABTRAM_lineage_presence.pdf'), p7, width = 4, height = 5)


# COCL2
lin.count.day0.for.week5.COCL2 <- metadata.day0 %>%
  select(assigned_lineage, lineage_in_week5_COCL2) %>%
  distinct() %>%
  group_by(lineage_in_week5_COCL2) %>% 
  summarize(count = n())

midpoint <- (max(metadata.day0$bias.COCL2) + min(metadata.day0$bias.COCL2)) / 2
midpoint <- 0.043111934
lin.count.day0.for.week5.COCL2.high.bias <- metadata.day0 %>% 
  filter(bias.COCL2 > midpoint) %>%
  select(assigned_lineage, lineage_in_week5_COCL2) %>%
  distinct() %>%
  group_by(lineage_in_week5_COCL2) %>%
  summarize(count = n())


# Compute percentages
lin.count.day0.for.week5.COCL2$fraction = lin.count.day0.for.week5.COCL2$count / sum(lin.count.day0.for.week5.COCL2$count)

# Compute the cumulative percentages (top of each rectangle)
lin.count.day0.for.week5.COCL2$ymax = cumsum(lin.count.day0.for.week5.COCL2$fraction)

# Compute the bottom of each rectangle
lin.count.day0.for.week5.COCL2$ymin = c(0, head(lin.count.day0.for.week5.COCL2$ymax, n=-1))

# Compute label position
lin.count.day0.for.week5.COCL2$labelPosition <- (lin.count.day0.for.week5.COCL2$ymax + lin.count.day0.for.week5.COCL2$ymin) / 2

lin.count.day0.for.week5.COCL2$label <- ifelse(lin.count.day0.for.week5.COCL2$lineage_in_week5_COCL2, 'Present', 'Absent')
# Make the plot
ggplot(lin.count.day0.for.week5.COCL2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=c("#3EDBF0", "#FFB433")) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()

# Compute percentages
lin.count.day0.for.week5.COCL2.high.bias$fraction = lin.count.day0.for.week5.COCL2.high.bias$count / sum(lin.count.day0.for.week5.COCL2.high.bias$count)

# Compute the cumulative percentages (top of each rectangle)
lin.count.day0.for.week5.COCL2.high.bias$ymax = cumsum(lin.count.day0.for.week5.COCL2.high.bias$fraction)

# Compute the bottom of each rectangle
lin.count.day0.for.week5.COCL2.high.bias$ymin = c(0, head(lin.count.day0.for.week5.COCL2.high.bias$ymax, n=-1))

# Compute label position
lin.count.day0.for.week5.COCL2.high.bias$labelPosition <- (lin.count.day0.for.week5.COCL2.high.bias$ymax + lin.count.day0.for.week5.COCL2.high.bias$ymin) / 2

lin.count.day0.for.week5.COCL2.high.bias$label <- ifelse(lin.count.day0.for.week5.COCL2.high.bias$lineage_in_week5_COCL2, 'Present', 'Absent')
# Make the plot
ggplot(lin.count.day0.for.week5.COCL2.high.bias, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=c("#3EDBF0", "#FFB433")) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()


# CIS
lin.count.day0.for.week5.CIS <- metadata.day0 %>%
  select(assigned_lineage, lineage_in_week5_CIS) %>%
  distinct() %>%
  group_by(lineage_in_week5_CIS) %>% 
  summarize(count = n())

midpoint <- (max(metadata.day0$bias.CIS) + min(metadata.day0$bias.CIS)) / 2
midpoint <- 0.5
lin.count.day0.for.week5.CIS.high.bias <- metadata.day0 %>% 
  filter(bias.CIS > midpoint) %>%
  select(assigned_lineage, lineage_in_week5_CIS) %>%
  distinct() %>%
  group_by(lineage_in_week5_CIS) %>%
  summarize(count = n())


# Compute percentages
lin.count.day0.for.week5.CIS$fraction = lin.count.day0.for.week5.CIS$count / sum(lin.count.day0.for.week5.CIS$count)

# Compute the cumulative percentages (top of each rectangle)
lin.count.day0.for.week5.CIS$ymax = cumsum(lin.count.day0.for.week5.CIS$fraction)

# Compute the bottom of each rectangle
lin.count.day0.for.week5.CIS$ymin = c(0, head(lin.count.day0.for.week5.CIS$ymax, n=-1))

# Compute label position
lin.count.day0.for.week5.CIS$labelPosition <- (lin.count.day0.for.week5.CIS$ymax + lin.count.day0.for.week5.CIS$ymin) / 2

lin.count.day0.for.week5.CIS$label <- ifelse(lin.count.day0.for.week5.CIS$lineage_in_week5_CIS, 'Present', 'Absent')
# Make the plot
ggplot(lin.count.day0.for.week5.CIS, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=c("#3EDBF0", "#FFB433")) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()

# Compute percentages
lin.count.day0.for.week5.CIS.high.bias$fraction = lin.count.day0.for.week5.CIS.high.bias$count / sum(lin.count.day0.for.week5.CIS.high.bias$count)

# Compute the cumulative percentages (top of each rectangle)
lin.count.day0.for.week5.CIS.high.bias$ymax = cumsum(lin.count.day0.for.week5.CIS.high.bias$fraction)

# Compute the bottom of each rectangle
lin.count.day0.for.week5.CIS.high.bias$ymin = c(0, head(lin.count.day0.for.week5.CIS.high.bias$ymax, n=-1))

# Compute label position
lin.count.day0.for.week5.CIS.high.bias$labelPosition <- (lin.count.day0.for.week5.CIS.high.bias$ymax + lin.count.day0.for.week5.CIS.high.bias$ymin) / 2

lin.count.day0.for.week5.CIS.high.bias$label <- ifelse(lin.count.day0.for.week5.CIS.high.bias$lineage_in_week5_CIS, 'Present', 'Absent')
# Make the plot
ggplot(lin.count.day0.for.week5.CIS.high.bias, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values=c("#3EDBF0", "#FFB433")) +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()

