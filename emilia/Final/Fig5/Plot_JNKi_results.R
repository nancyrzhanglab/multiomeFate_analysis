rm(list = ls())
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)

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

results_dir <- "/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Experiments/JNKi/"
figure_dir <- "/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig5/"

use.pal <- c('#4285F4', '#C00000')
names(use.pal) <- c('DMSO', 'JNKi')
# ==============================================================================
# Load data
# ==============================================================================
data <- read.csv(paste0(results_dir, "CJC011_Full_Timecourse_Plate_1.csv"))
metadat <- read.csv(paste0(results_dir, "CJC011_Full_Timecourse_Plate_1_metadata.csv"))

data2 <- read.csv(paste0(results_dir, "CJC011_Full_Timecourse_Plate_2.csv"))

# ==============================================================================
# Wrangle data
# ==============================================================================

# Data rep1
data <- as.data.frame(t(data))
colnames(data) <- data[1,]
data <- data[-1,]
data$Sample <- rownames(data)

data.m <- melt(data, id.vars = "Sample")
colnames(data.m) <- c('Sample', 'Elaspe.Time', 'Cell.Count')

data.m <- merge(data.m, metadat, by = "Sample")

# standardize time
data.m$Elaspe.Time.days <- as.numeric(as.character(data.m$Elaspe.Time))
data.m$Elaspe.Time.days <- data.m$Elaspe.Time.days / 24

# baseline cell count

count.baseline <- data.m %>% 
  filter(Elaspe.Time == '0')

count.baseline <- count.baseline[, c('Sample', 'Cell.Count')]

colnames(count.baseline) <- c('Sample', 'Cell.Count.d0')

data.m <- merge(data.m, count.baseline, by = "Sample")

data.m$Cell.Count <- as.numeric(data.m$Cell.Count)
data.m$Cell.Count.d0 <- as.numeric(data.m$Cell.Count.d0)
data.m$Cell.Count.log2FC <- log2(data.m$Cell.Count / data.m$Cell.Count.d0)

data.m.to.plot <- data.m[data.m$Elaspe.Time %in% c('288.9333333', '979.3'), ]
data.m.to.plot$Timepoint <- ifelse(data.m.to.plot$Elaspe.Time == '288.9333333', 'Day 12', 'Day 41')
data.m.to.plot$Rep <- 'Rep1'

# Data rep2
data2 <- as.data.frame(t(data2))
colnames(data2) <- data2[1,]
data2 <- data2[-1,]
data2$Sample <- rownames(data2)

data2.m <- melt(data2, id.vars = "Sample")
colnames(data2.m) <- c('Sample', 'Elaspe.Time', 'Cell.Count')
data2.m <- merge(data2.m, metadat, by = "Sample")

# standardize time
data2.m$Elaspe.Time.days <- as.numeric(as.character(data2.m$Elaspe.Time))
data2.m$Elaspe.Time.days <- data2.m$Elaspe.Time.days / 24

# baseline cell count
count.baseline2 <- data2.m %>% 
  filter(Elaspe.Time == '0')
count.baseline2 <- count.baseline2[, c('Sample', 'Cell.Count')]

count.baseline2 <- count.baseline2 %>%
  rename(Cell.Count.d0 = Cell.Count)

data2.m <- merge(data2.m, count.baseline2, by = "Sample")
data2.m$Cell.Count <- as.numeric(data2.m$Cell.Count)
data2.m$Cell.Count.d0 <- as.numeric(data2.m$Cell.Count.d0)
data2.m$Cell.Count.log2FC <- log2(data2.m$Cell.Count / data2.m$Cell.Count.d0)

data2.m.to.plot <- data2.m[data2.m$Elaspe.Time %in% c('267.5833333', '960'), ]
data2.m.to.plot$Timepoint <- ifelse(data2.m.to.plot$Elaspe.Time == '267.5833333', 'Day 12', 'Day 41')
data2.m.to.plot$Rep <- 'Rep2'

data.m.to.plot.all <- rbind(data.m.to.plot, data2.m.to.plot)

# ==============================================================================
# Plot
# ==============================================================================

ggplot(data.m, aes(x = Elaspe.Time.days, y = Cell.Count.log2FC, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # geom_line(aes(group = Treatment)) +
  geom_smooth(aes(group = Treatment), method = "gam", formula = y ~ s(x, k = 7, bs = "cs"), fill = "lightgrey") +
  geom_jitter(width = 0.35) +
  scale_color_manual(values = c('#B0B0B0', '#FF8400')) +
  labs(title = "CJC011_Full_Timecourse_Plate_1",
       x = "Time (days)",
       y = "log2FC Cell Count") +
  theme_Publication()

data.m.to.plot.all$Category <- paste0(data.m.to.plot$Timepoint, "-", data.m.to.plot$Treatment)

ggplot(data.m.to.plot.all, aes(x = Treatment, y = Cell.Count.log2FC)) +
  # geom_boxplot(aes(color = Treatment)) +
  geom_jitter(aes(color = Treatment), shape = 21, size = 5, width = 0.2, height = 0, stroke = 2) +
  facet_wrap(~Timepoint) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("DMSO", "JNKi")),
    label = "p.signif",
    label.y = 0.8  # Adjust label position
  ) +
  scale_color_manual(values = c('#808080', '#0F4C75')) +
  labs(x = '', y = 'log2FC cell count') +
  theme_Publication() +
  theme(legend.position = "none")

# ggsave(paste0(figure_dir, "Fig5G_boxplot_log2FC_cell_count_clean2.pdf"), width = 4.2, height = 3.3)
ggsave(paste0(figure_dir, "Fig5G_boxplot_log2FC_cell_count_clean2.pdf"), width = 5, height = 3.5)

