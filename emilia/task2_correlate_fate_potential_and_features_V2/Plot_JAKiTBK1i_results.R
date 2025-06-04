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

results_dir <- "/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Experiments/JAKi_TBK1i/"

# ==============================================================================
# Load data
# ==============================================================================
data <- read.csv(paste0(results_dir, "CJC020_Plate2_Counts.csv"))
metadat <- read.csv(paste0(results_dir, "CJC020_Full_Timecourse_Plate_1_metadata.csv"))

# ==============================================================================
# Wrangle data
# ==============================================================================
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
  filter(Elaspe.Time == '2.266666667')

count.baseline <- count.baseline[, c('Sample', 'Cell.Count')]

count.baseline <- count.baseline %>%
  rename(Cell.Count.d0 = Cell.Count)

data.m <- merge(data.m, count.baseline, by = "Sample")

data.m$Cell.Count <- as.numeric(data.m$Cell.Count)
data.m$Cell.Count.d0 <- as.numeric(data.m$Cell.Count.d0)
data.m$Cell.Count.log2FC <- log2(data.m$Cell.Count / data.m$Cell.Count.d0)

data.m.to.plot <- data.m[data.m$Elaspe.Time %in% c('314.2166667', '963.55'), ]
data.m.to.plot$Timepoint <- ifelse(data.m.to.plot$Elaspe.Time == '314.2166667', 'Day 13', 'Day 40')

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



# ===============================================================================
# Analysis 2
# ==============================================================================
data.m <- melt(data, id.vars = "Sample")

colnames(data.m) <- c('Sample', 'Elaspe.Time', 'Cell.Count')

data.m <- merge(data.m, metadat, by = "Sample")

# standardize time
data.m$Elaspe.Time.days <- as.numeric(as.character(data.m$Elaspe.Time))
data.m$Elaspe.Time.days <- data.m$Elaspe.Time.days / 24

# baseline cell count

count.baseline <- data.m %>% 
  filter(Elaspe.Time == '2.266666667')

count.baseline <- count.baseline[, c('Sample', 'Cell.Count')]

count.baseline <- count.baseline %>%
  rename(Cell.Count.d0 = Cell.Count)

data.m <- merge(data.m, count.baseline, by = "Sample")

data.day13 <- data.m[data.m$Elaspe.Time == '314.2166667', ]
data.day13$log2FC <- log2(data.day13$Cell.Count / data.day13$Cell.Count.d0)

# day13 cell count
count.d13 <- data.m %>% 
  filter(Elaspe.Time == '314.2166667')

count.d13 <- count.d13[, c('Sample', 'Cell.Count')]

count.d13 <- count.d13 %>%
  rename(Cell.Count.d13 = Cell.Count)

data.m <- merge(data.m, count.d13, by = "Sample")

data.day40 <- data.m[data.m$Elaspe.Time == '963.55', ]
data.day40$log2FC <- log2(data.day40$Cell.Count / data.day40$Cell.Count.d13)

# assemble
data.to_plot <- rbind(data.day13[, c('Sample', 'Elaspe.Time', 'Treatment', 'Elaspe.Time.days', 'log2FC')],
                      data.day40[, c('Sample', 'Elaspe.Time', 'Treatment', 'Elaspe.Time.days', 'log2FC')])

ggplot(data.to_plot, aes(x = Treatment, y = log2FC)) +
  geom_jitter() +
  stat_compare_means(comparisons = list(c('DMSO', 'JAKi_TBK1i')), method = 't.test') +
  facet_wrap(~Elaspe.Time.days, nrow = 1) 
