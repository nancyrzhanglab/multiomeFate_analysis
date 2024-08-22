library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/'

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


# ==============================================================================
# Read data
# ==============================================================================
anova.day0 <- read.csv(paste0(in_dir, 'day0_processed_RNA_ANOVA_pvals.csv'))
anova.day10CIS <- read.csv(paste0(in_dir, 'day10_CIS_processed_RNA_ANOVA_pvals.csv'))
anova.day10COCL2 <- read.csv(paste0(in_dir, 'day10_COCL2_processed_RNA_ANOVA_pvals.csv'))
anova.day10DABTRAM <- read.csv(paste0(in_dir, 'day10_DABTRAM_processed_RNA_ANOVA_pvals.csv'))
anova.week5CIS <- read.csv(paste0(in_dir, 'week5_CIS_processed_RNA_ANOVA_pvals.csv'))
anova.week5COCL2 <- read.csv(paste0(in_dir, 'week5_COCL2_processed_RNA_ANOVA_pvals.csv'))
anova.week5DABTRAM <- read.csv(paste0(in_dir, 'week5_DABTRAM_processed_RNA_ANOVA_pvals.csv'))

anova.day0$timepoint <- 'day0'
anova.day0$treatment <- 'baseline'

anova.day10CIS$timepoint <- 'day10'
anova.day10CIS$treatment <- 'CIS'

anova.day10COCL2$timepoint <- 'day10'
anova.day10COCL2$treatment <- 'COCL2'

anova.day10DABTRAM$timepoint <- 'day10'
anova.day10DABTRAM$treatment <- 'DABTRAM'

anova.week5CIS$timepoint <- 'week5'
anova.week5CIS$treatment <- 'CIS'

anova.week5COCL2$timepoint <- 'week5'
anova.week5COCL2$treatment <- 'COCL2'

anova.week5DABTRAM$timepoint <- 'week5'
anova.week5DABTRAM$treatment <- 'DABTRAM'


to_plot <- rbind(anova.day0, anova.day10CIS)
to_plot <- rbind(to_plot, anova.day10COCL2)
to_plot <- rbind(to_plot, anova.day10DABTRAM)
to_plot <- rbind(to_plot, anova.week5CIS)
to_plot <- rbind(to_plot, anova.week5COCL2)
to_plot <- rbind(to_plot, anova.week5DABTRAM)

to_plot$dataset <- paste0(to_plot$timepoint, '_', to_plot$treatment)
to_plot$neg_log10_pval <- (-1) * log10(to_plot$p_val)
# ==============================================================================
# Plotting
# ==============================================================================
order <-  c("day0_baseline", "day10_DABTRAM", "day10_COCL2", "day10_CIS", "week5_DABTRAM", "week5_COCL2", "week5_CIS")

ggplot(to_plot, aes(x = dataset, y = neg_log10_pval)) +
  geom_violin(scale = 'width', aes(fill = timepoint)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  # geom_jitter(width = 0.05, alpha = 0.01) +
  scale_fill_manual(values = c('gray', '#a1dfff', '#2ba4fc')) +
  # stat_summary(fun=median, geom="point", size=2, color="red") +
  xlab('') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

ggsave('~/Downloads/ANOVA_pvals.png', dpi = 300, width = 5, heigh = 4)











