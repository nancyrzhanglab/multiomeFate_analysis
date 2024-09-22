library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# in_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task1_ANOVA_lineage_specific_features/'
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/'
# ==============================================================================
# Read data
# ==============================================================================
anova.day0 <- read.csv(paste0(in_dir, 'day0_chromVAR_pvals.csv'))
anova.day10CIS <- read.csv(paste0(in_dir, 'day10_CIS_chromVAR_pvals.csv'))
anova.day10COCL2 <- read.csv(paste0(in_dir, 'day10_COCL2_chromVAR_pvals.csv'))
anova.day10DABTRAM <- read.csv(paste0(in_dir, 'day10_DABTRAM_chromVAR_pvals.csv'))
anova.week5CIS <- read.csv(paste0(in_dir, 'week5_CIS_chromVAR_pvals.csv'))
anova.week5COCL2 <- read.csv(paste0(in_dir, 'week5_COCL2_chromVAR_pvals.csv'))
anova.week5DABTRAM <- read.csv(paste0(in_dir, 'week5_DABTRAM_chromVAR_pvals.csv'))

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
dataset_colors <- c(day0_baseline = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

p = ggplot(to_plot, aes(x = factor(dataset, levels = order), y = neg_log10_pval)) +
  geom_violin(scale = 'width', aes(fill = dataset)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  # geom_jitter(width = 0.05, alpha = 0.01) +
  # stat_summary(fun=median, geom="point", size=2, color="red") +
  xlab('') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        legend.position = 'none')

ggsave(paste0(in_dir, 'ANOVA_pvals_chromVAR_violin_plot.png'), p, width = 6, height = 4.5, dpi = 300)

