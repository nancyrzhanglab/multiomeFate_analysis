library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# ==============================================================================
# Read data
# ==============================================================================
results_day0 <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day0_processed_RNA_ANOVA_pvals.csv")
results_day10_dabtram <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day10_DABTRAM_processed_RNA_ANOVA_pvals.csv")
results_week5_dabtram <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/week5_DABTRAM_processed_RNA_ANOVA_pvals.csv")

results_day10_cocl2 <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day10_COCL2_processed_RNA_ANOVA_pvals.csv")
results_week5_cocl2 <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/week5_COCL2_processed_RNA_ANOVA_pvals.csv")

results_day10_cis <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day10_CIS_processed_RNA_ANOVA_pvals.csv")
results_week5_cis <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/week5_CIS_processed_RNA_ANOVA_pvals.csv")


# visualize p-values
par(mfrow=c(3,1))
hist(results_day0$p_val, breaks=100)
hist(results_day10_dabtram$p_val, breaks=100)
hist(results_week5_dabtram$p_val, breaks=100)

par(mfrow=c(3,1))
hist(results_day0$p_val, breaks=100)
hist(results_day10_cocl2$p_val, breaks=100)
hist(results_week5_cocl2$p_val, breaks=100)

par(mfrow=c(3,1))
hist(results_day0$p_val, breaks=100)
hist(results_day10_cis$p_val, breaks=100)
hist(results_week5_cis$p_val, breaks=100)

results <- results[order(results$neg_log10_pval,decreasing=TRUE),]
results$order <- c(1: nrow(results))
ggplot(results) +
  geom_point(aes(x = order, y = neg_log10_pval))

results_day0$timepoint <- 'day0'
results_day10_dabtram$timepoint <- 'day10_DABTRAM'
results_week5_dabtram$timepoint <- 'week5_DABTRAM'

results_day10_cocl2$timepoint <- 'day10_COCL2'
results_week5_cocl2$timepoint <- 'week5_COCL2'

results_day10_cis$timepoint <- 'day10_CIS'
results_week5_cis$timepoint <- 'week5_CIS'

to_plot <- rbind(results_day0, results_day10_dabtram, results_week5_dabtram,
                 results_day10_cocl2, results_week5_cocl2, results_day10_cis, results_week5_cis)
to_plot$neg_log10_pval <- (-1) * log10(to_plot$p_val)
to_plot$neg_log10_pval <- ifelse(to_plot$neg_log10_pval > 150, 150, to_plot$neg_log10_pval)

order <- c('day0', 'day10_DABTRAM', 'week5_DABTRAM', 'day10_COCL2', 'week5_COCL2', 'day10_CIS', 'week5_CIS')
ggplot(to_plot) +
  geom_jitter(aes(x= factor(timepoint, levels=order), y=neg_log10_pval), alpha=0.2, size=1, width=0.2) +
  geom_boxplot(aes(x=factor(timepoint, levels=order), y=neg_log10_pval, color=timepoint), alpha=0.8, width=0.5, outlier.shape = NA) +
  scale_color_manual(values=c('#A9A9A9', '#72c3db', '#72c3db', '#72c3db', '#0174BE','#0174BE', '#0174BE' ))+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(to_plot) +
  geom_violin(aes(x=timepoint, y=neg_log10_pval), scale='width')
