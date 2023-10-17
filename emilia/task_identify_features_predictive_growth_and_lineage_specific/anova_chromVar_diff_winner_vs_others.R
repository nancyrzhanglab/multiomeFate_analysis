library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# ==============================================================================
# Read data
# ==============================================================================
results_day0 <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day0_ANOVA_ChromVar_pvals.csv")
results_day10_dabtram <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day10_DABTRAM_ANOVA_ChromVar_pvals.csv")
results_week5_dabtram <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/week5_DABTRAM_ANOVA_ChromVar_pvals.csv")

results_day10_cocl2 <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day10_COCL2_ANOVA_ChromVar_pvals.csv")
results_week5_cocl2 <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/week5_COCL2_ANOVA_ChromVar_pvals.csv")

results_day10_cis <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day10_CIS_ANOVA_ChromVar_pvals.csv")
results_week5_cis <- read.csv("~/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/week5_CIS_ANOVA_ChromVar_pvals.csv")


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
