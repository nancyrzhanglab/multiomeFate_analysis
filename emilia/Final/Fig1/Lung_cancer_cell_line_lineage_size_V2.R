rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(hdrcde)
library(ggdensity)
library(RColorBrewer)
library(ComplexHeatmap)

in_dir <- '/Users/emiliac/Dropbox/MultiomeFate/data/ShafferLab/robert_2024_05_03_HCC827_gDNA_Sequencing_SMS008/R_Scripts_and_Plots/Organized_Data/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task9_lung_clone_size_robert/'
figure_dir <- '/Users/emiliac/Library/CloudStorage/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig1/'

# ==============================================================================
# Read data
# ==============================================================================
data.without.spike <- read.csv(paste0(in_dir, 'Data_Without_Spike_Ins.csv'), row.names = 1)
data.spike_in <- read.csv(paste0(in_dir, 'Spike_Ins_With_Data.csv'), row.names = 1)
meta.data <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Lung_robert_metadata.csv')

# ==============================================================================
# Wrangle data
# ==============================================================================
meta.data$Sample <- paste0('Sample_', meta.data$Sample)
shared_samples <- intersect(meta.data$Sample, colnames(data.without.spike))

rownames(meta.data) <- meta.data$Sample
meta.data$Sample_Name <- paste0(meta.data$Sample, '_', meta.data$Treatment, '_', meta.data$Time)

meta.data <- meta.data[colnames(data.without.spike), ]
colnames(data.without.spike) <- meta.data$Sample_Name

meta.data1 <- meta.data[colnames(data.spike_in)[5:16], ]
colnames(data.spike_in) <- c(colnames(data.spike_in)[1:4], meta.data1$Sample_Name) 

# ==============================================================================
# Day10
# ==============================================================================

day10.spike <- data.spike_in[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_15_Cisplatin_day10',
                                 'Sample_16_Cisplatin_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
# day10.spike <- log10(day10.spike + 1 + 10e-5)
min <- min(day10.spike)
thres <- min / 5

day10 <- data.without.spike[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10', 'Sample_15_Cisplatin_day10',
                                'Sample_16_Cisplatin_day10', 'Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')]
# day10[day10 < thres] <- 0
day10$all <- rowSums(day10)
day10 <- day10[day10$all > 0, ]
day10 <- log10(day10+1)
day10 <- subset(day10, select = -c(all))
colnames(day10) <- c('Osi_rep1', 'Osi_rep2', 'Cisplatin_rep1', 'Cisplatin_rep2', 'Cocl2_rep1', 'CoCl2_rep2')

ggpairs(
  day10,
  upper = list(continuous = wrap("cor", size = 6, color = 'black')),
  lower = list(continuous = wrap("points", size = 1, color = 'pink'))   # Keep scatter plots in upper triangle
) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14)
  )
ggsave(paste0(figure_dir, 'Supp_Day10_comparison_lung.png'), width = 8, height = 8, dpi = 600)

# ==============================================================================
# End
# ==============================================================================
end.spike <- data.spike_in[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                               'Sample_38_Cisplatin_END')]
min <- min(end.spike)
thres <- min / 5

end <- data.without.spike[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                              'Sample_38_Cisplatin_END')]

# end[end < thres] <- 0
end$all <- rowSums(end)
end <- end[end$all > 0, ]
end <- log10(end+1)
end <- subset(end, select = -c(all))
colnames(end) <- c('Osi_rep1', 'Osi_rep2', 'CoCl2', 'Cisplatin')

ggpairs(
  end,
  upper = list(continuous = wrap("cor", size = 6, color = 'black')),
  lower = list(continuous = wrap("points", size = 1, color = '#A53860'))   # Keep scatter plots in upper triangle
) +
  theme_bw()
ggsave(paste0(figure_dir, 'Supp_End_comparison_lung.png'), width = 5, height = 5, dpi = 600)


# ==============================================================================
# Day10 vs End
# ==============================================================================
colnames(day10) <- paste0(colnames(day10), '_Day10')
day10$lineage <- rownames(day10)

colnames(end) <- paste0(colnames(end), '_End')
end$lineage <- rownames(end)
day10.end <- merge(day10, end, by = 'lineage', all = TRUE)

day10.end.CoCl2 <- day10.end[, c('Cocl2_rep1_Day10', 'CoCl2_rep2_Day10', 'CoCl2_End')] %>% drop_na()
day10.end.CoCl2$all <- rowSums(day10.end.CoCl2)
day10.end.CoCl2 <- day10.end.CoCl2[day10.end.CoCl2$all > 0, ]
day10.end.CoCl2 <- subset(day10.end.CoCl2, select = -c(all))

ggpairs(
  day10.end.CoCl2,
  upper = list(continuous = wrap("cor", size = 6, color = 'black')),
  lower = list(continuous = wrap("points", size = 1, color = 'black'))   # Keep scatter plots in upper triangle
) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14)
  )


day10.end.Cisplatin <- day10.end[, c('Cisplatin_day10_mean', 'Sample_38_Cisplatin_END')] %>% drop_na()
day10.end.Cisplatin$all <- rowSums(day10.end.Cisplatin)
day10.end.Cisplatin <- day10.end.Cisplatin[day10.end.Cisplatin$all > 0, ]
p9 <- ggplot(day10.end.Cisplatin, aes(x = Cisplatin_day10_mean, y = Sample_38_Cisplatin_END)) +
  geom_point(color = 'black') +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 6, size = 6) +
  xlab('Cisplatin day10 (mean)') +
  ylab('Cisplatin end (mean)') +
  theme_bw()

day10.end.Osi <- day10.end[, c('Osi_day10_mean', 'Osi_end_mean')] %>% drop_na()
day10.end.Osi$all <- rowSums(day10.end.Osi)
day10.end.Osi <- day10.end.Osi[day10.end.Osi$all > 0, ]
p10 <- ggplot(day10.end.Osi, aes(x = Osi_day10_mean, y = Osi_end_mean)) +
  geom_point(color = 'black') +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 6, size = 6) +
  xlab('Osi day10 (mean)') +
  ylab('Osi end (mean)') +
  theme_bw()

p11 <- grid.arrange(p8, p9, p10, nrow = 1)


# ==============================================================================
# Day 10 and End
# ==============================================================================

colnames(day10) <- paste0(colnames(day10), '_Day10')
day10$lineage <- rownames(day10)

colnames(end) <- paste0(colnames(end), '_End')
end$lineage <- rownames(end)

day10.end <- merge(day10, end, by = 'lineage', all = TRUE)
cor.df <- data.frame(matrix(nrow = ncol(day10.end)-1, ncol = ncol(day10.end)-1))
colnames(cor.df) <- colnames(day10.end)[-1]
rownames(cor.df) <- colnames(day10.end)[-1]

for(i in 2:ncol(day10.end)) {
  for (j in 2:ncol(day10.end)) {
    cor.df[i-1, j-1] <- cor(day10.end[, i], day10.end[, j], method = 'pearson', use = 'pairwise.complete.obs')
  }
}

# remove lower triangle
cor.df[lower.tri(cor.df)] <- NA
cor.df['Osi_rep1_Day10', 'Cisplatin_End'] <- NA
cor.df['Osi_rep2_Day10', 'Cisplatin_End'] <- NA
cor.df['Osi_rep1_Day10', 'CoCl2_End'] <- NA
cor.df['Osi_rep2_Day10', 'CoCl2_End'] <- NA
cor.df['Cisplatin_rep1_Day10', 'Osi_rep1_End'] <- NA
cor.df['Cisplatin_rep2_Day10', 'Osi_rep1_End'] <- NA
cor.df['Cisplatin_rep1_Day10', 'Osi_rep2_End'] <- NA
cor.df['Cisplatin_rep2_Day10', 'Osi_rep2_End'] <- NA
cor.df['Cisplatin_rep1_Day10', 'CoCl2_End'] <- NA
cor.df['Cisplatin_rep2_Day10', 'CoCl2_End'] <- NA
cor.df['Cocl2_rep1_Day10', 'Cisplatin_End'] <- NA
cor.df['CoCl2_rep2_Day10', 'Cisplatin_End'] <- NA
cor.df['Cocl2_rep1_Day10', 'Osi_rep1_End'] <- NA
cor.df['CoCl2_rep2_Day10', 'Osi_rep1_End'] <- NA
cor.df['Cocl2_rep1_Day10', 'Osi_rep2_End'] <- NA
cor.df['CoCl2_rep2_Day10', 'Osi_rep2_End'] <- NA

pdf(paste0(figure_dir, 'Supp_End_comparison_lung_all.pdf'), width = 8, height = 8)
pheatmap(as.matrix(cor.df),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         display_numbers = TRUE,
         main = 'Correlation of Day10 and End',
         cellwidth = 30,
         cellheight = 30,
         gaps_col = c(6), gaps_row = c(6),
         na_col = 'white',
         border_color = 'black'
)
dev.off()

# ggsave(paste0(figure_dir, 'Supp_End_comparison_lung_all.pdf'), width = 5, height = 5)


day10$Osi_day10_mean <- day10[, c('Osi_rep1_Day10', 'Osi_rep2_Day10')] %>% rowMeans()
day10$Cisplatin_day10_mean <- day10[, c('Cisplatin_rep1_Day10', 'Cisplatin_rep2_Day10')] %>% rowMeans()
day10$CoCl2_day10_mean <- day10[, c('Cocl2_rep1_Day10', 'CoCl2_rep2_Day10')] %>% rowMeans()

end$Osi_end_mean <- end[, c('Osi_rep1_End', 'Osi_rep2_End')] %>% rowMeans()

day10.end <- merge(day10, end, by = 'lineage', all = TRUE)
day10.end <- day10.end[, c('lineage', 'Osi_day10_mean', 'Cisplatin_day10_mean', 'CoCl2_day10_mean', 
                           'Osi_end_mean', 'Cisplatin_End', 'CoCl2_End')]

cor.df <- data.frame(matrix(nrow = ncol(day10.end)-1, ncol = ncol(day10.end)-1))
colnames(cor.df) <- colnames(day10.end)[-1]
rownames(cor.df) <- colnames(day10.end)[-1]
for(i in 2:ncol(day10.end)) {
  for (j in 2:ncol(day10.end)) {
    cor.df[i-1, j-1] <- cor(day10.end[, i], day10.end[, j], method = 'pearson', use = 'pairwise.complete.obs')
  }
}
# remove lower triangle
cor.df[lower.tri(cor.df)] <- NA
cor.df['Osi_day10_mean', 'Cisplatin_End'] <- NA
cor.df['Osi_day10_mean', 'CoCl2_End'] <- NA
cor.df['Cisplatin_day10_mean', 'Osi_end_mean'] <- NA
cor.df['Cisplatin_day10_mean', 'CoCl2_End'] <- NA
cor.df['CoCl2_day10_mean', 'Osi_end_mean'] <- NA
cor.df['CoCl2_day10_mean', 'Cisplatin_End'] <- NA

pdf(paste0(figure_dir, 'Supp_End_comparison_lung_all_mean.pdf'), width = 5, height = 5)
pheatmap(as.matrix(cor.df),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         display_numbers = TRUE,
         main = 'Correlation of Day10 and End',
         cellwidth = 30,
         cellheight = 30,
         gaps_col = c(3), gaps_row = c(3),
         breaks = seq(0, 1, by = 0.1),
         na_col = 'white',
         border_color = 'black'
)
dev.off()
