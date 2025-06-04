rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)
library(hdrcde)
library(ggdensity)
library(RColorBrewer)

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
day10$all <- rowSums(day10)
day10 <- day10[day10$all > 0, ]
# day10[day10 < thres] <- 0
day10 <- log10(day10+1)
day10 <- subset(day10, select = -c(all))



# ggpairs(day10, columns = 1:6,
#         lower = list(continuous=wrap("points"))) +
#   scale_color_manual(values = c('gray', 'red')) +
#   theme_bw()
# 
# 
# day10$Osi_day10_mean <- day10[, c('Sample_11_Osi_day10', 'Sample_12_Osi_day10')] %>% rowMeans()
# day10$Cisplatin_day10_mean <- day10[, c('Sample_15_Cisplatin_day10', 'Sample_16_Cisplatin_day10')] %>% rowMeans()
# day10$CoCl2_day10_mean <- day10[, c('Sample_5_CoCl2_day10', 'Sample_6_CoCl2_day10')] %>% rowMeans()
# 
# ggpairs(day10[, c('Osi_day10_mean', 'Cisplatin_day10_mean', 'CoCl2_day10_mean')],
#         upper = list(continuous = wrap("cor", method = "spearman")),
#         lower = list(continuous= wrap("points",size = 0.5))) +
#   theme_bw()

# plot a contour plot
day10.density <- day10[day10$CoCl2_day10_mean > 0, ]
day10.density <- day10.density[day10.density$Cisplatin_day10_mean > 0, ]
p1 <- ggplot(day10, aes(x = CoCl2_day10_mean, y = Cisplatin_day10_mean)) +
  geom_point(color = 'pink') +
  # geom_hdr(data = day10.density, aes(x = CoCl2_day10_mean, y = Cisplatin_day10_mean, fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45, size = 6) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  xlab('CoCl2 day10 (mean)') +
  ylab('Cisplatin day10 (mean)') +
  theme_bw()

day10.density <- day10[day10$CoCl2_day10_mean > 0, ]
day10.density <- day10.density[day10.density$Osi_day10_mean > 0, ]
p2 <- ggplot(day10, aes(x = CoCl2_day10_mean, y = Osi_day10_mean)) +
  geom_point(color = 'pink') +
  # geom_hdr(data = day10.density, aes(x = CoCl2_day10_mean, y = Osi_day10_mean, fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45, size = 6) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  xlab('CoCl2 day10 (mean)') +
  ylab('Osi day10 (mean)') +
  theme_bw()

day10.density <- day10[day10$Cisplatin_day10_mean > 0, ]
day10.density <- day10.density[day10.density$Osi_day10_mean > 0, ]
p3 <- ggplot(day10, aes(x = Cisplatin_day10_mean, y = Osi_day10_mean)) +
  geom_point(color = 'pink') +
  # geom_hdr(data = day10.density, aes(x = Cisplatin_day10_mean, y = Osi_day10_mean, fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45, size = 6) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  xlab('Cisplatin day10 (mean)') +
  ylab('Osi day10 (mean)') +
  theme_bw()

p4 <- grid.arrange(p1, p2, p3, nrow = 1)
ggsave(paste0(figure_dir, 'Supp_Day10_comparison_lung.png'), p4, width = 9, height = 3, dpi = 600)

# ==============================================================================
# End
# ==============================================================================
end.spike <- data.spike_in[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                               'Sample_38_Cisplatin_END')]
min <- min(end.spike)
thres <- min / 5
# end.spike <- log10(end.spike + 1)

end <- data.without.spike[, c('Sample_40_Osi_END', 'Sample_41_Osi_END', 'Sample_43_CoCl2_END',
                              'Sample_38_Cisplatin_END')]
end$all <- rowSums(end)
end <- end[end$all > 0, ]
end[end < thres] <- 0
end <- log10(end+1)

end <- subset(end, select = -c(all))



end$is.spike <- 'non-spike'
end.spike$is.spike <- 'spike'

end <- rbind(end, end.spike)

end$Osi_end_mean <- end[, c('Sample_40_Osi_END', 'Sample_41_Osi_END')] %>% rowMeans()

ggpairs(end[, c('Osi_end_mean', 'Sample_43_CoCl2_END', 'Sample_38_Cisplatin_END')],
        upper = list(continuous = wrap("cor", method = "pearson")),
        lower = list(continuous=wrap("points",size = 0.5))) +
  theme_bw()


p4 <- ggplot(end, aes(x = Sample_43_CoCl2_END, y = Sample_38_Cisplatin_END)) +
  geom_point(color = '#A53860') +
  # geom_hdr(data = day10.density, aes(x = CoCl2_day10_mean, y = Cisplatin_day10_mean, fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 6, size = 6) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  xlab('CoCl2 end (mean)') +
  ylab('Cisplatin end (mean)') +
  theme_bw()

p5 <- ggplot(end, aes(x = Osi_end_mean, y = Sample_38_Cisplatin_END)) +
  geom_point(color = '#A53860') +
  # geom_hdr(data = day10.density, aes(x = CoCl2_day10_mean, y = Cisplatin_day10_mean, fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 6, size = 6) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  xlab('Osi end (mean)') +
  ylab('Cisplatin end (mean)') +
  theme_bw()

p6 <- ggplot(end, aes(x = Osi_end_mean, y = Sample_43_CoCl2_END)) +
  geom_point(color = '#A53860') +
  # geom_hdr(data = day10.density, aes(x = CoCl2_day10_mean, y = Cisplatin_day10_mean, fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 4, size = 6) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  xlab('Osi end (mean)') +
  ylab('CoCl2 end (mean)') +
  theme_bw()

p7 <- grid.arrange(p4, p5, p6, nrow = 1)
ggsave(paste0(figure_dir, 'Supp_End_comparison_lung.png'), p7, width = 9, height = 3, dpi = 600)

# ==============================================================================
# Day10 vs End
# ==============================================================================

day10$lineage <- rownames(day10)
end$lineage <- rownames(end)
day10.end <- merge(day10, end, by = 'lineage', all = TRUE)

day10.end.CoCl2 <- day10.end[, c('CoCl2_day10_mean', 'Sample_43_CoCl2_END')] %>% drop_na()
day10.end.CoCl2$all <- rowSums(day10.end.CoCl2)
day10.end.CoCl2 <- day10.end.CoCl2[day10.end.CoCl2$all > 0, ]

p8 <- ggplot(day10.end.CoCl2, aes(x = CoCl2_day10_mean, y = Sample_43_CoCl2_END)) +
  geom_point(color = 'black') +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 6, size = 6) +
  xlab('CoCl2 day10 (mean)') +
  ylab('CoCl2 end (mean)') +
  theme_bw()

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
ggsave(paste0(figure_dir, 'Supp_Day10_End_comparison_lung.png'), p11, width = 9, height = 3, dpi = 600)
