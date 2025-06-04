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

data.spike_in$category <- paste0('spike_in_', data.spike_in$Num_Of_Cells)

data.spike_in <- data.spike_in[, c('category', colnames(data.without.spike))]
data.without.spike$category <- 'non_spike_in'

data.all <- rbind(data.spike_in, data.without.spike)

# calculate library size per sample
library.size <- colSums(data.all[, 2: ncol(data.all)])
names(library.size) <- colnames(data.all)[2: ncol(data.all)]

# divide the data by library size
data.all.freq <- sapply(seq(2, ncol(data.all)), function(x) data.all[x] / library.size[x-1])
# combine the data
data.all.freq <- as.data.frame(data.all.freq)
rownames(data.all.freq) <- rownames(data.all)
data.all.freq$category <- data.all$category


data.spike.freq <- data.all.freq[data.all.freq$category != 'non_spike_in', ]
data.nonspike.freq <- data.all.freq[data.all.freq$category == 'non_spike_in', ]
ggplot(data.all.freq, aes(x = log2(Sample_11_Osi_day10 + 10e-6), y = log2(Sample_12_Osi_day10+ 10e-6))) +
  geom_point() +
  geom_point(data = data.spike.freq, aes(x = log2(Sample_11_Osi_day10 + 10e-6), y = log2(Sample_12_Osi_day10+ 10e-6)), color = 'red') +
  stat_cor(method = 'spearman') +
  theme_classic()

ggplot(data.all.freq, aes(x = log2(Sample_11_Osi_day10 + 10e-6), y = log2(Sample_40_Osi_END+ 10e-6))) +
  geom_point() +
  stat_cor(method = 'spearman') +
  theme_classic()

data.nonspike.freq.toplot <- data.nonspike.freq[, 1: ncol(data.nonspike.freq)-1]
data.nonspike.freq.toplot <- log2(data.nonspike.freq.toplot + 10e-6)

png('~/Downloads/freq_compare.png', width = 3000, height = 3000, res = 200)
ggpairs(data.nonspike.freq.toplot)
dev.off()

data.nonspike.freq.toplot <- data.nonspike.freq.toplot[, c("Sample_11_Osi_day10", "Sample_12_Osi_day10", 
                                                           "Sample_5_CoCl2_day10", "Sample_6_CoCl2_day10",
                                                           "Sample_15_Cisplatin_day10", "Sample_16_Cisplatin_day10",
                                                           "Sample_17_DMSO_day10", "Sample_18_DMSO_day10",
                                                           "Sample_40_Osi_END", "Sample_41_Osi_END",
                                                           "Sample_43_CoCl2_END", "Sample_38_Cisplatin_END")]
cor.df <- data.frame(matrix(nrow = ncol(data.nonspike.freq.toplot), ncol = ncol(data.nonspike.freq.toplot)))
colnames(cor.df) <- colnames(data.nonspike.freq.toplot)
rownames(cor.df) <- colnames(data.nonspike.freq.toplot)
for(i in 1:ncol(data.nonspike.freq.toplot)) {
  for (j in 1:ncol(data.nonspike.freq.toplot)) {
    cor.df[i, j] <- cor(data.nonspike.freq.toplot[, i], 
                        data.nonspike.freq.toplot[, j], method = 'pearson', use = 'pairwise.complete.obs')
  }
}

# remove lower triangle
cor.df[lower.tri(cor.df)] <- NA
cor.df['Sample_11_Osi_day10', 'Sample_38_Cisplatin_END'] <- NA
cor.df['Sample_11_Osi_day10', 'Sample_43_CoCl2_END'] <- NA
cor.df['Sample_12_Osi_day10', 'Sample_38_Cisplatin_END'] <- NA
cor.df['Sample_12_Osi_day10', 'Sample_43_CoCl2_END'] <- NA
cor.df['Sample_15_Cisplatin_day10', 'Sample_40_Osi_END'] <- NA
cor.df['Sample_15_Cisplatin_day10', 'Sample_41_Osi_END'] <- NA
cor.df['Sample_16_Cisplatin_day10', 'Sample_40_Osi_END'] <- NA
cor.df['Sample_16_Cisplatin_day10', 'Sample_41_Osi_END'] <- NA
cor.df['Sample_15_Cisplatin_day10', 'Sample_43_CoCl2_END'] <- NA
cor.df['Sample_16_Cisplatin_day10', 'Sample_43_CoCl2_END'] <- NA
cor.df['Sample_5_CoCl2_day10', 'Sample_40_Osi_END'] <- NA
cor.df['Sample_5_CoCl2_day10', 'Sample_41_Osi_END'] <- NA
cor.df['Sample_5_CoCl2_day10', 'Sample_38_Cisplatin_END'] <- NA
cor.df['Sample_6_CoCl2_day10', 'Sample_40_Osi_END'] <- NA
cor.df['Sample_6_CoCl2_day10', 'Sample_41_Osi_END'] <- NA
cor.df['Sample_6_CoCl2_day10', 'Sample_38_Cisplatin_END'] <- NA

pdf(paste0(figure_dir, 'Supp_gDNA_barcode_freq_lung_correlation.pdf'), width = 10, height = 10)
pheatmap(as.matrix(cor.df),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         display_numbers = TRUE,
         cellwidth = 30,
         cellheight = 30,
         gaps_col = c(8), gaps_row = c(8),
         breaks = seq(0, 1, by = 0.1),
         na_col = 'white',
         border_color = 'black'
)
dev.off()
