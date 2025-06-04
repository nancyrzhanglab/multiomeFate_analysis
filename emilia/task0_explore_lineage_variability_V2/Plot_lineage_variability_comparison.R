library(tidyverse)
library(ggplot2)
library(ggpubr)

result_dir <- '/Users/emiliac/Library/CloudStorage/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/'
figure_dir <- '/Users/emiliac/Library/CloudStorage/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig3/'
dataset_colors <- c(Day0 = "gray",
                    Day0_shuffle = "gray",
                    Day10_CIS = "#FBD08C",
                    Day10_CIS_shuffle = "#FBD08C",
                    Day10_COCL2 = "#6DC49C",
                    Day10_COCL2_shuffle = "#6DC49C",
                    Day10_DABTRAM = "#9D85BE",
                    Day10_DABTRAM_shuffle = "#9D85BE",
                    Week5_CIS = "#C96D29",
                    Week5_CIS_shuffle = "#C96D29",
                    Week5_COCL2 = "#0F8241",
                    Week5_COCL2_shuffle = "#0F8241",
                    Week5_DABTRAM = "#623594",
                    Week5_DABTRAM_shuffle = "#623594",
                    D2 = "#AFDDFF",
                    D2_shuffle = "#AFDDFF",
                    D4 = "#60B5FF",
                    D4_shuffle = "#60B5FF",
                    D6 = "#0072FF",
                    D6_shuffle = "#0072FF")
# =============================================================================
# Load data
# =============================================================================
hematopo_d2 <- read.csv(paste0(result_dir, 'Weinreb_Larry_Hematopoiesis/lineage_variability_Weinreb_Larry_Hematopoiesis_d2.csv'), row.names = 1)
hematopo_d4 <- read.csv(paste0(result_dir, 'Weinreb_Larry_Hematopoiesis/lineage_variability_Weinreb_Larry_Hematopoiesis_d4.csv'), row.names = 1)
hematopo_d6 <- read.csv(paste0(result_dir, 'Weinreb_Larry_Hematopoiesis/lineage_variability_Weinreb_Larry_Hematopoiesis_d6.csv'), row.names = 1)

hematopo_d2.shuffle <- read.csv(paste0(result_dir, 'Weinreb_Larry_Hematopoiesis/lineage_variability_shuffled_Weinreb_Larry_Hematopoiesis_d2.csv'), row.names = 1)
hematopo_d4.shuffle <- read.csv(paste0(result_dir, 'Weinreb_Larry_Hematopoiesis/lineage_variability_shuffled_Weinreb_Larry_Hematopoiesis_d4.csv'), row.names = 1)
hematopo_d6.shuffle <- read.csv(paste0(result_dir, 'Weinreb_Larry_Hematopoiesis/lineage_variability_shuffled_Weinreb_Larry_Hematopoiesis_d6.csv'), row.names = 1)

day0 <- read.csv(paste0(result_dir, 'day0/lineage_variability_day0_saver_sample.csv'), row.names = 1)
day0.shuffle <- read.csv(paste0(result_dir, 'day0/lineage_variability_shuffledday0_saver_sample.csv'), row.names = 1)

day10_dabtram <- read.csv(paste0(result_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'), row.names = 1)
day10_dabtram.shuffle <- read.csv(paste0(result_dir, 'day10_DABTRAM/lineage_variability_shuffledday10_DABTRAM_saver_sample.csv'), row.names = 1)

day10_cocl2 <- read.csv(paste0(result_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'), row.names = 1)
day10_cocl2.shuffle <- read.csv(paste0(result_dir, 'day10_COCL2/lineage_variability_shuffledday10_COCL2_saver_sample.csv'), row.names = 1)

day10_cis <- read.csv(paste0(result_dir, 'day10_CIS/lineage_variability_day10_CIS_saver_sample.csv'), row.names = 1)
day10_cis.shuffle <- read.csv(paste0(result_dir, 'day10_CIS/lineage_variability_shuffledday10_CIS_saver_sample.csv'), row.names = 1)

week5_dabtram <- read.csv(paste0(result_dir, 'week5_DABTRAM/lineage_variability_week5_DABTRAM_saver_sample.csv'), row.names = 1)
week5_dabtram.shuffle <- read.csv(paste0(result_dir, 'week5_DABTRAM/lineage_variability_shuffledweek5_DABTRAM_saver_sample.csv'), row.names = 1)

week5_cocl2 <- read.csv(paste0(result_dir, 'week5_COCL2/lineage_variability_week5_COCL2_saver_sample.csv'), row.names = 1)
week5_cocl2.shuffle <- read.csv(paste0(result_dir, 'week5_COCL2/lineage_variability_shuffledweek5_COCL2_saver_sample.csv'), row.names = 1)

week5_cis <- read.csv(paste0(result_dir, 'week5_CIS/lineage_variability_week5_CIS_saver_sample.csv'), row.names = 1)
week5_cis.shuffle <- read.csv(paste0(result_dir, 'week5_CIS/lineage_variability_shuffledweek5_CIS_saver_sample.csv'), row.names = 1)

# =============================================================================
# Wrangle
# =============================================================================

hematopo_d2$category <- 'D2'
hematopo_d4$category <- 'D4'
hematopo_d6$category <- 'D6'
hematopo_d2.shuffle$category <- 'D2_shuffle'
hematopo_d4.shuffle$category <- 'D4_shuffle'
hematopo_d6.shuffle$category <- 'D6_shuffle'

day0$category <- 'Day0'
day0.shuffle$category <- 'Day0_shuffle'
day10_dabtram$category <- 'Day10_DABTRAM'
day10_dabtram.shuffle$category <- 'Day10_DABTRAM_shuffle'
day10_cocl2$category <- 'Day10_COCL2'
day10_cocl2.shuffle$category <- 'Day10_COCL2_shuffle'
day10_cis$category <- 'Day10_CIS'
day10_cis.shuffle$category <- 'Day10_CIS_shuffle'
week5_dabtram$category <- 'Week5_DABTRAM'
week5_dabtram.shuffle$category <- 'Week5_DABTRAM_shuffle'
week5_cocl2$category <- 'Week5_COCL2'
week5_cocl2.shuffle$category <- 'Week5_COCL2_shuffle'
week5_cis$category <- 'Week5_CIS'
week5_cis.shuffle$category <- 'Week5_CIS_shuffle'

# Combine all data
compare.df <- rbind(hematopo_d2, hematopo_d2.shuffle, 
                    hematopo_d4, hematopo_d4.shuffle, 
                    hematopo_d6, hematopo_d6.shuffle,
                   day0, day0.shuffle,
                   day10_dabtram, day10_dabtram.shuffle,
                   day10_cocl2, day10_cocl2.shuffle,
                   day10_cis, day10_cis.shuffle,
                   week5_dabtram, week5_dabtram.shuffle,
                   week5_cocl2, week5_cocl2.shuffle,
                   week5_cis, week5_cis.shuffle)

compare.df$category <- factor(compare.df$category, 
                              levels = c('D2', 'D2_shuffle',
                                         'D4', 'D4_shuffle',
                                         'D6', 'D6_shuffle',
                                         'Day0', 'Day0_shuffle',
                                         'Day10_DABTRAM', 'Day10_DABTRAM_shuffle',
                                         'Day10_COCL2', 'Day10_COCL2_shuffle',
                                         'Day10_CIS', 'Day10_CIS_shuffle',
                                         'Week5_DABTRAM', 'Week5_DABTRAM_shuffle',
                                         'Week5_COCL2', 'Week5_COCL2_shuffle',
                                         'Week5_CIS', 'Week5_CIS_shuffle'))

compare.df.hema <- compare.df[compare.df$category %in% c('D2', 'D2_shuffle',
                                                         'D4', 'D4_shuffle',
                                                         'D6', 'D6_shuffle'), ]
compare.df.cancer <- compare.df[!compare.df$category %in% c('D2', 'D2_shuffle',
                                                            'D4', 'D4_shuffle',
                                                            'D6', 'D6_shuffle'), ]

p1 <- ggplot(compare.df.hema, aes(x = category, y = normalized_avg_eud_dist_by_shuffle)) +
  geom_violin(scale = 'width', aes(fill = category)) +
  geom_boxplot(outlier.size = 0.5, outlier.shape = NA, width = 0.3)+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
  scale_fill_manual(values = dataset_colors) +
  ylab('Clonal variability') +
  ylim(min(compare.df$normalized_avg_eud_dist_by_shuffle), 
       max(compare.df$normalized_avg_eud_dist_by_shuffle)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = 'bottom')

p2 <- ggplot(compare.df.cancer, aes(x = category, y = normalized_avg_eud_dist_by_shuffle)) +
  geom_violin(scale = 'width', aes(fill = category)) +
  geom_boxplot(outlier.size = 0.5, outlier.shape = NA, width = 0.3)+
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
  scale_fill_manual(values = dataset_colors) +
  ylab('Clonal variability') +
  ylim(min(compare.df$normalized_avg_eud_dist_by_shuffle), 
       max(compare.df$normalized_avg_eud_dist_by_shuffle)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = 'bottom')

ggarrange(p1, p2, widths = c(7, 14), align = 'hv')

ggsave(paste0(figure_dir, 'Supp_lineage_variability_all_samples.pdf'), 
       width = 8, height = 5)
        