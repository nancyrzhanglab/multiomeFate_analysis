rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggplot2)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE


# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data@misc <- all_data_fatepotential

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# =============================================================================
# Calculate current size
# =============================================================================

metadat <- all_data@meta.data

lineage.size <- metadat %>% 
  group_by(dataset, assigned_lineage) %>% 
  summarise(n_cells = n()) %>% 
  arrange(desc(n_cells))

lineage.size.wide <- spread(lineage.size, key = dataset, value = n_cells)

# =============================================================================
# Get predicted lineage size
# =============================================================================

fp_w5_DABTRAM <- as.data.frame(all_data@misc[['fatepotential_DABTRAM_d10_w5']][['lineage_imputed_count']])
colnames(fp_w5_DABTRAM) <- c('lineage_imputed_count.DABTRAM')

fp_w5_COCL2 <- as.data.frame(all_data@misc[['fatepotential_COCL2_d10_w5']][['lineage_imputed_count']])
colnames(fp_w5_COCL2) <- c('lineage_imputed_count.COCL2')

fp_w5_CIS <- as.data.frame(all_data@misc[['fatepotential_CIS_d10_w5']][['lineage_imputed_count']])
colnames(fp_w5_CIS) <- c('lineage_imputed_count.CIS')

fp_w5_DABTRAM$assigned_lineage <- rownames(fp_w5_DABTRAM)
fp_w5_COCL2$assigned_lineage <- rownames(fp_w5_COCL2)
fp_w5_CIS$assigned_lineage <- rownames(fp_w5_CIS)

fp_w5 <- merge(fp_w5_DABTRAM, fp_w5_COCL2, by = 'assigned_lineage', all = TRUE)
fp_w5 <- merge(fp_w5, fp_w5_CIS, by = 'assigned_lineage', all = TRUE)


# =============================================================================
# Combine
# =============================================================================

comp_df <- merge(lineage.size.wide, fp_w5, by = 'assigned_lineage', all = TRUE)

comp_df.DABTRAM <- comp_df[, c('assigned_lineage', 'day10_DABTRAM', 'lineage_imputed_count.DABTRAM')]
comp_df.COCL2 <- comp_df[, c('assigned_lineage', 'day10_COCL2', 'lineage_imputed_count.COCL2')]
comp_df.CIS <- comp_df[, c('assigned_lineage', 'day10_CIS', 'lineage_imputed_count.CIS')]

comp_df.DABTRAM.actual <- comp_df[, c('assigned_lineage', 'day10_DABTRAM', 'week5_DABTRAM')]
comp_df.COCL2.actual <- comp_df[, c('assigned_lineage', 'day10_COCL2', 'week5_COCL2')]
comp_df.CIS.actual <- comp_df[, c('assigned_lineage', 'day10_CIS', 'week5_CIS')]

# =============================================================================
# Plot
# =============================================================================

p1 <- ggplot(comp_df.DABTRAM, aes(x = log10(day10_DABTRAM + 1), y = log10(lineage_imputed_count.DABTRAM+1))) +
  geom_point() +
  stat_cor() +
  labs(x = 'Observed lineage size', y = 'Predicted lineage size') +
  ggtitle('DABTRAM') +
  theme_bw()

p2 <- ggplot(comp_df.COCL2, aes(x = log10(day10_COCL2 + 1), y = log10(lineage_imputed_count.COCL2+1))) +
  geom_point() +
  stat_cor() +
  labs(x = 'Observed lineage size', y = 'Predicted lineage size') +
  ggtitle('COCL2') +
  theme_bw()

p3 <- ggplot(comp_df.CIS, aes(x = log10(day10_CIS + 1), y = log10(lineage_imputed_count.CIS+1))) +
  geom_point() +
  stat_cor() +
  labs(x = 'Observed lineage size', y = 'Predicted lineage size') +
  ggtitle('CIS') +
  theme_bw()


p4 <- ggplot(comp_df.DABTRAM.actual, aes(x = log10(day10_DABTRAM + 1), y = log10(week5_DABTRAM+1))) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day10 lineage size', y = 'Week5 lineage size') +
  ggtitle('DABTRAM') +
  theme_bw()

p5 <- ggplot(comp_df.COCL2.actual, aes(x = log10(day10_COCL2 + 1), y = log10(week5_COCL2+1))) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day10 lineage size', y = 'Week5 lineage size') +
  ggtitle('COCL2') +
  theme_bw()

p6 <- ggplot(comp_df.CIS.actual, aes(x = log10(day10_CIS + 1), y = log10(week5_CIS+1))) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day10 lineage size', y = 'Week5 lineage size') +
  ggtitle('CIS') +
  theme_bw()

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
