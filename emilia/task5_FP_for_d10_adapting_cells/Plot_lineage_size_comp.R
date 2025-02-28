rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)
library(gridExtra)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))


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
# calculate day0 and week5 lineage size
# =============================================================================
metadat <- all_data@meta.data

lineage.size <- table(metadat$dataset, metadat$assigned_lineage)
lineage.size <- as.data.frame(lineage.size)
colnames(lineage.size) <- c('dataset', 'assigned_lineage', 'n_cells')
lineage.size$n_cells <- log10(lineage.size$n_cells + 1)
lineage.size <- spread(lineage.size, key = dataset, value = n_cells)

lineage.size.dabtram <- lineage.size[, c('assigned_lineage', 'day0', 'week5_DABTRAM')]
lineage.size.cocl2 <- lineage.size[, c('assigned_lineage', 'day0', 'week5_COCL2')]
lineage.size.cis <- lineage.size[, c('assigned_lineage', 'day0', 'week5_CIS')]

# =============================================================================
# plot
# =============================================================================

# DABTRAM
n_lineage_dabtram_both <- lineage.size.dabtram %>% 
  filter(day0 != 0 & week5_DABTRAM != 0) %>%
  nrow()

n_lineage_dabtram_either <- lineage.size.dabtram %>% 
  filter(day0 != 0 | week5_DABTRAM != 0) %>%
  nrow()

lineage.size.dabtram <- lineage.size.dabtram %>% 
  filter(day0 != 0 | week5_DABTRAM != 0)

p1 <- ggplot(lineage.size.dabtram, aes(x = day0, y = week5_DABTRAM)) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day 0 lineage size', 
       y = 'Week 5 lineage size', 
       title = paste0('DABTRAM: \n',
                      '# lin in both = ', n_lineage_dabtram_both,
                      ' out of # lin in either = ', n_lineage_dabtram_either)) +
  theme_bw()
p1

# COCL2
n_lineage_cocl2_both <- lineage.size.cocl2 %>% 
  filter(day0 != 0 & week5_COCL2 != 0) %>%
  nrow()

n_lineage_cocl2_either <- lineage.size.cocl2 %>% 
  filter(day0 != 0 | week5_COCL2 != 0) %>%
  nrow()

lineage.size.cocl2 <- lineage.size.cocl2 %>% 
  filter(day0 != 0 | week5_COCL2 != 0)
p2 <- ggplot(lineage.size.cocl2, aes(x = day0, y = week5_COCL2)) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day 0 lineage size', 
       y = 'Week 5 lineage size', 
       title = paste0('COCL2: \n',
                      '# lin in both = ', n_lineage_cocl2_both,
                      ' out of # lin in either = ', n_lineage_cocl2_either)) +
  theme_bw()
p2

# CIS
n_lineage_cis_both <- lineage.size.cis %>% 
  filter(day0 != 0 & week5_CIS != 0) %>%
  nrow()

n_lineage_cis_either <- lineage.size.cis %>% 
  filter(day0 != 0 | week5_CIS != 0) %>%
  nrow()

lineage.size.cis <- lineage.size.cis %>% 
  filter(day0 != 0 | week5_CIS != 0)
p3 <- ggplot(lineage.size.cis, aes(x = day0, y = week5_CIS)) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day 0 lineage size', 
       y = 'Week 5 lineage size', 
       title = paste0('CIS: \n',
                      '# lin in both = ', n_lineage_cis_both,
                      ' out of # lin in either = ', n_lineage_cis_either)) +
  theme_bw()
p3

grid.arrange(p1, p2, p3, ncol = 3)


# =============================================================================
# read final fits
# =============================================================================
final_fit_dabtram <- readRDS('~/Downloads/final_fit_d0_w5_DABTRAM.rds')
lineage_imputed_count.dabtram <- as.data.frame(final_fit_dabtram$lineage_imputed_count)
colnames(lineage_imputed_count.dabtram) <- 'd0_predicted_w5_lineage_size'
lineage_imputed_count.dabtram$assigned_lineage <- rownames(lineage_imputed_count.dabtram)

final_fit_cocl2 <- readRDS('~/Downloads/final_fit_d0_w5_COCL2_scaled.rds')
lineage_imputed_count.cocl2 <- as.data.frame(final_fit_cocl2$lineage_imputed_count)
colnames(lineage_imputed_count.cocl2) <- 'd0_predicted_w5_lineage_size'
lineage_imputed_count.cocl2$assigned_lineage <- rownames(lineage_imputed_count.cocl2)

final_fit_cis <- readRDS('~/Downloads/final_fit_d0_w5_CIS.rds')
lineage_imputed_count.cis <- as.data.frame(final_fit_cis$lineage_imputed_count)
colnames(lineage_imputed_count.cis) <- 'd0_predicted_w5_lineage_size'
lineage_imputed_count.cis$assigned_lineage <- rownames(lineage_imputed_count.cis)

lineage.size.dabtram <- merge(lineage.size.dabtram, lineage_imputed_count.dabtram, by = 'assigned_lineage', all = T)
lineage.size.cocl2 <- merge(lineage.size.cocl2, lineage_imputed_count.cocl2, by = 'assigned_lineage', all = T)
lineage.size.cis <- merge(lineage.size.cis, lineage_imputed_count.cis, by = 'assigned_lineage', all = T)

# =============================================================================
# plot
# =============================================================================
lineage.size.dabtram.plot <- lineage.size.dabtram %>% 
  filter(day0 != 0 & week5_DABTRAM != 0)

p4 <- ggplot(lineage.size.dabtram, aes(x = week5_DABTRAM, y = d0_predicted_w5_lineage_size)) +
  geom_point() +
  stat_cor() +
  labs(x = 'Day 0 lineage size', 
       y = 'Week 5 lineage size', 
       title = paste0('DABTRAM: \n') +
  theme_bw()
p4
