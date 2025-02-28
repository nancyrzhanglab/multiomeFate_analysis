rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

treatment <- 'DABTRAM'

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
# Wrangle data
# =============================================================================

# Extract day10 future size
all_data_day10 <- subset(all_data, dataset == paste0('day10_', treatment))
metadat.day10 <- all_data_day10@meta.data
metadat.day10$cell_id <- rownames(metadat.day10)

fp_name <- paste0('fatepotential_', treatment, '_d10_w5')

fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['lineage_imputed_count']])
colnames(fp_d10_w5_treatment) <- c('lineage_imputed_count')
fp_d10_w5_treatment$assigned_lineage <- rownames(fp_d10_w5_treatment)

# Calculate day10 current

lineage_size.day10 <- metadat.day10 %>% 
  group_by(assigned_lineage) %>%
  summarise(lineage_size = n())

lineage_size.day10 <- merge(lineage_size.day10, fp_d10_w5_treatment, by = 'assigned_lineage', all.x = TRUE)

ggplot(lineage_size.day10, aes(x = log10(lineage_size + 1), y = log10(lineage_imputed_count+1))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  stat_cor() +
  labs(x = 'Current lineage size', y = 'Future lineage size') +
  theme_bw()
