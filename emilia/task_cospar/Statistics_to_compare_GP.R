library(ggplot2)
library(tidyverse)
library(dplyr)

# ==============================================================================
# Read data
# ==============================================================================

data_dir <- "~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/"
treatment_day10 <- 'day10_DABTRAM'
load(paste0(data_dir, treatment_day10, "/Writeup6n_DABTRAM_day10_lineage-imputation_stepdown_concise-postprocessed.RData"))
metadat_day10 <- read.csv(paste0(data_dir, treatment_day10, '/', treatment_day10, '_meta.csv'), row.names = 1)

cell_imputed_count <- as.data.frame(cell_imputed_count[!is.na(cell_imputed_count)])
colnames(cell_imputed_count) <- c('Growth_Potential')
cell_imputed_count$cell_id <- row.names(cell_imputed_count)

lineage_imputed_count <- as.data.frame(lineage_imputed_count)
lineage_imputed_count$assigned_lineage <- row.names(lineage_imputed_count)

cospar <- read.csv('~/Downloads/Day10_DABTRAM_GP_CoSpar.csv')

# ==============================================================================
# Plot
# ==============================================================================
to_plot <- merge(cospar, lineage_imputed_count, by = 'assigned_lineage')

ggplot(to_plot, aes(x = lineage_imputed_count, y = `Growth.potential`)) +
  geom_point() +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_smooth(method = 'lm') +
  theme_bw() +
  xlab('Growth potential') +
  ylab('CoSpar growth potential metric')



