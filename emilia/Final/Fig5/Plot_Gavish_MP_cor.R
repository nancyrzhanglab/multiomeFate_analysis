rm(list = ls())

library(tidyverse)
library(ggplot2)
library(pheatmap)

remove_unassigned_cells <- TRUE
# =============================================================================
# reading data
# =============================================================================
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

gavish.mp.cor.d0 <- read.csv(paste0(results_dir, 'GavishMP_UCell_cor_d0_d10.csv'))
gavish.mp.cor.d10 <- read.csv(paste0(results_dir, 'GavishMP_UCell_cor_d10_w5.csv'))

mps <- c('Cell Cycle - G1/S.rho', 'Cell Cycle - G2/M.rho', 'Cell Cycle HMG-rich.rho', 'Chromatin.rho',
         'EMT-I.rho', 'EMT-II.rho', 'EMT-III.rho', 'EMT-IV.rho', 'Epithelial Senescence.rho', 
         'Hypoxia.rho', 'Interferon/MHC-II (I).rho', 'Interferon/MHC-II (II).rho', 'MYC.rho', 
         'Proteasomal degradation.rho', 'Protein maturation.rho', 'Respiration.rho', 'Secreted I.rho',
         'Secreted II.rho', 'Skin-pigmentation.rho', 'Stress (in vitro).rho', 'Stress.rho',
         'Translation initiation.rho', 'Unfolded protein response.rho')

gavish.mp.cor.d0 <- gavish.mp.cor.d0[gavish.mp.cor.d0$MetaProgram %in% mps,]
gavish.mp.cor.d10 <- gavish.mp.cor.d10[gavish.mp.cor.d10$MetaProgram %in% mps,]

# =============================================================================
# Plot
# =============================================================================
to_plot <- merge(gavish.mp.cor.d0, gavish.mp.cor.d10, by = c('MetaProgram'))
rownames(to_plot) <- to_plot$MetaProgram
to_plot <- to_plot[,-c(1)]

pheatmap(to_plot, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T, 
         border_color = 'black', cellwidth = 20, cellheight = 20,
         gaps_col = c(3, 6),
         fontsize = 8, fontsize_row = 8, fontsize_col = 8)
         

