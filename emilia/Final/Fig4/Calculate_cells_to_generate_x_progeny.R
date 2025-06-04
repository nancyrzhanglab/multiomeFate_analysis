rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'

remove_unassigned_cells <- TRUE

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

# =============================================================================
# Wrangle
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# DABTRAM
dabtram_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(dabtram_d10_w5) <- 'DABTRAM_d10_w5'
dabtram_d10_w5$cell_id <- rownames(dabtram_d10_w5)

dabtram_d10_w5 <- merge(dabtram_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
dabtram_d10_w5$progeny_size <- 10**dabtram_d10_w5$DABTRAM_d10_w5

dabtram_d10_w5_total <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["lineage_imputed_count"]])
colnames(dabtram_d10_w5_total) <- 'DABTRAM_d10_w5_total'
dabtram_d10_w5_total$assigned_lineage <- rownames(dabtram_d10_w5_total)


dabtram_d10_w5 <- merge(dabtram_d10_w5, dabtram_d10_w5_total, by = 'assigned_lineage')
dabtram_d10_w5$proportion <- dabtram_d10_w5$progeny_size / dabtram_d10_w5$DABTRAM_d10_w5_total * 100


cum.prop.dabtram <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cum.prop.dabtram) <- c('lineage', 'n_cells_greater_50', 'n_cells_total')
lineages <- unique(dabtram_d10_w5$assigned_lineage)
thres = 0.5
for(lin in lineages){
  df <- dabtram_d10_w5[dabtram_d10_w5$assigned_lineage == lin, ]
  
  progeny_sizes <- df$progeny_size
  
  # get progeny
  progeny_sizes=sort(progeny_sizes, decreasing=TRUE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  
  # find the index of the first value that is greater than or equal to 0.5
  index = which(sizeprop >= thres)[1]
  
  cum.prop.dabtram <- rbind(cum.prop.dabtram, data.frame(lineage = lin, n_cells_greater_50 = index, n_cells_total = length(progeny_sizes)))
  
}

cum.prop.dabtram$fraction <- cum.prop.dabtram$n_cells_greater_50 / cum.prop.dabtram$n_cells_total

cum.prop.dabtram.large <- cum.prop.dabtram[cum.prop.dabtram$n_cells_total > 10, ]
mean(cum.prop.dabtram.large$fraction)

# COCL2
cocl2_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_COCL2_d10_w5"]][["cell_imputed_score"]])
colnames(cocl2_d10_w5) <- 'COCL2_d10_w5'
cocl2_d10_w5$cell_id <- rownames(cocl2_d10_w5)
cocl2_d10_w5 <- merge(cocl2_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cocl2_d10_w5$progeny_size <- 10**cocl2_d10_w5$COCL2_d10_w5
cocl2_d10_w5_total <- as.data.frame(all_data_fatepotential[["fatepotential_COCL2_d10_w5"]][["lineage_imputed_count"]])
colnames(cocl2_d10_w5_total) <- 'COCL2_d10_w5_total'
cocl2_d10_w5_total$assigned_lineage <- rownames(cocl2_d10_w5_total)
cocl2_d10_w5 <- merge(cocl2_d10_w5, cocl2_d10_w5_total, by = 'assigned_lineage')

cocl2_d10_w5$proportion <- cocl2_d10_w5$progeny_size / cocl2_d10_w5$COCL2_d10_w5_total * 100

cum.prop.cocl2 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cum.prop.cocl2) <- c('lineage', 'n_cells_greater_50', 'n_cells_total')
lineages <- unique(cocl2_d10_w5$assigned_lineage)
thres = 0.5
for(lin in lineages){
  df <- cocl2_d10_w5[cocl2_d10_w5$assigned_lineage == lin, ]
  
  progeny_sizes <- df$progeny_size
  
  # get progeny
  progeny_sizes=sort(progeny_sizes, decreasing=TRUE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  
  # find the index of the first value that is greater than or equal to 0.5
  index = which(sizeprop >= thres)[1]
  
  cum.prop.cocl2 <- rbind(cum.prop.cocl2, data.frame(lineage = lin, n_cells_greater_50 = index, n_cells_total = length(progeny_sizes)))
  
  
  
}

cum.prop.cocl2$fraction <- cum.prop.cocl2$n_cells_greater_50 / cum.prop.cocl2$n_cells_total

cum.prop.cocl2.large <- cum.prop.cocl2[cum.prop.cocl2$n_cells_total > 10, ]
mean(cum.prop.cocl2.large$fraction)

# CIS
cis_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_CIS_d10_w5"]][["cell_imputed_score"]])
colnames(cis_d10_w5) <- 'CIS_d10_w5'
cis_d10_w5$cell_id <- rownames(cis_d10_w5)
cis_d10_w5 <- merge(cis_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cis_d10_w5$progeny_size <- 10**cis_d10_w5$CIS_d10_w5
cis_d10_w5_total <- as.data.frame(all_data_fatepotential[["fatepotential_CIS_d10_w5"]][["lineage_imputed_count"]])
colnames(cis_d10_w5_total) <- 'CIS_d10_w5_total'
cis_d10_w5_total$assigned_lineage <- rownames(cis_d10_w5_total)
cis_d10_w5 <- merge(cis_d10_w5, cis_d10_w5_total, by = 'assigned_lineage')
cis_d10_w5$proportion <- cis_d10_w5$progeny_size / cis_d10_w5$CIS_d10_w5_total * 100


cum.prop.cis <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cum.prop.cis) <- c('lineage', 'n_cells_greater_50', 'n_cells_total')
lineages <- unique(cis_d10_w5$assigned_lineage)
thres = 0.5
for(lin in lineages){
  df <- cis_d10_w5[cis_d10_w5$assigned_lineage == lin, ]
  
  progeny_sizes <- df$progeny_size
  
  # get progeny
  progeny_sizes=sort(progeny_sizes, decreasing=TRUE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  
  # find the index of the first value that is greater than or equal to 0.5
  index = which(sizeprop >= thres)[1]
  
  cum.prop.cis <- rbind(cum.prop.cis, data.frame(lineage = lin, n_cells_greater_50 = index, n_cells_total = length(progeny_sizes)))
  
}
cum.prop.cis$fraction <- cum.prop.cis$n_cells_greater_50 / cum.prop.cis$n_cells_total

cum.prop.cis.large <- cum.prop.cis[cum.prop.cis$n_cells_total > 10, ]
mean(cum.prop.cis.large$fraction)

# plot

cum.prop.cis.large$dataset <- 'week5_CIS'
cum.prop.dabtram.large$dataset <- 'week5_DABTRAM'
cum.prop.cocl2.large$dataset <- 'week5_COCL2'

df <- rbind(cum.prop.cis.large, cum.prop.cocl2.large, cum.prop.dabtram.large)
df$dataset <- factor(df$dataset, levels = c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS'))

ggplot(df, aes(x = dataset, y = fraction)) +
  geom_violin(aes(fill = dataset), scale = 'width', alpha = 0.6) +
  geom_jitter(size = 1, width = 0.15) +
  geom_boxplot(width = 0.3, alpha = 0.8) +
  scale_fill_manual(values = dataset_colors) +
  ylab('Fraction of cells with >50% progeny size') +
  ylim(0, 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')
ggsave(paste0(figure_dir, 'Supp_Fraction_Progeny_From_X_cells.pdf'), width = 3, height = 4.5)
