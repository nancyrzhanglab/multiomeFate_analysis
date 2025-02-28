library(multiomeFate)
library(Seurat)
library(dplyr)

all_data = multiomeFate:::data_loader(which_files = c("fatepotential"))

# ==============================================================================
# Calculate lineage size on week5
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_barcode <- rownames(metadat)
metadat.wk5 <- metadat[metadat$dataset %in% c('week5_CIS', 'week5_COCL2', 'week5_DABTRAM')]

lineage_size_wk5 <- metadat.wk5 %>% 
  group_by(dataset, assigned_lineage) %>% 
  summarise(n = n())

# ==============================================================================
# Extracting fate potential scores d0 -> d10 DABTRAM
# ==============================================================================

cur_time = 'd0'
fut_time = 'd10'
treatment = 'DABTRAM'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp_d0_d10 <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
colnames(fp_d0_d10) <- c(fp_name)
fp_d0_d10$cell_barcode <- rownames(fp_d0_d10)

fp_d0_d10 <- merge(fp_d0_d10, metadat[, c('cell_barcode', 'assigned_lineage')], by = 'cell_barcode')
fp_d0_d10 <- fp_d0_d10[order(fp_d0_d10[[fp_name]], decreasing = TRUE), ]
rownames(fp_d0_d10) <- seq(1, nrow(fp_d0_d10))
fp_d0_d10_top <- fp_d0_d10[seq(1, as.integer(nrow(fp_d0_d10)) * 0.05), ]

wk5_dabtram <- lineage_size_wk5[lineage_size_wk5$dataset == 'week5_DABTRAM', ]
wk5_dabtram <- wk5_dabtram[order(wk5_dabtram$n, decreasing = TRUE), ]
rownames(wk5_dabtram) <- seq(1, nrow(wk5_dabtram))
wk5_dabtram$is_top_d0 <- ifelse(wk5_dabtram$assigned_lineage %in% fp_d0_d10_top$assigned_lineage, 'Top_d0', 'Not_top_d0')

n_top_d0_in_wk5 <- wk5_dabtram %>% 
  group_by(is_top_d0) %>% 
  summarise(n_top_d0_by_lineage = n_distinct(assigned_lineage))

metadat.wk5 <- metadat.wk5[metadat.wk5$dataset == 'week5_DABTRAM', ]
metadat.wk5$is_top_d0 <- ifelse(metadat.wk5$assigned_lineage %in% fp_d0_d10_top$assigned_lineage, 'Top_d0', 'Not_top_d0')
n_cell_in_wk5_from_lin_top_d0 <- metadat.wk5 %>% 
  group_by(is_top_d0) %>% 
  summarise(n_top_d0_by_cell_num = n())


