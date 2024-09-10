rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)

# =============================================================================
# Load data
# =============================================================================

TIME = 'day10' # 'day10', or 'day0'
TREATMENT = 'DABTRAM' # 'COCL2', 'DABTRAM', or 'CIS'
MODALITY = 'saver' # or 'peakvi' 
SAMPLE_NAME = paste0(TIME, '_', TREATMENT)
output_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task1_ANOVA_lineage_specific_features/"

load('~/nzhanglab/project/Multiome_fate/out/kevin/Writeup10a/Writeup10a_data_empty.RData')

load('~/nzhanglab/project/Multiome_fate/out/kevin/Writeup10a/Writeup10a_data_saver.RData')
all_data[["Saver"]] <- all_data_saver
all_data[["Saver.pca"]] <- all_data_saver_pca
all_data[["Saver.umap"]] <- all_data_saver_umap

all_data$keep <- !is.na(all_data$assigned_lineage)
if(any(!all_data$keep)){
  print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
  all_data <- subset(all_data, keep == TRUE)
}

# =============================================================================
# Wrangle data
# =============================================================================

all_data_use = subset(all_data, dataset == SAMPLE_NAME)

Seurat::DefaultAssay(all_data_use) <- "Saver"
all_data_use <- NormalizeData(all_data_use)
all_data_use <- FindVariableFeatures(all_data_use)
all_data_use <- ScaleData(all_data_use)

data <- all_data_use@assays[["Saver"]]@scale.data
metadat <- all_data_use@meta.data
metadat$cell_ids <- row.names(metadat)

# ==============================================================================
# Match lineage information
# ==============================================================================
lineage_info <- metadat[, c('assigned_lineage', 'cell_ids')]

data <- as.data.frame(t(data))
data$cell_ids <- rownames(data)
data <- merge(data, lineage_info, by = 'cell_ids')

# remove data from lineages that have fewer than 4 cells
lineage_size <- lineage_info %>% 
  group_by(assigned_lineage) %>% 
  summarise(size = n())
lineage_size$larger.than.4 <- ifelse(lineage_size$size >= 4, 'YES', 'NO')
lineage_size <- lineage_size[lineage_size$larger.than.4 == 'YES', ]

data <- data[data$assigned_lineage %in% lineage_size$assigned_lineage, ]

# ==============================================================================
# Perform ANOVA
# ==============================================================================

result_df <- data.frame(matrix(ncol=3, nrow=0))
colnames(result_df) <- c('feature', 'F_val', 'p_val')

for (f in features) {
  gene_df <- data[, c('cell_ids', 'assigned_lineage', f)]
  colnames(gene_df) <- c('cell_ids', 'assigned_lineage', 'feature')
  res.aov <- aov(feature ~ assigned_lineage, data = gene_df)
  summary.aov <- summary(res.aov)
  f_val <- summary.aov[[1]][["F value"]][1]
  p_val <- summary.aov[[1]][["Pr(>F)"]][1]
  
  feature_result_df <- data.frame(matrix(ncol=3, nrow=0))
  colnames(feature_result_df) <- c('feature', 'F_val', 'p_val')
  feature_result_df[1, ] <- list(f, f_val, p_val)
  
  result_df <- rbind(result_df, feature_result_df)
}

write.csv(result_df, paste0(output_dir, SAMPLE_NAME, '_RNA_pvals.csv'), row.names=FALSE)


