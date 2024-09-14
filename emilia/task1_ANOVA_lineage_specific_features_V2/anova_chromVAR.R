rm(list=ls())
library(Seurat)
library(dplyr)

# =============================================================================
# Load data
# =============================================================================

TIME = 'day0' # 'day10', or 'day0'
TREATMENT = 'COCL2' # 'COCL2', 'DABTRAM', or 'CIS'
# SAMPLE_NAME = paste0(TIME, '_', TREATMENT)
SAMPLE_NAME = TIME
output_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task1_ANOVA_lineage_specific_features_V2/"

all_data <- multiomeFate::data_loader(which_files = c("chromvar"))

# =============================================================================
# Wrangle data
# =============================================================================

all_data_use = subset(all_data, dataset == SAMPLE_NAME)

data <- all_data[[paste0("chromVar.", SAMPLE_NAME)]]@data
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
features <- colnames(data)[seq(2, ncol(data)-1)]
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

write.csv(result_df, paste0(output_dir, SAMPLE_NAME, '_chromVAR_pvals.csv'), row.names=FALSE)


