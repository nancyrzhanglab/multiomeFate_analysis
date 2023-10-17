library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(matrixStats)

args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]
# ==============================================================================
# Read data
# ==============================================================================

dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/data/'
out_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/outs/ANOVA/'
sc_obj <- readRDS(paste0(dir, file_name, '.rds'))

data <- sc_obj@assays[["Saver"]]@scale.data
metadat <- sc_obj@meta.data
metadat$cell_ids <- row.names(metadat)

# ==============================================================================
# Subset data
# ==============================================================================
# Subset for testing
features <- row.names(data)
data <- data[features, ]

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

write.csv(result_df, paste0(out_dir, file_name, '_RNA_pvals.csv'), row.names=FALSE)
# ggplot(data = result_df) +
#   geom_histogram(aes(x=p_val), bins=100)
# sig <- result_df[result_df$p_val < 0.05 / nrow(result_df), ]

# ==============================================================================
# Plot expression levels
# ==============================================================================
# data_aov_sig <- data[, c('cell_ids', sig$feature)]
# data_aov_sig <- melt(data_aov_sig, id.vars=c('cell_ids'))
# colnames(data_aov_sig) <- c('cell_ids', 'gene', 'scale.data')
# ggplot(data = data_aov_sig) +
#   geom_violin(aes(x=gene, y=scale.data), scale = "width") +
#   geom_jitter(aes(x=gene, y=scale.data), size = 0.5, alpha = 0.3) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

