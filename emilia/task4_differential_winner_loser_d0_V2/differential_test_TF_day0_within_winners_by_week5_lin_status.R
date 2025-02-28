library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
remove_unassigned_cells <- TRUE

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))

all_data@misc <- all_data_fatepotential
chromVAR_day0 <- all_data_chromVar_day0

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

## Read features to test
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_TFs.RData')

jun <- colnames(chromVAR_day0)[grepl('JUN', colnames(chromVAR_day0))]
fos <- colnames(chromVAR_day0)[grepl('FOS', colnames(chromVAR_day0))]
mitf <- colnames(chromVAR_day0)[grepl('MITF', colnames(chromVAR_day0))]
sox10 <- colnames(chromVAR_day0)[grepl('SOX10', colnames(chromVAR_day0))]
snai <- colnames(chromVAR_day0)[grepl('SNAI', colnames(chromVAR_day0))]
tead <- colnames(chromVAR_day0)[grepl('TEAD', colnames(chromVAR_day0))]

keyTFs <- unique(c(jun, fos, mitf, sox10, snai, tead))
# ==============================================================================
# Get day0 winner cells
# ==============================================================================

# get fate potential
cur_time = 'd0'
fut_time = 'd10'
treatment = 'DABTRAM'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

# get winner cells
fp_median <- median(fp[[fp_name]])
fp_winner <- fp[fp[[fp_name]] > fp_median, ]

# get saver
metadat <- all_data@meta.data
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
metadat.day0$cell_id <- rownames(metadat.day0)

chromVAR_day0 <- t(all_data_chromVar_day0@data)
chromVAR_day0 <- chromVAR_day0[metadat.day0$cell_id, ]

chromVAR_day0_winner <- chromVAR_day0[fp_winner$cell_id, ]

# ==============================================================================
# Get lineages survived to Week5 DABTRAM
# ==============================================================================
lin.week5 <- metadat[metadat$dataset == 'week5_DABTRAM', ]$assigned_lineage
lin.day0 <- metadat[metadat$dataset == 'day0', ]$assigned_lineage

lin.day0.week5 <- intersect(lin.day0, lin.week5)

fp_winner <- merge(fp_winner, metadat.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
fp_winner$LinInWeek5 <- ifelse(fp_winner$assigned_lineage %in% lin.day0.week5, 'Yes', 'No')

table(fp_winner$LinInWeek5)

# ==============================================================================
# Differential tests
# ==============================================================================
columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

# features <- lineage_specific_adapation_TFs[[paste0(tolower(treatment), "_", fut_time, "_saver_cor_vec_top25")]][[1]]
features <- colnames(chromVAR_day0_winner)

day0_winner_lin_in_week5 <- fp_winner[fp_winner$LinInWeek5 == 'Yes', 'cell_id']
day0_winner_lin_Notin_week5 <- fp_winner[fp_winner$LinInWeek5 == 'No', 'cell_id']

for (f in features) {
  
  feature_day0_winner_lin_in_week5 <- chromVAR_day0_winner[day0_winner_lin_in_week5, f]
  feature_day0_winner_lin_Notin_week5 <- chromVAR_day0_winner[day0_winner_lin_Notin_week5, f]
  
  feature_day0_winner_lin_in_week5 <- feature_day0_winner_lin_in_week5[!is.na(feature_day0_winner_lin_in_week5)]
  feature_day0_winner_lin_Notin_week5 <- feature_day0_winner_lin_Notin_week5[!is.na(feature_day0_winner_lin_Notin_week5)]
  
  variance <- var(feature_day0_winner_lin_in_week5) + var(feature_day0_winner_lin_Notin_week5)
  
  if(variance == 0) {
    next
  }
  res <- t.test(feature_day0_winner_lin_in_week5,
                feature_day0_winner_lin_Notin_week5,
                alternative = 'two.sided')
  
  t_statistics <- res[["statistic"]][["t"]]
  t_test_p_val <- res[["p.value"]] 
  
  t_test_results[nrow(t_test_results) + 1, ] <- c(
    f, 
    mean(feature_day0_winner_lin_in_week5), 
    mean(feature_day0_winner_lin_Notin_week5), 
    t_statistics, 
    t_test_p_val
  )
}

t_test_results$p_val <- as.numeric(t_test_results$p_val)
t_test_results$p_val_adj <- p.adjust(t_test_results$p_val, method = 'BH')
t_test_results$neg_log10_p_val <- -log10(t_test_results$p_val)

t_test_results$mean_winner_lin_in_week5 <- as.numeric(t_test_results$mean_winner_lin_in_week5)
t_test_results$mean_winner_lin_notin_week5 <- as.numeric(t_test_results$mean_winner_lin_notin_week5)
t_test_results$mean_diff <- t_test_results$mean_winner_lin_in_week5 - t_test_results$mean_winner_lin_notin_week5
# write.csv(t_test_results, paste0('~/Downloads/t_test.csv'), row.names = FALSE)

t_test_results.key.TF <- t_test_results[t_test_results$feature %in% keyTFs, ]
t_test_results.key.TF <- t_test_results.key.TF[!(grepl('var.2', t_test_results.key.TF$feature)), ]

# t_test_results.key.TF <- t_test_results[t_test_results$feature %in% keyTFs, ]
p_val_thres <- t_test_results[t_test_results$p_val_adj < 0.05, ]
if(nrow(p_val_thres) == 0) {
  p_val_thres <- max(t_test_results$neg_log10_p_val)
}else {
  p_val_thres <- min(p_val_thres$neg_log10_p_val)
}



ggplot(t_test_results, aes(x = mean_diff, 
                           y = neg_log10_p_val)) +
  geom_point(color = 'gray') +
  geom_point(data = t_test_results.key.TF, aes(x = mean_diff, y = neg_log10_p_val), color = 'red') +
  ggrepel::geom_text_repel(data = t_test_results.key.TF, 
                           aes(label = feature), 
                           box.padding = 0.5, max.overlaps = 20) +
  geom_hline(yintercept = p_val_thres, linetype = 'dashed') +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = 'Day 0 (DABTRAM)') +
  theme_Publication() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))
ggsave('~/Downloads/volcano_d0_DABTRAM_TF.png', width = 5, height = 5, dpi = 300)


