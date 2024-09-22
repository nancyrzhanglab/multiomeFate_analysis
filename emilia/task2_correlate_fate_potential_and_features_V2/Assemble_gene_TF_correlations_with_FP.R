library(tidyverse)

# GENE EXPRESSION
# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/saver_cor_vec.RData')

dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']] 

dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']] 

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_genes_saver.RData')
dabtram_d10_saver_cor_vec_top25 <- lineage_specific_adapation_TFs[['dabtram_d10_saver_cor_vec_top25']]
cocl2_d10_saver_cor_vec_top25 <- lineage_specific_adapation_TFs[['cocl2_d10_saver_cor_vec_top25']]
cis_d10_saver_cor_vec_top25 <- lineage_specific_adapation_TFs[['cis_d10_saver_cor_vec_top25']]

d0_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day0_SAVER_RNA_processed_ANOVA_pvals.csv')
dabtram_d10_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day10_DABTRAM_SAVER_RNA_processed_ANOVA_pvals.csv')
dabtram_wk5_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/week5_DABTRAM_SAVER_RNA_processed_ANOVA_pvals.csv')

cocl2_d10_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day10_COCL2_SAVER_RNA_processed_ANOVA_pvals.csv')
cocl2_wk5_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/week5_COCL2_SAVER_RNA_processed_ANOVA_pvals.csv')

cis_d10_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day10_CIS_SAVER_RNA_processed_ANOVA_pvals.csv')
cis_wk5_saver_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/week5_CIS_SAVER_RNA_processed_ANOVA_pvals.csv')

# load t-test results
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/differential_winner_loser_saver_lineage_specific_adaptation_gene_DABTRAM_t_test_results.RData')
dabtram_ttest <- t_test_results
colnames(dabtram_ttest) <- c('gene', 'mean_winner', 'mean_other', 't_statistics', 'p_val')
dabtram_ttest$mean_winner <- as.numeric(dabtram_ttest$mean_winner)
dabtram_ttest$mean_other <- as.numeric(dabtram_ttest$mean_other)
dabtram_ttest$mean_diff <- dabtram_ttest$mean_winner - dabtram_ttest$mean_other
dabtram_ttest$log2_fold_change <- log2(dabtram_ttest$mean_winner / dabtram_ttest$mean_other)
dabtram_ttest$t_test_p_val_adj <- p.adjust(dabtram_ttest$p_val, method = 'BH')

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/differential_winner_loser_saver_lineage_specific_adaptation_gene_COCL2_t_test_results.RData')
cocl2_ttest <- t_test_results
colnames(cocl2_ttest) <- c('gene', 'mean_winner', 'mean_other', 't_statistics', 'p_val')
cocl2_ttest$mean_winner <- as.numeric(cocl2_ttest$mean_winner)
cocl2_ttest$mean_other <- as.numeric(cocl2_ttest$mean_other)
cocl2_ttest$mean_diff <- cocl2_ttest$mean_winner - cocl2_ttest$mean_other
cocl2_ttest$log2_fold_change <- log2(cocl2_ttest$mean_winner / cocl2_ttest$mean_other)
cocl2_ttest$t_test_p_val_adj <- p.adjust(cocl2_ttest$p_val, method = 'BH')

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/differential_winner_loser_saver_lineage_specific_adaptation_gene_CIS_t_test_results.RData')
cis_ttest <- t_test_results
colnames(cis_ttest) <- c('gene', 'mean_winner', 'mean_other', 't_statistics', 'p_val')
cis_ttest$mean_winner <- as.numeric(cis_ttest$mean_winner)
cis_ttest$mean_other <- as.numeric(cis_ttest$mean_other)
cis_ttest$mean_diff <- cis_ttest$mean_winner - cis_ttest$mean_other
cis_ttest$log2_fold_change <- log2(cis_ttest$mean_winner / cis_ttest$mean_other)
cis_ttest$t_test_p_val_adj <- p.adjust(cis_ttest$p_val, method = 'BH')

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_saver_cor_vec) <- paste0(colnames(dabtram_d0_saver_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_saver_cor_vec) <- paste0(colnames(cocl2_d0_saver_cor_vec), '.COCL2_d0')
colnames(cis_d0_saver_cor_vec) <- paste0(colnames(cis_d0_saver_cor_vec), '.CIS_d0')

colnames(dabtram_d10_saver_cor_vec) <- paste0(colnames(dabtram_d10_saver_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_saver_cor_vec) <- paste0(colnames(cocl2_d10_saver_cor_vec), '.COCL2_d10')
colnames(cis_d10_saver_cor_vec) <- paste0(colnames(cis_d10_saver_cor_vec), '.CIS_d10')

colnames(dabtram_d10_saver_anova) <- c('gene', 'ANOVA_F_val.DABTRAM_d10', 'ANOVA_p_val.DABTRAM_d10')
colnames(cocl2_d10_saver_anova) <- c('gene', 'ANOVA_F_val.COCL2_d10', 'ANOVA_p_val.COCL2_d10')
colnames(cis_d10_saver_anova) <- c('gene', 'ANOVA_F_val.CIS_d10', 'ANOVA_p_val.CIS_d10')

dabtram_d10_saver_anova$ANOVA_p_val_adj.DABTRAM_d10 <- p.adjust(dabtram_d10_saver_anova$ANOVA_p_val.DABTRAM_d10, method = 'BH')
cocl2_d10_saver_anova$ANOVA_p_val_adj.COCL2_d10 <- p.adjust(cocl2_d10_saver_anova$ANOVA_p_val.COCL2_d10, method = 'BH')
cis_d10_saver_anova$ANOVA_p_val_adj.CIS_d10 <- p.adjust(cis_d10_saver_anova$ANOVA_p_val.CIS_d10, method = 'BH')

# Join d0 and d10 data

# Dabtram
dabtram_saver_cor_vec <- merge(dabtram_d0_saver_cor_vec, dabtram_d10_saver_cor_vec, by = 'row.names') %>% 
  rename('gene' = 'Row.names')
dabtram_saver_cor_vec <- merge(dabtram_saver_cor_vec, dabtram_d10_saver_anova[ , c('gene', 'ANOVA_p_val_adj.DABTRAM_d10')], by = 'gene')
dabtram_saver_cor_vec$adaptation_gene <- ifelse(dabtram_saver_cor_vec$gene %in% dabtram_d10_saver_cor_vec_top25$feature, 'yes', 'no')

# Cocl2
cocl2_saver_cor_vec <- merge(cocl2_d0_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names') %>% 
  rename('gene' = 'Row.names')
cocl2_saver_cor_vec <- merge(cocl2_saver_cor_vec, cocl2_d10_saver_anova[ , c('gene', 'ANOVA_p_val_adj.COCL2_d10')], by = 'gene')
cocl2_saver_cor_vec$adaptation_gene <- ifelse(cocl2_saver_cor_vec$gene %in% cocl2_d10_saver_cor_vec_top25$feature, 'yes', 'no')

# Cis
cis_saver_cor_vec <- merge(cis_d0_saver_cor_vec, cis_d10_saver_cor_vec, by = 'row.names') %>% 
  rename('gene' = 'Row.names')
cis_saver_cor_vec <- merge(cis_saver_cor_vec, cis_d10_saver_anova[ , c('gene', 'ANOVA_p_val_adj.CIS_d10')], by = 'gene')
cis_saver_cor_vec$adaptation_gene <- ifelse(cis_saver_cor_vec$gene %in% cis_d10_saver_cor_vec_top25$feature, 'yes', 'no')

dabtram_saver_cor_vec <- merge(dabtram_saver_cor_vec, dabtram_ttest[ , c('gene', 'mean_diff', 'log2_fold_change', 't_test_p_val_adj')], by = 'gene', all.x = TRUE)
cocl2_saver_cor_vec <- merge(cocl2_saver_cor_vec, cocl2_ttest[ , c('gene', 'mean_diff', 'log2_fold_change', 't_test_p_val_adj')], by = 'gene', all.x = TRUE)
cis_saver_cor_vec <- merge(cis_saver_cor_vec, cis_ttest[ , c('gene', 'mean_diff', 'log2_fold_change', 't_test_p_val_adj')], by = 'gene', all.x = TRUE)

write.csv(dabtram_saver_cor_vec, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_saver_cor_d0_d10.csv', row.names = F)
write.csv(cocl2_saver_cor_vec, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_saver_cor_d0_d10.csv', row.names = F)
write.csv(cis_saver_cor_vec, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_saver_cor_d0_d10.csv', row.names = F)

# CHROMVAR
# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d0_chromVAR_cor_vec']] 
cocl2_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d0_chromVAR_cor_vec']] 
cis_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d0_chromVAR_cor_vec']] 

dabtram_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d10_chromVAR_cor_vec']] 
cocl2_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d10_chromVAR_cor_vec']] 
cis_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d10_chromVAR_cor_vec']] 

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_TFs.RData')
dabtram_d10_chromVAR_cor_vec_top25 <- lineage_specific_adapation_TFs[['dabtram_d10_chromVAR_cor_vec_top25']]
cocl2_d10_chromVAR_cor_vec_top25 <- lineage_specific_adapation_TFs[['cocl2_d10_chromVAR_cor_vec_top25']]
cis_d10_chromVAR_cor_vec_top25 <- lineage_specific_adapation_TFs[['cis_d10_chromVAR_cor_vec_top25']]

dabtram_d10_chromVAR_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day10_DABTRAM_chromVAR_pvals.csv')

cocl2_d10_chromVAR_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day10_COCL2_chromVAR_pvals.csv')

cis_d10_chromVAR_anova <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/day10_CIS_chromVAR_pvals.csv')

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_chromVAR_cor_vec) <- paste0(colnames(dabtram_d0_chromVAR_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_chromVAR_cor_vec) <- paste0(colnames(cocl2_d0_chromVAR_cor_vec), '.COCL2_d0')
colnames(cis_d0_chromVAR_cor_vec) <- paste0(colnames(cis_d0_chromVAR_cor_vec), '.CIS_d0')

colnames(dabtram_d10_chromVAR_cor_vec) <- paste0(colnames(dabtram_d10_chromVAR_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_chromVAR_cor_vec) <- paste0(colnames(cocl2_d10_chromVAR_cor_vec), '.COCL2_d10')
colnames(cis_d10_chromVAR_cor_vec) <- paste0(colnames(cis_d10_chromVAR_cor_vec), '.CIS_d10')

colnames(dabtram_d10_chromVAR_anova) <- c('TF', 'ANOVA_F_val.DABTRAM_d10', 'ANOVA_p_val.DABTRAM_d10')
colnames(cocl2_d10_chromVAR_anova) <- c('TF', 'ANOVA_F_val.COCL2_d10', 'ANOVA_p_val.COCL2_d10')
colnames(cis_d10_chromVAR_anova) <- c('TF', 'ANOVA_F_val.CIS_d10', 'ANOVA_p_val.CIS_d10')

dabtram_d10_chromVAR_anova$ANOVA_p_val_adj.DABTRAM_d10 <- p.adjust(dabtram_d10_chromVAR_anova$ANOVA_p_val.DABTRAM_d10, method = 'BH')
cocl2_d10_chromVAR_anova$ANOVA_p_val_adj.COCL2_d10 <- p.adjust(cocl2_d10_chromVAR_anova$ANOVA_p_val.COCL2_d10, method = 'BH')
cis_d10_chromVAR_anova$ANOVA_p_val_adj.CIS_d10 <- p.adjust(cis_d10_chromVAR_anova$ANOVA_p_val.CIS_d10, method = 'BH')

# Join d0 and d10 data

# Dabtram
dabtram_chromVAR_cor_vec <- merge(dabtram_d0_chromVAR_cor_vec, dabtram_d10_chromVAR_cor_vec, by = 'row.names') %>% 
  rename('TF' = 'Row.names')
dabtram_chromVAR_cor_vec <- merge(dabtram_chromVAR_cor_vec, dabtram_d10_chromVAR_anova[, c('TF', 'ANOVA_p_val_adj.DABTRAM_d10')], by = 'TF')
dabtram_chromVAR_cor_vec$adaptation_TF <- ifelse(dabtram_chromVAR_cor_vec$TF %in% dabtram_d10_chromVAR_cor_vec_top25$feature, 'yes', 'no')

# Cocl2
cocl2_chromVAR_cor_vec <- merge(cocl2_d0_chromVAR_cor_vec, cocl2_d10_chromVAR_cor_vec, by = 'row.names') %>% 
  rename('TF' = 'Row.names')
cocl2_chromVAR_cor_vec <- merge(cocl2_chromVAR_cor_vec, cocl2_d10_chromVAR_anova[, c('TF', 'ANOVA_p_val_adj.COCL2_d10')], by = 'TF')
cocl2_chromVAR_cor_vec$adaptation_TF <- ifelse(cocl2_chromVAR_cor_vec$TF %in% cocl2_d10_chromVAR_cor_vec_top25$feature, 'yes', 'no')

# Cis
cis_chromVAR_cor_vec <- merge(cis_d0_chromVAR_cor_vec, cis_d10_chromVAR_cor_vec, by = 'row.names') %>% 
  rename('TF' = 'Row.names')
cis_chromVAR_cor_vec <- merge(cis_chromVAR_cor_vec, cis_d10_chromVAR_anova[, c('TF', 'ANOVA_p_val_adj.CIS_d10')], by = 'TF')
cis_chromVAR_cor_vec$adaptation_TF <- ifelse(cis_chromVAR_cor_vec$TF %in% cis_d10_chromVAR_cor_vec_top25$feature, 'yes', 'no')


# load t-test results
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/differential_winner_loser_chromVAR_lineage_specific_adapation_TF_DABTRAM_t_test_results.RData')
dabtram_ttest <- t_test_results
colnames(dabtram_ttest) <- c('TF', 'mean_winner', 'mean_other', 't_statistics', 'p_val')
dabtram_ttest$mean_winner <- as.numeric(dabtram_ttest$mean_winner)
dabtram_ttest$mean_other <- as.numeric(dabtram_ttest$mean_other)
dabtram_ttest$mean_diff <- dabtram_ttest$mean_winner - dabtram_ttest$mean_other
dabtram_ttest$t_test_p_val_adj <- p.adjust(dabtram_ttest$p_val, method = 'BH')

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/differential_winner_loser_chromVAR_lineage_specific_adaptation_TF_COCL2_t_test_results.RData')
cocl2_ttest <- t_test_results
colnames(cocl2_ttest) <- c('TF', 'mean_winner', 'mean_other', 't_statistics', 'p_val')
cocl2_ttest$mean_winner <- as.numeric(cocl2_ttest$mean_winner)
cocl2_ttest$mean_other <- as.numeric(cocl2_ttest$mean_other)
cocl2_ttest$mean_diff <- cocl2_ttest$mean_winner - cocl2_ttest$mean_other
cocl2_ttest$t_test_p_val_adj <- p.adjust(cocl2_ttest$p_val, method = 'BH')

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/differential_winner_loser_chromVAR_lineage_specific_adaptation_TF_CIS_t_test_results.RData')
cis_ttest <- t_test_results
colnames(cis_ttest) <- c('TF', 'mean_winner', 'mean_other', 't_statistics', 'p_val')
cis_ttest$mean_winner <- as.numeric(cis_ttest$mean_winner)
cis_ttest$mean_other <- as.numeric(cis_ttest$mean_other)
cis_ttest$mean_diff <- cis_ttest$mean_winner - cis_ttest$mean_other
cis_ttest$t_test_p_val_adj <- p.adjust(cis_ttest$p_val, method = 'BH')

dabtram_chromVAR_cor_vec <- merge(dabtram_chromVAR_cor_vec, dabtram_ttest[, c('TF', 'mean_diff', 't_test_p_val_adj')], by = 'TF', all.x = TRUE)
cocl2_chromVAR_cor_vec <- merge(cocl2_chromVAR_cor_vec, cocl2_ttest[, c('TF', 'mean_diff', 't_test_p_val_adj')], by = 'TF', all.x = TRUE)
cis_chromVAR_cor_vec <- merge(cis_chromVAR_cor_vec, cis_ttest[, c('TF', 'mean_diff', 't_test_p_val_adj')], by = 'TF', all.x = TRUE)

# write.csv(dabtram_chromVAR_cor_vec, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_chromVAR_cor_d0_d10.csv', row.names = F)
# write.csv(cocl2_chromVAR_cor_vec, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_chromVAR_cor_d0_d10.csv', row.names = F)
# write.csv(cis_chromVAR_cor_vec, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_chromVAR_cor_d0_d10.csv', row.names = F)
