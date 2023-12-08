library(tidyverse)
library(data.table)

dir <- '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
future_treatment <- 'day10_CIS'
keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
  DABTRAM = sort(c("AXL", "EGFP", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
  COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
                 "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7", "CAV1"))
)

# ==============================================================================
# Read data
# ==============================================================================

# growth potential
load(paste0(dir, future_treatment, '/Writeup6n_CIS_day10_lineage-imputation_stepdown_concise-postprocessed.RData'))
growth_potential_use <- cell_imputed_count
growth_potential_use <- growth_potential_use[!is.na(growth_potential_use)]

# adaptation genes
adaptation_genes <- read.csv(paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", future_treatment, "_adaptation_genes.csv"))

# day0 data
file_name <- 'day0_processed.rds'
sc_obj <- readRDS(paste0(dir, 'Raw_and_Processed/', file_name))
data <- sc_obj@assays[["Saver"]]@scale.data
metadata_day0 <- read.csv(paste0(dir, 'day0/day0_meta.csv'), row.names = 1)
metadata_day0$cell_id <- rownames(metadata_day0)

# day10 data
metadata_day10 <- read.csv(paste0(dir, future_treatment, '/', future_treatment, '_meta.csv'), row.names = 1)
metadata_day10$cell_id <- rownames(metadata_day10)

# ==============================================================================
# Wrangle data
# ==============================================================================

# =====================================
# find winning lineages
# =====================================
winners <- as.data.frame(growth_potential_use[growth_potential_use > 0])
others <- as.data.frame(growth_potential_use[growth_potential_use <= 0])

colnames(winners) <- c('growth_potential')
winners$category  <- 'winners'
winners$cell_id <- row.names(winners)
winners <- merge(winners, metadata_day10[, c('assigned_lineage', 'cell_id')], by='cell_id')

metadata_day10$winning_lineage <- ifelse(metadata_day10$assigned_lineage %in% winners$assigned_lineage, 'winning_lineage', 'other_lineage')

day10_lineage_status <- metadata_day10[, c('assigned_lineage', 'winning_lineage')]
row.names(day10_lineage_status) <- NULL
day10_lineage_status <- unique( day10_lineage_status )

day10_winners <- day10_lineage_status[day10_lineage_status$winning_lineage == 'winning_lineage', ]
day10_others <- day10_lineage_status[day10_lineage_status$winning_lineage == 'other_lineage', ]

# =====================================
# pick out cells at day0 in winner vs other lineages 
# =====================================
metadata_day0 <- metadata_day0[metadata_day0$assigned_lineage %in% day10_lineage_status$assigned_lineage, ]

metadata_day0$day10_lineage_status <- ifelse(metadata_day0$assigned_lineage %in% day10_winners$assigned_lineage, 'day10_winning_lineage', 'NA')
metadata_day0$day10_lineage_status <- ifelse(metadata_day0$assigned_lineage %in% day10_others$assigned_lineage, 'day10_other_lineage', metadata_day0$day10_lineage_status)

# metadata_day0_day10_lineage_status_summary <- metadata_day0 %>%
#   group_by(day10_lineage_status) %>%
#   summarise(size = n_distinct(cell_id))


winners_df <- metadata_day0[metadata_day0$day10_lineage_status == 'day10_winning_lineage', ]
others_df <- metadata_day0[metadata_day0$day10_lineage_status == 'day10_other_lineage', ]

data_winners <- as.data.frame(data[, winners_df$cell_id])
data_others <- as.data.frame(data[, others_df$cell_id])

# =====================================
# calculate the mean expression levels
# =====================================
data_winners$winners_mean <- rowMeans(data_winners, dims = 1)
data_others$others_mean <- rowMeans(data_others, dims = 1)

data_winners$feature <- rownames(data_winners)
data_others$feature <- rownames(data_others)

data_winners <- data_winners[, c('feature', 'winners_mean')]
data_others <- data_others[, c('feature','others_mean')]


# =====================================
# looking at raw data
# =====================================
data_winners_all <- as.data.frame(data[, winners_df$cell_id])
data_others_all <- as.data.frame(data[, others_df$cell_id])

data_winners_all <- as.data.frame(t(data_winners_all))
data_others_all <- as.data.frame(t(data_others_all))

data_winners_all$category <- 'winners'
data_others_all$category <- 'others'

data_all <- as.data.frame(rbind(data_winners_all, data_others_all))

adaptation_genes <- adaptation_genes[adaptation_genes$gene %in% colnames(data_all), ]
# adaptation_genes$abs_correlation <- abs(adaptation_genes$correlation)
# adaptation_genes$corr_dir <- ifelse(adaptation_genes$correlation > 0, 'Pos', 'None')
# adaptation_genes$corr_dir <- ifelse(adaptation_genes$correlation < 0, 'Neg', adaptation_genes$corr_dir)
# adaptation_genes <- adaptation_genes %>%
#   arrange(desc(abs_correlation))
# top_adaptation_genes <- head(adaptation_genes, 20)
# bottom_adaptation_genes <- tail(adaptation_genes, 5)
# other_genes <- adaptation_genes %>% 
#   filter(! gene %in% top_adaptation_genes$gene) %>% 
#   filter(! gene %in% bottom_adaptation_genes$gene) %>% 
#   sample_n(10)
# 
# top_adaptation_genes$gene_category <- 'Top'
# bottom_adaptation_genes$gene_category <- 'Bottom'
# other_genes$gene_category <- 'Others'

# genes_to_check <- rbind(top_adaptation_genes, bottom_adaptation_genes, other_genes)
genes_to_check <- adaptation_genes

data_to_check <- data_all[, c(genes_to_check$gene, 'category')]
data_to_check$cell_id <- row.names(data_to_check)
data_to_check <- melt(data_to_check, id_vars=c('cell_id', 'category'))
colnames(data_to_check) <- c('cell_lineage_category', 'cell_id', 'gene', 'RNA_SAVER')

data_to_check <- merge(data_to_check, genes_to_check, by='gene')

# =====================================
# calculate the difference in expression level between winner and other cells
# =====================================
data_mean <- merge(data_winners, data_others, by='feature')
data_mean$winner_minus_others <- data_mean$winners_mean - data_mean$others_mean

results <- merge(adaptation_genes, data_mean, by.x = c('gene'), by.y =c('feature'), all=TRUE)
results$corr_dir <- ifelse(results$correlation > 0, 'Pos', 'None')
results$corr_dir <- ifelse(results$correlation < 0, 'Neg', results$corr_dir)

results[is.na(results)] <- 1
results$corr_dir <- ifelse(results$corr_dir == 1, 'Others', results$corr_dir)
# results <- results[results$winner_minus_others != 1, ]
# results <- results[results$correlation != 1, ]

genes_of_interest <- c('BACE2', 'ACTB', 'FN1', 'BAAT', 'MYO1D', 'PSD3', 'ANXA2')
ggplot(results,aes(x = correlation, y = winner_minus_others)) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(results, gene %in% genes_of_interest),
                           ggplot2::aes(label = gene),
                           box.padding = ggplot2::unit(2, 'lines'),
                           point.padding = ggplot2::unit(0.2, 'lines'),
                           color = 'red',
                           max.overlaps = 80) +
  xlim(c(-0.6, 0.6)) +
  ylim(c(-0.3, 0.3))

# =====================================
# Performing t-test
# =====================================
columns <- c('gene', 'mean_winner', 'mean_other', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

for (g in genes_to_check$gene) {
  g_winner <- data_to_check[data_to_check$gene == g & data_to_check$cell_lineage_category == 'winners', ]
  g_other <- data_to_check[data_to_check$gene == g & data_to_check$cell_lineage_category == 'others', ]
  g_t_test <- t.test(g_winner$RNA_SAVER, g_other$RNA_SAVER, alternative = 'two.sided')
  t_statistics <- g_t_test[["statistic"]][["t"]]
  t_test_p_val <- g_t_test[["p.value"]] 
  one_gene <- data.frame(matrix(ncol = 5, nrow = 1))
  colnames(one_gene) <- columns
  one_gene$t_statistic <- t_statistics
  one_gene$p_val <- t_test_p_val
  one_gene$mean_winner <- mean(g_winner$RNA_SAVER)
  one_gene$mean_other <- mean(g_other$RNA_SAVER)
  one_gene$gene <- g
  
  t_test_results <- rbind(t_test_results, one_gene)
}

t_test_results$p_val_adjust <- p.adjust(t_test_results$p_val, 'BH')
t_test_results$neg_log10_pval <- (-1) * log10(t_test_results$p_val)
# =====================================
# annotate genes
# =====================================
t_test_results$gene_set <- 'Others'
t_test_results$gene_set <- ifelse(t_test_results$gene %in% keygenes['jackpot'][[1]], 'jackpot', t_test_results$gene_set)
t_test_results$gene_set <- ifelse(t_test_results$gene %in% keygenes['DABTRAM'][[1]], 'DABTRAM', t_test_results$gene_set)
t_test_results$gene_set <- ifelse(t_test_results$gene %in% keygenes['COCL2'][[1]], 'COCL2', t_test_results$gene_set)
t_test_results$winner_minus_others <- t_test_results$mean_winner - t_test_results$mean_other

# ==============================================================================
#  Plotting
# ==============================================================================

t_test_results <- read_csv(paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", future_treatment, "_adaptation_genes_diff_t_tests_on_day0_cells.csv"))
t_test_results_sig <- t_test_results[t_test_results$p_val_adjust < 0.05, ]
results$sig_by_t_test <- ifelse(results$gene %in% t_test_results_sig$gene, 'Yes', 'No')
results_sig <- results[results$sig_by_t_test == 'Yes', ]

ggplot(results, aes(x = factor(corr_dir, levels = c('Pos', 'Neg', 'Others')), y=winner_minus_others)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha=0.5, size=1, aes(color=sig_by_t_test)) +
  scale_color_manual(values=c("#A9A9A9", "#ff0000")) +
  xlab('Correlation direction') +
  ylim(c(-0.4, 0.4)) +
  theme_bw()

ggsave(paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", future_treatment, "_adaptation_genes_diff_t_tests_on_day0_cells.png"))

ggplot(data_to_check, aes(x = reorder(gene, -abs_correlation), y = RNA_SAVER, fill=cell_lineage_category, color=corr_dir)) +
  geom_boxplot(lwd=0.8) +
  scale_fill_manual(values=c("#A7C7E7", "#ff0000")) +
  scale_color_manual(values=c("#A9A9A9", "#50C878")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust = 0.2, hjust=0.2)) +
  coord_cartesian(ylim = c(-1, 1))

# =====================================
# plot volcano plots
# =====================================
t_test_results_keygenes <- t_test_results[t_test_results$gene_set != 'Others', ]
ggplot() +
  geom_point(data=t_test_results, aes(x=winner_minus_others, y=neg_log10_pval), color='gray', size=1) +
  geom_point(data=t_test_results_keygenes, aes(x=winner_minus_others, y=neg_log10_pval), color='red', size=1) +
  geom_hline(yintercept = 2.6, linetype='dashed') +
  xlim(c(-0.5, 0.5)) +
  theme_minimal()

hist(t_test_results$p_val, breaks=100)

# write.csv(t_test_results, paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", future_treatment, "_adaptation_genes_diff_t_tests_on_day0_cells.csv"), row.names = FALSE)
