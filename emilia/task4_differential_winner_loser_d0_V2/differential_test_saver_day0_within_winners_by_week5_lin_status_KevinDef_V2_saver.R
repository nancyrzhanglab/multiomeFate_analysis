rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
  DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
  COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
                 "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
                 "CAV1"))
)
keygenes <- unlist(keygenes)

treatment <- 'DABTRAM'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
# load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))

all_data@misc <- all_data_fatepotential
# all_data[['chromVar_day0']] <- all_data_chromVar_day0
all_data[['saver']] <- all_data_saver

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

all_data_day10 <- subset(all_data, dataset == paste0('day10_', treatment))
all_data_day0 <- subset(all_data, dataset == 'day0')

# do the analysis 
tab_mat <- table(all_data_day10$assigned_lineage, all_data_day10$dataset)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# =============================================================================
# Lineage specific genes
# =============================================================================
anova_in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/'
anova_saver_dabtram_d10 <- read.csv(paste0(anova_in_dir, 'day10_DABTRAM_SAVER_RNA_processed_ANOVA_pvals.csv'))
# Take lineage specific saver features
anova_saver_dabtram_d10$p_val_adj <- p.adjust(anova_saver_dabtram_d10$p_val, 'BH')
anova_saver_dabtram_d10_sig <- anova_saver_dabtram_d10[anova_saver_dabtram_d10$p_val_adj < 0.05, ]


# =============================================================================
# Define winner and loser
# =============================================================================
def = 'cell'
# Extract the fate potential scores for each cell
fp_name <- paste0('fatepotential_', treatment, '_d10_w5')

if(def == 'cell') {
  fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['cell_imputed_score']])
  colnames(fp_d10_w5_treatment) <- c('cell_imputed_score')
  
  fp_d10_w5_treatment$cell_id <- rownames(fp_d10_w5_treatment)
  fp_d10_w5_treatment <- merge(fp_d10_w5_treatment, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
  
  # score each lineage by number of cells with non-negative day10 score
  lineage_score <- fp_d10_w5_treatment %>% 
    group_by(assigned_lineage) %>%
    summarise(n_cells = n(),
              n_pos = sum(cell_imputed_score > 0))
  
  # score each lineage by mean day10 score
  lineage_score2 <- fp_d10_w5_treatment %>% 
    group_by(assigned_lineage) %>%
    summarise(n_cells = n(),
              mean_score = mean(cell_imputed_score))
  
  # define winners and losers
  lineage_score_all <- merge(lineage_score, lineage_score2, by = c('assigned_lineage', 'n_cells'))
  lineage_score_all$winner <- ifelse(lineage_score_all$n_pos > 0 & 
                                     lineage_score_all$mean_score > 0, 'Winning', NA)
  lineage_score_all$winner <- ifelse(lineage_score_all$mean_score <= 0, 'Losing', lineage_score_all$winner)
  
} else if (def == 'lineage') {
  fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['lineage_imputed_count']])
  colnames(fp_d10_w5_treatment) <- c('lineage_imputed_count')
  fp_d10_w5_treatment$assigned_lineage <- rownames(fp_d10_w5_treatment)
  
  lineage_score_all <- fp_d10_w5_treatment
  lineage_score_all$winner <- ifelse(lineage_score_all$lineage_imputed_count > 1, 'Winning', NA)
  lineage_score_all$winner <- ifelse(lineage_score_all$lineage_imputed_count < 1, 'Losing', lineage_score_all$winner)
}



table(lineage_score_all$winner)

# define winners and losers in day0
metadat.day0 <- all_data_day0@meta.data
metadat.day0$isWinningLineage <- ifelse(metadat.day0$assigned_lineage %in% lineage_score_all$assigned_lineage[lineage_score_all$winner == 'Winning'], 'Winning', NA)
metadat.day0$isWinningLineage <- ifelse(metadat.day0$assigned_lineage %in% lineage_score_all$assigned_lineage[lineage_score_all$winner == 'Losing'], 'Losing', metadat.day0$isWinningLineage)

# metadat.day0$isWinningLineage <- ifelse(metadat.day0$assigned_lineage %in% lineage_score_all$assigned_lineage[lineage_score_all$winner == 'Winning'], 'Winning', 'Losing')

table(metadat.day0$isWinningLineage)
all_data_day0 <- AddMetaData(all_data_day0, metadat.day0)

metadat.day10 <- all_data_day10@meta.data
nrow(metadat.day0[metadat.day0$assigned_lineage %in% metadat.day10$assigned_lineage, ])
# =============================================================================
# Differential tests
# =============================================================================

Seurat::Idents(all_data_day0) <- "isWinningLineage"

Seurat::DefaultAssay(all_data_day0) <- "saver"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data_day0,
  ident.1 = "Winning",
  ident.2 = "Losing",
  only.pos = FALSE,
  feature = anova_saver_dabtram_d10_sig$feature,
  min.pct = 0,
  logfc.threshold = -Inf,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  # fc.name = 'avg_log2FC',
  verbose = F
)

summary(de_res)

p_val_thres <- de_res[de_res$p_val_adj < 0.05, ]
if(nrow(p_val_thres) == 0) {
  p_val_thres <- max(de_res$neg_log10_p_val)
}else {
  p_val_thres <- min(de_res$neg_log10_p_val)
}

de_res$neg_log10_p_val <- -log10(de_res$p_val)

de_res.keygenes <- de_res[keygenes, ]
ggplot(de_res) +
  geom_point(aes(x = avg_diff, y = neg_log10_p_val)) +
  geom_point(data = de_res.keygenes, aes(x = avg_diff, y = neg_log10_p_val), color = 'red') +
  geom_hline(yintercept = p_val_thres, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # coord_cartesian(xlim = c(-0.1, 0.1)) +
  theme_minimal()

# =============================================================================
# t test
# =============================================================================
columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns


saver_day0 <- all_data_day0@assays$saver@scale.data
saver_day0 <- t(saver_day0)
features <- colnames(saver_day0) # anova_saver_dabtram_d10_sig$feature


metadat.day0$cell_id <- rownames(metadat.day0)
day0_winner_lin_in_week5 <- metadat.day0 %>% filter(isWinningLineage == 'Winning') %>% pull(cell_id)
day0_winner_lin_Notin_week5 <- metadat.day0 %>% filter(isWinningLineage == 'Losing') %>% pull(cell_id)

for (f in features) {
  
  feature_day0_winner_lin_in_week5 <- saver_day0[day0_winner_lin_in_week5, f]
  feature_day0_winner_lin_Notin_week5 <- saver_day0[day0_winner_lin_Notin_week5, f]
  
  feature_day0_winner_lin_in_week5 <- feature_day0_winner_lin_in_week5[!is.na(feature_day0_winner_lin_in_week5)]
  feature_day0_winner_lin_Notin_week5 <- feature_day0_winner_lin_Notin_week5[!is.na(feature_day0_winner_lin_Notin_week5)]
  
  variance <- var(feature_day0_winner_lin_in_week5) + var(feature_day0_winner_lin_Notin_week5)
  
  if(diff(range(feature_day0_winner_lin_in_week5)) <= 1e-6 || diff(range(feature_day0_winner_lin_Notin_week5)) <= 1e-6 ) {
    next
  }
  
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


p_val_thres <- t_test_results[t_test_results$p_val_adj < 0.05, ]
if(nrow(p_val_thres) == 0) {
  p_val_thres <- max(t_test_results$neg_log10_p_val)
}else {
  p_val_thres <- min(p_val_thres$neg_log10_p_val)
}

t_test_results$log2FC <- log2(t_test_results$mean_winner_lin_in_week5) - log2(t_test_results$mean_winner_lin_notin_week5)
t_test_results.keygenes <- t_test_results[t_test_results$feature %in% keygenes, ]
ggplot(t_test_results, aes(x = mean_diff, 
                           y = neg_log10_p_val)) +
  geom_point(color = 'gray') +
  geom_point(data = t_test_results.keygenes, aes(x = mean_diff, y = neg_log10_p_val), color = 'red') +
  ggrepel::geom_text_repel(data = t_test_results.keygenes,
                           aes(label = feature),
                           box.padding = 0.5, max.overlaps = 20) +
  geom_hline(yintercept = p_val_thres, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = 'Day 0 (DABTRAM)') +
  # theme_Publication() +
  # coord_cartesian(xlim = c(-0.15, 0.15), ylim = c(0, 4.5)) +
  # coord_cartesian(ylim = c(0, 4.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))


t_test_results.sig <- t_test_results %>% 
  filter(p_val_adj < 0.05 &
         mean_diff > 0)

t_test_results %>% 
  filter(p_val_adj < 0.05 ) %>% 
  nrow()


###=== Try Kevin's code ===###
pvalue_list <- lapply(all_data_day0[["saver"]]@var.features, function(gene){
  x_vec <- all_data_day0[["saver"]]@scale.data[gene, which(all_data_day0$isWinningLineage == "Winning")]
  y_vec <- all_data_day0[["saver"]]@scale.data[gene, which(all_data_day0$isWinningLineage == "Losing")]
  
  if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec) - mean(y_vec), pvalue = test_res$p.value)
})
names(pvalue_list) <- all_data_day0[["saver"]]@var.features
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 5)]] <- "3"
labeling_vec[intersect(keygenes, names(labeling_vec))] <- "2"
table(labeling_vec)


df <- data.frame(difference = diff_vec,
                 log10pval = logpval_vec,
                 name = names(logpval_vec),
                 labeling = labeling_vec)

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = difference, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red", "blue"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2", "3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                               color = "red", linewidth=2)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 DE based on lineage's mean growth potential at Day10\n(Top 25% to Bottom 25%)")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)")
p1 <- p1 + Seurat::NoLegend()


de_res$neg_log10_p_val <- -log10(de_res$p_val)
de_res.keygenes <- de_res[keygenes, ]
de_res.keygenes$feature <- rownames(de_res.keygenes)
ggplot(de_res, aes(x = avg_diff, 
                   y = neg_log10_p_val)) +
  geom_point(color = 'gray') +
  geom_point(data = de_res.keygenes, aes(x = avg_diff, y = neg_log10_p_val), color = 'red') +
  ggrepel::geom_text_repel(data = de_res.keygenes,
                           aes(label = feature),
                           box.padding = 0.5, max.overlaps = 20) +
  geom_hline(yintercept = p_val_thres, linetype = 'dashed') +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = 'Day 0 (DABTRAM)') +
  theme_bw()

de_res$feature <- rownames(de_res)
comp <- merge(de_res, t_test_results, by = 'feature')
ggplot(comp, aes(x = neg_log10_p_val.x, y = neg_log10_p_val.y)) +
  geom_point()


