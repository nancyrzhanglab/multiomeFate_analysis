rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

treatment <- 'DABTRAM'
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
# load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))

all_data@misc <- all_data_fatepotential
all_data[['chromVar_day0']] <- all_data_chromVar_day0

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
# Define winner and loser
# =============================================================================

# Extract the fate potential scores for each cell
fp_name <- paste0('fatepotential_', treatment, '_d10_w5')
# fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['cell_imputed_score']])
# colnames(fp_d10_w5_treatment) <- c('cell_imputed_score')

fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['lineage_imputed_count']])
colnames(fp_d10_w5_treatment) <- c('lineage_imputed_count')
fp_d10_w5_treatment$assigned_lineage <- rownames(fp_d10_w5_treatment)

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
                                   lineage_score_all$mean_score >= -1, 'Winning', 'Losing')

# lineage_score_all <- fp_d10_w5_treatment
# lineage_score_all$winner <- ifelse(lineage_score_all$lineage_imputed_count > 1, 'Winning', 'Losing')

table(lineage_score_all$winner)

# define winners and losers in day0
metadat.day0 <- all_data_day0@meta.data
metadat.day0$isWinningLineage <- ifelse(metadat.day0$assigned_lineage %in% lineage_score_all$assigned_lineage[lineage_score_all$winner == 'Winning'], 'Winning', NA)
metadat.day0$isWinningLineage <- ifelse(metadat.day0$assigned_lineage %in% lineage_score_all$assigned_lineage[lineage_score_all$winner == 'Losing'], 'Losing', metadat.day0$isWinningLineage)

table(metadat.day0$isWinningLineage)
all_data_day0 <- AddMetaData(all_data_day0, metadat.day0)
# =============================================================================
# Differential tests
# =============================================================================

Seurat::Idents(all_data_day0) <- "isWinningLineage"

Seurat::DefaultAssay(all_data_day0) <- "chromVar_day0"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data_day0,
  ident.1 = "Winning",
  ident.2 = "Losing",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = -Inf,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)

columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns


chromVAR_day0 <- all_data_day0@assays$chromVar_day0@data
chromVAR_day0 <- t(chromVAR_day0)
features <- colnames(chromVAR_day0)


metadat.day0$cell_id <- rownames(metadat.day0)
day0_winner_lin_in_week5 <- metadat.day0 %>% filter(isWinningLineage == 'Winning') %>% pull(cell_id)
day0_winner_lin_Notin_week5 <- metadat.day0 %>% filter(isWinningLineage == 'Losing') %>% pull(cell_id)

for (f in features) {
  
  feature_day0_winner_lin_in_week5 <- chromVAR_day0[day0_winner_lin_in_week5, f]
  feature_day0_winner_lin_Notin_week5 <- chromVAR_day0[day0_winner_lin_Notin_week5, f]
  
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
  # theme_Publication() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))

de_res$neg_log10_p_val <- -log10(de_res$p_val)
ggplot(de_res, aes(x = avg_diff, 
                           y = neg_log10_p_val)) +
  geom_point(color = 'gray') +
  # geom_point(data = t_test_results.key.TF, aes(x = mean_diff, y = neg_log10_p_val), color = 'red') +
  # ggrepel::geom_text_repel(data = t_test_results.key.TF, 
  #                          aes(label = feature), 
  #                          box.padding = 0.5, max.overlaps = 20) +
  geom_hline(yintercept = p_val_thres, linetype = 'dashed') +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = 'Day 0 (DABTRAM)') +
  theme_bw()

de_res$feature <- rownames(de_res)
comp <- merge(de_res, t_test_results, by = 'feature')
ggplot(comp, aes(x = neg_log10_p_val.x, y = neg_log10_p_val.y)) +
  geom_point()


