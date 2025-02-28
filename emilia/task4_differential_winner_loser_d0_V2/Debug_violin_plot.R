rm(list=ls())
library(Seurat)
library(Signac)
library(tidyverse)


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
# Define winner and loser
# =============================================================================
def = 'linaege'
# Extract the fate potential scores for each cell
fp_name <- paste0('fatepotential_', treatment, '_d10_w5')
fp_name2 <- paste0('fatepotential_', treatment, '_d0_d10')

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
  lineage_score_all$winner <- ifelse(lineage_score_all$mean_score <= -0.8440993, 'Losing', lineage_score_all$winner)
  
} else if (def == 'lineage') {
  fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['lineage_imputed_count']])
  colnames(fp_d10_w5_treatment) <- c('lineage_imputed_count')
  fp_d10_w5_treatment$assigned_lineage <- rownames(fp_d10_w5_treatment)
  
  lineage_score_all <- fp_d10_w5_treatment
  lineage_score_all$winner <- ifelse(lineage_score_all$lineage_imputed_count > 1, 'Winning', NA)
  lineage_score_all$winner <- ifelse(lineage_score_all$lineage_imputed_count < 0.5, 'Losing', lineage_score_all$winner)
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
# t test
# =============================================================================
columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns


saver_day0 <- all_data_day0@assays$saver@scale.data
saver_day0 <- t(saver_day0)
# features <- anova_saver_dabtram_d10_sig$feature

metadat.day0$cell_id <- rownames(metadat.day0)
day0_winner_lin_in_week5 <- metadat.day0 %>% filter(isWinningLineage == 'Winning') %>% pull(cell_id)
day0_winner_lin_Notin_week5 <- metadat.day0 %>% filter(isWinningLineage == 'Losing') %>% pull(cell_id)

gene <- 'FN1'

feature_day0_winner_lin_in_week5 <- saver_day0[day0_winner_lin_in_week5, gene]
feature_day0_winner_lin_Notin_week5 <- saver_day0[day0_winner_lin_Notin_week5, gene]

winner <- as.data.frame(feature_day0_winner_lin_in_week5)
loser <- as.data.frame(feature_day0_winner_lin_Notin_week5)

colnames(winner) <- gene
colnames(loser) <- gene

winner$category <- 'winner_day0'
loser$category <- 'loser_day0'
comp_df <- rbind(winner, loser)
comp_df$cell_id <- rownames(comp_df)

fp_d0_d10_treatment <- as.data.frame(all_data_day10@misc[[fp_name2]][['cell_imputed_score']])
colnames(fp_d0_d10_treatment) <- c('cell_imputed_score_d0_d10')
fp_d0_d10_treatment$cell_id <- rownames(fp_d0_d10_treatment)

comp_df <- merge(comp_df, fp_d0_d10_treatment, by = 'cell_id')

t.test(feature_day0_winner_lin_in_week5, feature_day0_winner_lin_Notin_week5)

comp_df <- comp_df[comp_df$FN1 < 1, ]
ggplot(comp_df, aes(x = category, y = FN1)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim = quantile(comp_df$FN1, c(0.05, 0.95))) +
  theme_bw() +
  labs(title = 'CAV1 expression in day0', x = 'Lineage', y = 'CAV1 expression')



ggplot(comp_df, aes(x = cell_imputed_score_d0_d10, y = FN1)) +
  geom_point() +
  facet_wrap(~category) +
  theme_bw()

comp_df %>% 
  filter(category == 'loser_day0') %>%
  pull(cell_imputed_score_d0_d10) %>% 
  hist(breaks = 50)
