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
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_wnn.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[["wnn.umap"]] <- all_data_wnn

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
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_genes_saver.RData')

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
metadat.day0_winner <- metadat.day0[fp_winner$cell_id, ]

saver_day0 <- t(all_data[["Saver"]]@data)
saver_day0 <- saver_day0[metadat.day0$cell_id, ]

saver_day0_winner <- saver_day0[fp_winner$cell_id, ]

# ==============================================================================
# Get lineages survived to Week5
# ==============================================================================
# get predicted week5 lineage status from day10 cells
fp_name_d10_w5 = paste0("fatepotential_", treatment, '_', fut_time, '_w5')
imputedLinSize.day10ToWeek5 <- as.data.frame(all_data_fatepotential[[fp_name_d10_w5]][["lineage_imputed_count"]])
colnames(imputedLinSize.day10ToWeek5) <- 'lineage_imputed_count'
imputedLinSize.day10ToWeek5$assigned_lineage <- rownames(imputedLinSize.day10ToWeek5)

# get estimated to exist in week5
lin.in.week5 <- imputedLinSize.day10ToWeek5 %>% 
  filter(lineage_imputed_count >= 1) %>% 
  select(assigned_lineage)
lin.not.in.week5 <- imputedLinSize.day10ToWeek5 %>% 
  filter(lineage_imputed_count < 1) %>% 
  select(assigned_lineage)

# get lineages where we have day0 data
lin.in.week5 <- unique(intersect(lin.in.week5$assigned_lineage, metadat.day0_winner$assigned_lineage))
lin.not.in.week5 <- unique(intersect(lin.not.in.week5$assigned_lineage, metadat.day0_winner$assigned_lineage))

metadat.day0_winner$is.lin.in.week5 <- ifelse(metadat.day0_winner$assigned_lineage %in% lin.in.week5, 'lin.in.week5', NA)
metadat.day0_winner$is.lin.in.week5 <- ifelse(!metadat.day0_winner$assigned_lineage %in% lin.in.week5, 'lin.not.in.week5', metadat.day0$is.lin.in.week5)
# metadat.day0 <- metadat.day0 %>% drop_na(is.lin.in.week5)

# ==============================================================================
# Differential tests
# ==============================================================================
columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

features <- colnames(saver_day0_winner)

day0_winner_lin_in_week5 <- metadat.day0_winner %>% 
  filter(is.lin.in.week5 == 'lin.in.week5') %>% 
  select(cell_id)

# day0_winner_lin_Notin_week5 <- metadat.day0 %>% 
#   filter(is.lin.in.week5 == 'lin.not.in.week5') %>% 
#   filter(cell_id %in% rownames(saver_day0_winner)) %>% 
#   select(cell_id)

day0_winner_lin_Notin_week5 <- metadat.day0_winner %>% 
  filter(is.lin.in.week5 != 'lin.in.week5') %>% 
  filter(cell_id %in% rownames(saver_day0_winner)) %>% 
  select(cell_id)

for (f in features) {
  
  feature_day0_winner_lin_in_week5 <- saver_day0_winner[day0_winner_lin_in_week5$cell_id, f]
  feature_day0_winner_lin_Notin_week5 <- saver_day0_winner[day0_winner_lin_Notin_week5$cell_id, f]
  
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
t_test_results$mean_diff <- log2(t_test_results$mean_winner_lin_in_week5) - log2(t_test_results$mean_winner_lin_notin_week5)
t_test_results$t_statistic <- as.numeric(t_test_results$t_statistic)


t_test_results.key.gene <- t_test_results[t_test_results$feature %in% keygenes, ]
p_val_thres <- t_test_results[t_test_results$p_val_adj < 0.05, ]
p_val_thres <- min(p_val_thres$neg_log10_p_val)

ggplot(t_test_results, aes(x = mean_diff, 
                           y = neg_log10_p_val)) +
  geom_point(color = 'gray') +
  geom_point(data = t_test_results.key.gene, aes(x = mean_diff, y = neg_log10_p_val), color = 'red') +
  ggrepel::geom_text_repel(data = t_test_results.key.gene, aes(label = feature), box.padding = 0.5) +
  geom_hline(yintercept = p_val_thres, linetype = 'dashed') +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = paste0('Day 0 (', treatment, ')')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))
