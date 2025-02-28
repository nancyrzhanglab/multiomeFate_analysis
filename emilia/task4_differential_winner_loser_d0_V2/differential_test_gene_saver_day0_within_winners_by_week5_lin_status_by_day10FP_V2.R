library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggthemes)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
remove_unassigned_cells <- TRUE

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

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

treatment <- 'DABTRAM'
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# ==============================================================================
# Get imputed day 10 status per cell
# ==============================================================================
# imputed progeny size day 0 to day10 
fp_name <- paste0('fatepotential_', treatment, '_d0_d10')
imputedCell.d0_d10 <- as.data.frame(all_data_fatepotential[[fp_name]][["cell_imputed_score"]])
colnames(imputedCell.d0_d10) <- paste0('cell_imputed_score_', treatment, '_d0_d10')
imputedCell.d0_d10$cell_id <- rownames(imputedCell.d0_d10)

imputedCell.d0_d10 <- merge(imputedCell.d0_d10, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')


# imputed lineage size day 10 to week5 
fp_name_2 <- paste0('fatepotential_', treatment, '_d10_w5')
imputedLineage.d10_w5 <- as.data.frame(all_data_fatepotential[[fp_name_2]][["lineage_imputed_count"]])
colnames(imputedLineage.d10_w5) <- paste0('lineage_imputed_count_', treatment, '_d10_w5')
imputedLineage.d10_w5$assigned_lineage <- rownames(imputedLineage.d10_w5)


# combine day0 and week5 through day 10
combine.d0.d10.w5 <- merge(imputedCell.d0_d10, imputedLineage.d10_w5, by = 'assigned_lineage', all.x = T)
# combine.d0.d10.w5 <- combine.d0.d10.w5 %>% drop_na()


fp_median <- median(combine.d0.d10.w5[[paste0('cell_imputed_score_', treatment, '_d0_d10')]])
combine.d0.d10.w5$winner_day0 <- ifelse(combine.d0.d10.w5[[paste0('cell_imputed_score_', treatment, '_d0_d10')]] > fp_median, 'Yes', 'No')
combine.d0.d10.w5$winnerLin_week5 <- ifelse(combine.d0.d10.w5[[paste0('lineage_imputed_count_', treatment, '_d10_w5')]] > 1, 'Yes', 'No')
combine.d0.d10.w5$winnerLin_week5 <- ifelse(is.na(combine.d0.d10.w5$winnerLin_week5), 'No', combine.d0.d10.w5$winnerLin_week5)
combine.d0.d10.w5 %>% 
  filter(winner_day0 == 'Yes',
         winnerLin_week5 == 'Yes') %>% 
  nrow()

combine.d0.d10.w5 %>% 
  filter(winner_day0 == 'Yes',
         winnerLin_week5 != 'Yes') %>% 
  nrow()

day0_primed <- combine.d0.d10.w5 %>% 
  filter(winner_day0 == 'Yes',
         winnerLin_week5 == 'Yes') %>% 
  pull('cell_id')

day0_non.primed <- combine.d0.d10.w5 %>% 
  filter(winner_day0 == 'Yes',
         winnerLin_week5 != 'Yes') %>% 
  pull('cell_id')

# ==============================================================================
# Get day0 winner cells
# ==============================================================================

# get saver
metadat <- all_data@meta.data
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
metadat.day0$cell_id <- rownames(metadat.day0)

saver_day0 <- t(all_data[["Saver"]]@data)
saver_day0 <- saver_day0[metadat.day0$cell_id, ]

saver_day0_primed <- saver_day0[day0_primed, ]
saver_day0_non.primed <- saver_day0[day0_non.primed, ]

# ==============================================================================
# Differential tests
# ==============================================================================
columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

features <- colnames(saver_day0)


for (f in features) {
  
  feature_day0_winner_lin_in_week5 <- saver_day0_primed[, f]
  feature_day0_winner_lin_Notin_week5 <- saver_day0_non.primed[, f]
  
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
  geom_vline(xintercept = 0, linetype = 'dashed') +
  coord_cartesian(xlim = c(-0.3, 2)) +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = paste0('Day 0 (', treatment, ')')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))

