library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
early_t <- 'day0'
late_t <- 'week5_CIS'
early_tp <- readRDS(paste0(in_dir, 'Raw_and_Processed/', early_t, '_processed.rds'))
early_metadat <- early_tp@meta.data
late_tp <- readRDS(paste0(in_dir, 'Raw_and_Processed/',late_t, '_processed.rds'))
late_metadat <- late_tp@meta.data


# ==============================================================================
# Wrangling
# ==============================================================================

late_lineage_size <- late_metadat %>% 
  group_by(assigned_lineage) %>% 
  summarise(n_cells_late = n())
# week5_lineage_size$isWinner <- ifelse(week5_lineage_size$n_cells_wk5 > 20, 'Winner', 'Non-winner')

early_metadat$nFeature_ATAC_norm <- early_metadat$nFeature_ATAC / early_metadat$nCount_ATAC
early_metadat$nFeature_RNA_norm <- early_metadat$nFeature_RNA / early_metadat$nCount_RNA

early_metadat_lineage_feature <- early_metadat %>% 
  group_by(assigned_lineage) %>% 
  summarise(nFeature_ATAC_mean = mean(nFeature_ATAC),
            nFeature_ATAC_norm_mean = mean(nFeature_ATAC_norm),
            nFeature_RNA_mean = mean(nFeature_RNA),
            nFeature_RNA_norm_mean = mean(nFeature_RNA_norm))

early_metadat_lineage_feature <- merge(early_metadat_lineage_feature, late_lineage_size, by = 'assigned_lineage', all = TRUE)
early_metadat_lineage_feature[is.na(early_metadat_lineage_feature)] <- 0
early_metadat_lineage_feature$log10_nCells_late <- log10(early_metadat_lineage_feature$n_cells_late + 1)

early_metadat <- merge(early_metadat, late_lineage_size, by = 'assigned_lineage', all.x = TRUE)
early_metadat[is.na(early_metadat)] <- 0
early_metadat$log10_nCells_late <- log10(early_metadat$n_cells_late + 1)
early_metadat$isWinner <- ifelse(early_metadat$n_cells_late > 20, 'Winner', 'Non-winner')

# ==============================================================================
# Plotting
# ==============================================================================


res <- cor.test(early_metadat_lineage_feature$nFeature_RNA_mean, y = early_metadat_lineage_feature$log10_nCells_late, method = 'spearman')
rho <- res[["estimate"]][["rho"]]
p_val <- res[["p.value"]]
ggplot(early_metadat_lineage_feature, aes(x = nFeature_ATAC_mean, y = log10_nCells_late)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  ggtitle(paste0('r = ', rho, '; p = ', p_val))


p1 <- ggplot(early_metadat, aes(x = isWinner, y = nFeature_RNA)) +
  geom_jitter(alpha=0.5, width = 0.1) +
  geom_boxplot(alpha=0.9, width = 0.4, outlier.shape = NA) +
  theme_bw()
p2 <- ggplot(early_metadat, aes(x = isWinner, y = nCount_RNA)) +
  geom_jitter(alpha=0.5, width = 0.1) +
  geom_boxplot(alpha=0.9, width = 0.4, outlier.shape = NA) +
  theme_bw()
p3 <- ggplot(early_metadat, aes(x = isWinner, y = nFeature_RNA_norm)) +
  geom_jitter(alpha=0.5, width = 0.1) +
  geom_boxplot(alpha=0.9, width = 0.4, outlier.shape = NA) +
  theme_bw()
p4 <- ggplot(early_metadat, aes(x = isWinner, y = nFeature_ATAC)) +
  geom_jitter(alpha=0.5, width = 0.1) +
  geom_boxplot(alpha=0.9, width = 0.4, outlier.shape = NA) +
  theme_bw()
p5 <- ggplot(early_metadat, aes(x = isWinner, y = nCount_ATAC)) +
  geom_jitter(alpha=0.5, width = 0.1) +
  geom_boxplot(alpha=0.9, width = 0.4, outlier.shape = NA) +
  theme_bw()
p6 <- ggplot(early_metadat, aes(x = isWinner, y = nFeature_ATAC_norm)) +
  geom_jitter(alpha=0.5, width = 0.1) +
  geom_boxplot(alpha=0.9, width = 0.4, outlier.shape = NA) +
  theme_bw()
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

features <- c('nFeature_RNA', 'nCount_RNA', 'nFeature_RNA_norm', 'nFeature_ATAC', 'nCount_ATAC', 'nFeature_ATAC_norm')
for (f in features) {
  winner <- early_metadat[early_metadat$isWinner == 'Winner', ][[f]]
  loser <- early_metadat[early_metadat$isWinner == 'Non-winner', ][[f]]
  res <- wilcox.test(winner, loser) 
  p_val <- res[["p.value"]]
  print(paste0(f, ' ', p_val))
}


cor.test(day10_metadat$nCount_RNA_log, day10_metadat$log10_nCells_wk5, method = 'spearman')


par(mfrow = c(2, 2))
hist(early_tp$nCount_RNA, breaks=50, main = paste0(early_t, ' nCount_RNA'))
abline(v = median(early_tp$nCount_RNA), col = 'red', lwd = 3) 

hist(early_tp$nFeature_RNA, breaks=50,  main = paste0(early_t, ' nFeature_RNA'))
abline(v = median(early_tp$nFeature_RNA), col = 'red', lwd = 3) 

hist(early_tp$nCount_ATAC, breaks=50, main = paste0(early_t, ' nCount_ATAC'))
abline(v = median(early_tp$nCount_ATAC), col = 'red', lwd = 3) 

hist(early_tp$nFeature_ATAC, breaks=50,  main = paste0(early_t, ' nFeature_ATAC'))
abline(v = median(early_tp$nFeature_ATAC), col = 'red', lwd = 3) 


hist(log10(late_lineage_size$n_cells_late), breaks=20, main = paste0('week5_CIS', ' lineage_size'))


# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/day0_processed.rds'))
day0_meta <- sc@meta.data

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/day10_DABTRAM_processed.rds'))
day10_DABTRAM_meta <- sc@meta.data

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/day10_COCL2_processed.rds'))
day10_COCL2_meta <- sc@meta.data

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/day10_CIS_processed.rds'))
day10_CIS_meta <- sc@meta.data

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/week5_DABTRAM_processed.rds'))
week5_DABTRAM_meta <- sc@meta.data

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/week5_COCL2_processed.rds'))
week5_COCL2_meta <- sc@meta.data

sc <- readRDS(paste0(in_dir, 'Raw_and_Processed/week5_CIS_processed.rds'))
week5_CIS_meta <- sc@meta.data

full_meta <- rbind(day0_meta, day10_DABTRAM_meta)
full_meta <- rbind(full_meta, day10_COCL2_meta)
full_meta <- rbind(full_meta, day10_CIS_meta)
full_meta <- rbind(full_meta, week5_DABTRAM_meta)
full_meta <- rbind(full_meta, week5_COCL2_meta)
full_meta <- rbind(full_meta, week5_CIS_meta)

ggplot(full_meta, aes(x = dataset, y = nCount_RNA)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  scale_y_continuous(breaks=seq(0,max(full_meta$nCount_RNA),by=10000)) +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  coord_cartesian(ylim=c(0, 50000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.3))

ggplot(full_meta, aes(x = dataset, y = nCount_ATAC)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  scale_y_continuous(breaks=seq(0,max(full_meta$nCount_ATAC),by=10000)) +
  coord_cartesian(ylim=c(0, 50000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=0.3))

