library(ggplot2)
library(tidyverse)
library(dplyr)

# ==============================================================================
# Read data
# ==============================================================================

data_dir <- "~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/"
treatment_day10 <- 'day10_CIS'
load(paste0(data_dir, 'chromVar_day0_data.RData'))
load(paste0(data_dir, treatment_day10, "/Writeup6n_CIS_day10_lineage-imputation_stepdown_concise-postprocessed.RData"))
metadat_day0 <- read.csv(paste0(data_dir, 'day0/day0_meta.csv'), row.names = 1)
metadat_day10 <- read.csv(paste0(data_dir, treatment_day10, '/', treatment_day10, '_meta.csv'), row.names = 1)

cell_imputed_count <- as.data.frame(cell_imputed_count[!is.na(cell_imputed_count)])
colnames(cell_imputed_count) <- c('Growth_Potential')
cell_imputed_count$cell_id <- row.names(cell_imputed_count)

lineage_imputed_count <- as.data.frame(lineage_imputed_count)
lineage_imputed_count$assigned_lineage <- row.names(lineage_imputed_count)

motifs <- read.csv('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/motif_info.csv')

chromvar_results_day0 <- as.data.frame(chromvar_results_day0)
chromvar_results_day0$motif_code <- rownames(chromvar_results_day0)
chromvar_results_day0 <- merge(chromvar_results_day0, motifs, by = 'motif_code')
rownames(chromvar_results_day0) <- chromvar_results_day0$motif_names
chromvar_results_day0 <- subset(chromvar_results_day0, select= -c(motif_code, motif_names))
chromvar_results_day0 <- as.data.frame(t(chromvar_results_day0))
chromvar_results_day0$cell_id <- row.names(chromvar_results_day0)

adaptation_motifs <- read.csv('~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_DABTRAM_adaptation_motifs.csv')

# ==============================================================================
# Get lineage barcode info
# ==============================================================================
metadat_day10$cell_id <- row.names(metadat_day10)
metadat_day10 <- metadat_day10[, c('cell_id', 'assigned_lineage')]
cell_imputed_count <- merge(cell_imputed_count, metadat_day10, by='cell_id')
cell_imputed_count_summary <- cell_imputed_count %>% 
  group_by(assigned_lineage) %>% 
  summarise(mean_GP = mean(Growth_Potential))

metadat_day0$cell_id <- row.names(metadat_day0)
metadat_day0 <- metadat_day0[, c('cell_id', 'assigned_lineage')]

chromvar_results_day0 <- merge(chromvar_results_day0, metadat_day0, by = 'cell_id')

# ==============================================================================
# Calculate mean
# ==============================================================================
chromvar_results_day0 <- subset(chromvar_results_day0, select = -c(cell_id))
mean_chromvar <- chromvar_results_day0 %>% 
  group_by(assigned_lineage) %>%
  summarise(across(everything(), mean))

mean_chromvar <- subset(mean_chromvar, select = c('assigned_lineage', adaptation_motifs$motif_names))

# ==============================================================================
# Plotting
# ==============================================================================
# results <- merge(mean_chromvar, cell_imputed_count_summary, by = 'assigned_lineage')
results <- merge(mean_chromvar, lineage_imputed_count, by = 'assigned_lineage')
results$log10_lineage_imputed_count <- log10(results$lineage_imputed_count)

# results_sm <- results[, c('assigned_lineage', 'FOS', 'mean_GP')]
results_sm <- results[, c('assigned_lineage', 'FOSB::JUNB', 'lineage_imputed_count')]
results_sm$log10_lineage_imputed_count <- log10(results_sm$lineage_imputed_count)

colnames(results_sm) <- c("assigned_lineage", "FOSB_JUNB", "lineage_imputed_count", "log10_lineage_imputed_count")
ggplot(results_sm,aes(x = log10_lineage_imputed_count, y = FOSB_JUNB)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()

# ==============================================================================
# Correlate with GP from day10 to week5
# ==============================================================================
cor_vec <- sapply(seq(2, ncol(results)-2), function(j){
  res <- stats::cor.test(results$log10_lineage_imputed_count, results[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cor_vec <- as.data.frame(t(cor_vec))
colnames(cor_vec) <- c("correlation", "p.value")
rownames(cor_vec) <- colnames(results)[seq(2, ncol(results)-2)]
cor_vec$p_val_adjust <- p.adjust(cor_vec$p.value, 'BH')
cor_vec$feature <- rownames(cor_vec)

# ==============================================================================
# Plotting
# ==============================================================================
cor_vec <- cor_vec[order(-cor_vec$correlation),]
cor_vec$isSig <- ifelse(cor_vec$p_val_adjust < 0.05, 'Yes', 'No')
cor_vec$order <- seq(1: nrow(cor_vec))
features_to_label <- head(cor_vec, 5)$feature
features_to_label <- c(features_to_label, tail(cor_vec, 5)$feature)


ggplot(cor_vec, aes(x = order, y = correlation, color=isSig)) +
  geom_point() +
  scale_color_manual(values=c('gray', 'red')) + 
  ggrepel::geom_text_repel(data = subset(cor_vec, feature %in% features_to_label),
                           ggplot2::aes(label = feature),
                           box.padding = ggplot2::unit(1.2, 'lines'),
                           point.padding = ggplot2::unit(0.1, 'lines'),
                           color = 'blue',
                           max.overlaps = 80) +
  ylab("Correlation between mean(chromVAR) in day0\nand their lineage's GP at d10 for wk5") +
  theme_bw()


