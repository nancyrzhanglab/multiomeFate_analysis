rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

keyTFs <- c('JUN', 'FOS::JUN', 'FOSB::JUN', 'FOSL1::JUN', 'FOSL2::JUN', 'JUN::JUNB', 
            'FOS::JUNB', 'FOSB::JUNB', 'FOSL1::JUNB', 'FOSL2::JUNB', 'FOS::JUND',
            'FOSL1::JUND', 'FOSL2::JUND', 'JUNB', 'JUND',
            'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 
            'SOX10', 'MITF', 
            'SNAI1', 'SNAI2', 'SNAI3')

treatment <- 'COCL2'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))

all_data[['chromVAR_day0']] <- all_data_chromVar_day0

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}



# final_fit <- readRDS('~/Downloads/final_fit_d0_w5_COCL2_1.rds')
# fp_adapting <- as.data.frame(final_fit[["cell_imputed_score"]])
fp_adapting <- readRDS('~/Downloads/cell_imputed_score2_d0_d10_predicted_w5_w_d10_size.COCL2_1.rds')
fp_adapting <- as.data.frame(fp_adapting)
colnames(fp_adapting) <- 'FP_adapting'
fp_adapting$cell_id <- rownames(fp_adapting)

# =============================================================================
# subset to day0
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')

chromVAR.day0 <- all_data_day0@assays[["chromVAR_day0"]]@data

cor_vec <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cor_vec) <- c('TF', 'cor', 'p_val')
features <- rownames(chromVAR.day0)

for(tf in features) {
  tf.act <- as.data.frame(chromVAR.day0[tf, ])
  colnames(tf.act) <- 'TF'
  tf.act$cell_id <- rownames(tf.act)
  
  tf.act <- merge(tf.act, fp_adapting, by = 'cell_id')
  res <- cor.test(tf.act$FP_adapting, tf.act$TF)  
  rho <- res$estimate
  p_val <- res$p.value
  cor_vec <- rbind(cor_vec, c(tf, rho, p_val))
}

colnames(cor_vec) <- c('TF', 'cor', 'p_val')
cor_vec$cor <- as.numeric(cor_vec$cor)
cor_vec$p_val <- as.numeric(cor_vec$p_val)
hist(cor_vec$cor)
cor_vec <- cor_vec %>% arrange(desc(cor))
cor_vec$order <- rownames(cor_vec)
cor_vec$order <- as.numeric(cor_vec$order)
cor_vec$keyTF <- ifelse(cor_vec$TF %in% keyTFs, 'Key TF', 'Other TF')
cor_vec$p_val_adj <- p.adjust(cor_vec$p_val, method = 'BH')
cor_vec$is_significant <- cor_vec$p_val_adj < 0.05
cor_vec <- cor_vec %>% drop_na()

cor_upper <- cor_vec %>% filter(cor > 0 & is_significant) %>% arrange(desc(cor)) %>% tail(1) %>% pull(cor)
cor_lower <- cor_vec %>% filter(cor < 0 & is_significant) %>% arrange(cor) %>% tail(1) %>% pull(cor)

ggplot(cor_vec, aes(x = order, y = cor)) +
  geom_point() +
  geom_point(data = subset(cor_vec, keyTF == 'Key TF'), aes(color = 'key TF')) +
  ggrepel::geom_text_repel(data = subset(cor_vec, keyTF == 'Key TF'), aes(label = TF), nudge_y = 0.05) +
  geom_hline(yintercept = cor_upper, linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = cor_lower, linetype = 'dashed', color = 'black') +
  labs(title = paste0('Correlation of TF chromVAR with fate potential (', treatment, ')'),
       x = 'TF', y = 'Correlation') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



write.csv(cor_vec, paste0('~/Downloads/cor_vec_TF_', treatment, '_1.csv'), row.names = FALSE)

# =============================================================================
# Compare with day0 cor with d0 to d10 FP
# =============================================================================
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
load(paste0(result_dir, 'chromVAR_cor_vec.RData'))
d0_chromVAR_cor_vec <- chromVAR_cor_vec[[paste0(tolower(treatment), '_d0_chromVAR_cor_vec')]]
d0_chromVAR_cor_vec$TF <- rownames(d0_chromVAR_cor_vec)
colnames(d0_chromVAR_cor_vec) <- c('correlation.d0_d10_FP', 'p_value.d0_d10_FP', 'TF')

d10_chromVAR_cor_vec <- chromVAR_cor_vec[[paste0(tolower(treatment), '_d10_chromVAR_cor_vec')]]
d10_chromVAR_cor_vec$TF <- rownames(d10_chromVAR_cor_vec)
colnames(d10_chromVAR_cor_vec) <- c('correlation.d10_w5_FP', 'p_value.d10_w5_FP', 'TF')

comp_df <- merge(d0_chromVAR_cor_vec, cor_vec[, c('TF', 'cor')],  by = 'TF')
comp_df <- merge(d10_chromVAR_cor_vec, comp_df,  by = 'TF')
comp_df$keyTF <- ifelse(comp_df$TF %in% keyTFs, 'Key TF', 'Other TF')

ggplot(comp_df, aes(x = correlation.d0_d10_FP, y = cor)) +
  geom_point() +
  geom_point(data = subset(comp_df, keyTF == 'Key TF'), aes(color = 'key TF')) +
  ggrepel::geom_text_repel(data = subset(comp_df, keyTF == 'Key TF'), aes(label = TF)) +
  stat_cor() +
  # geom_smooth(method = 'lm') +
  theme_minimal()


comp_df.dis <- comp_df[which(comp_df$correlation.d10_w5_FP < 0 & comp_df$cor > 0), ]

