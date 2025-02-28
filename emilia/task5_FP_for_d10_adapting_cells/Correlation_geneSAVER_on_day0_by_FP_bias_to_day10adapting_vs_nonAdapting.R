rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

treatment <- 'DABTRAM'

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# keygenes <- list(
#   jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
#                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
#                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
#                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
#                    "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
#   DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
#                    "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
#                    "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
#                    "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
#                    "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44"))
# )

keygenes <- read.csv('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/ISG.RS.txt', header = FALSE)

keygenes <- read.csv('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Persistent IFN ISG Groups/Memory ISGs Human.csv', header = FALSE)
keygenes <- keygenes$V1

# 
# keygenes <- list(COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                                 "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                                 "CAV1")))

# keygenes <- list(CIS = sort(c("YY1AP1", "LGALS3", "MCF2L", "TIMM50", "AC207130.1",
#                               "SLC25A6", "EIF3L", "CTSD", "NQO1", "HNMT", "ZFYVE16",
#                               "PHACTR1", "TNFRSF14", "RAI14", "TRPM1", "HIST1H1C",
#                               "HIST2H2AC", "SPARC", "TRIM63", "TUBA1B", "HIST1H1A",
#                               "HIST1H1D", "PYCARD", "FSTL1", "DCT", "CTSK", "HIST1H4C",
#                               "GDF15", "HIST1H1B")))
keygenes <- unlist(keygenes)
# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

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


nonAdaptingFP <- readRDS(paste0(out_dir, 'final_fit_d0_d10_Nonadpating_thres_0_', treatment, '.rds'))
nonAdaptingFP <- nonAdaptingFP[["cell_imputed_score"]]

adaptingFP <- readRDS(paste0(out_dir, 'final_fit_d0_d10_adpating_thres_0_', treatment, '.rds'))
adaptingFP <- adaptingFP[["cell_imputed_score"]]

# =============================================================================
# Check
# =============================================================================
adaptingFP <- as.data.frame(adaptingFP)
nonAdaptingFP <- as.data.frame(nonAdaptingFP)

df <- merge(adaptingFP, nonAdaptingFP, by = 'row.names')
df$bias <- 10**(df$adaptingFP) / (10**(df$adaptingFP) + 10**(df$nonAdaptingFP))
hist(df$bias, breaks = 50)

ggplot(df, aes(x = nonAdaptingFP, y = bias)) +
  geom_jitter(width = 0.1) +
  # geom_abline(intercept = 0, slope = 1, col = 'red') +
  stat_cor(method = 'spearman') +
  theme_minimal()

colnames(df)[1] <- 'cell_id' 

write.csv(df, paste0(out_dir, 'adapting_bias_thres_0_', treatment, '.csv'), row.names = FALSE)
# =============================================================================
# subset to day0
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')

saver.day0 <- all_data_day0@assays[["saver"]]@scale.data

# =============================================================================
# correlation
# =============================================================================

cor_vec <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cor_vec) <- c('gene', 'cor', 'p_val')
features <- rownames(saver.day0)
for(g in features) {
  gene.exp <- as.data.frame(saver.day0[g, ])
  colnames(gene.exp) <- 'gene'
  gene.exp$cell_id <- rownames(gene.exp)
  
  gene.exp <- merge(gene.exp, df, by = 'cell_id')
  res <- cor.test(gene.exp$bias, gene.exp$gene)  
  rho <- res$estimate
  p_val <- res$p.value
  cor_vec <- rbind(cor_vec, c(g, rho, p_val))
}
colnames(cor_vec) <- c('gene', 'cor', 'p_val')
cor_vec$cor <- as.numeric(cor_vec$cor)
cor_vec$p_val <- as.numeric(cor_vec$p_val)
hist(cor_vec$cor)
cor_vec <- cor_vec %>% arrange(desc(cor))
cor_vec$order <- rownames(cor_vec)
cor_vec$order <- as.numeric(cor_vec$order)
cor_vec$keygene <- ifelse(cor_vec$gene %in% keygenes, 'Key gene', 'Other gene')
cor_vec$p_val_adj <- p.adjust(cor_vec$p_val, method = 'BH')
cor_vec$is_significant <- cor_vec$p_val_adj < 0.05
cor_vec <- cor_vec %>% drop_na()

cor_upper <- cor_vec %>% filter(cor > 0 & is_significant) %>% arrange(desc(cor)) %>% tail(1) %>% pull(cor)
cor_lower <- cor_vec %>% filter(cor < 0 & is_significant) %>% arrange(cor) %>% tail(1) %>% pull(cor)

ggplot(cor_vec, aes(x = order, y = cor)) +
  geom_point() +
  geom_point(data = subset(cor_vec, keygene == 'Key gene'), aes(color = 'key gene')) +
  ggrepel::geom_text_repel(data = subset(cor_vec, keygene == 'Key gene'), aes(label = gene), nudge_y = 0.1) +
  geom_hline(yintercept = cor_upper, linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = cor_lower, linetype = 'dashed', color = 'black') +
  labs(title = paste0('Correlation of gene expression with fate bias (', treatment, ')'),
       x = 'Gene', y = 'Correlation') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor_vec <- cor_vec[, c('gene', 'cor', 'p_val')]
save(date_of_run, session_info,
     cor_vec, 
     file = paste0(out_dir, 'geneSaver_on_day0_cor_vec_', treatment, '.RData'))
