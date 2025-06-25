rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)

library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
fatebias_out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'

figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V4/Fig5/'

remove_unassigned_cells <- TRUE

# =============================================================================
# Read data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVAR_day0.RData'))


all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
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
all_data.day0 <- subset(all_data, dataset == 'day0')
all_data.day10_DABTRAM <- subset(all_data, dataset == 'day10_DABTRAM')

lin.category.df <- read.csv('~/Downloads/lin_category_DABTRAM.csv')
# =============================================================================
# Classify priming and plasticity at day 10
# =============================================================================
metadat.day10_DABTRAM <- all_data.day10_DABTRAM@meta.data

fp.summary <- metadat.day10_DABTRAM %>% 
  group_by(assigned_lineage) %>%
  summarise(
    mean_fatepotential = mean(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    sd_fatepotential = sd(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    n_cells = n()
  )

fp.summary.large <- fp.summary[fp.summary$n_cells > 1, ]
ggplot(fp.summary.large, aes(x = mean_fatepotential, y = sd_fatepotential)) +
  geom_point() +
  theme_bw()

fp.summary$category <- ifelse((fp.summary$mean_fatepotential > 0) & 
                              (fp.summary$sd_fatepotential < median(fp.summary$sd_fatepotential, na.rm = T)), 
                              'priming', NA)
fp.summary$category <- ifelse((fp.summary$mean_fatepotential > 0.5) & 
                              (fp.summary$n_cells == 1),
                              'priming', fp.summary$category)

fp.summary$category <- ifelse((fp.summary$sd_fatepotential > median(fp.summary$sd_fatepotential, na.rm = T)) &
                              (fp.summary$mean_fatepotential > -1) &
                              # (fp.summary$n_cells > 5) &
                              (fp.summary$n_cells < 50) &
                              is.na(fp.summary$category),
                              'plastic', fp.summary$category )

# lin.priming <- fp.summary %>% 
#   filter(category == 'priming') %>%
#   pull(assigned_lineage)
# lin.plastic <- fp.summary %>%
#   filter(category == 'plastic') %>%
#   pull(assigned_lineage)

lin.priming <- lin.category.df %>%
  filter(category == 'priming') %>%
  pull(assigned_lineage)
lin.plastic <- lin.category.df %>%
  filter(category == 'plastic') %>%
  pull(assigned_lineage)


saver.mat <- all_data.day0[["Saver"]]@data
saver.mat <- t(saver.mat)

metadat.day0 <- all_data.day0@meta.data
priming_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.priming]
plastic_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.plastic]

p <- ncol(saver.mat)

wilcox_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = saver.mat[priming_progenitors,j],
    y = saver.mat[plastic_progenitors,j]
  )
  logfc <- log2(mean(saver.mat[priming_progenitors,j])) - log2(mean(saver.mat[plastic_progenitors,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results <- t(wilcox_results)

gene_df <- as.data.frame(wilcox_results)
rownames(gene_df) <- colnames(saver.mat)
colnames(gene_df) <- c("logfc", "p.value")

gene_df$padj <- p.adjust(gene_df$p.value, method = "BH")
gene_df$neglog10_pval <- -log10(gene_df$p.value)

ggplot(gene_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = padj < 0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change (Adapting vs Non-Adapting)",
       y = "-log10(p-value)",
       title = "Differential Expression of Adapting Progenitors") +
  theme_minimal() +
  theme(legend.position = "none")

gene_df$gene <- rownames(gene_df)


# chromVar
cv.mat <- all_data.day0[["chromVar_day0"]]@data
cv.mat <- t(cv.mat)

metadat.day0 <- all_data.day0@meta.data
priming_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.priming]
plastic_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.plastic]

p <- ncol(cv.mat)

wilcox_results_chromVar <- sapply(1:p, function(j){
  # print(j)
  var_j = var(cv.mat[c(priming_progenitors, plastic_progenitors), j])
  if(var_j == 0 | is.na(var_j)) {
    print(paste("Skipping feature", j, "due to zero variance."))
    return(c(logfc = NA, p.value = NA))
  }
  
  tmp <- stats::wilcox.test(
    x = cv.mat[priming_progenitors,j],
    y = cv.mat[plastic_progenitors,j]
  )
  logfc <- mean(cv.mat[priming_progenitors,j]) - mean(cv.mat[plastic_progenitors,j])
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})


wilcox_results_chromVar <- t(wilcox_results_chromVar)

tf_df <- as.data.frame(wilcox_results_chromVar)
rownames(tf_df) <- colnames(cv.mat)
colnames(tf_df) <- c("logfc", "p.value")

tf_df$padj <- p.adjust(tf_df$p.value, method = "BH")
tf_df$neglog10_pval <- -log10(tf_df$p.value)
tf_df$tf <- rownames(tf_df)

ggplot(tf_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = padj < 0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change (Adapting vs Non-Adapting)",
       y = "-log10(p-value)",
       title = "Differential Expression of Adapting Progenitors") +
  theme_minimal() +
  theme(legend.position = "none")


gavish.mp <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/saver.GavishMP.UCellScores.csv')
rownames(gavish.mp) <- gavish.mp$cell_id

metaprograms <- colnames(gavish.mp)
metaprograms <- metaprograms[!metaprograms %in% c("cell_id", 
                                                  "fatepotential_DABTRAM_d0_d10", "fatepotential_DABTRAM_d10_w5",
                                                  "fatepotential_COCL2_d0_d10", "fatepotential_COCL2_d10_w5",
                                                  "fatepotential_CIS_d0_d10", "fatepotential_CIS_d10_w5",
                                                  "ISG.Mem")]

metaprograms <- metaprograms[!metaprograms %in% c("MES..glioma.", "Cilia", "Astrocytes", "NPC.Glioma", "Oligo.normal",
                                                  "Oligo.Progenitor", "NPC.OPC", "PDAC.classical",
                                                  "Alveolar", "RBCs", "Platelet.activation", "Hemato.related.I",
                                                  "IG", "Hemato.related.II", "Glutathione", "PDAC.related", "Unassigned", "AXL_Program", "MITF_Program")]

wilcox_results_gavish_mp <- sapply(metaprograms, function(j){
  tmp <- stats::wilcox.test(
    x = gavish.mp[priming_progenitors,j],
    y = gavish.mp[plastic_progenitors,j]
  )
  logfc <- log2(mean(gavish.mp[priming_progenitors,j])) - log2(mean(gavish.mp[plastic_progenitors,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results_gavish_mp <- t(wilcox_results_gavish_mp)
wilcox_results_gavish_mp_df <- as.data.frame(wilcox_results_gavish_mp)
rownames(wilcox_results_gavish_mp_df) <- metaprograms
colnames(wilcox_results_gavish_mp_df) <- c("logfc", "p.value")

wilcox_results_gavish_mp_df$neglog10_pval <- -log10(wilcox_results_gavish_mp_df$p.value)
wilcox_results_gavish_mp_df$logfc <- ifelse(wilcox_results_gavish_mp_df$logfc < -0.1, -0.1, wilcox_results_gavish_mp_df$logfc)

ggplot(wilcox_results_gavish_mp_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = p.value < 0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change (Plastic <-> Priming)",
       y = "-log10(p-value)",
       title = "Differential Expression of Metaprograms") +
  theme_bw() +
  theme(legend.position = "none")

wilcox_results_gavish_mp_df$metaprogram <- rownames(wilcox_results_gavish_mp_df)
wilcox_results_gavish_mp_df <- wilcox_results_gavish_mp_df[order(wilcox_results_gavish_mp_df$logfc),]
wilcox_results_gavish_mp_df <- wilcox_results_gavish_mp_df[!is.na(wilcox_results_gavish_mp_df$logfc),]

rownames(wilcox_results_gavish_mp_df) <- NULL
ggplot(wilcox_results_gavish_mp_df, aes(y = reorder(metaprogram, logfc),  x = logfc)) +
  geom_bar(stat = "identity", aes(fill = p.value < 0.05)) +
  scale_fill_manual(values = c("grey", "red")) +
  ylab('') +
  xlab('Log2 Fold Change (Plastic <-> Priming)') +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

fp_d0_d10_dabtram <- all_data_fatepotential[["fatepotential_DABTRAM_d0_d10"]][["cell_imputed_score"]]
fp_d0_d10_dabtram <- as.data.frame(fp_d0_d10_dabtram)
fp_d0_d10_dabtram$cell_id <- rownames(fp_d0_d10_dabtram)

metadat.day0$cell_id <- rownames(metadat.day0)
fp_d0_d10_dabtram <- merge(fp_d0_d10_dabtram, metadat.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

fp_d0_d10_dabtram$category <- ifelse(fp_d0_d10_dabtram$assigned_lineage %in% lin.priming, 'Priming', NA)
fp_d0_d10_dabtram$category <- ifelse(fp_d0_d10_dabtram$assigned_lineage %in% lin.plastic, 'Plastic', fp_d0_d10_dabtram$category)

ggplot(fp_d0_d10_dabtram, aes(x = category, y = fp_d0_d10_dabtram)) +
  geom_violin(scale = 'width') +
  geom_jitter(width = 0.1) +
  geom_boxplot(width = 0.1, fill = 'white', alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("Priming" = "blue", "Plastic" = "red")) +
  labs(x = "Cateogry",
       y = "Fate Potential D0 to D10",
       title = "Fate Potential Distribution") +
  theme_minimal() +
  theme(legend.position = "top")

fp_d10_w5_dabtram <- all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]]
fp_d10_w5_dabtram <- as.data.frame(fp_d10_w5_dabtram)
fp_d10_w5_dabtram$cell_id <- rownames(fp_d10_w5_dabtram)

metadat.day10_DABTRAM$cell_id <- rownames(metadat.day10_DABTRAM)
fp_d10_w5_dabtram <- merge(fp_d10_w5_dabtram, metadat.day10_DABTRAM[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

fp_d10_w5_dabtram$category <- ifelse(fp_d10_w5_dabtram$assigned_lineage %in% lin.priming, 'Priming', NA)
fp_d10_w5_dabtram$category <- ifelse(fp_d10_w5_dabtram$assigned_lineage %in% lin.plastic, 'Plastic', fp_d10_w5_dabtram$category)

ggplot(fp_d10_w5_dabtram, aes(x = category, y = fp_d10_w5_dabtram)) +
  geom_violin(scale = 'width') +
  geom_jitter(width = 0.1) +
  geom_boxplot(width = 0.1, fill = 'white', alpha = 0.7, outlier.shape = NA) +
  labs(x = "Cateogry",
       y = "Fate Potential D10 to W5",
       title = "Fate Potential Distribution") +
  theme_minimal() +
  theme(legend.position = "top")
