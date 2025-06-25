rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)

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
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_wnn.RData'))
load(paste0(data_dir, 'Writeup10a_data_peakVI_DABTRAM.RData'))


all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[['chromVar_day0']] <- all_data_chromVar_day0
all_data[['chromVar_day10_DABTRAM']] <- all_data_chromVar_day10_DABTRAM
all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap
all_data@reductions[['wnn']] <- all_data_wnn
all_data[["peakVI.DABTRAM"]] <- all_data_peakVI_DABTRAM
all_data[["pVI.DABTRAM.umap"]] <- all_data_pVI_DABTRAM_umap

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

adaptation_bias_d0 <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/adapting_bias_thres_0_DABTRAM.csv')
adapting_progenitor <- adaptation_bias_d0[adaptation_bias_d0$bias > median(adaptation_bias_d0$bias, na.rm = TRUE), 'cell_id']
non_adapting_progenitor <- adaptation_bias_d0[adaptation_bias_d0$bias <= median(adaptation_bias_d0$bias, na.rm = TRUE), 'cell_id']

# umap
ft.umap.DABTRAM <- all_data@reductions[[paste0("ft.DABTRAM.umap")]]@cell.embeddings
ft.umap.DABTRAM <- as.data.frame(ft.umap.DABTRAM)
ft.umap.DABTRAM$cell_id <- rownames(ft.umap.DABTRAM)

#wnn
wnn.umap <- all_data@reductions[['wnn']]@cell.embeddings
wnn.umap <- as.data.frame(wnn.umap)
wnn.umap$cell_id <- rownames(wnn.umap)

# peakVI
pVI.umap.DABTRAM <- all_data@reductions[['pVI.DABTRAM.umap']]@cell.embeddings
pVI.umap.DABTRAM <- as.data.frame(pVI.umap.DABTRAM)
pVI.umap.DABTRAM$cell_id <- rownames(pVI.umap.DABTRAM)

# =============================================================================
# Classify priming and plasticity at day 10
# =============================================================================
metadat.day0 <- all_data.day0@meta.data
metadat.day10_DABTRAM <- all_data.day10_DABTRAM@meta.data

linsize.d0 <- metadat.day0 %>% 
  group_by(assigned_lineage) %>%
  summarise(n_cells.day0 = n())

fp.summary <- metadat.day10_DABTRAM %>% 
  group_by(assigned_lineage) %>%
  summarise(
    min_fatepotential = min(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    max_fatepotential = max(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    mean_fatepotential = mean(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    sd_fatepotential = sd(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    n_cells = n()
  )


fp.summary$category <- ifelse((fp.summary$min_fatepotential >= 0), 
                               'priming', NA)

fp.summary$category <- ifelse((fp.summary$min_fatepotential < 0) & 
                              (fp.summary$max_fatepotential > 0) &
                              (fp.summary$sd_fatepotential > quantile(fp.summary$sd_fatepotential, 0.5, na.rm = TRUE)), 
                              'plasticity', fp.summary$category)



lin.priming <- fp.summary %>%
  filter(category == 'priming') %>%
  pull(assigned_lineage)

lin.plasticity <- fp.summary %>%
  filter(category == 'plasticity') %>%
  pull(assigned_lineage)

# intersect(lin.priming, lin.plasticity)


saver.mat <- all_data.day0[["Saver"]]@data
saver.mat <- t(saver.mat)

metadat.day0 <- all_data.day0@meta.data
priming_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.priming]
plastic_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.plasticity]

priming_progenitors <- intersect(priming_progenitors, adapting_progenitor)
plastic_progenitors <- intersect(plastic_progenitors, adapting_progenitor)


ggplot(ft.umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 0.5, color = 'gray') +
  geom_point(data = ft.umap.DABTRAM[priming_progenitors, ], 
             color = 'blue', size = 0.5) +
  geom_point(data = ft.umap.DABTRAM[plastic_progenitors, ], 
             color = 'red', size = 0.5) +
  labs(title = "Priming Progenitors at Day 0") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(wnn.umap, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(size = 0.5, color = 'gray') +
  geom_point(data = wnn.umap[priming_progenitors, ], 
             color = 'blue', size = 0.5) +
  geom_point(data = wnn.umap[plastic_progenitors, ], 
             color = 'red', size = 0.5) +
  labs(title = "Priming Progenitors at Day 0") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(pVI.umap.DABTRAM, aes(x = pvidabtramumap_1, y = pvidabtramumap_2)) +
  geom_point(size = 0.5, color = 'gray') +
  geom_point(data = pVI.umap.DABTRAM[priming_progenitors, ], 
             color = 'blue', size = 0.5) +
  geom_point(data = pVI.umap.DABTRAM[plastic_progenitors, ], 
             color = 'red', size = 0.5) +
  labs(title = "Priming Progenitors at Day 0") +
  theme_minimal() +
  theme(legend.position = "none")


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
gene_df$gene <- rownames(gene_df)

gene_df$padj <- p.adjust(gene_df$p.value, method = "BH")
gene_df$neglog10_pval <- -log10(gene_df$p.value)

jackpot <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
             "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
             "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
             "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
             "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")

isg.rs <- c('IFI27', 'IRF7', 'USP18', 'BST2', 'CXCL10', 'DDX60',
            'HERC6', 'HLA-B', 'HLA-G', 'IFI35', 'IFI44', 'IFI44L',
            'IFIT1', 'IFIT3', 'ISG15', 'LGALS3BP', 'LY6E', 'MX1',
            'MX2', 'OAS3', 'OASL', 'PLSCR1', 'STAT1', 'TRIM14',
            'HSD17B1','OAS1','CA2','CCNA1','CXCL1','GALC','IFI6',
            'IFITM1','LAMP3','MCL1','ROBO1','SLC6A15','THBS1','TIMP3')

gene_df$category <- ifelse(gene_df$gene %in% jackpot, 'Jackpot', 'Other')
gene_df$category <- ifelse(gene_df$gene %in% isg.rs, 'ISG.RS', gene_df$category)
gene_df$category <- factor(gene_df$category, 
                           levels = c('Jackpot', 'ISG.RS', 'Other'))
gene_df <- gene_df[order(gene_df$category, decreasing = T), ]

ggplot(gene_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = c('Jackpot' = 'red', 'ISG.RS' = 'blue', 'Other' = "grey")) +
  labs(x = "Log2 Fold Change (Adapting vs Non-Adapting)",
       y = "-log10(p-value)",
       title = "Differential Expression of Adapting Progenitors") +
  theme_minimal() +
  theme(legend.position = "none")

gene_df$gene <- rownames(gene_df)


saver.mat.d10 <- all_data.day10_DABTRAM[["Saver"]]@data
saver.mat.d10 <- t(saver.mat.d10)
wilcox_results_d10 <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = saver.mat.d10[priming_current,j],
    y = saver.mat.d10[plastic_current,j]
  )
  logfc <- log2(mean(saver.mat.d10[priming_current,j])) - log2(mean(saver.mat.d10[plastic_current,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results_d10 <- t(wilcox_results_d10)

gene_df_d10 <- as.data.frame(wilcox_results_d10)
rownames(gene_df_d10) <- colnames(saver.mat.d10)
colnames(gene_df_d10) <- c("logfc", "p.value")
gene_df$gene <- rownames(gene_df)

gene_df$padj <- p.adjust(gene_df$p.value, method = "BH")
gene_df$neglog10_pval <- -log10(gene_df$p.value)


# chromVar
cv.mat <- all_data.day0[["chromVar_day0"]]@data
cv.mat <- t(cv.mat)

metadat.day0 <- all_data.day0@meta.data
priming_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.priming]
plastic_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.plasticity]




# priming_progenitors <- intersect(priming_progenitors, adapting_progenitor)
# plastic_progenitors <- intersect(plastic_progenitors, adapting_progenitor)

rest <- non_adapting_progenitor

p <- ncol(cv.mat)

wilcox_results_chromVar <- sapply(1:p, function(j){
  # print(j)
  var_j = var(cv.mat[c(priming_progenitors, plastic_progenitors), j], na.rm = T)
  if(var_j == 0 | is.na(var_j)) {
    print(paste("Skipping feature", j, "due to zero variance."))
    return(c(logfc = NA, p.value = NA))
  }
  
  tmp <- stats::wilcox.test(
    x = cv.mat[priming_progenitors,j],
    y = cv.mat[plastic_progenitors,j]
  )
  logfc <- median(cv.mat[priming_progenitors,j], na.rm = T) - median(cv.mat[plastic_progenitors,j], na.rm = T)
  
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

ap1 <- grep('FOS|JUN', tf_df$tf, value = TRUE)
ap1 <- ap1[!grepl('var.2', ap1)]  # Exclude FOSB and JUNB

ifn <- c('STAT1', 'STAT3', 'IRF3', 'STAT1::STAT2')
emt <- c('SNAI1', 'SNAI2', 'SNAI3', 'SNAI4', 'ZEB1', 'TWIST1')

tf_df$category <- ifelse(tf_df$tf %in% ap1, 'AP1', 'Other')
tf_df$category <- ifelse(tf_df$tf %in% ifn, 'IFN', tf_df$category)
tf_df$category <- ifelse(tf_df$tf %in% emt, 'EMT', tf_df$category)
tf_df$category <- factor(tf_df$category, levels = c('AP1', 'IFN', 'EMT', 'Other'))
tf_df <- tf_df[order(tf_df$category, decreasing = T), ]

ggplot(tf_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  ggrepel::geom_text_repel(data = subset(tf_df, category == 'AP1'), aes(label = tf),  color = 'red', 
                       size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df, category == 'IFN'), aes(label = tf),  color = 'blue', 
                           size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df, category == 'EMT'), aes(label = tf),  color = '#059212', 
                           size = 3, show.legend = FALSE) +
  scale_color_manual(values = c('Other' = "grey", 'AP1' = "red", 'IFN' = 'blue', 'EMT' = '#059212')) +
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change (Plastic <-> Priming)",
       y = "-log10(p-value)",
       title = "Differential TF") +
  theme_minimal()

x = cv.mat[priming_progenitors, 'FOSL1']
y = cv.mat[plastic_progenitors,'FOSL1']
wilcox.test(x, y)

df <- data_frame(c(rep('Priming Progenitors', length(x)),
                   rep('Plastic progenitors', length(y))),
                 value = c(x, y))
colnames(df) <- c('category', 'FOSL1')

ggplot(df, aes(x = category, y = FOSL1)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", size=2, color="red") +
  labs(x = "Category",
       y = "FOSL1 Expression",
       title = "FOSL1 Expression in Adapting Progenitors") +
  theme_minimal() +
  theme(legend.position = "none")


metadat.day10_DABTRAM <- all_data.day10_DABTRAM@meta.data
priming_current <- rownames(metadat.day10_DABTRAM)[metadat.day10_DABTRAM$assigned_lineage %in% lin.priming]
plastic_current <- rownames(metadat.day10_DABTRAM)[metadat.day10_DABTRAM$assigned_lineage %in% lin.plasticity]

cv.mat.d10 <- all_data.day10_DABTRAM[["chromVar_day10_DABTRAM"]]@data
cv.mat.d10 <- t(cv.mat.d10)


wilcox_results_chromVar_d10 <- sapply(1:p, function(j){
  # print(j)
  var_j = var(cv.mat.d10[c(priming_current, plastic_current), j], na.rm = T)
  if(var_j == 0 | is.na(var_j)) {
    print(paste("Skipping feature", j, "due to zero variance."))
    return(c(logfc = NA, p.value = NA))
  }
  
  tmp <- stats::wilcox.test(
    x = cv.mat.d10[priming_current,j],
    y = cv.mat.d10[plastic_current,j]
  )
  logfc <- median(cv.mat.d10[priming_current,j], na.rm = T) - median(cv.mat.d10[plastic_current,j], na.rm = T)
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})


wilcox_results_chromVar_d10 <- t(wilcox_results_chromVar_d10)

tf_df_d10 <- as.data.frame(wilcox_results_chromVar_d10)
rownames(tf_df_d10) <- colnames(cv.mat.d10)
colnames(tf_df_d10) <- c("logfc", "p.value")

tf_df_d10$padj <- p.adjust(tf_df_d10$p.value, method = "BH")
tf_df_d10$neglog10_pval <- -log10(tf_df_d10$p.value)

tf_df_d10$tf <- rownames(tf_df_d10)

tf_df_d10$category <- ifelse(tf_df_d10$tf %in% ap1, 'AP1', 'Other')
tf_df_d10$category <- ifelse(tf_df_d10$tf %in% ifn, 'IFN', tf_df_d10$category)
tf_df_d10$category <- ifelse(tf_df_d10$tf %in% emt, 'EMT', tf_df_d10$category)
tf_df_d10$category <- factor(tf_df_d10$category, levels = c('AP1', 'IFN', 'EMT', 'Other'))
tf_df_d10 <- tf_df_d10[order(tf_df_d10$category, decreasing = T), ]

ggplot(tf_df_d10, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'AP1'), aes(label = tf),  color = 'red', 
                           size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'IFN'), aes(label = tf),  color = 'blue', 
                           size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'EMT'), aes(label = tf),  color = '#059212', 
                           size = 3, show.legend = FALSE) +
  scale_color_manual(values = c('Other' = "grey", 'AP1' = "red", 'IFN' = 'blue', 'EMT' = '#059212')) +
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change (Plastic <-> Priming)",
       y = "-log10(p-value)",
       title = "Differential TF") +
  theme_minimal()



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




wilcox_results_gavish_mp_d10 <- sapply(metaprograms, function(j){
  tmp <- stats::wilcox.test(
    x = gavish.mp[priming_current,j],
    y = gavish.mp[plastic_current,j]
  )
  logfc <- log2(mean(gavish.mp[priming_current,j])) - log2(mean(gavish.mp[plastic_current,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results_gavish_mp_d10 <- t(wilcox_results_gavish_mp_d10)
wilcox_results_gavish_mp_d10_df <- as.data.frame(wilcox_results_gavish_mp_d10)
rownames(wilcox_results_gavish_mp_d10_df) <- metaprograms
colnames(wilcox_results_gavish_mp_d10_df) <- c("logfc", "p.value")


