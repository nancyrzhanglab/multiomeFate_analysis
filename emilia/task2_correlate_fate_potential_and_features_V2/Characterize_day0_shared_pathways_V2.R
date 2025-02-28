library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(circlize)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

# ==============================================================================
# Read correlation results
# ==============================================================================
load(paste0(result_dir, 'saver_cor_vec.RData'))

# ==============================================================================
# Get and sort table
# ==============================================================================

getAndSortTable <- function(name) {
  cor_vec <- as.data.frame(saver_cor_vec[[name]]) %>% drop_na()
  cor_vec$p.value <- as.numeric(cor_vec$p.value)
  cor_vec$p_adj <- p.adjust(cor_vec$p.value, method = 'BH')
  cor_vec <- cor_vec[cor_vec$p_adj < 0.05, ]
  
  cor_vec$gene <- rownames(cor_vec)
  cor_vec$correlation <- as.numeric(cor_vec$correlation)
  cor_vec <- cor_vec[order(cor_vec$correlation, decreasing = TRUE),]
  
  return(cor_vec)
}

day0_dabtram <- getAndSortTable('dabtram_d0_saver_cor_vec')
day0_cocl2 <- getAndSortTable('cocl2_d0_saver_cor_vec')
day0_cis <- getAndSortTable('cis_d0_saver_cor_vec')

day10_dabtram <- getAndSortTable('dabtram_d10_saver_cor_vec')
day10_cocl2 <- getAndSortTable('cocl2_d10_saver_cor_vec')
day10_cis <- getAndSortTable('cis_d10_saver_cor_vec')

# ==============================================================================
# Geneset enrichment analysis (POSITIVE)
# ==============================================================================

# day0 dabtram positive
day0_dabtram.positive <- day0_dabtram[day0_dabtram$correlation > 0, ]
day0_dabtram.positive <- day0_dabtram.positive[!grepl('\\.', day0_dabtram.positive$gene), ]
eg.day0_dabtram = bitr(day0_dabtram.positive$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x.day0_dabtram <- enrichPathway(gene=eg.day0_dabtram$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x.day0_dabtram_res <- x.day0_dabtram@result
x.day0_dabtram_res <- x.day0_dabtram_res[x.day0_dabtram_res$p.adjust < 0.05, ]
x.day0_dabtram_res <- x.day0_dabtram_res[, c('ID', 'Description', 'p.adjust', 'GeneRatio')]
colnames(x.day0_dabtram_res) <- c('ID', 'Description', 'p.adjust.day0_DABTRAM', 'GeneRatio.day0_DABTRAM')

# day0 cocl2 positive
day0_cocl2.positive <- day0_cocl2[day0_cocl2$correlation > 0.08, ]
day0_cocl2.positive <- day0_cocl2.positive[!grepl('\\.', day0_cocl2.positive$gene), ]
eg.day0_cocl2 = bitr(day0_cocl2.positive$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x.day0_cocl2 <- enrichPathway(gene=eg.day0_cocl2$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x.day0_cocl2_res <- x.day0_cocl2@result
x.day0_cocl2_res <- x.day0_cocl2_res[x.day0_cocl2_res$p.adjust < 0.05, ]
x.day0_cocl2_res <- x.day0_cocl2_res[, c('ID', 'Description', 'p.adjust', 'GeneRatio')]
colnames(x.day0_cocl2_res) <- c('ID', 'Description', 'p.adjust.day0_COCL2', 'GeneRatio.day0_COCL2')

# day0 cis positive
day0_cis.positive <- day0_cis[day0_cis$correlation > 0.08619867, ]
day0_cis.positive <- day0_cis.positive[!grepl('\\.', day0_cis.positive$gene), ]
eg.day0_cis = bitr(day0_cis.positive$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x.day0_cis <- enrichPathway(gene=eg.day0_cis$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x.day0_cis_res <- x.day0_cis@result
x.day0_cis_res <- x.day0_cis_res[x.day0_cis_res$p.adjust < 0.05, ]
x.day0_cis_res <- x.day0_cis_res[, c('ID', 'Description', 'p.adjust', 'GeneRatio')]
colnames(x.day0_cis_res) <- c('ID', 'Description', 'p.adjust.day0_CIS', 'GeneRatio.day0_CIS')


# find common pathways
x.day0 <- merge(x.day0_dabtram_res, x.day0_cocl2_res, by = c('ID', 'Description'), all = TRUE)
x.day0 <- merge(x.day0, x.day0_cis_res, by = c('ID', 'Description'), all = TRUE)

# count the number of NAs in table
x.day0$GeneRatio.day0_DABTRAM[is.na(x.day0$GeneRatio.day0_DABTRAM)] <- 0
x.day0$GeneRatio.day0_COCL2[is.na(x.day0$GeneRatio.day0_COCL2)] <- 0
x.day0$GeneRatio.day0_CIS[is.na(x.day0$GeneRatio.day0_CIS)] <- 0
x.day0$n_NA <- rowSums(is.na(x.day0))
x.day0 <- x.day0[x.day0$n_NA < 2, ]
x.day0 <- x.day0[order(x.day0$n_NA), ]
write.csv(x.day0, '~/Downloads/day0_positive_pathways.csv', row.names = FALSE)

x.day0 <- read.csv(paste0(result_dir, 'day0_positive_pathways.csv'))
x.day0 <- x.day0[x.day0$Plot != 'No', ]

x.day0$p.adjust.day0_DABTRAM <- (-1) * log10(x.day0$p.adjust.day0_DABTRAM)
x.day0$p.adjust.day0_COCL2 <- (-1) * log10(x.day0$p.adjust.day0_COCL2)
x.day0$p.adjust.day0_CIS <- (-1) * log10(x.day0$p.adjust.day0_CIS)

# plot
rownames(x.day0) <- x.day0$Description
colnames(x.day0) <- c('ID', 'Description', 'day0_DABTRAM', 'GeneRatio.day0_DABTRAM', 'day0_COCL2', 'GeneRatio.day0_COCL2', 'day0_CIS', 'GeneRatio.day0_CIS', 'n_NA')
pheatmap::pheatmap(as.matrix(x.day0[, c('day0_DABTRAM', 'day0_COCL2', 'day0_CIS')]), 
                  show_rownames = T, show_colnames = T, cluster_rows = FALSE, cluster_cols = FALSE, 
                  scale = 'none', border_color = NA, fontsize = 8, cellwidth = 18, cellheight = 18,
                  color = colorRampPalette(brewer.pal(5, "RdPu"))(6), main = 'Day0 Positive Pathways',
                  filename = paste0(result_dir, 'day0_positive_pathways.png'))



plot(x.day0_cis_res)

length(intersect(x.day0_dabtram_res$ID, x.day0_cocl2_res$ID))
length(intersect(x.day0_dabtram_res$ID, x.day0_cis_res$ID))
length(intersect(x.day0_cocl2_res$ID, x.day0_cis_res$ID))




# ==============================================================================
# Geneset enrichment analysis (NEGATIVE) - DOES NOT WORK
# ==============================================================================

# day0 dabtram negative
day0_dabtram.negative <- day0_dabtram[day0_dabtram$correlation < -0.21578313, ]
day0_dabtram.negative <- day0_dabtram.negative[!grepl('\\.', day0_dabtram.negative$gene), ]
eg.day0_dabtram = bitr(day0_dabtram.negative$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x.day0_dabtram <- enrichPathway(gene=eg.day0_dabtram$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
x.day0_dabtram_res <- x.day0_dabtram@result
x.day0_dabtram_res <- x.day0_dabtram_res[x.day0_dabtram_res$p.adjust < 0.05, ]

# # day0 cocl2 negative
# day0_cocl2.positive <- day0_cocl2[day0_cocl2$correlation > 0, ]
# eg.day0_cocl2 = bitr(day0_cocl2.positive$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# x.day0_cocl2 <- enrichPathway(gene=eg.day0_cocl2$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
# x.day0_cocl2_res <- x.day0_cocl2@result
# x.day0_cocl2_res <- x.day0_cocl2_res[x.day0_cocl2_res$p.adjust < 0.05, ]
# 
# # day0 cis positive
# day0_cis.positive <- day0_cis[day0_cis$correlation > 0, ]
# eg.day0_cis = bitr(day0_cis.positive$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# x.day0_cis <- enrichPathway(gene=eg.day0_cis$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
# x.day0_cis_res <- x.day0_cis@result
# x.day0_cis_res <- x.day0_cis_res[x.day0_cis_res$p.adjust < 0.05, ]

