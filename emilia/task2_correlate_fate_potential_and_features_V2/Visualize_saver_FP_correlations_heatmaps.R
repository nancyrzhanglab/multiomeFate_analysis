library(tidyverse)
library(ggplot2)
library(GGally)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/saver_cor_vec.RData')

dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']] 

dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']] 

dabtram_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_saver_cor_SMS.csv')
cocl2_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_saver_cor_SMS.csv')
cis_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_saver_cor_SMS.csv')

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_saver_cor_vec) <- paste0(colnames(dabtram_d0_saver_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_saver_cor_vec) <- paste0(colnames(cocl2_d0_saver_cor_vec), '.COCL2_d0')
colnames(cis_d0_saver_cor_vec) <- paste0(colnames(cis_d0_saver_cor_vec), '.CIS_d0')

colnames(dabtram_d10_saver_cor_vec) <- paste0(colnames(dabtram_d10_saver_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_saver_cor_vec) <- paste0(colnames(cocl2_d10_saver_cor_vec), '.COCL2_d10')
colnames(cis_d10_saver_cor_vec) <- paste0(colnames(cis_d10_saver_cor_vec), '.CIS_d10')

# Day0
d0_cor <- merge(dabtram_d0_saver_cor_vec, cocl2_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

d0_cor <- merge(d0_cor, cis_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

# Day10
d10_cor <- merge(dabtram_d10_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

# Filter for genes of interest
dabtram_annotations <- dabtram_annotations[dabtram_annotations$Gene.of.interest == 'Yes', ]
cocl2_annotations <- cocl2_annotations[cocl2_annotations$Gene.of.interest == 'Yes', ]
cis_annotations <- cis_annotations[cis_annotations$Gene.of.interest == 'Yes', ]

annotations <- rbind(dabtram_annotations[, c('gene', 'Pathway.')], 
                     cocl2_annotations[, c('gene', 'Pathway.')],
                     cis_annotations[, c('gene', 'Pathway.')]) %>% distinct()

# group annotations by gene and merge pathway information
annotations <- annotations %>% group_by(gene) %>% summarise(Pathway. = paste(Pathway., collapse = ', '))


annotations <- annotations[order(annotations$Pathway.), ]
annotations$order <- seq(1, nrow(annotations))

d0_cor_anno <- d0_cor[rownames(d0_cor) %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]
d10_cor_anno <- d10_cor[rownames(d10_cor) %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]


pheatmap::pheatmap(d0_cor_anno[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0')])

d0_cor_anno <- merge(d0_cor_anno, annotations, by.x = 'row.names', by.y = 'gene')
d0_cor_anno <- d0_cor_anno[order(d0_cor_anno$order), ]
rownames(d0_cor_anno) <- d0_cor_anno$Row.names

pdf("~/Downloads/heatmap_output.pdf", width = 8, height = 10)
Heatmap(d0_cor_anno[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0')],
        col = colorRamp2(seq(-1,1,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE, row_split = d0_cor_anno$`Pathway.`, show_row_names = T, row_title_rot = 0,
        border = TRUE)
dev.off()


d10_cor_anno <- merge(d10_cor_anno, annotations, by.x = 'row.names', by.y = 'gene')
d10_cor_anno <- d10_cor_anno[order(d10_cor_anno$order), ]
rownames(d10_cor_anno) <- d10_cor_anno$Row.names

pdf("~/Downloads/heatmap_output_d10.pdf", width = 8, height = 10)
Heatmap(d10_cor_anno[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        col = colorRamp2(seq(-1,1,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE, row_split = d10_cor_anno$`Pathway.`, show_row_names = T, row_title_rot = 0,
        border = TRUE)
dev.off()

d0_cor_anno1 <- d0_cor_anno[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'Pathway.')]
d0_cor_anno1$Gene <- rownames(d0_cor_anno1)
d0_cor_anno1$Day <- 'Day0'
# colnames(d0_cor_anno1) <- c('DABTRAM', 'COCL2', 'CIS', 'Pathway', 'Day')

d10_cor_anno1 <- d10_cor_anno[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10', 'Pathway.')]
d10_cor_anno1$Day <- 'Day10'
d10_cor_anno1$Gene <- rownames(d10_cor_anno1)
colnames(d10_cor_anno1) <- c('DABTRAM', 'COCL2', 'CIS', 'Pathway', 'Day')

d0_d10 <- merge(d0_cor_anno1, d10_cor_anno1, by = c('Pathway.', 'Gene'))

write.csv(d0_d10, '~/Downloads/d0_d10_correlation.csv', row.names = F)
d0_d10 <- read.csv('~/Downloads/d0_d10_correlation.csv')
pdf("~/Downloads/heatmap_output_d0_d10.pdf", width = 7, height = 11)
Heatmap(d0_d10[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE, row_split = d0_d10$`Pathway.`, show_row_names = T, row_title_rot = 0, 
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()

pdf("~/Downloads/heatmap_output_d0_d10_2.pdf", width = 7, height = 11)
Heatmap(d0_d10[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = T, show_row_names = F,
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()

d0_d10_2 <- merge(d0_cor, d10_cor, by = 'row.names') %>% drop_na()
pdf("~/Downloads/heatmap_output_d0_d10_2.pdf", width = 4, height = 10)
Heatmap(d0_d10_2[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        col = colorRamp2(seq(-1,1, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = T, show_row_names = F, row_title_rot = 0,
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()


write.csv(d0_d10, '~/Downloads/d0_d10_correlation_mod.csv', row.names = F)
d0_d10 <- read.csv('~/Downloads/d0_d10_correlation_mod.csv')
pathway_order <- unique(d0_d10$Pathway.)
d0_d10$Pathway. <- factor(d0_d10$Pathway., levels = pathway_order)
pdf("~/Downloads/heatmap_output_d0_d10_3.pdf", width = 8, height = 13)
Heatmap(d0_d10[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE, row_split = d0_d10$`Pathway.`, show_row_names = T, row_title_rot = 0, 
        cluster_columns = F, row_labels = d0_d10$Gene,
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()


# ==============================================================================
# overlap
# ==============================================================================
library(VennDiagram)

topN <- 0.10
dabtram_d0_saver_cor_vec <- dabtram_d0_saver_cor_vec[order(dabtram_d0_saver_cor_vec$correlation.DABTRAM_d0, decreasing = T), ]
dabtram_d0_saver_cor_vec.topN <- dabtram_d0_saver_cor_vec[1:round(nrow(dabtram_d0_saver_cor_vec)*topN), ]
dabtram_d0_saver_cor_vec.topN <- dabtram_d0_saver_cor_vec.topN[dabtram_d0_saver_cor_vec.topN$correlation.DABTRAM_d0 > 0, ]

cocl2_d0_saver_cor_vec <- cocl2_d0_saver_cor_vec[order(cocl2_d0_saver_cor_vec$correlation.COCL2_d0, decreasing = T), ]
cocl2_d0_saver_cor_vec.topN <- cocl2_d0_saver_cor_vec[1:round(nrow(cocl2_d0_saver_cor_vec)*topN), ]
cocl2_d0_saver_cor_vec.topN <- cocl2_d0_saver_cor_vec.topN[cocl2_d0_saver_cor_vec.topN$correlation.COCL2_d0 > 0, ]

cis_d0_saver_cor_vec <- cis_d0_saver_cor_vec[order(cis_d0_saver_cor_vec$correlation.CIS_d0, decreasing = T), ]
cis_d0_saver_cor_vec.topN <- cis_d0_saver_cor_vec[1:round(nrow(cis_d0_saver_cor_vec)*topN), ]
cis_d0_saver_cor_vec.topN <- cis_d0_saver_cor_vec.topN[cis_d0_saver_cor_vec.topN$correlation.CIS_d0 > 0, ]

venn.diagram(list(dabtram_d0 = rownames(dabtram_d0_saver_cor_vec.topN), 
                  cocl2_d0 = rownames(cocl2_d0_saver_cor_vec.topN),
                  cis_d0 = rownames(cis_d0_saver_cor_vec.topN)),
             fill = c("red", "green", 'blue'),
             alpha = c(0.5, 0.5, 0.5), cex = 2, cat.fontface = 4,lty =2,
             filename = "~/Downloads/test.png", output = TRUE)

library(eulerr)

# day0
#specify values to use in venn diagram
fit <- euler(c('DABTRAM' = length(rownames(dabtram_d0_saver_cor_vec.topN)), 
               'COCL2' = length(rownames(cocl2_d0_saver_cor_vec.topN)),
               'CIS' = length(rownames(cis_d0_saver_cor_vec.topN)),
               'DABTRAM&COCL2' = length(intersect(rownames(dabtram_d0_saver_cor_vec.topN), 
                                                  rownames(cocl2_d0_saver_cor_vec.topN))),
               'DABTRAM&CIS' = length(intersect(rownames(dabtram_d0_saver_cor_vec.topN), 
                                                rownames(cis_d0_saver_cor_vec.topN))),
               'COCL2&CIS' = length(intersect(rownames(cocl2_d0_saver_cor_vec.topN),
                                              rownames(cis_d0_saver_cor_vec.topN))),
               'DABTRAM&COCL2&CIS' = length(intersect(intersect(rownames(dabtram_d0_saver_cor_vec.topN), 
                                                                rownames(cocl2_d0_saver_cor_vec.topN)),
                                                      rownames(cis_d0_saver_cor_vec.topN)))))

#create venn diagram
plot(fit,
     fills = c("#9D85BE", "#6DC49C", "#FBD08C"),
     fontsize = 32,
     quantities = list(fontsize = 32))

# Day 10

topN <- 0.10
dabtram_d10_saver_cor_vec <- dabtram_d10_saver_cor_vec[order(dabtram_d10_saver_cor_vec$correlation.DABTRAM_d10, decreasing = T), ]
dabtram_d10_saver_cor_vec.topN <- dabtram_d10_saver_cor_vec[1:round(nrow(dabtram_d10_saver_cor_vec)*topN), ]
dabtram_d10_saver_cor_vec.topN <- dabtram_d10_saver_cor_vec.topN[dabtram_d10_saver_cor_vec.topN$correlation.DABTRAM_d10 > 0, ]

cocl2_d10_saver_cor_vec <- cocl2_d10_saver_cor_vec[order(cocl2_d10_saver_cor_vec$correlation.COCL2_d10, decreasing = T), ]
cocl2_d10_saver_cor_vec.topN <- cocl2_d10_saver_cor_vec[1:round(nrow(cocl2_d10_saver_cor_vec)*topN), ]
cocl2_d10_saver_cor_vec.topN <- cocl2_d10_saver_cor_vec.topN[cocl2_d10_saver_cor_vec.topN$correlation.COCL2_d10 > 0, ]

cis_d10_saver_cor_vec <- cis_d10_saver_cor_vec[order(cis_d10_saver_cor_vec$correlation.CIS_d10, decreasing = T), ]
cis_d10_saver_cor_vec.topN <- cis_d10_saver_cor_vec[1:round(nrow(cis_d10_saver_cor_vec)*topN), ]
cis_d10_saver_cor_vec.topN <- cis_d10_saver_cor_vec.topN[cis_d10_saver_cor_vec.topN$correlation.CIS_d10 > 0, ]

#specify values to use in venn diagram
fit <- euler(c('DABTRAM' = length(rownames(dabtram_d10_saver_cor_vec.topN)), 
               'COCL2' = length(rownames(cocl2_d10_saver_cor_vec.topN)),
               'CIS' = length(rownames(cis_d10_saver_cor_vec.topN)),
               'DABTRAM&COCL2' = length(intersect(rownames(dabtram_d10_saver_cor_vec.topN), 
                                                  rownames(cocl2_d10_saver_cor_vec.topN))),
               'DABTRAM&CIS' = length(intersect(rownames(dabtram_d10_saver_cor_vec.topN), 
                                                rownames(cis_d10_saver_cor_vec.topN))),
               'COCL2&CIS' = length(intersect(rownames(cocl2_d10_saver_cor_vec.topN),
                                              rownames(cis_d10_saver_cor_vec.topN))),
               'DABTRAM&COCL2&CIS' = length(intersect(intersect(rownames(dabtram_d10_saver_cor_vec.topN), 
                                                                rownames(cocl2_d10_saver_cor_vec.topN)),
                                                      rownames(cis_d10_saver_cor_vec.topN)))))

#create venn diagram
plot(fit,
     fills = c("#623594", "#0F8241", "#C96D29"),
     fontsize = 32,
     quantities = list(fontsize = 32))

