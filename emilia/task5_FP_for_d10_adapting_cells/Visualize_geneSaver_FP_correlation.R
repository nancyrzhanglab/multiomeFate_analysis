rm(list = ls())

library(ggplot2)
library(GGally)
library(hdrcde)
library(ggdensity)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(gridExtra)

theme_Publication<- function(base_size=14, base_family="sans") {
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

# =============================================================================
# Read data
# =============================================================================
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- cor_vec

load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_COCL2.RData'))
cor_vec.COCL2 <- cor_vec

load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_CIS.RData'))
cor_vec.CIS <- cor_vec

ref_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Markers/'
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
                 "CAV1")),
  CIS = sort(c("YY1AP1", "LGALS3", "MCF2L", "TIMM50", "AC207130.1",
               "SLC25A6", "EIF3L", "CTSD", "NQO1", "HNMT", "ZFYVE16",
               "PHACTR1", "TNFRSF14", "RAI14", "TRPM1", "HIST1H1C",
               "HIST2H2AC", "SPARC", "TRIM63", "TUBA1B", "HIST1H1A",
               "HIST1H1D", "PYCARD", "FSTL1", "DCT", "CTSK", "HIST1H4C",
               "GDF15", "HIST1H1B"))
)
keygenes.list <- unlist(keygenes)

# res_genes <- c("WNT5A", "AXL", "EGFR", "PDGFRB", "JUN", "NGFR", "PCNA")
res_genes <- c("WNT5A", "AXL", "EGFR", "JUN", "NGFR", "PCNA")
# =============================================================================
# Wrangle
# =============================================================================
colnames(cor_vec.DABTRAM) <- c('gene', 'cor.DATBRAM', 'p_val.DABTRAM')
colnames(cor_vec.COCL2) <- c('gene', 'cor.COCL2', 'p_val.COCL2')
colnames(cor_vec.CIS) <- c('gene', 'cor.CIS', 'p_val.CIS')

comp_df <- merge(cor_vec.DABTRAM, cor_vec.COCL2, by = 'gene')
comp_df <- merge(comp_df, cor_vec.CIS, by = 'gene')

ggpairs(comp_df[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')], 
        lower = list(continuous = 'points'))

comp_df <- as.data.frame(comp_df)

# comp_df$annotate <- ifelse(comp_df$gene %in% keygenes[['jackpot']], 'Jackpot', NA)
# comp_df$annotate[comp_df$gene %in% keygenes[['DABTRAM']]] <- 'DABTRAM'
# comp_df$annotate[comp_df$gene %in% keygenes[['CIS']]] <- 'CIS'
# comp_df$annotate[comp_df$gene %in% keygenes[['COCL2']]] <- 'COCL2'
# 
# comp_df$annotate <- factor(comp_df$annotate, levels = c('Jackpot', 'DABTRAM', 'COCL2', 'CIS'))

# write.csv(comp_df, paste0(out_dir, 'correlation_comparison.csv'), row.names = F)

# comp_df.toplot <- comp_df[!is.na(comp_df$annotate),]
# comp_df.toplot <- comp_df.toplot[order(comp_df.toplot$annotate),]
# rownames(comp_df.toplot) <- comp_df.toplot$gene

# =============================================================================
# Plot pair-wise comparison
# =============================================================================

x_min <- min(comp_df[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')])
x_max <- max(comp_df[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')])
p1 <- ggplot(comp_df, aes(x = cor.DATBRAM, y = cor.COCL2)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(x_min, x_max)) +
  ylim(c(x_min, x_max)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p2 <- ggplot(comp_df, aes(x = cor.CIS, y = cor.COCL2)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(x_min, x_max)) +
  ylim(c(x_min, x_max)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p3 <- ggplot(comp_df, aes(x = cor.DATBRAM, y = cor.CIS)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(x_min, x_max)) +
  ylim(c(x_min, x_max)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p4 <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave(paste0(out_dir, 'geneSaver_correlation_day0_panel.pdf'), 
       p4, 
       width = 8, height = 3, dpi = 300)

# =============================================================================
# Plot heatmap
# =============================================================================

dabtram_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_saver_cor_SMS.csv')
cocl2_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_saver_cor_SMS.csv')
cis_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_saver_cor_SMS.csv')
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

comp_df.toplot <- comp_df[comp_df$gene %in% annotations$gene,]
rownames(comp_df.toplot) <- comp_df.toplot$gene
comp_df.toplot <- comp_df.toplot[annotations$gene,]
# =============================================================================
# Plotting
# =============================================================================

# write.csv(comp_df.toplot, paste0(out_dir, 'keygenes_comparisons.csv'), row.names = T)
# write.csv(annotations, paste0(out_dir, 'keygenes_annotations.csv'), row.names = T)

comp_df.toplot <- read.csv(paste0(out_dir, 'keygenes_comparisons.csv'))
rownames(comp_df.toplot) <- comp_df.toplot$gene
annotations <- read.csv(paste0(out_dir, 'keygenes_annotations.csv'))
pdf("~/Downloads/heatmap_output_geneSaver_cor_d0.pdf", width = 7, height = 11)
Heatmap(comp_df.toplot[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')],
        col = colorRamp2(seq(-1,1,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE, cluster_columns = FALSE, row_split = annotations$Pathway., show_row_names = T, row_title_rot = 0,
        row_dend_width = unit(0, "mm"),
        border = TRUE)
dev.off()

# DABTRAM
comp_df <- comp_df[order(comp_df$cor.DATBRAM),]
comp_df$order.DABTRAM <- 1:nrow(comp_df)
ggplot(comp_df, aes(x = order.DABTRAM, y = cor.DATBRAM)) +
  geom_point() +
  geom_point(data = subset(comp_df, annotate == 'DABTRAM'), aes(color = 'key gene')) +
  ggrepel::geom_text_repel(data = subset(comp_df, annotate == 'DABTRAM'), aes(label = gene), 
                           nudge_x = 0.1, nudge_y = 0.1) +
  labs(title = paste0('Correlation of gene expression with fate bias (DABTRAM)'),
       x = 'Rank', y = 'Correlation') +
  theme_Publication()
ggsave(paste0(out_dir, 'correlation_DABTRAM.pdf'), width = 6, height = 4)

# COCL2
comp_df <- comp_df[order(comp_df$cor.COCL2),]
comp_df$order.COCL2 <- 1:nrow(comp_df)
ggplot(comp_df, aes(x = order.COCL2, y = cor.COCL2)) +
  geom_point() +
  geom_point(data = subset(comp_df, annotate == 'COCL2'), aes(color = 'key gene')) +
  ggrepel::geom_text_repel(data = subset(comp_df, annotate == 'COCL2'), aes(label = gene), 
                           nudge_x = 0.1, nudge_y = 0.1) +
  labs(title = paste0('Correlation of gene expression with fate bias (COCL2)'),
       x = 'Rank', y = 'Correlation') +
  theme_Publication()
ggsave(paste0(out_dir, 'correlation_COCL2.pdf'), width = 6, height = 4)

# CIS
comp_df <- comp_df[order(comp_df$cor.CIS),]
comp_df$order.CIS <- 1:nrow(comp_df)
ggplot(comp_df, aes(x = order.CIS, y = cor.CIS)) +
  geom_point() +
  geom_point(data = subset(comp_df, annotate == 'CIS'), aes(color = 'key gene')) +
  ggrepel::geom_text_repel(data = subset(comp_df, annotate == 'CIS'), aes(label = gene), 
                           nudge_x = 0.1, nudge_y = 0.1) +
  labs(title = paste0('Correlation of gene expression with fate bias (CIS)'),
       x = 'Rank', y = 'Correlation') +
  theme_Publication()
ggsave(paste0(out_dir, 'correlation_CIS.pdf'), width = 6, height = 4)




# DABTRAM day0
comp_df <- comp_df %>% arrange(cor.DATBRAM)
comp_df$order.DABTRAM <- 1:nrow(comp_df)
p1 <- ggplot(comp_df, aes(x = order.DABTRAM, y = cor.DATBRAM)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = comp_df[comp_df$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, gene %in% res_genes), 
                           max.overlaps = Inf,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# COCL2 day0
comp_df <- comp_df %>% arrange(cor.COCL2)
comp_df$order.COCL2 <- 1:nrow(comp_df)
p2 <- ggplot(comp_df, aes(x = order.COCL2, y = cor.COCL2)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = comp_df[comp_df$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, gene %in% res_genes), 
                           max.overlaps = Inf,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# CIS day0
comp_df <- comp_df %>% arrange(cor.CIS)
comp_df$order.CIS <- 1:nrow(comp_df)
p3 <- ggplot(comp_df, aes(x = order.CIS, y = cor.CIS)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = comp_df[comp_df$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, gene %in% res_genes), 
                           max.overlaps = Inf,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

p4 <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave(filename = paste0(figure_dir, 'Fig6.D0.RNA.pdf'), p4, width = 8, height = 2.5)







# DABTRAM
comp_df <- comp_df[order(comp_df$cor.DATBRAM),]
hist(comp_df$cor.DATBRAM, breaks = 20, main = 'DABTRAM', xlab = 'Correlation')

# Plot
keygenes.dabtram <- comp_df[comp_df$gene %in% keygenes[['DABTRAM']], 'cor.DATBRAM']
keygenes.cocl2 <- comp_df[comp_df$gene %in% keygenes[['COCL2']], 'cor.COCL2']
keygenes.cis <- comp_df[comp_df$gene %in% keygenes[['CIS']], 'cor.CIS']

p1 <- ggplot(comp_df, aes(x = cor.DATBRAM)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.dabtram, color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, gene %in% keygenes[['DABTRAM']]), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-0.5, 0.65) +
  labs(ylab = 'Density', xlab = 'Correlation (DABTRAM)')

p2 <- ggplot(comp_df, aes(x = cor.COCL2)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.cocl2, color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, gene %in% keygenes[['COCL2']]), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-0.5, 0.65) +
  labs(ylab = 'Density', xlab = 'Correlation (COCL2)')


p3 <- ggplot(comp_df, aes(x = cor.CIS)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.cis, color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, gene %in% keygenes[['CIS']]), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-0.5, 0.65) +
  labs(ylab = 'Density', xlab = 'Correlation (CIS)')


p4 <- grid.arrange(p1, p2, p3, ncol = 1)
ggsave(paste0(out_dir, 'correlation_density.pdf'), plot = p4, width = 6, height = 4)

# ggplot(comp_df, aes(x = cor.DATBRAM)) +
#   stat_density(aes(y = 1, fill = ..density..), geom = "tile", adjust = 1, n = 200) +
#   scale_fill_gradient(low = "#F5F5F5", high = "#696969") +
#   geom_vline(xintercept = keygenes.dabtram, color = "darkred", alpha = 0.8, size = 0.8) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.title.x = element_blank(),
#     panel.grid = element_blank()
#   ) +
#   labs(fill = "Density", 
#        title = "Correlation between gene expression and fate bias (DABTRAM)")


# =============================================================================
# Read data
# =============================================================================
fb.DABTRAM <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
fb.COCL2 <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_COCL2.csv'))
fb.CIS <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_CIS.csv'))

colnames(fb.DABTRAM) <- c('cell_id', 'adaptingFP.DABTRAM', 'nonAdaptingFP.DABTRAM', 'bias.DABTRAM')
colnames(fb.COCL2) <- c('cell_id', 'adaptingFP.COCL2', 'nonAdaptingFP.COCL2', 'bias.COCL2')
colnames(fb.CIS) <- c('cell_id', 'adaptingFP.CIS', 'nonAdaptingFP.CIS', 'bias.CIS')

comp_df <- merge(fb.DABTRAM, fb.COCL2, by = 'cell_id')
comp_df <- merge(comp_df, fb.CIS, by = 'cell_id')

ggpairs(comp_df[, c('bias.DABTRAM', 'bias.COCL2', 'bias.CIS')], 
        lower = list(continuous = 'points'))
