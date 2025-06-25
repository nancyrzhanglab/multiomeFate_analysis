rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)

library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis

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
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
fatebias_out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task10_priming_vs_plasticity/'

figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V5/Figure5/'

remove_unassigned_cells <- TRUE

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

stress <- c('ATF3', 'EGR1', 'FOS', 'PPP1R15A', 'JUN', 'DUSP1','JUNB',
            'FOSB','IER2','GADD45B','BTG2','ZFP36','DNAJB1','NR4A1',
            'SERTAD1','DDIT3','KLF6','NFKBIA','HSPA1A','HSPA1B','CYR61',
            'MAFF','NFKBIZ','IER3', 'CCNL1','MCL1','RHOB','SOCS3','IRF1',
            'CDKN1A','HES1','ID2','KLF4','KLF10','CXCL2','ID1','TOB1',
            'TRIB1','HSPH1','RASD1','SAT1','DDIT4','EGR2','HERPUD1',
            'PLK2','DNAJA1','HSPA6','PMAIP1','CITED2','ID3')

# =============================================================================
# Read data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVAR_day0.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_week5_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[['chromVar_day0']] <- all_data_chromVar_day0
all_data[['chromVar_day10_DABTRAM']] <- all_data_chromVar_day10_DABTRAM
all_data[['chromVar_week5_DABTRAM']] <- all_data_chromVar_week5_DABTRAM


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
all_data.week5_DABTRAM <- subset(all_data, dataset == 'week5_DABTRAM')

saver.mat.d0 <- all_data.day0@assays[["Saver"]]@data
saver.mat.d0 <- t(saver.mat.d0)

saver.mat.d10 <- all_data.day10_DABTRAM@assays[["Saver"]]@data
saver.mat.d10 <- t(saver.mat.d10)

cv.mat.d0 <- all_data.day0@assays[["chromVar_day0"]]@data
cv.mat.d0 <- t(cv.mat.d0)

cv.mat.d10 <- all_data.day10_DABTRAM@assays[["chromVar_day10_DABTRAM"]]@data
cv.mat.d10 <- t(cv.mat.d10)

cv.mat.w5 <- all_data.week5_DABTRAM@assays[["chromVar_week5_DABTRAM"]]@data
cv.mat.w5 <- t(cv.mat.w5)

adapting_fate_bias <- read.csv(paste0(fatebias_out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
high.adapting.fate.d0 <- adapting_fate_bias[adapting_fate_bias$bias > median(adapting_fate_bias$bias), ]
low.adapting.fate.d0 <- adapting_fate_bias[adapting_fate_bias$bias <= median(adapting_fate_bias$bias), ]

# =============================================================================
# Classify priming and plasticity at day 10
# =============================================================================
metadat.day0 <- all_data.day0@meta.data
metadat.day10_DABTRAM <- all_data.day10_DABTRAM@meta.data
metadat.day10_DABTRAM$cell_id <- rownames(metadat.day10_DABTRAM)

metadat.week5_DABTRAM <- all_data.week5_DABTRAM@meta.data

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


fp.summary$category <- ifelse((fp.summary$min_fatepotential >= 0) &
                              (fp.summary$n_cells >= 2),
                              'priming', NA)

fp.summary$category <- ifelse((fp.summary$min_fatepotential < 0) &
                                (fp.summary$max_fatepotential > 0) &
                                (fp.summary$sd_fatepotential > quantile(fp.summary$sd_fatepotential, 0.5, na.rm = TRUE)),
                              'plasticity', fp.summary$category)


# fp.summary$category <- ifelse((fp.summary$min_fatepotential >= 0),
#                               'priming', NA)
# 
# fp.summary$category <- ifelse((fp.summary$min_fatepotential < 0) &
#                                 (fp.summary$max_fatepotential > 0) &
#                                 (fp.summary$sd_fatepotential > quantile(fp.summary$sd_fatepotential, 0.75, na.rm = TRUE)),
#                               'plasticity', fp.summary$category)


lin.priming <- fp.summary %>%
  filter(category == 'priming') %>%
  pull(assigned_lineage)

lin.plasticity <- fp.summary %>%
  filter(category == 'plasticity') %>%
  pull(assigned_lineage)


high.fp.cells <- metadat.day10_DABTRAM %>%
  filter(fatepotential_DABTRAM_d10_w5 > 0) %>%
  pull(cell_id)

low.fp.cells <- metadat.day10_DABTRAM %>%
  filter(fatepotential_DABTRAM_d10_w5 < median(metadat.day10_DABTRAM$fatepotential_DABTRAM_d10_w5)) %>%
  pull(cell_id)


# ==============================================================================
# Differential tests on current cells
# ==============================================================================
priming_current <- rownames(metadat.day10_DABTRAM)[metadat.day10_DABTRAM$assigned_lineage %in% lin.priming]
plastic_current <- rownames(metadat.day10_DABTRAM)[metadat.day10_DABTRAM$assigned_lineage %in% lin.plasticity]

priming_current <- intersect(priming_current, high.fp.cells)
plastic_current <- intersect(plastic_current, high.fp.cells)




# gene
p <- ncol(saver.mat.d10)

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
gene_df_d10$gene <- rownames(gene_df_d10)
gene_df_d10$padj <- p.adjust(gene_df_d10$p.value, method = "BH")
gene_df_d10$neglog10_pval <- -log10(gene_df_d10$p.value)


ap1_related <- grep('JUN|FOS', gene_df_d10$gene, value = TRUE)

gene_df_d10$category <- ifelse(gene_df_d10$gene %in% jackpot, 'Jackpot', 'Other')
gene_df_d10$category <- ifelse(gene_df_d10$gene %in% isg.rs, 'ISG.RS', gene_df_d10$category)
gene_df_d10$category <- ifelse(gene_df_d10$gene %in% stress, 'Stress', gene_df_d10$category)
gene_df_d10$category <- factor(gene_df_d10$category, 
                              levels = c('Jackpot', 'ISG.RS', 'Stress', 'Other'))
gene_df_d10 <- gene_df_d10[order(gene_df_d10$category, decreasing = T), ]

thres <- gene_df_d10 %>% 
  filter(padj < 0.05) %>%
  summarise(sig_thres = min(neglog10_pval, na.rm = TRUE)) 

ggplot(gene_df_d10, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = c('Jackpot' = 'red', 'ISG.RS' = 'blue', 'Stress' = 'darkgreen', 'Other' = "lightgrey")) +
  # ggrepel::geom_text_repel(data = subset(gene_df_d10, category == 'Jackpot'), aes(label = gene),  color = 'red',
  #                          size = 3) +
  # ggrepel::geom_text_repel(data = subset(gene_df_d10, category == 'ISG.RS'), aes(label = gene),  color = 'blue',
  #                          size = 3) +
  # ggrepel::geom_text_repel(data = subset(gene_df_d10, category == 'Stress'), aes(label = gene),  color = 'darkgreen',
  #                          size = 3) +
  ggrepel::geom_text_repel(data = subset(gene_df_d10, category != 'Other'), aes(label = gene, color = category),
                           size = 3) +
  geom_hline(yintercept = thres$sig_thres, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change (Plastic <-> Priming)",
       y = "-log10(p-value)",
       title = "Differential expression of d10 clones") +
  xlim(-2, 2) +
  theme_Publication() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

write.csv(gene_df_d10, 
          file = paste0(output_dir, "DABTRAM_d10_priming_vs_plastic.csv"),
          row.names = FALSE)

ggsave(paste0(figure_dir, "Supp_DABTRAM_d10_DE_priming_vs_plastic.pdf"),
       width = 6, height = 5)

# ==============================================================================
# GSEA DABTRAM d10 priming vs plastic
# ==============================================================================

# Read signatures

ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(ref_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
threeca <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))

gene_df_d10 <- gene_df_d10[order(gene_df_d10$logfc, decreasing = TRUE), ]
gene_df_d10$gene <- rownames(gene_df_d10)

gsea_input.DABTRAM.d10 <- gene_df_d10$logfc
names(gsea_input.DABTRAM.d10) <- gene_df_d10$gene

set.seed(123)
GSEA_res.DABTRAM <- GSEA(geneList = gsea_input.DABTRAM.d10, 
                         TERM2GENE = threeca, 
                         pvalueCutoff = 0.2,
                         seed = T,
                         verbose = F)
GSEA_res.DABTRAM.df <- as_tibble(GSEA_res.DABTRAM@result)
GSEA_res.DABTRAM.df <- GSEA_res.DABTRAM.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]

GSEA_res.DABTRAM.df <- GSEA_res.DABTRAM.df[GSEA_res.DABTRAM.df$qvalue < 0.05, ]
GSEA_res.DABTRAM.df <- GSEA_res.DABTRAM.df[GSEA_res.DABTRAM.df$p.adjust < 0.05, ]

toplot.df <- GSEA_res.DABTRAM.df
toplot.df <- toplot.df[!toplot.df$ID %in% c('GAVISH_3CA_MALIGNANT_METAPROGRAM_26_NPC_GLIOMA', 'GAVISH_3CA_MALIGNANT_METAPROGRAM_34_PLATELET_ACTIVATION', 
                                      'GAVISH_3CA_MALIGNANT_METAPROGRAM_16_MES_GLIOMA', 'GAVISH_3CA_MALIGNANT_METAPROGRAM_28_OLIGO_NORMAL',
                                      'GAVISH_3CA_MALIGNANT_METAPROGRAM_40_PDAC_RELATED', 'GAVISH_3CA_MALIGNANT_METAPROGRAM_25_ASTROCYTES', 
                                      'GAVISH_3CA_MALIGNANT_METAPROGRAM_29_NPC_OPC'), ]
toplot.df$ID <- gsub('GAVISH_3CA_MALIGNANT_METAPROGRAM_', '', toplot.df$ID)
toplot.df$ID <- paste0('MP', toplot.df$ID)

toplot.df <- toplot.df[order(toplot.df$NES, decreasing = F), ]

toplot.df$ID <- gsub('_', ' ', toplot.df$ID)
toplot.df$ID <- factor(toplot.df$ID, levels = unique(toplot.df$ID))


ggplot(toplot.df, aes(x = NES, y = ID)) +
  geom_bar(stat = 'identity', color = 'lightgray') +
  theme_Publication()

# tf

p <- ncol(cv.mat.d10)

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
  
  logfc <- mean(cv.mat.d10[priming_current,j], na.rm = T) - mean(cv.mat.d10[plastic_current,j], na.rm = T)
  
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

ap1 <- grep('FOS|JUN', colnames(cv.mat.d10), value = TRUE)
ap1 <- ap1[!grepl('var.2', ap1)]  

ifn <- c('STAT1', 'STAT3', 'IRF3', 'STAT1::STAT2')
emt <- c('SNAI1', 'SNAI2', 'SNAI3', 'SNAI4', 'ZEB1', 'TWIST1')
misc <- c('SOX2', 'SOX10', 'CTCF')

tf_df_d10$category <- ifelse(tf_df_d10$tf %in% ap1, 'AP1', 'Other')
tf_df_d10$category <- ifelse(tf_df_d10$tf %in% ifn, 'IFN', tf_df_d10$category)
tf_df_d10$category <- ifelse(tf_df_d10$tf %in% emt, 'EMT', tf_df_d10$category)
tf_df_d10$category <- ifelse(tf_df_d10$tf %in% misc, 'SOX/CTCF', tf_df_d10$category)
tf_df_d10$category <- factor(tf_df_d10$category, levels = c('AP1', 'IFN', 'EMT', 'SOX/CTCF',  'Other'))
tf_df_d10 <- tf_df_d10[order(tf_df_d10$category, decreasing = T), ]

thres <- tf_df_d10 %>% 
  filter(padj < 0.05) %>%
  summarise(sig_thres = min(neglog10_pval, na.rm = TRUE)) 

tf_df_d10$logfc.plot <- ifelse(tf_df_d10$logfc > 2, 2,tf_df_d10$logfc)

ggplot(tf_df_d10, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'AP1'), aes(label = tf),  color = '#B22222', 
                           size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'IFN'), aes(label = tf),  color = 'blue', 
                           size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'EMT'), aes(label = tf),  color = '#059212', 
                           size = 3, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10, category == 'SOX/CTCF'), aes(label = tf),  color = 'orange', 
                           size = 3, show.legend = FALSE) +
  geom_hline(yintercept = thres$sig_thres, linetype = "dashed", color = "black") +
  scale_color_manual(values = c('Other' = "grey", 'AP1' = "red", 'IFN' = 'blue', 'EMT' = '#059212', 'SOX/CTCF' = 'orange')) +
  labs(x = "Log2 Fold Change (Plastic <-> Priming)",
       y = "-log10(p-value)",
       title = "Differential TF") +
  theme_Publication()

tf_df_d10.sig.up <- tf_df_d10[tf_df_d10$padj < 0.05 & tf_df_d10$logfc > 0, ]
tf_df_d10.sig.up2 <- tf_df_d10[tf_df_d10$tf %in% c('NFE2L1', 'MAF::NFE2', 'MAFK'), ]

tf_df_d10.sig.down <- tf_df_d10[tf_df_d10$padj < 0.05 & tf_df_d10$logfc < 0, ]

ggplot(tf_df_d10, aes(x = logfc, y = neglog10_pval)) +
  geom_point(color = 'grey') +
  geom_point(data = subset(tf_df_d10.sig.up, category %in% c('AP1', 'IFN')), color = '#82A4A8', size = 2) +
  geom_point(data = tf_df_d10.sig.up2, color = '#82A4A8', size = 2) +
  geom_point(data = subset(tf_df_d10.sig.down, category %in% c('SOX/CTCF', 'EMT')), color = '#C6983E', size = 2) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10.sig.up, category %in% c('AP1', 'IFN')), 
                           aes(label = tf),  color = '#82A4A8', 
                           size = 5, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = tf_df_d10.sig.up2, aes(label = tf),  color = '#82A4A8', 
                           size = 5, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(tf_df_d10.sig.down, category %in% c('SOX/CTCF', 'EMT')), aes(label = tf),  color = '#C6983E', 
                           size = 5, show.legend = FALSE) +
  geom_hline(yintercept = thres$sig_thres, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change (Plastic <-> Priming)",
       y = "-log10(p-value)") +
  theme_Publication()

ggsave('~/Downloads/diff_tf_d10.pdf', width = 6, height = 4)


cv.mat.d10.ap1 <- cv.mat.d10[, ap1]
cv.mat.d10.ap1 <- as.data.frame(cv.mat.d10.ap1)
colnames(cv.mat.d10.ap1) <- paste0("chromVar_", colnames(cv.mat.d10.ap1))
cv.mat.d10.ap1$cell_id <- rownames(cv.mat.d10.ap1)

gavish.mp <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/saver.GavishMP.UCellScores.csv')
rownames(gavish.mp) <- gavish.mp$cell_id
gavish.mp <- gavish.mp[rownames(gavish.mp) %in% rownames(cv.mat.d10.ap1), ]


comp.df <- merge(gavish.mp, cv.mat.d10.ap1, by = 'cell_id', all.x = TRUE)
comp.df <- comp.df[comp.df$cell_id %in% high.fp.cells, ]
comp.df$category <- ifelse(comp.df$cell_id %in% priming_current, 'Priming', 'Other')
comp.df$category <- ifelse(comp.df$cell_id %in% plastic_current, 'Plastic', comp.df$category )
ggplot(comp.df, aes(x = Stress, y = ISG.RS)) +
  geom_point(aes(color = category)) +
  stat_cor()

for(i in ap1) {
  res <- cor.test(comp.df[[paste0('chromVar_', i)]], comp.df$ISG.RS, method = "spearman", use = "pairwise.complete.obs")
  print(res$estimate)
}

ggplot(comp.df, aes(x = category, y = ISG.RS)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "chromVar FOSL2 Score") +
  theme_minimal()


# ==============================================================================
# Time course
# ==============================================================================
priming_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.priming]
plastic_progenitors <- rownames(metadat.day0)[metadat.day0$assigned_lineage %in% lin.plasticity]

cv.mat.d0.ap1 <- cv.mat.d0[, c(ap1, 'STAT1')]
cv.mat.d0.ap1 <- as.data.frame(cv.mat.d0.ap1)

cv.mat.d0.ap1$cell_id <- rownames(cv.mat.d0.ap1)
cv.mat.d0.ap1$category <- ifelse(cv.mat.d0.ap1$cell_id %in% priming_progenitors, 'Priming', 'Other')
cv.mat.d0.ap1$category <- ifelse(cv.mat.d0.ap1$cell_id %in% plastic_progenitors, 'Plastic', cv.mat.d0.ap1$category )
cv.mat.d0.ap1$category <- ifelse(cv.mat.d0.ap1$cell_id %in% low.adapting.fate.d0$cell_id, 'non-w5-progenitors', cv.mat.d0.ap1$category)
cv.mat.d0.ap1 <- cv.mat.d0.ap1[cv.mat.d0.ap1$category %in% c('Priming', 'Plastic', 'non-w5-progenitors'), ]
cv.mat.d0.ap1$time <- 'D0'

ggplot(cv.mat.d0.ap1, aes(x = category, y = log10(`STAT1`))) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  stat_compare_means(comparisons = list('non-w5-progenitors vs priming' = c('non-w5-progenitors', 'Priming'),
                                        'non-w5-progenitors vs plastic' = c('non-w5-progenitors', 'Plastic'),
                                        'Priming vs Plastic' = c('Priming', 'Plastic')),
                     label = "p.signif", method = "wilcox.test") +
  theme_Publication()

cv.mat.d10.ap1 <- cv.mat.d10[, c(ap1, 'STAT1')]
cv.mat.d10.ap1 <- as.data.frame(cv.mat.d10.ap1)

cv.mat.d10.ap1$cell_id <- rownames(cv.mat.d10.ap1)
cv.mat.d10.ap1$category <- ifelse(cv.mat.d10.ap1$cell_id %in% priming_current, 'Priming', 'Other')
cv.mat.d10.ap1$category <- ifelse(cv.mat.d10.ap1$cell_id %in% plastic_current, 'Plastic', cv.mat.d10.ap1$category )
cv.mat.d10.ap1$category <- ifelse(cv.mat.d10.ap1$cell_id %in% low.fp.cells, 'non-w5-progenitors', cv.mat.d10.ap1$category)
cv.mat.d10.ap1 <- cv.mat.d10.ap1[cv.mat.d10.ap1$category %in% c('Priming', 'Plastic', 'non-w5-progenitors'), ]

ggplot(cv.mat.d10.ap1, aes(x = category, y = log10(`STAT1`))) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  stat_compare_means(comparisons = list('non-w5-progenitors vs priming' = c('non-w5-progenitors', 'Priming'),
                                        'non-w5-progenitors vs plastic' = c('non-w5-progenitors', 'Plastic'),
                                        'Priming vs Plastic' = c('Priming', 'Plastic')),
                     label = "p.signif", method = "wilcox.test") +
  theme_Publication()

cv.mat.d10.ap1$time <- 'D10'


cv.mat.plot <- rbind(cv.mat.d0.ap1, cv.mat.d10.ap1)
# cv.mat.plot <- rbind(cv.mat.plot, cv.mat.w5.ap1)
ggplot(cv.mat.plot, aes(x = category, y = `JUNB`)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_violin(scale = 'width', alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.5) +
  facet_wrap(.~ time, scale = 'free_x') +
  stat_compare_means(comparisons = list('non-w5-progenitors vs priming' = c('non-w5-progenitors', 'Priming'),
                                        'non-w5-progenitors vs plastic' = c('non-w5-progenitors', 'Plastic'),
                                        'Priming vs Plastic' = c('Priming', 'Plastic')),
                     label = "p.signif", method = "wilcox.test") +
  xlab('') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




cv.mat.w5.ap1 <- cv.mat.w5[, c(ap1)]
cv.mat.w5.ap1 <- as.data.frame(cv.mat.w5.ap1)
cv.mat.w5.ap1 <- cv.mat.w5.ap1[rownames(metadat.week5_DABTRAM), ]
cv.mat.w5.ap1$cell_id <- rownames(cv.mat.w5.ap1)

priming_w5 <- rownames(metadat.week5_DABTRAM)[metadat.week5_DABTRAM$assigned_lineage %in% lin.priming]
plastic_w5 <- rownames(metadat.week5_DABTRAM)[metadat.week5_DABTRAM$assigned_lineage %in% lin.plasticity]

cv.mat.w5.ap1$category <- ifelse(cv.mat.w5.ap1$cell_id %in% priming_w5, 'Priming', 'Other')
cv.mat.w5.ap1$category <- ifelse(cv.mat.w5.ap1$cell_id %in% plastic_w5, 'Plastic', cv.mat.w5.ap1$category )
cv.mat.w5.ap1$time <- 'W5'

ggplot(cv.mat.w5.ap1, aes(x = category, y = FOSL1)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  stat_compare_means(comparisons = list('Priming vs Plastic' = c('Priming', 'Plastic')),
                     label = "p.signif", method = "wilcox.test") +
  theme_Publication()

ft_umap <- all_data_ft_DABTRAM_umap@cell.embeddings
ft_umap <- as.data.frame(ft_umap)
ft_umap$cell_id <- rownames(ft_umap)
ft_umap$category <- ifelse(ft_umap$cell_id %in% priming_w5, 'Priming', 'Other')
ft_umap$category <- ifelse(ft_umap$cell_id %in% plastic_w5, 'Plastic', ft_umap$category )
ft_umap$category <- factor(ft_umap$category, levels = c('Priming', 'Plastic', 'Other'))
ft_umap <- ft_umap[order(ft_umap$category, decreasing = T),]
                           
ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = category)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c('Priming' = 'blue', 'Plastic' = 'red', 'Other' = '#999999')) +
  theme_Publication() +
  labs(color = '') +
  theme(legend.position = 'bottom')

ft_umap.w5 <- ft_umap[ft_umap$cell_id %in% rownames(metadat.week5_DABTRAM), ]

library(hdrcde)
library(ggdensity)
library(RColorBrewer)

p.plastic <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = 'lightgray') +
  geom_hdr(data = subset(ft_umap, category == 'Plastic'), 
           aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.9, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  ggtitle('Plastic clones') +
  theme_Publication()

p.priming <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = 'lightgray') +
  geom_hdr(data = subset(ft_umap, category == 'Priming'),
           aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.9, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  ggtitle('Priming clones') +
  theme_Publication() 

p.rest <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = 'lightgray') +
  geom_hdr(data = subset(ft_umap.w5, category == 'Other'),
           aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.9, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  ggtitle('Other clones') +
  theme_Publication() 

ggarrange(p.priming, p.plastic,
          ncol = 2, nrow = 1, 
          common.legend = T, legend = 'right')
ggsave('~/Downloads/priming_plasticity_umap_hdr.pdf', width = 6, height = 2.8)





# plot day 10
ft_umap <- all_data_ft_DABTRAM_umap@cell.embeddings
ft_umap <- as.data.frame(ft_umap)
ft_umap$cell_id <- rownames(ft_umap)
ft_umap$category <- ifelse(ft_umap$cell_id %in% priming_current, 'Priming', 'Other')
ft_umap$category <- ifelse(ft_umap$cell_id %in% plastic_current, 'Plastic', ft_umap$category )
ft_umap$category <- factor(ft_umap$category, levels = c('Priming', 'Plastic', 'Other'))
ft_umap <- ft_umap[order(ft_umap$category, decreasing = T),]


p.priming.d10 <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = category)) +
  geom_point(aes(size = category)) +
  scale_color_manual(values = c('Priming' = '#82A4A8', 'Plastic' = '#E8E8E8', 'Other' = '#E8E8E8')) +
  scale_size_manual(values =  c('Priming' = 0.8, 'Plastic' = 0.1, 'Other' = 0.1), guide = 'none') +
  theme_Publication() +
  labs(color = '', title = 'Cells in priming clones (high growth potential)') +
  theme(legend.position = 'right')

ft_umap$category <- factor(ft_umap$category, levels = c('Plastic','Priming', 'Other'))
ft_umap <- ft_umap[order(ft_umap$category, decreasing = T),]
p.plastic.d10 <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = category)) +
  geom_point(aes(size = category)) +
  scale_color_manual(values = c('Priming' = '#E8E8E8', 'Plastic' = '#C6983E', 'Other' = '#E8E8E8')) +
  scale_size_manual(values =  c('Priming' = 0.1, 'Plastic' = 0.8, 'Other' = 0.1), guide = 'none') +
  theme_Publication() +
  labs(color = '', title = 'Cells in plastic clones (high growth potential)') +
  theme(legend.position = 'right')

ggarrange(p.priming.d10, p.plastic.d10,
          ncol = 1, nrow = 2, legend = 'right')

ggsave(paste0(figure_dir, 'Supp_priming_plasticity_umap_d10.pdf'), width = 4, height = 6)

ft_umap.w5 <- ft_umap[ft_umap$cell_id %in% rownames(metadat.week5_DABTRAM), ]


p.plastic <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = 'lightgray') +
  geom_hdr(data = subset(ft_umap, category == 'Plastic'), 
           aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.9, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  ggtitle('Plastic clones') +
  theme_Publication()

p.priming <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = 'lightgray') +
  geom_hdr(data = subset(ft_umap, category == 'Priming'),
           aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.9, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  ggtitle('Priming clones') +
  theme_Publication() 

p.rest <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = 'lightgray') +
  geom_hdr(data = subset(ft_umap.w5, category == 'Other'),
           aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.9, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "Purples")) +
  ggtitle('Other clones') +
  theme_Publication() 