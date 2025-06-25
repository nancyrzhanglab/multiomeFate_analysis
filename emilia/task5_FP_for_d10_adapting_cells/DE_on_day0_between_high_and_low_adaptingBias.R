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

figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V5/Figure5/'

remove_unassigned_cells <- TRUE

jackpot = c("SOX10", 'MITF', 'FN1', 'AXL', 'EGFR', 'NT5E',
            'C1S', 'FRZB', 'SERPINB2', 'SERPINE1', 'NGFR',
            'SERPINE2', 'NDRG1', 'FEZF1', 'EGR3', 'VGF',
            'WNT5A', 'POSTN', 'PDGFRB', 'NRG1', 'VEGFC', 'FOSL1',
            'RUNX2', 'LOXL2', 'JUN', 'PDGFRC', 'CD44', 'ID3')

isg.rs = c('IFI27', 'IRF7','USP18', 'BST2', 'CXCL10', 'DDX60',
           'HERC6', 'HLA-B', 'HLA-G', 'IFI35','IFI44','IFI44L',
           'IFIT1', 'IFIT3', 'ISG15', 'LGALS3BP', 'LY6E', 'MX1',
           'MX2', 'OAS3', 'OASL', 'PLSCR1', 'STAT1', 'TRIM14',
           'HSD17B1', 'OAS1', 'CA2', 'CCNA1', 'CXCL1', 'GALC',
           'IFI6', 'IFITM1', 'LAMP3', 'MCL1', 'ROBO1', 'SLC6A15',
           'THBS1', 'TIMP3', 'DDX58', 'IFIH1')

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

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVAR_day0.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[["chromVar_day0"]] <- all_data_chromVar_day0

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

# fate bias from d0 to week5 adapting
df.bias.DABTRAM <- read.csv(paste0(fatebias_out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))

# =============================================================================
# Wrangle data
# =============================================================================
rownames(df.bias.DABTRAM) <- df.bias.DABTRAM$cell_id
df.bias.DABTRAM$fate_bias_class <- ifelse(df.bias.DABTRAM$bias > median(df.bias.DABTRAM$bias), 'Adapting_Biased', 'Non_Adapting_Biased')

all_data.day0 <- AddMetaData(all_data.day0, 
                             metadata = df.bias.DABTRAM[,c('bias', 'fate_bias_class')])

adapting_progenitors <- which(all_data.day0$fate_bias_class == 'Adapting_Biased')
nonadapting_progenitors <- which(all_data.day0$fate_bias_class == 'Non_Adapting_Biased')

saver.mat <- all_data.day0[["Saver"]]@data
saver.mat <- t(saver.mat)

p <- ncol(saver.mat)

wilcox_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = saver.mat[adapting_progenitors,j],
    y = saver.mat[nonadapting_progenitors,j]
  )
  logfc <- log2(mean(saver.mat[adapting_progenitors,j])) - log2(mean(saver.mat[nonadapting_progenitors,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results <- t(wilcox_results)

gene_df <- as.data.frame(wilcox_results)
rownames(gene_df) <- colnames(saver.mat)
colnames(gene_df) <- c("logfc", "p.value")

gene_df$padj <- p.adjust(gene_df$p.value, method = "BH")
gene_df$neglog10_pval <- -log10(gene_df$p.value)

gene_df <- read.csv('~/Downloads/de.csv', row.names = 1, stringsAsFactors = FALSE)
gene_df$gene <- rownames(gene_df)
gene_df$category <- ifelse(gene_df$gene %in% jackpot, 'Jackpot', 'Other')
gene_df$category <- ifelse(gene_df$gene %in% isg.rs, 'ISG.RS', gene_df$category)
thres <- min(gene_df$neglog10_pval[gene_df$padj < 0.05], na.rm = T)
gene_df <- gene_df[order(gene_df$category, decreasing = T), ]
gene_df$neglog10_pval <- ifelse(gene_df$neglog10_pval > 40, 40, gene_df$neglog10_pval)

ggplot(gene_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = c("blue", "red", "#E8E8E8")) +
  # ggrepel::geom_text_repel(data = subset(gene_df, category == 'Jackpot'), color = 'red',
  #                 aes(label = gene), size = 3, max.overlaps = 10) +
  # ggrepel::geom_text_repel(data = subset(gene_df, category == 'ISG.RS'), color = 'blue',
  #                          aes(label = gene), size = 3, max.overlaps = 10) +
  ggrepel::geom_text_repel(data = subset(gene_df, category != 'Other'),
                           aes(label = gene, color = category), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = thres, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change\n(Week 5 estimate low <->Week 5 estimate high cells at d0)",
       y = "-log10(p-value)") +
  theme_Publication() +
  theme(legend.position = "right", legend.direction = 'vertical')

ggsave(paste0(figure_dir, 'Supp_day0_DE_selected_genes.pdf'), width = 6, height = 5)
write.csv(gene_df, paste0(output_dir, 'day0_DE_selected_genes.csv'))

cv.mat <- all_data.day0[["chromVar_day0"]]@data
cv.mat <- t(cv.mat)

p <- ncol(cv.mat)

wilcox_results_tf <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = cv.mat[adapting_progenitors,j],
    y = cv.mat[nonadapting_progenitors,j]
  )
  logfc <- mean(cv.mat[adapting_progenitors,j], na.rm = T) - mean(cv.mat[nonadapting_progenitors,j], na.rm = T)
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results_tf <- t(wilcox_results_tf)

tf_df <- as.data.frame(wilcox_results_tf)
rownames(tf_df) <- colnames(cv.mat)
colnames(tf_df) <- c("logfc", "p.value")

tf_df$padj <- p.adjust(tf_df$p.value, method = "BH")
tf_df$neglog10_pval <- -log10(tf_df$p.value)
tf_df$tf <- rownames(tf_df)

tfs <- rownames(tf_df)
ap1 <- tfs[grepl("FOS|JUN", tf_df$tf, ignore.case = TRUE)]
ap1 <- ap1[!grepl('var', ap1)]

isg <- c('STAT1', 'IRF3', 'STAT3')

tf_df$category <- ifelse(tf_df$tf %in% ap1, 'AP1', 'Other')
tf_df$category <- ifelse(tf_df$tf %in% isg, 'IFN', tf_df$category)

thres <- min(tf_df$neglog10_pval[tf_df$padj < 0.05], na.rm = T)
tf_df <- tf_df[order(tf_df$category, decreasing = T), ]
# tf_df$neglog10_pval <- ifelse(tf_df$neglog10_pval > 40, 40, tf_df$neglog10_pval)
ggplot(tf_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = c("red", "blue", "gray")) +
  ggrepel::geom_text_repel(data = subset(tf_df, category != 'Other'),
                           aes(label = tf, color = category), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = thres, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change (Adapting vs Non-Adapting Progenitors at d0)",
       y = "-log10(p-value)") +
  theme_Publication() +
  theme(legend.position = "none")
ggsave(paste0(figure_dir, 'Supp_day0_DE_selected_TFs.pdf'), width = 6, height = 5)
write.csv(tf_df, paste0(output_dir, 'day0_DE_selected_TFs.csv'))

  
# ==============================================================================
# Read signatures
# ==============================================================================
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(ref_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
threeca <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))

# ==============================================================================
# GSEA
# ==============================================================================

getAndSortTable <- function(df) {
  df <- as.data.frame(df) %>% drop_na()
  df <- df[df$padj < 0.05, ]
  df$gene <- rownames(df)
  
  df$logfc <- as.numeric(df$logfc)
  df <- df[order(df$logfc, decreasing = TRUE),]
  
  return(df)
}

gene_df <- getAndSortTable(gene_df)

# GSEA DABTRAM
gsea_input.DABTRAM <- gene_df$logfc
names(gsea_input.DABTRAM) <- gene_df$gene

set.seed(123)
GSEA_res.DABTRAM <- GSEA(geneList = gsea_input.DABTRAM, 
                               TERM2GENE = threeca, 
                               pvalueCutoff = 0.2,
                               seed = T,
                               verbose = F)
GSEA_res.DABTRAM.df <- as_tibble(GSEA_res.DABTRAM@result)
GSEA_res.DABTRAM.df <- GSEA_res.DABTRAM.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.DABTRAM.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM', 'setSize')
