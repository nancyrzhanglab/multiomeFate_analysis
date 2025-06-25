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

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# =============================================================================
# Wrangle data
# =============================================================================
day0 <- which(all_data$dataset == 'day0')
week5_dabtram <- which(all_data$dataset == 'week5_DABTRAM')

saver.mat <- all_data[["Saver"]]@data
saver.mat <- t(saver.mat)

p <- ncol(saver.mat)

wilcox_results_d0_w5 <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = saver.mat[day0,j],
    y = saver.mat[week5_dabtram,j]
  )
  logfc <- log2(mean(saver.mat[week5_dabtram,j])) - log2(mean(saver.mat[day0,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results_d0_w5 <- t(wilcox_results_d0_w5)

gene_df_d0_w5 <- as.data.frame(wilcox_results_d0_w5)
rownames(gene_df_d0_w5) <- colnames(saver.mat)
colnames(gene_df_d0_w5) <- c("logfc", "p.value")
# gene_df_d0_w5$logfc <- (-1) * gene_df_d0_w5$logfc

gene_df_d0_w5$padj <- p.adjust(gene_df_d0_w5$p.value, method = "BH")
gene_df_d0_w5$neglog10_pval <- -log10(gene_df_d0_w5$p.value)

ggplot(gene_df_d0_w5, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = padj < 0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change (Adapting vs Non-Adapting)",
       y = "-log10(p-value)",
       title = "Differential Expression of Adapting Progenitors") +
  theme_minimal() +
  theme(legend.position = "none")

# write.csv(gene_df_d0_w5, '~/Downloads/de_d0_w5.csv')

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

gene_df_d0_w5 <- getAndSortTable(gene_df_d0_w5)

# GSEA DABTRAM
gsea_input.DABTRAM <- gene_df_d0_w5$logfc
names(gsea_input.DABTRAM) <- gene_df_d0_w5$gene

set.seed(123)
GSEA_res.DABTRAM <- GSEA(geneList = gsea_input.DABTRAM, 
                         TERM2GENE = threeca, 
                         pvalueCutoff = 0.2,
                         seed = T,
                         verbose = F)
GSEA_res.DABTRAM.df <- as_tibble(GSEA_res.DABTRAM@result)
GSEA_res.DABTRAM.df <- GSEA_res.DABTRAM.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.DABTRAM.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM', 'setSize')
