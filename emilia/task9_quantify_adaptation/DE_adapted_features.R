rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(hdrcde)
library(ggdensity)
library(RColorBrewer)
library(circlize)
library(ggExtra)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
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

fp_gini <- read.csv(paste0(output_dir, 'fate_potential_gini_index.csv'))

ada_index.DABTRAM <- read.csv(paste0(output_dir, 'adaptation_index_DABTRAM_day10_week5.csv'))

ada_index.DABTRAM.minmax <- read.csv(paste0(output_dir, 'adaptation_index_DABTRAM_day10_week5_MinMax.csv'))

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVAR_day0.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

all_data.day10_DABTRAM <- subset(all_data, dataset == 'day10_DABTRAM')
all_data.week5_DABTRAM <- subset(all_data, dataset == 'week5_DABTRAM')

# =============================================================================
# Wrangle
# =============================================================================

ada_index.DABTRAM$adaptation.index <- (ada_index.DABTRAM$dist / mean(ada_index.DABTRAM.minmax$dist.min) ) 

lin.high_adapt <- ada_index.DABTRAM[ada_index.DABTRAM$adaptation.index > quantile(ada_index.DABTRAM$adaptation.index, 0.75), 'lineage']


adapting_progenitors <- which(all_data.day10_DABTRAM$assigned_lineage %in% lin.high_adapt)
adapted_cells <- which(all_data.week5_DABTRAM$assigned_lineage %in% lin.high_adapt)

umap <- all_data_ft_DABTRAM_umap@cell.embeddings
umap <- as.data.frame(umap)
umap$cell_id <- rownames(umap)

metadat.day10_DABTRAM <- all_data.day10_DABTRAM@meta.data
metadat.week5_DABTRAM <- all_data.week5_DABTRAM@meta.data

adapting_progenitors_cell_id <- rownames(metadat.day10_DABTRAM[metadat.day10_DABTRAM$assigned_lineage %in% lin.high_adapt, ])
adapted_cells_cell_id <- rownames(metadat.week5_DABTRAM[metadat.week5_DABTRAM$assigned_lineage %in% lin.high_adapt, ])


umap$category <- ifelse(umap$cell_id %in% adapting_progenitors_cell_id, 'adapting_progenitors_d10', 'Other')
umap$category <- ifelse(umap$cell_id %in% adapted_cells_cell_id, 'adapted_cells_w5', umap$category)
umap <- umap[order(umap$category, decreasing = T), ]

ggplot(umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = category)) +
  geom_point(aes(size = category)) +
  scale_color_manual(values = c('adapting_progenitors_d10' = '#9D85BE', 'adapted_cells_w5' = '#623594', 'Other' = 'lightgray')) +
  scale_size_manual(values = c('adapting_progenitors_d10' = 1, 'adapted_cells_w5' = 1, 'Other' = 0.1)) +
  labs(title = "UMAP of Adapting Progenitors and Adapted Cells",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(paste0(figure_dir, 'Supp_UMAP_adapting_progenitors_adapted_cells_DABTRAM.pdf'), width = 3, height = 3.5)



saver.mat.day10 <- all_data.day10_DABTRAM[["Saver"]]@data
saver.mat.day10 <- t(saver.mat.day10)

saver.mat.week5 <- all_data.week5_DABTRAM[["Saver"]]@data
saver.mat.week5 <- t(saver.mat.week5)

p <- ncol(saver.mat.day10)

wilcox_results <- sapply(1:p, function(j){
  tmp <- stats::wilcox.test(
    x = saver.mat.week5[adapted_cells,j],
    y = saver.mat.day10[adapting_progenitors,j]
  )
  logfc <- log2(mean(saver.mat.week5[adapted_cells,j])) - log2(mean(saver.mat.day10[adapting_progenitors,j]))
  
  c(logfc = logfc,
    p.value = tmp$p.value)
})
wilcox_results <- t(wilcox_results)

gene_df <- as.data.frame(wilcox_results)
rownames(gene_df) <- colnames(saver.mat.day10)
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

write.csv(gene_df, paste0(output_dir, 'day0_DE_adapted_genes_D10_W5_DABTRAM.csv'))



# gene_df <- read.csv('~/Downloads/de_d10_w5_adapting_cells.csv', row.names = 1)

gene_df$gene <- rownames(gene_df)
gene_df$category <- ifelse(gene_df$gene %in% jackpot, 'Jackpot', 'Other')
gene_df$category <- ifelse(gene_df$gene %in% isg.rs, 'ISG.RS', gene_df$category)
thres <- min(gene_df$neglog10_pval[gene_df$padj < 0.05], na.rm = T)
gene_df <- gene_df[order(gene_df$category, decreasing = T), ]

ggplot(gene_df, aes(x = logfc, y = neglog10_pval)) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = c("blue", "red", "#E8E8E8")) +
  ggrepel::geom_text_repel(data = subset(gene_df, category != 'Other'),
                           aes(label = gene, color = category), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = thres, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change (Adapting clones at d10 vs Adapted clones at w5)",
       y = "-log10(p-value)") +
  theme_Publication() +
  theme(legend.position = "none")

ggsave(paste0(figure_dir, 'Supp_DE_adapted_genes_D10_W5_DABTRAM.pdf'), width = 5, height = 5)
