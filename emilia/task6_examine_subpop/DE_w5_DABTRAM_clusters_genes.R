rm(list=ls())
library(Seurat)
library(clusterProfiler)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

keygenes <- list(
  # jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
  #                  "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
  #                  "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
  #                  "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
  #                  "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
  DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44"))
)

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


# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

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

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# ==============================================================================
# Get high fp cells from d10 and cells in w5
# ==============================================================================

# d0
metadat.d0 <- metadat %>% filter(dataset == 'day0')

# high fp d10
fp.dabtram.d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(fp.dabtram.d10_w5) <- c('fatepotential_DABTRAM_d10_w5')
fp.dabtram.d10_w5$cell_id <- rownames(fp.dabtram.d10_w5)
high_fp_d10 <- fp.dabtram.d10_w5 %>% filter(fatepotential_DABTRAM_d10_w5 > 0)

high_fp_d10 <- merge(high_fp_d10, metadat, by=c('fatepotential_DABTRAM_d10_w5', 'cell_id'))

low_fp_d10 <- fp.dabtram.d10_w5 %>% filter(fatepotential_DABTRAM_d10_w5 < 0)
low_fp_d10 <- merge(low_fp_d10, metadat, by=c('fatepotential_DABTRAM_d10_w5', 'cell_id'))

# week5 
metadat.week5_DABTRAM <- metadat %>% filter(dataset == 'week5_DABTRAM')

metadat.DABTRAM <- read.csv(paste0(data_dir, 'metadat.all_data_dabtram_recluster.csv'), row.names = 1)
metadat.week5_DABTRAM <- metadat.DABTRAM[metadat.DABTRAM$dataset == 'week5_DABTRAM', ]
metadat.week5_DABTRAM.clust0 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 0)
metadat.week5_DABTRAM.clust3 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 3)

# ==============================================================================
# DE
# ==============================================================================

all_data.week5_DABTRAM <- subset(all_data, dataset == 'week5_DABTRAM')
# metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM', ]

pvalue_list <- lapply(all_data.week5_DABTRAM[["saver"]]@var.features, function(gene){
  x_vec <- all_data.week5_DABTRAM[["saver"]]@scale.data[gene, rownames(metadat.week5_DABTRAM.clust0)]
  y_vec <- all_data.week5_DABTRAM[["saver"]]@scale.data[gene, rownames(metadat.week5_DABTRAM.clust3)]
  
  if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec) - mean(y_vec), pvalue = test_res$p.value)
})
names(pvalue_list) <- all_data.week5_DABTRAM[["saver"]]@var.features
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})


df <- data.frame(difference = diff_vec,
                 pvalue_adj =  pvalue_vec,
                 log10pval = logpval_vec,
                 name = names(logpval_vec))


# ==============================================================================
# Plot
# ==============================================================================
df$log10pval <- ifelse(df$log10pval > 100, 100, df$log10pval)
ggplot(df, aes(x = difference, y = log10pval)) +
  geom_point() +
  geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                               color = "gray", linewidth=2) +
  geom_vline(xintercept=0, linetype="dashed",
                               color = "gray", linewidth=2) +
  geom_point(data = subset(df, name %in% keygenes), color = "red", size = 3) +
  ggrepel::geom_text_repel(data = subset(df, name %in% keygenes),
                           aes(label = name),
                           box.padding = unit(0.3, 'lines'),
                           point.padding = unit(1.6, 'lines'),
                           max.overlaps = 50, color = "red") +
  ylim(0, 110) +
  theme_Publication() 

# ==============================================================================
# GSEA
# ==============================================================================
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
threeCA <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))
hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))


# GSEA
df <- df[order(df$difference, decreasing = TRUE),]

gsea_input <- df$difference
names(gsea_input) <- df$name

set.seed(123)
GSEA_res_dabtram <- GSEA(geneList = gsea_input, 
                         TERM2GENE = hallmark, 
                         pvalueCutoff = 0.1,
                         seed = T,
                         verbose = F)
GSEA_res_dabtram.df <- as_tibble(GSEA_res_dabtram@result)
GSEA_res_dabtram.df <- GSEA_res_dabtram.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')

# ==============================================================================
# Plot
# ==============================================================================

ggplot(GSEA_res_dabtram.df, aes(x = NES.DABTRAM, y = reorder(ID, NES.DABTRAM))) +
  geom_bar(stat = "identity", width = 0.5, fill = 'gray') +
  theme_Publication()

