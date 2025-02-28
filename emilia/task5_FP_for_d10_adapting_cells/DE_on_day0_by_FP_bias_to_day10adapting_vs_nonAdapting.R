rm(list = ls())

library(Seurat)
library(tidyverse)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

treatment <- 'DABTRAM'

remove_unassigned_cells <- TRUE

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
# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

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


nonAdaptingFP <- readRDS(paste0(out_dir, 'final_fit_d0_d10_Nonadpating_thres_0_', treatment, '.rds'))
nonAdaptingFP <- nonAdaptingFP[["cell_imputed_score"]]

adaptingFP <- readRDS(paste0(out_dir, 'final_fit_d0_d10_adpating_thres_0_', treatment, '.rds'))
adaptingFP <- adaptingFP[["cell_imputed_score"]]

# =============================================================================
# Define winner and loser
# =============================================================================
adaptingFP <- as.data.frame(adaptingFP)
nonAdaptingFP <- as.data.frame(nonAdaptingFP)

fate_bias <- merge(adaptingFP, nonAdaptingFP, by = 'row.names')
fate_bias$bias <- 10**(fate_bias$adaptingFP) / (10**(fate_bias$adaptingFP) + 10**(fate_bias$nonAdaptingFP))
quantile(fate_bias$bias)
hist(fate_bias$bias, breaks = 100)

# fate_bias$isWinning <- ifelse(fate_bias$bias >= 0.043111934, 'adapting', NA)
# fate_bias$isWinning <- ifelse(fate_bias$bias <= 0.024694730, 'nonAdapting', df$isWinning)
fate_bias$isWinning <- ifelse(fate_bias$bias >= 0.5, 'adapting', 'nonAdapting')
table(fate_bias$isWinning)
colnames(fate_bias)[1] <- 'cell_id' 

cell_id.adapting <- fate_bias$cell_id[fate_bias$isWinning == 'adapting']
cell_id.nonAdapting <- fate_bias$cell_id[fate_bias$isWinning == 'nonAdapting']

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

metadat$isWinning <- ifelse(metadat$cell_id %in% cell_id.adapting, 'adapting', NA)
metadat$isWinning <- ifelse(metadat$cell_id %in% cell_id.nonAdapting, 'nonAdapting', metadat$isWinning)

all_data <- AddMetaData(all_data, metadat)

# =============================================================================
# subset to day0
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')
metadat_day0 <- all_data_day0@meta.data
metadat_day0$cell_id <- rownames(metadat_day0)
table(metadat_day0$isWinning)

saver.day0 <- all_data_day0@assays[["saver"]]@scale.data

pvalue_list <- lapply(all_data_day0[["saver"]]@var.features, function(gene){
  x_vec <- all_data_day0[["saver"]]@scale.data[gene, which(all_data_day0$isWinning == "adapting")]
  y_vec <- all_data_day0[["saver"]]@scale.data[gene, which(all_data_day0$isWinning == "nonAdapting")]
  
  if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
  
  test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
  list(stat = mean(x_vec) - mean(y_vec), pvalue = test_res$p.value)
})
names(pvalue_list) <- all_data_day0[["saver"]]@var.features
pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
pvalue_vec <- sapply(pvalue_list, function(x){x$pvalue})
names(pvalue_vec) <- names(pvalue_list)
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

diff_vec <- sapply(pvalue_list, function(x){x$stat})
logpval_vec <- sapply(pvalue_list, function(x){-log10(x$pvalue)})

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 5)]] <- "3"
labeling_vec[intersect(keygenes, names(labeling_vec))] <- "2"
table(labeling_vec)


df <- data.frame(difference = diff_vec,
                 log10pval = logpval_vec,
                 name = names(logpval_vec),
                 labeling = labeling_vec)

df$labeling <- factor(df$labeling, levels = c("1", "3", "2"))
df <- df[order(df$labeling), ]
df$log10pval <- ifelse(df$log10pval > 20, 20, df$log10pval)
df.top <- df[order(df$difference, decreasing = TRUE), ][1:10, ]
df.bottom <- df[order(df$difference, decreasing = FALSE), ][1:10, ]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = difference, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::geom_point(data = subset(df, labeling %in% c("2", "3")), 
                               ggplot2::aes(color = labeling), size = 2)
p1 <- p1 + ggplot2::scale_color_manual(values = c("1" = "black", "2" = "red", "3" = "blue"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                               color = "gray", linewidth=2)
p1 <- p1 + ggplot2::geom_vline(xintercept=0, linetype="dashed",
                               color = "gray", linewidth=2)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 DE based on fate bias for week5 adapting")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)")
p1 <- p1 + ggplot2::theme_bw()
p1 <- p1 + Seurat::NoLegend()
p1

ggplot(df, aes(x = difference, y = log10pval)) +
  geom_point(color = 'gray') +
  geom_point(data = df.top, color = 'red', size = 2) +
  geom_point(data = df.bottom, color = 'blue', size = 2) +
  geom_point(data = df[df$labeling == "2" & !df$name %in% df.top$name, ], color = '#7C00FE', size = 2) +
  ggrepel::geom_text_repel(data = df.top, aes(label = name), color = 'red',
                           box.padding = 0.5, point.padding = 1.6, max.overlaps = 50) +
  ggrepel::geom_text_repel(data = df.bottom, aes(label = name), color = 'blue',
                           box.padding = 0.5, point.padding = 1.6, max.overlaps = 50) +
  geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
             color = "black", linewidth=2) +
  geom_vline(xintercept=0, linetype="dashed",
             color = "black", linewidth=2) +
  ggplot2::ggtitle(paste0("DABTRAM Day0 DE based on fate bias for week5 adapting")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)") +
  theme_bw()

x_vec <- as.data.frame(all_data_day0[["saver"]]@scale.data['FN1', ])
colnames(x_vec) <- 'FN1_saver'
x_vec$cell_id <- rownames(x_vec)
metadat_day0 <- merge(metadat_day0, x_vec, by = 'cell_id')
metadat_day0 <- merge(metadat_day0, fate_bias, by = c('cell_id', 'isWinning'))

ggplot(metadat_day0, aes(x = isWinning, y = FN1_saver)) +
  geom_violin(aes(fill = isWinning), scale = 'width') +
  geom_jitter(width = 0.1) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  theme_bw()

ggplot(metadat_day0, aes(x = bias, y = NQO1_saver)) +
  geom_point(aes(color = isWinning))
