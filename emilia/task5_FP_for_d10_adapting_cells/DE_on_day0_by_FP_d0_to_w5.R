rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE
# 
# keygenes <- list(
#   jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
#                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
#                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
#                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
#                    "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
#   DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
#                    "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
#                    "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
#                    "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
#                    "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
#   COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                  "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                  "CAV1"))
# )

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

# keygenes <- list(
#   jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
#                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
#                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
#                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1")),
#   COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                  "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                  "CAV1"))
# )

# keygenes <- list(
#   jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
#                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
#                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
#                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
#                    "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"))
# )
keygenes <- unlist(keygenes)

treatment <- 'DABTRAM'

# =============================================================================
# Read data
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

final_fit <- readRDS('~/Downloads/final_fit_d0_w5_DABTRAM_scaled.rds')
fp_adapting <- as.data.frame(final_fit[["cell_imputed_score"]])

# fp_adapting <- readRDS('~/Downloads/cell_imputed_score2_d0_d10_predicted_w5_w_d10_size.DABTRAM.rds')
# fp_adapting <- as.data.frame(fp_adapting)
colnames(fp_adapting) <- c('fp_d0_w5')
fp_adapting$cell_id <- rownames(fp_adapting)
# =============================================================================
# subset to day0
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')
metadat_day0 <- all_data_day0@meta.data
metadat_day0$cell_id <- rownames(metadat_day0)

# =============================================================================
# Define winner and loser
# =============================================================================

quantile(fp_adapting$fp_d0_w5)
thres_up <- quantile(fp_adapting$fp_d0_w5, 0.75)
thres_bottom <- quantile(fp_adapting$fp_d0_w5, 0.25)

metadat_day0 <- merge(metadat_day0, fp_adapting, by = 'cell_id')

metadat_day0$isWinning <- ifelse(metadat_day0$fp_d0_w5 >= thres_up, 'Winning', NA)
metadat_day0$isWinning <- ifelse(metadat_day0$fp_d0_w5 <  thres_bottom, 'Losing', metadat_day0$isWinning)

all_data_day0 <- AddMetaData(all_data_day0, metadat_day0)
table(metadat_day0$isWinning)

saver.day0 <- all_data_day0@assays[["saver"]]@scale.data

pvalue_list <- lapply(all_data_day0[["saver"]]@var.features, function(gene){
  x_vec <- all_data_day0[["saver"]]@scale.data[gene, which(all_data_day0$isWinning == "Winning")]
  y_vec <- all_data_day0[["saver"]]@scale.data[gene, which(all_data_day0$isWinning == "Losing")]
  
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

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = difference, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::geom_point(data = subset(df, labeling %in% c("2", "3")), 
                               ggplot2::aes(color = labeling), size = 2)
p1 <- p1 + ggplot2::scale_color_manual(values = c("1" = "black", "2" = "red", "3" = "blue"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2", "3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                               color = "gray", linewidth=2)
p1 <- p1 + ggplot2::geom_vline(xintercept=0, linetype="dashed", 
                               color = "gray", linewidth=2)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 DE based on presence or absence on Week5")) +
  ggplot2::xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + ggplot2::ylab("T test p-value (-Log10)")
p1 <- p1 + ggplot2::theme_bw()
p1 <- p1 + Seurat::NoLegend()
p1


x_vec <- all_data_day0[["saver"]]@scale.data['CBX5', which(all_data_day0$isWinning == "Winning")]
y_vec <- all_data_day0[["saver"]]@scale.data['CBX5', which(all_data_day0$isWinning == "Losing")]
test_res <- stats::wilcox.test(x = x_vec, y = y_vec)
test_res

x_vec <- as.data.frame(all_data_day0[["saver"]]@scale.data['VGF', ])
colnames(x_vec) <- 'VGF_saver'
x_vec$cell_id <- rownames(x_vec)
metadat_day0 <- merge(metadat_day0, x_vec, by = 'cell_id')

ggplot(metadat_day0, aes(x = VGF_saver, y = fp_d0_w5)) +
  geom_point(aes(color = isWinning)) +
  geom_smooth(method = "lm") +
  stat_cor() +
  theme_bw()


x_vec <- as.data.frame(all_data_day0[["saver"]]@scale.data['FOSL1', ])
colnames(x_vec) <- 'FOSL1_saver'
x_vec$cell_id <- rownames(x_vec)
metadat_day0 <- merge(metadat_day0, x_vec, by = 'cell_id')

ggplot(metadat_day0, aes(x = isWinning, y = FOSL1_saver)) +
  geom_violin(aes(color = isWinning)) +
  geom_jitter() +
  geom_boxplot(width = 0.1) +
  geom_smooth(method = "lm") +
  stat_cor() +
  theme_bw()
