rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)
source("../Writeup5a/color_palette.R")

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_chromatinAct_fasttopics_CIS.RData")
ft_CIS_cact <- topic_res
load("../../../../out/kevin/Writeup6b/Writeup6b_chromatinAct_fasttopics_COCL2.RData")
ft_COCL2_cact <- topic_res
load("../../../../out/kevin/Writeup6b/Writeup6b_chromatinAct_fasttopics_DABTRAM.RData")
ft_DABTRAM_cact <- topic_res

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

##################

ft_CIS_cact_load <- ft_CIS_cact$F
ft_COCl2_cact_load <- ft_COCL2_cact$F
ft_DT_cact_load <- ft_DABTRAM_cact$F

ft_CIS_rna_load <- all_data[["fasttopic_CIS"]]@feature.loadings
ft_COCl2_rna_load <- all_data[["fasttopic_COCL2"]]@feature.loadings
ft_DT_rna_load <- all_data[["fasttopic_DABTRAM"]]@feature.loadings

tmp <- table(c(rownames(ft_CIS_cact_load), rownames(ft_COCl2_cact_load), rownames(ft_DT_cact_load),
               rownames(ft_CIS_rna_load), rownames(ft_COCl2_rna_load), rownames(ft_DT_rna_load)))
gene_vec <- sort(names(tmp)[tmp == 6])

cor_CIS_CIS <- cor(sqrt(ft_CIS_rna_load[gene_vec,]), sqrt(ft_CIS_cact_load[gene_vec,]))
cor_CIS_COCL2 <- cor(sqrt(ft_CIS_rna_load[gene_vec,]), sqrt(ft_COCl2_cact_load[gene_vec,]))
cor_CIS_DT <- cor(sqrt(ft_CIS_rna_load[gene_vec,]), sqrt(ft_DT_cact_load[gene_vec,]))

cor_COCL2_CIS <- cor(sqrt(ft_COCl2_rna_load[gene_vec,]), sqrt(ft_CIS_cact_load[gene_vec,]))
cor_COCL2_COCL2 <- cor(sqrt(ft_COCl2_rna_load[gene_vec,]), sqrt(ft_COCl2_cact_load[gene_vec,]))
cor_COCL2_DT <- cor(sqrt(ft_COCl2_rna_load[gene_vec,]), sqrt(ft_DT_cact_load[gene_vec,]))

cor_DT_CIS <- cor(sqrt(ft_DT_rna_load[gene_vec,]), sqrt(ft_CIS_cact_load[gene_vec,]))
cor_DT_COCL2 <- cor(sqrt(ft_DT_rna_load[gene_vec,]), sqrt(ft_COCl2_cact_load[gene_vec,]))
cor_DT_DT <- cor(sqrt(ft_DT_rna_load[gene_vec,]), sqrt(ft_DT_cact_load[gene_vec,]))

png(paste0("../../../../out/figures/Writeup6b/Writeup6b_RNA_chromatinAct_ftTopic_correlation.png"),
    width = 3000, height = 3000, units = "px", res = 300)
par(mfrow=c(3,3))
par(mar=c(5,5,1,1))
fields::image.plot(cor_CIS_CIS, xlab="CIS (RNA)", ylab="CIS (Chr)", cex.lab=2)
fields::image.plot(cor_CIS_COCL2, xlab="CIS (RNA)", ylab="COCL2 (Chr)", cex.lab=2)
fields::image.plot(cor_CIS_DT, xlab="CIS (RNA)", ylab="DABTRAM (Chr)", cex.lab=2)

fields::image.plot(cor_COCL2_CIS, xlab="COCL2 (RNA)", ylab="CIS (Chr)", cex.lab=2)
fields::image.plot(cor_COCL2_COCL2, xlab="COCL2 (RNA)", ylab="COCL2 (Chr)", cex.lab=2)
fields::image.plot(cor_COCL2_DT, xlab="COCL2 (RNA)", ylab="DABTRAM (Chr)", cex.lab=2)

fields::image.plot(cor_DT_CIS, xlab="DABTRAM (RNA)", ylab="CIS (Chr)", cex.lab=2)
fields::image.plot(cor_DT_COCL2, xlab="DABTRAM (RNA)", ylab="COCL2 (Chr)", cex.lab=2)
fields::image.plot(cor_DT_DT, xlab="DABTRAM (RNA)", ylab="DABTRAM (Chr)", cex.lab=2)
graphics.off()

#####################

rna_idx <- 10; cact_idx <- 18
# rna_idx <- 26; cact_idx <- 7

png(paste0("../../../../out/figures/Writeup6b/Writeup6b_RNA", rna_idx, "_chromatinAct", cact_idx, "_ftTopic_DABTRAM.png"),
    width = 2000, height = 2000, units = "px", res = 300)
plot(x = sqrt(ft_DT_rna_load[gene_vec,rna_idx]),
     y = sqrt(ft_DT_cact_load[gene_vec,cact_idx]),
     main = paste0("DABTRAM: Cor=", round(stats::cor(sqrt(ft_DT_rna_load[gene_vec,10]), sqrt(ft_DT_cact_load[gene_vec,18])),2)),
     xlab = paste0("FT (RNA, ", rna_idx, ", Genes, Sqrt)"),
     ylab = paste0("FT (Chromatin act, ", cact_idx, ", Genes, Sqrt)"),
     pch = 16)
graphics.off()

cell_names1 <- colnames(all_data)[which(all_data$dataset == "day0")]
cell_names2 <- colnames(all_data)[which(all_data$dataset == "day10_DABTRAM")]
cell_names3 <- colnames(all_data)[which(all_data$dataset == "week5_DABTRAM")]
cell_name_vec <- c(cell_names1, cell_names2, cell_names3)
dataset_vec <- c(rep("day0", length(cell_names1)),
                 rep("day10_DABTRAM", length(cell_names2)),
                 rep("week5_DABTRAM", length(cell_names3)))
names(dataset_vec) <- cell_name_vec

rna_vec <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[cell_name_vec,rna_idx]
chr_vec <- ft_DABTRAM_cact$L[cell_name_vec,cact_idx]

set.seed(10)
idx <- sample(1:length(rna_vec))

png(paste0("../../../../out/figures/Writeup6b/Writeup6b_RNA", rna_idx, "_chromatinAct", cact_idx, "_ftTopic_DABTRAM_cell.png"),
    width = 2000, height = 2000, units = "px", res = 300)
plot(x = log10(rna_vec)[idx],
     y = log10(chr_vec)[idx],
     main = paste0("DABTRAM: Cor=", round(stats::cor(log10(rna_vec), log10(chr_vec)),2)),
     xlab = paste0("FT (RNA, ", rna_idx, ", Cell, Log10)"),
     ylab = paste0("FT (Chromatin act, ", cact_idx, ", Cell, Log10)"),
     col = col_palette[dataset_vec][idx],
     pch = 16)
graphics.off()

#####################

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
# lineage_name <- rownames(tab_mat)[which.max(tab_mat[,"week5_DABTRAM"])]
lineage_name <- "Lin112532"
cell_names <- colnames(all_data)[all_data$assigned_lineage == lineage_name]

png(paste0("../../../../out/figures/Writeup6b/Writeup6b_RNA", rna_idx, "_chromatinAct", cact_idx, "_ftTopic_DABTRAM_cell_", lineage_name, ".png"),
    width = 2000, height = 2000, units = "px", res = 300)
plot(x = log10(rev(rna_vec[cell_names])),
     y = log10(rev(chr_vec[cell_names])),
     main = paste0("DABTRAM: Cor=", round(stats::cor(log10(rna_vec), log10(chr_vec)),2)),
     xlab = paste0("FT (RNA, ", rna_idx, ", Cell, Log10)"),
     ylab = paste0("FT (Chromatin act, ", cact_idx, ", Cell, Log10)"),
     col = col_palette[rev(dataset_vec[cell_names])],
     pch = 16)
graphics.off()


