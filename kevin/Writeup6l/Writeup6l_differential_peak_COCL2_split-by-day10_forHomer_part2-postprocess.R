rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(GGally)
library(ggplot2)

load("../../../../out/kevin/Writeup6l/Writeup6l_differential_peak_COCL2_split-by-day10_part2.RData")

# first, positives:
padjust_vec <- pos_motif_de$p.adjust
log10p_vec <- -log10(pos_motif_de$pvalue)
fold_vec <- pos_motif_de$fold.enrichment

motif_idx <- sort(unique(c(grep("JUN", pos_motif_de$motif.name),
                           grep("FOS", pos_motif_de$motif.name),
                           grep("TEAD", pos_motif_de$motif.name))))

labeling_vec <- rep(0, length(log10p_vec))
labeling_vec[which(padjust_vec <= 0.05)] <- 2
labeling_vec[intersect(which(padjust_vec <= 0.05), order(log10p_vec, decreasing = T)[1:20])] <- 3
labeling_vec[motif_idx] <- 1

df <- data.frame(fold = fold_vec,
                 log10pval = log10p_vec,
                 name = pos_motif_de$motif.name,
                 labeling = labeling_vec)
df[,"labeling"] <- as.factor(df[,"labeling"])
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == 0),  which(df[,"labeling"] == 2), which(df[,"labeling"] == 1), which(df[,"labeling"] == 3)),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = fold, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "dodgerblue3",  "red", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("1","3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::geom_hline(yintercept = min(log10p_vec[padjust_vec <= 0.05]), linetype = "dashed", color = "red", size = 1.5) 
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2 (Day0 winners vs Day0 losers, split by Day10):\nPositive peaks")) 
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6l/Writeup6l_differential_peak_COCL2_split-by-day10_positive-peaks.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

#########################

padjust_vec <- neg_motif_de$p.adjust
log10p_vec <- -log10(neg_motif_de$pvalue)
fold_vec <- neg_motif_de$fold.enrichment

motif_idx <- sort(unique(c(grep("JUN", neg_motif_de$motif.name),
                           grep("FOS", neg_motif_de$motif.name),
                           grep("TEAD", neg_motif_de$motif.name))))

labeling_vec <- rep(0, length(log10p_vec))
labeling_vec[which(padjust_vec <= 0.05)] <- 2
labeling_vec[intersect(which(padjust_vec <= 0.05), order(log10p_vec, decreasing = T)[1:20])] <- 3
labeling_vec[motif_idx] <- 1

df <- data.frame(fold = fold_vec,
                 log10pval = log10p_vec,
                 name = neg_motif_de$motif.name,
                 labeling = labeling_vec)
df[,"labeling"] <- as.factor(df[,"labeling"])
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == 0),  which(df[,"labeling"] == 2), which(df[,"labeling"] == 1), which(df[,"labeling"] == 3)),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = fold, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "dodgerblue3",  "red", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("1","3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::geom_hline(yintercept = min(log10p_vec[padjust_vec <= 0.05]), linetype = "dashed", color = "red", size = 1.5) 
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2 (Day0 winners vs Day0 losers, split by Day10):\nNegative peaks")) 
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6l/Writeup6l_differential_peak_COCL2_split-by-day10_negative-peaks.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

