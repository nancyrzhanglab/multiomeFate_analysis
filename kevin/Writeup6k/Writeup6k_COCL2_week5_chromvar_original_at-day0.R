rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

# do the chromvar analysis 
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)


# let's define winners and losers
# score each lineage by its mean day10 score
lineage_score <- tab_mat[,"week5_COCL2"]
names(lineage_score) <- rownames(tab_mat)
quantile(lineage_score, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score > 0))

cell_winning_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score >= 10)]),
                              which(all_data$dataset == "day0"))
cell_winning_names <- colnames(all_data)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score == 0)]),
                             which(all_data$dataset == "day0"))
cell_losing_names <- colnames(all_data)[cell_losing_idx]
length(cell_winning_names); length(cell_losing_names)
keep_vec <- rep(NA, ncol(all_data))
names(keep_vec) <- colnames(all_data)
keep_vec[cell_winning_names] <- "winning_day0"
keep_vec[cell_losing_names] <- "losing_day0"

all_data$ident <- keep_vec
Seurat::Idents(all_data) <- "ident"


Seurat::DefaultAssay(all_data) <- "chromvar"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = "winning_day0",
  ident.2 = "losing_day0",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)

motifs <- rownames(de_res)
data.use <- Signac::GetMotifData(object = all_data,
                                 assay = "ATAC",
                                 slot = "pwm")
data.use <- data.use[motifs]
names(data.use) <- Signac::GetMotifData(
  object = all_data,
  assay = "ATAC",
  slot = "motif.names"
)[motifs]
name_vec <- names(data.use) 
tmp <- de_res; rownames(tmp) <- names(data.use); tmp[1:50,]


diff_vec <- de_res$avg_diff
names(diff_vec) <- name_vec
pvalue_vec <- de_res$p_val
names(pvalue_vec) <- name_vec
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")
logpval_vec <- -log10(de_res$p_val)
names(logpval_vec) <- name_vec

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 1.5)]] <- "3"
labeling_vec[grep("JUN", name_vec)] <- "2"
labeling_vec[grep("FOS", name_vec)] <- "2"
labeling_vec[grep("NFE2", name_vec)] <- "2"
table(labeling_vec)

df <- data.frame(difference = diff_vec,
                 log10pval = logpval_vec,
                 name = name_vec,
                 labeling = labeling_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "1"), which(df[,"labeling"] %in% c("2", "3"))),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = difference, 
                                       y = log10pval))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red", "blue"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling %in% c("2", "3")),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
if(any(pvalue_vec <= 0.05)) {
  p1 <- p1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                                 color = "red", linewidth=2)
}
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 Motif, based on Week5 lineage size\n(Winner-Loser)")) +
  ggplot2::xlab("Mean difference in chromVar score: (Winner-Loser)") + ggplot2::ylab("DE p-value (-Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_chromVar_COCL2_day0-win-vs-lose_week5_original_DE_volcano.png"),
                p1, device = "png", width = 10, height = 10, units = "in")
