rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_postprocessed.RData")
all_data_original <- all_data
load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

# do the chromvar analysis 
tab_mat <- table(all_data_original$assigned_lineage, all_data_original$dataset)

#################

# let's define winners and losers
# score each lineage by its mean day10 score
lineage_score <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data_original)[which(all_data_original$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  # quantile(all_data_original$imputed_count[cell_names], probs = 0.75, na.rm = T)
  length(which(all_data_original$imputed_count[cell_names] >= 0))
})
names(lineage_score) <- rownames(tab_mat)
quantile(lineage_score, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score > 0))

lineage_score2 <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data_original)[which(all_data_original$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  mean(all_data_original$imputed_count[cell_names], na.rm = T)
})
names(lineage_score2) <- rownames(tab_mat)
quantile(lineage_score2, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score2 > 0))

cell_winning_idx <- intersect(which(all_data_original$assigned_lineage %in% names(lineage_score2)[intersect(
  which(lineage_score > 0),
  which(lineage_score2 >= -1)
)]),
which(all_data_original$dataset == "day0"))
cell_winning_names <- colnames(all_data_original)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data_original$assigned_lineage %in% names(lineage_score2)[which(lineage_score2 <= -1)]),
                             which(all_data_original$dataset == "day0"))
cell_losing_names <- colnames(all_data_original)[cell_losing_idx]
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

names(data.use) <- sapply(1:length(data.use), function(i){
  paste0(names(data.use)[i], ", ", rownames(de_res)[i], 
         "\n-Log10Pval=", round(-log10(de_res[i,"p_val"]), 2),
         ", value=", round(de_res[i,"avg_diff"], 2))
})
data.use <- data.use[1:28]

plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
plot1 <- plot1 + ggplot2::theme_bw()
width <- 15; height <- 10

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_chromVar_DABTRAM_day0-win-vs-lose_day10-lineage-imputation_stepdown.png"),
                plot1, device = "png", width = width, height = height, units = "in")

#######################

diff_vec <- de_res$avg_diff
names(diff_vec) <- name_vec
pvalue_vec <- de_res$p_val
names(pvalue_vec) <- name_vec
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")

logpval_vec <- -log10(de_res$p_val)
names(logpval_vec) <- name_vec

labeling_vec <- rep("1", length(diff_vec))
names(labeling_vec) <- names(diff_vec)
labeling_vec[names(logpval_vec)[which(logpval_vec >= 2)]] <- "3"
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
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day0 Motif, based on Day10 lineage growth potential\n(Winner-Loser)")) +
  ggplot2::xlab("Mean difference in chromVar score: (Winner-Loser)") + ggplot2::ylab("DE p-value (-Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_chromVar_DABTRAM_day0-win-vs-lose_day10-lineage-imputation_stepdown_DE_volcano.png"),
                p1, device = "png", width = 10, height = 10, units = "in")
