rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

keep_idx <- rep(NA, ncol(all_data))
keep_idx[intersect(
  intersect(
    which(all_data$assigned_posterior >= 0.5),
    which(!is.na(all_data$assigned_lineage))
  ),
  which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))
)] <- TRUE
all_data$keep <- keep_idx
all_data <- subset(all_data, keep == TRUE)

# do the chromvar analysis 
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = T),]

#################

# let's define winners and losers
# score each lineage by its mean day10 score
lineage_score <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data)[intersect(
    which(all_data$assigned_lineage == lineage_name),
    which(all_data$dataset == "day10_COCL2")
  )]
  length(cell_names)
})
names(lineage_score) <- rownames(tab_mat)
quantile(lineage_score, na.rm = T, probs = seq(0,1,length.out=11))
length(which(lineage_score > 0))

cell_winning_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score >= 10)]),
                              which(all_data$dataset == "day0"))
cell_winning_names <- colnames(all_data)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score <= 0)]),
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
tmp <- de_res; rownames(tmp) <- names(data.use); tmp[1:50,]

names(data.use) <- sapply(1:length(data.use), function(i){
  paste0(names(data.use)[i], ", ", rownames(de_res)[i], 
         "\n-Log10AdjPval=", round(-log10(de_res[i,"p_val_adj"]), 2),
         ", value=", round(de_res[i,"avg_diff"], 2))
})
data.use <- data.use[1:28]

plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
plot1 <- plot1 + ggplot2::theme_bw()
width <- 15; height <- 10

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_chromVar_COCL2_day0-win-vs-lose_day10_original.png"),
                plot1, device = "png", width = width, height = height, units = "in")
