rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")
all_data_dabtram <- all_data
load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

# do the chromvar analysis 
tab_mat <- table(all_data_dabtram$assigned_lineage, all_data_dabtram$dataset)

##################################

# let's define winners and losers
# score each lineage by its mean day10 score
lineage_score <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(all_data_dabtram)[which(all_data_dabtram$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  mean(all_data_dabtram$pseudotime[cell_names], na.rm = T)
})
names(lineage_score) <- rownames(tab_mat)
quantile(lineage_score, na.rm = T)

cell_winning_idx <- intersect(which(all_data_dabtram$assigned_lineage %in% names(lineage_score)[which(lineage_score > 0.5)]),
                              which(all_data_dabtram$dataset == "day0"))
cell_winning_names <- colnames(all_data_dabtram)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data_dabtram$assigned_lineage %in% names(lineage_score)[which(lineage_score < 0.3)]),
                             which(all_data_dabtram$dataset == "day0"))
cell_losing_names <- colnames(all_data_dabtram)[cell_losing_idx]
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

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_chromVar_DABTRAM_day0-win-vs-lose.png"),
                plot1, device = "png", width = width, height = height, units = "in")

##################################

cell_winning_idx <- intersect(which(all_data_dabtram$assigned_lineage %in% names(lineage_score)[which(lineage_score > 0.5)]),
                              which(all_data_dabtram$dataset == "day0"))
cell_winning_names <- colnames(all_data_dabtram)[cell_winning_idx]

dead_lineages <- rownames(tab_mat)[intersect(
  which(tab_mat[,"day10_DABTRAM"] == 0),
  which(tab_mat[,"week5_DABTRAM"] == 0)
)]
cell_losing_idx <- intersect(which(all_data_dabtram$assigned_lineage %in% dead_lineages),
                              which(all_data_dabtram$dataset == "day0"))
cell_losing_names <- colnames(all_data_dabtram)[cell_losing_idx]

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
         "\n-Log10pval=", round(-log10(de_res[i,"p_val"]), 2),
         ", value=", round(de_res[i,"avg_diff"], 2))
})
data.use <- data.use[1:28]

plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
plot1 <- plot1 + ggplot2::theme_bw()
width <- 15; height <- 10

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_chromVar_DABTRAM_day0-win-vs-death.png"),
                plot1, device = "png", width = width, height = height, units = "in")


