library(Seurat)
library(Signac)
library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/data/subsampled/'
day10 <- readRDS(paste0(in_dir, 'day10_DABTRAM_processed_sm.rds'))
DefaultAssay(day10) <- 'ATAC'

tp_early <- 'day10'
treatment <- "DABTRAM"
kevin_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup6r/'
load(paste0(kevin_dir, "Writeup6r_", treatment, "_", tp_early, "_lineage-imputation_postprocess.RData"))
dabtram_d10_imputed <- cell_imputed_score[!is.na(cell_imputed_score)]


# ==============================================================================
# Wrangle
# ==============================================================================
day10@assays[["ATAC"]]@fragments[[4]]@path  <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_DABTRAM/outs/atac_fragments.tsv.gz'
day10@assays[["ATAC"]]@fragments[[1]] <- NULL
day10@assays[["ATAC"]]@fragments[[1]] <- NULL
day10@assays[["ATAC"]]@fragments[[1]] <- NULL
day10@assays[["ATAC"]]@fragments[[2]] <- NULL
day10@assays[["ATAC"]]@fragments[[2]] <- NULL
day10@assays[["ATAC"]]@fragments[[2]] <- NULL

# ==============================================================================
# Label winner vs loser
# ==============================================================================
metadat <- day10@meta.data
metadat$cell_barcode <- rownames(metadat)

dabtram_d10_imputed_df <- data.frame(matrix(nrow=length(dabtram_d10_imputed), ncol = 0))
dabtram_d10_imputed_df$cell_barcode <- names(dabtram_d10_imputed)
dabtram_d10_imputed_df$cell_growth_potential <- dabtram_d10_imputed

metadat <- merge(metadat, dabtram_d10_imputed_df, by = 'cell_barcode')
rownames(metadat) <- metadat$cell_barcode

metadat$isWinner <- ifelse(metadat$cell_growth_potential > 0, 'Winner', 'Loser')

day10 <- AddMetaData(day10, metadat)

# ==============================================================================
# Plot coverage
# ==============================================================================
Idents(day10) <- "isWinner"
cov_plot <- CoveragePlot(
  object = day10,
  region = "EREG",
  idents = c("Winner", "Loser"),
  annotation = TRUE,
  peaks = TRUE
)

ggsave('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/cov_plot_EREG.png', width = 7, height = 4)

cov_plot <- CoveragePlot(
  object = day10,
  region = "chr4-74354000-74820000",
  idents = c("Winner", "Loser"),
  annotation = TRUE,
  peaks = TRUE
)
cov_plot

ggsave('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/cov_plot_BTC.png', width = 7, height = 4)

cov_plot <- CoveragePlot(
  object = day10,
  region = "EGF",
  idents = c("Winner", "Loser"),
  annotation = TRUE,
  peaks = TRUE
)
cov_plot

ggsave('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/cov_plot_EGF.png', width = 7, height = 4)

