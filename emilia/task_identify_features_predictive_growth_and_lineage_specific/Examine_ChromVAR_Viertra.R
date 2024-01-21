library(Seurat)
library(ggplot2)
library(ggpubr)

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'

# ==============================================================================
# Read data
# ==============================================================================

sc <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day0_processed.rds")
metadat <- sc@meta.data
metadat$cell_barcode <- rownames(metadat)
dev <- readRDS(paste0(in_dir, "mat.motifs_day0_Viertra_annotations.rds"))
dev_z <- as.data.frame(dev@assays@data@listData[["z"]])

"AC0201|TEAD|TEA", "AC0241|NFE2L/JUN|bZIP"

motifs <- rownames(dev_z)[grep("FOS|JUN|TEAD", rownames(dev_z))]

# ==============================================================================
# Extract data
# ==============================================================================
dev_z <- as.data.frame(t(dev_z))
dev_z <- dev_z[, motifs]
dev_z$cell_barcode <- rownames(dev_z)

dev_z <- merge(dev_z, metadat[, c('cell_barcode', 'nCount_ATAC')], by = 'cell_barcode')

# ==============================================================================
# Plotting
# ==============================================================================
ggplot(dev_z, aes(x = nCount_ATAC, y = `AC0508|FOSL|bZIP`)) +
  geom_point() +
  stat_cor(method="pearson")


