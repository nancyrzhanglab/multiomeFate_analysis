rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(reshape2)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/output/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
metadat <- seurat_object@meta.data
metadat$cell_barcode <- rownames(metadat)

umap <- as.data.frame(seurat_object@reductions[["umap"]]@cell.embeddings)
umap$cell_barcode <- rownames(umap)

metadat <- merge(metadat, umap, by = 'cell_barcode')

# ==============================================================================
# Lineage size
# ==============================================================================

total_cell_count <- metadat %>% 
  group_by(time_point) %>% 
  summarise(total_cell = n())

lineage_count <- metadat %>% 
  group_by(time_point, lineage_barcode) %>% 
  summarise(n_cell = n())
lineage_count <- lineage_count %>% drop_na()

lineage_count <- merge(lineage_count, total_cell_count, by = 'time_point')
lineage_count$freq <- lineage_count$n_cell / lineage_count$total_cell

lineage_count_use <- lineage_count[, c('lineage_barcode', 'time_point', 'freq')]
lineage_count_w <- dcast(lineage_count_use, lineage_barcode ~ time_point)

lineage_count_w[is.na(lineage_count_w)] <- 0
rownames(lineage_count_w) <- lineage_count_w$lineage_barcode
lineage_count_w <- subset(lineage_count_w, select = -c(lineage_barcode))

lineage_count_w$`0` <- log10(lineage_count_w$`0` + 10**-5)
lineage_count_w$`3` <- log10(lineage_count_w$`3`+ 10**-5)
lineage_count_w$`7` <- log10(lineage_count_w$`7`+ 10**-5)
lineage_count_w$`14` <- log10(lineage_count_w$`14`+ 10**-5)
ggpairs(lineage_count_w)
