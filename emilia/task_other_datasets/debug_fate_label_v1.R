rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(reshape2)
library(Seurat)
# library(multiomeFate)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/output/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
seurat_object <- readRDS(paste0(in_dir, 'pc9_time_course.rds'))
Idents(seurat_object) <- 'majority_fate'
DimPlot(seurat_object)
# ggsave('~/Downloads/watermelon_pc9_UMAP_majority_fate.png', dpi = 300, width = 4, height = 3)

metadat <- seurat_object@meta.data
metadat <- metadat %>% drop_na()
metadat$cell_barcode <- rownames(metadat)

supp.data <- read.csv('~/Downloads/watermelon_fig2c.csv')
colnames(supp.data) <- c('lineage_barcode', '0', '3', '7', '14','majority_fate')
supp.data.m <- melt(supp.data[, c('lineage_barcode', '0', '3', '7', '14')])

lineage_size <- metadat %>% 
  group_by(lineage_barcode, time_point) %>% 
  summarise(num_cells = n())

comp <- merge(lineage_size, supp.data.m, by.x = c('lineage_barcode', 'time_point'), by.y = c('lineage_barcode', 'variable'))

ggplot(comp, aes(x = num_cells, y = value)) +
  geom_point() +
  xlab('Lineage size (we calculated)') +
  ylab('Lineage size from source data 2c') +
  facet_wrap(. ~ time_point, scale = 'free')

# ==============================================================================
# Compare majority lineage fate information
# ==============================================================================

majority.fate.meta <- metadat[, c('lineage_barcode', 'majority_fate')]
rownames(majority.fate.meta) <- seq(1: nrow(majority.fate.meta))
majority.fate.meta <- majority.fate.meta %>% distinct()
colnames(majority.fate.meta) <- c('lineage_barcode', 'majority_fate_meta_data')

majority.fate.supp <- supp.data[, c('lineage_barcode', 'majority_fate')]
rownames(majority.fate.supp) <- seq(1: nrow(majority.fate.supp))
majority.fate.supp <- majority.fate.supp %>% distinct()
colnames(majority.fate.supp) <- c('lineage_barcode', 'majority_fate_supp_data')

comp_df <- merge(majority.fate.meta, majority.fate.supp, by = 'lineage_barcode')
comp_df$consistent <- ifelse(comp_df$majority_fate_meta_data == comp_df$majority_fate_supp_data, 'Yes', 'No')


# ==============================================================================
# Compare lineage size
# ==============================================================================

ggplot(supp.data, aes(x = majority_fate, y = `X14`)) +
  geom_violin(scale = 'width') +
  geom_jitter(width = 0.2, alpha = 0.2) +
  # coord_cartesian(ylim=c(0, 10))+
  stat_summary(fun=median, geom="point", size=2, color="red") +
  theme_bw()

supp.data$X14_log <- log10(supp.data$X14)
ggplot(supp.data, aes(x = majority_fate, y = `X14_log`)) +
  geom_violin(scale = 'width') +
  geom_jitter(width = 0.2, alpha = 0.2) +
  # coord_cartesian(ylim=c(0, 10))+
  stat_summary(fun=median, geom="point", size=2, color="red") +
  theme_bw()


supp.data2e <- read.csv('~/Downloads/watermelon_fig2e.csv')
ggplot(supp.data2e) +
  geom_jitter(aes(x = clone_size_by_rep_rep1, y = clone_size_by_rep_rep2))
supp.data2e$consistent_fate <- ifelse(supp.data2e$majority_fate_rep_rep1 == supp.data2e$majority_fate_rep_rep2, 'Yes', 'No')
ggplot(supp.data2e, aes(x = clone_size_by_rep_rep1, y = clone_size_by_rep_rep2)) +
  geom_jitter(aes(color = majority_fate_rep_rep1)) +
  ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  facet_wrap(. ~ consistent_fate)

supp.data2e.consistent <- supp.data2e[supp.data2e$consistent_fate == 'Yes', ]

supp.data2e$clone_size_by_rep_rep1_log <- log10(supp.data2e$clone_size_by_rep_rep1)
supp.data2e$clone_size_by_rep_rep2_log <- log10(supp.data2e$clone_size_by_rep_rep2)

ggplot(supp.data2e, aes(x = majority_fate_rep_rep1, y = clone_size_by_rep_rep1)) +
  geom_violin() +
  geom_jitter() +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  theme_bw()

ggplot(supp.data2e.consistent, aes(x = majority_fate_rep_rep2, y = clone_size_by_rep_rep2_log)) +
  geom_violin() +
  geom_jitter() +
  stat_summary(fun=median, geom="point", size=2, color="red")

ggplot(supp.data2e.consistent, aes(x = majority_fate_rep_rep1, y = clone_size_by_rep_rep1_log)) +
  geom_violin() +
  geom_jitter() +
  stat_summary(fun=median, geom="point", size=2, color="red") 


supp.data2c.e <- merge(supp.data, supp.data2e, by = 'lineage_barcode')
# ==============================================================================
# Plot UMAP
# ==============================================================================

umap <- seurat_object@reductions[["umap"]]@cell.embeddings
umap <- as.data.frame(umap)
umap$cell_barcode <- rownames(umap)

rep2 <- metadat[metadat$sample_name %in% c("0", "3_rep1", "3_rep2", "7_rep1", "7_rep2")]


