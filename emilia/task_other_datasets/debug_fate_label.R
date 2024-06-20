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
metadat <- metadat %>% drop_na()
metadat$cell_barcode <- rownames(metadat)

umap <- seurat_object@reductions[["umap"]]@cell.embeddings
umap <- as.data.frame(umap)
umap$cell_barcode <- rownames(umap)

Idents(seurat_object) <- 'sample_type'
DimPlot(seurat_object, split.by = 'sample_name')
DimPlot(seurat_object, split.by = 'majority_fate')

# ==============================================================================
# Relabel majority lineage fate on D14
# ==============================================================================
metadat_D14 <- metadat[metadat$time_point == 14, ]
metadat_D14_fate_count <- metadat_D14 %>% 
  group_by(lineage_barcode, sample_type) %>% 
  summarise(n = n())
metadat_D14_majority_fate <- metadat_D14 %>% 
  group_by(lineage_barcode, majority_fate) %>% 
  summarise(n = n())
metadat_D14_majority_fate_w <- dcast(metadat_D14_majority_fate, lineage_barcode ~ majority_fate)

d14_lineage_size <- metadat_D14 %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_total = n())

metadat_D14_fate_count <- merge(metadat_D14_fate_count, d14_lineage_size, by = 'lineage_barcode')
metadat_D14_fate_count$fraction <- metadat_D14_fate_count$n / metadat_D14_fate_count$n_total

metadat_D14_fate_count <- subset(metadat_D14_fate_count, select = -c(n, n_total))
metadat_D14_fate_count_w <- dcast(metadat_D14_fate_count, lineage_barcode ~ sample_type)
metadat_D14_fate_count_w[is.na(metadat_D14_fate_count_w)] <- 0
metadat_D14_fate_count_w$frac_max <- pmax(metadat_D14_fate_count_w$`14_high`, metadat_D14_fate_count_w$`14_med`, metadat_D14_fate_count_w$`14_low`)

for( row in 1:nrow(metadat_D14_fate_count_w)) {
  frac_high <- metadat_D14_fate_count_w[row, '14_high']
  frac_med <- metadat_D14_fate_count_w[row, '14_med']
  frac_low <- metadat_D14_fate_count_w[row, '14_low']
  frac_max <- metadat_D14_fate_count_w[row, 'frac_max']
  
  isHighMax <- FALSE
  isMedMax <- FALSE
  isLowMax <- FALSE
  maxFate <- ''
  if (frac_high == frac_max) {
    isHighMax <- TRUE
    maxFate <- paste0(maxFate, '_high')
  }
  if(frac_med == frac_max) {
    isMedMax <- TRUE
    maxFate <- paste0(maxFate, '_med')
  }
  if(frac_low == frac_max) {
    isLowMax <- TRUE
    maxFate <- paste0(maxFate, '_low')
  }
  
  numMajorityFate <- sum(c(isHighMax, isMedMax, isLowMax))
  metadat_D14_fate_count_w[row, 'numMajorityFate'] <- numMajorityFate
  metadat_D14_fate_count_w[row, 'MajorityFate.RE'] <- maxFate
}

metadat_D14_fate_count_w <- merge(metadat_D14_fate_count_w, metadat_D14[, c('majority_fate', 'lineage_barcode')], by = 'lineage_barcode') %>% distinct()
metadat_D14_fate_count_w.single <- metadat_D14_fate_count_w[metadat_D14_fate_count_w$numMajorityFate == 1, ]

metadat_D14 <- merge(metadat_D14, metadat_D14_fate_count_w.single[, c('majority_fate', 'lineage_barcode', 'MajorityFate.RE')], by = 'lineage_barcode')
ggplot(metadat_D14, aes(x = umap_1, y = umap_2, color = sample_type)) +
  geom_point() +
  scale_color_manual(values = c('#C80036', '#3572EF', '#EF9C66')) +
  geom_density_2d(color = 'black') +
  facet_wrap(. ~ factor(MajorityFate, levels = c('_high', '_med', '_low'))) +
  theme_bw()
ggplot(metadat_D14, aes(x = umap_1, y = umap_2, color = sample_type)) +
  geom_point() +
  scale_color_manual(values = c('#C80036', '#3572EF', '#EF9C66')) +
  geom_density_2d(color = 'black') +
  facet_wrap(. ~ majority_fate) +
  theme_bw()



metadat_D14_fate_count <- merge(metadat_D14_fate_count, metadat_D14_fate_count_w[, c('lineage_barcode', 'MajorityFate.RE')], by = 'lineage_barcode') %>% distinct()
metadat_D14_fate_count <- merge(metadat_D14_fate_count, d14_lineage_size, by = 'lineage_barcode')
metadat_D14_fate_count <- metadat_D14_fate_count[metadat_D14_fate_count$MajorityFate %in% c('_high', '_med', '_low'), ]
ggplot(metadat_D14_fate_count, aes(x = MajorityFate.RE, y = n_total)) +
  geom_violin() + 
  geom_jitter(width = 0.1) +
  coord_cartesian(ylim=c(0, 20)) +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  theme_bw()

large_lin <- d14_lineage_size[d14_lineage_size$n_total > 10, ]$lineage_barcode
metadat_D14_fate_count.large <- metadat_D14_fate_count[metadat_D14_fate_count$lineage_barcode %in% large_lin, ]
ggplot(metadat_D14_fate_count.large, aes(y = lineage_barcode, x = fraction, fill = sample_type)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c('#C80036', '#3572EF', '#EF9C66')) +
  facet_wrap(. ~ MajorityFate, ncol = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ==============================================================================
# Visualize original
# ==============================================================================
metadat_D14 <- merge(metadat_D14, d14_lineage_size, by = 'lineage_barcode')
metadat_D14 <- merge(metadat_D14, umap, by = 'cell_barcode')
metadat_D14$n_total_group <- ifelse(metadat_D14$n_total < 6, '0-5', 'N/A')
metadat_D14$n_total_group <- ifelse(metadat_D14$n_total > 5 & metadat_D14$n_total < 11, '5-10', metadat_D14$n_total_group)
metadat_D14$n_total_group <- ifelse(metadat_D14$n_total > 10 & metadat_D14$n_total < 21, '10-20', metadat_D14$n_total_group)
metadat_D14$n_total_group <- ifelse(metadat_D14$n_total > 20, '>20', metadat_D14$n_total_group)


ggplot(metadat_D14, aes(x = umap_1, y = umap_2, color = n_total_group)) +
  geom_point() +
  facet_wrap(. ~ n_total_group) +
  geom_density_2d(color = 'black') +
  scale_color_manual(values = c('#0570b0', '#bdc9e1','#2b8cbe', '#74a9cf'))

metadat_D14.small <- metadat_D14[, c('lineage_barcode', 'time_point', 'n_total', 'sample_type', 'majority_fate')]

# ==============================================================================
# Visualize original
# ==============================================================================
metadat_D14_cycling <- metadat[metadat$time_point == 14, ]
metadat_D14_cycling_size <- metadat_D14_cycling %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_cycling = n())


metadat_D14.test <- merge(metadat_D14, metadat_D14_cycling_size, by = 'lineage_barcode')
metadat_D14.test$n_total_group <- ifelse(metadat_D14.test$n_cycling < 6, '0-5', 'N/A')
metadat_D14.test$n_total_group <- ifelse(metadat_D14.test$n_cycling > 5 & metadat_D14.test$n_cycling < 11, '5-10', metadat_D14.test$n_total_group)
metadat_D14.test$n_total_group <- ifelse(metadat_D14.test$n_cycling > 10 & metadat_D14.test$n_cycling < 21, '10-20', metadat_D14.test$n_total_group)
metadat_D14.test$n_total_group <- ifelse(metadat_D14.test$n_cycling > 20, '>20', metadat_D14.test$n_total_group)

ggplot(metadat_D14.test, aes(x = umap_1, y = umap_2, color = n_total_group)) +
  geom_point() +
  facet_grid(. ~ n_total_group) +
  geom_density_2d(color = 'black') +
  scale_color_manual(values = c('#0570b0', '#bdc9e1','#2b8cbe', '#74a9cf')) +
  theme_bw()


