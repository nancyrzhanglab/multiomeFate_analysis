rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
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
# Basic info
# ==============================================================================
metadat$lineage_barcode_present <- ifelse(is.na(metadat$lineage_barcode), 'No', 'Yes')
total_cells <- metadat %>% 
  group_by(time_point) %>% 
  summarise(n = n())

lineage_size <- metadat %>% 
  group_by(lineage_barcode, time_point) %>% 
  summarise(num_cells = n())
lineage_size_d14 <- lineage_size[lineage_size$time_point == 14, ]

lineage_size <- merge(lineage_size, total_cells, by = c('time_point'))

metadat <- merge(metadat, lineage_size, by = c('time_point', 'lineage_barcode'))
metadat <- metadat %>% drop_na()
metadat$log10_num_cells <- log10(metadat$num_cells + 1)
metadat$num_cells_cap <- ifelse(metadat$num_cells > 20, 20, metadat$num_cells)
metadat$num_cells_cap <- ifelse(metadat$num_cells_cap > 10 & metadat$num_cells < 20, 10, metadat$num_cells_cap)
metadat$num_cells_cap <- ifelse(metadat$num_cells_cap < 10, 5, metadat$num_cells_cap)
metadat$num_cells_cap <- as.factor(metadat$num_cells_cap )
ggplot(metadat, aes(x = umap_1, y = umap_2, after_stat(level))) +
  # geom_point(aes(x = umap_1, y = umap_2, color = num_cells_cap), size = 0.5) +
  geom_density_2d_filled() +
  # scale_color_gradient(low = 'blue', high = 'red') +
  facet_grid(. ~ num_cells_cap) +
  theme_bw()

ggplot(metadat, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = log10_num_cells), size = 0.5) +
  scale_color_gradient(low = 'blue', high = 'red') +
  facet_grid(. ~ num_cells_cap) +
  theme_bw()
  

order <- c("0", "3_rep1", "3_rep2", "7_rep1", "7_rep2", "14_rep1_low", "14_rep2_low", "14_rep1_med", "14_rep2_med", "14_rep1_high", "14_rep2_high")
ggplot(metadat, aes(x = factor(sample_name, levels = order), y = log10_num_cells)) +
  geom_violin(scale = 'width') +
  geom_jitter(alpha = 0.1, width = 0.1) +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  # facet_wrap(. ~ time_point) +
  theme_bw()

fate <- metadat[, c('lineage_barcode', 'majority_fate')] %>% distinct()
lineage_size <- merge(lineage_size, fate, by = 'lineage_barcode')
lineage_size$num_cells_log10 <- log10(lineage_size$num_cells)
ggplot(lineage_size, aes(x = majority_fate, y = num_cells_log10)) +
  geom_violin(scale = 'width') +
  geom_jitter(alpha = 0.1, width = 0.1) +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  facet_wrap(. ~ time_point) +
  theme_bw()

sample_type_size <- metadat %>% 
  group_by(lineage_barcode, sample_type) %>% 
  summarise(n_cells = n())
sample_type_size$num_cells_log10 <- log10(sample_type_size$n_cells)

order <- c('0', '3', '7', '14_low', '14_med', '14_high')
ggplot(sample_type_size, aes(x = factor(sample_type, levels = order), y = num_cells_log10)) +
  geom_violin(scale = 'width') +
  geom_jitter(alpha = 0.1, width = 0.1) +
  # ylim(0, 20) +
  stat_summary(fun=mean, geom="point", size=2, color="red") +
  theme_bw()

# ==============================================================================
# D14 status
# ==============================================================================
metadat_d14 <- metadat[metadat$time_point == 14, ]
sample_type_size <- metadat_d14 %>% 
  group_by(sample_type) %>% 
  summarise(n_cells = n())

lineage_fate_d14 <- metadat_d14 %>% 
  group_by(lineage_barcode, sample_type, majority_fate) %>% 
  summarise(n_cells = n())
lineage_size_d14 <- metadat_d14 %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_total = n())

lineage_fate_d14 <- merge(lineage_fate_d14, lineage_size_d14, by = 'lineage_barcode')
lineage_fate_d14$freq <- lineage_fate_d14$n_cells / lineage_fate_d14$n_total
lineage_fate_d14$log2_freq <- log2(lineage_fate_d14$freq)
lineage_fate_d14$log10_freq <- log10(lineage_fate_d14$freq)
lineage_fate_d14.large <- lineage_fate_d14[lineage_fate_d14$n_total > 10, ]
lineage_fate_d14.large <- lineage_fate_d14.large[order(lineage_fate_d14.large$freq),]

ggplot(lineage_fate_d14.large, aes(y = lineage_barcode, x = freq, fill = sample_type)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c('#C80036', '#3572EF', '#EF9C66')) +
  facet_wrap(. ~ majority_fate, ncol = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave('~/Downloads/lineage_bias.png', dpi = 300)

plot(ecdf(lineage_size_d14$n_total))

quantile(lineage_size_d14$n_total, probs = seq(0, 1, 0.1))
