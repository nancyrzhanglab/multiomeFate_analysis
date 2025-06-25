rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
fatebias_out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'

figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V4/Fig5/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

dataset_colors <- c(day0 = "#696969",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")
# =============================================================================
# Read data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}
all_data.day0 <- subset(all_data, dataset == 'day0')
all_data.day10_DABTRAM <- subset(all_data, dataset == 'day10_DABTRAM')

fp.d0_d10 <- as.data.frame(all_data_fatepotential[[paste0("fatepotential_", treatment, "_d0_d10")]][["cell_imputed_score"]])
fp.d10_w5 <- as.data.frame(all_data_fatepotential[[paste0("fatepotential_", treatment, "_d10_w5")]][["cell_imputed_score"]])

# lineage data
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.day0 <- metadat[metadat$dataset == 'day0',]
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM',]
metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM',]

# umap
ft.umap.DABTRAM <- all_data@reductions[[paste0("ft.DABTRAM.umap")]]@cell.embeddings
ft.umap.DABTRAM <- as.data.frame(ft.umap.DABTRAM)
ft.umap.DABTRAM$cell_id <- rownames(ft.umap.DABTRAM)

ft.umap.DABTRAM <- merge(ft.umap.DABTRAM, metadat[, c('cell_id', 'assigned_lineage', 'dataset')], by = 'cell_id')

colnames(fp.d10_w5) <- c('fp.d10_w5')
fp.d10_w5$cell_id <- rownames(fp.d10_w5)

ada_index.plot <- read.csv(paste0(output_dir, 'adaptation_index_norm_by_NoAdapt.csv'))
ada_index.plot.DABTRAM <- ada_index.plot[ada_index.plot$dataset == 'DABTRAM_d10_w5', ]

# lineage size
lineage_size.day0 <- metadat.day0 %>%
  group_by(assigned_lineage) %>%
  summarise(n_cells.d0 = n()) %>%
  ungroup()

ada_index.plot.DABTRAM <- merge(ada_index.plot.DABTRAM, lineage_size.day0, by.x = 'lineage', by.y = 'assigned_lineage', all.x = TRUE)
# =============================================================================
# Classify priming and plasticity at day 10
# =============================================================================
metadat.day10_DABTRAM <- all_data.day10_DABTRAM@meta.data

fp.summary <- metadat.day10_DABTRAM %>% 
  group_by(assigned_lineage) %>%
  summarise(
    min_fatepotential = min(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    mean_fatepotential = mean(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    median_fatepotential = median(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    max_fatepotential = max(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    sd_fatepotential = sd(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    range_fatepotential = max(fatepotential_DABTRAM_d10_w5, na.rm = TRUE) - min(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
    n_cells = n()
  )

ggplot(fp.summary, aes(x = min_fatepotential, y = max_fatepotential)) +
  geom_point()

fp.summary$is.in.d0 <- ifelse(fp.summary$assigned_lineage %in% all_data.day0$assigned_lineage, TRUE, FALSE)
fp.summary <- merge(fp.summary, lineage_size.day0, by = 'assigned_lineage', all.x = TRUE)

ada_index.plot.DABTRAM$is.in.d0 <- ifelse(ada_index.plot.DABTRAM$lineage %in% all_data.day0$assigned_lineage, TRUE, FALSE)

ft.umap.DABTRAM.d10 <-  merge(ft.umap.DABTRAM, fp.d10_w5, by = 'cell_id')
ft.umap.DABTRAM.d10 <- ft.umap.DABTRAM.d10[order(ft.umap.DABTRAM.d10$fp.d10_w5, decreasing = FALSE), ]
ft.umap.DABTRAM.d0.w5 <- ft.umap.DABTRAM[ft.umap.DABTRAM$dataset %in% c('day0', 'week5_DABTRAM'), ]




max_val <- stats::quantile(ft.umap.DABTRAM.d10$fp.d10_w5, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(ft.umap.DABTRAM.d10$fp.d10_w5, 
                           probs = 0.01, 
                           na.rm = TRUE)
ft.umap.DABTRAM.d10$fp.d10_w5_scaled <- scales::rescale(
  pmax(pmin(ft.umap.DABTRAM.d10$fp.d10_w5, 
            max_val), 
       min_val),
  to = c(-1, 1)
)


lin.to.plot <- 'Lin61686'
ggplot(ft.umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = '#E8E8E8', size = 0.5) +
  geom_point(data = ft.umap.DABTRAM.d0.w5[ft.umap.DABTRAM.d0.w5$assigned_lineage == lin.to.plot, ], 
             aes(color = dataset), size = 5, shape = 18) +
  geom_point(data = ft.umap.DABTRAM.d10[ft.umap.DABTRAM.d10$assigned_lineage == lin.to.plot, ],
             aes(fill = fp.d10_w5_scaled), shape = 21, size = 5, stroke = 0) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'bisque', midpoint = 0, limits = c(min_val, max_val)) +
  xlab('') +
  ylab('') +
  theme_bw() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')

metadat.day10_DABTRAM[metadat.day10_DABTRAM$assigned_lineage == lin.to.plot, 'fatepotential_DABTRAM_d10_w5']

fp.summary.priming <- fp.summary[(fp.summary$min_fatepotential > 0) & (fp.summary$is.in.d0 == T), ]
sum(fp.summary.priming$n_cells.d0)
fp.summary.priming$category <- 'priming'

fp.summary.plastic <- fp.summary %>% 
  filter(min_fatepotential < 0 & 
         max_fatepotential > 0.5 &
         is.in.d0 == TRUE)
sum(fp.summary.plastic$n_cells.d0, na.rm = T)
fp.summary.plastic$category <- 'plastic'

lin.category.df <- rbind(
  fp.summary.plastic[, c('assigned_lineage', 'category')],
  fp.summary.priming[, c('assigned_lineage', 'category')]
)

write.csv(lin.category.df, '~/Downloads/lin_category_DABTRAM_1.csv', row.names = FALSE)

lin.to.plot <- fp.summary.priming$assigned_lineage

ggplot(ft.umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = '#E8E8E8', size = 0.5) +
  geom_point(data = ft.umap.DABTRAM.d0.w5[ft.umap.DABTRAM.d0.w5$assigned_lineage %in% lin.to.plot, ], 
             aes(color = dataset), size = 5, shape = 18) +
  # geom_point(data = ft.umap.DABTRAM.d10[ft.umap.DABTRAM.d10$assigned_lineage == lin.to.plot, ],
  #            aes(fill = fp.d10_w5_scaled), shape = 21, size = 5, stroke = 0) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'bisque', midpoint = 0, limits = c(min_val, max_val)) +
  xlab('') +
  ylab('') +
  theme_bw() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')
