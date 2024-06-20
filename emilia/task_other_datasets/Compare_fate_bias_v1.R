library(ggpubr)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
metadat <- seurat_object@meta.data
metadat$cell_barcode <- rownames(metadat)

d14_fate <- 'high' # high = non-cycling, median = mid, low = cycling
current <- 0
future <- 14
load(paste0(out_dir, 'PC9_time_course_d', current, '_to_d14_', d14_fate, '.RData'))
non_cycling_fate <- final_fit[["cell_imputed_score"]]
non_cycling_fate <- 10**(non_cycling_fate)

d14_fate <- 'mid_low' # high = non-cycling, median = mid, low = cycling
current <- 0
future <- 14
load(paste0(out_dir, 'PC9_time_course_d', current, '_to_d14_', d14_fate, '.RData'))
cycling_fate <- final_fit[["cell_imputed_score"]]
cycling_fate <- 10**(cycling_fate)


hist(non_cycling_fate, right = FALSE, breaks = 100)
hist(cycling_fate, right = FALSE, xlim = c(1.3, 1.32), breaks = 100)

non_cycling_fate_df <- as.data.frame(non_cycling_fate)
non_cycling_fate_df$cell_barcode <- rownames(non_cycling_fate_df)

cycling_fate_df <- as.data.frame(cycling_fate)
cycling_fate_df$cell_barcode <- rownames(cycling_fate_df)

fate_potential_df <- merge(cycling_fate_df, non_cycling_fate_df, by = 'cell_barcode')

ggplot(fate_potential_df, aes(x = cycling_fate, y = non_cycling_fate)) +
  geom_point(size = 0.5) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  theme_bw()


# ==============================================================================
# D14 status (observed)
# ==============================================================================
metadat_d14 <- metadat[metadat$time_point == 14, ]
metadat_d14$lineage_fate <- ifelse(metadat_d14$sample_type == '14_high', 'non-cycling', 'cycling')
sample_type_size <- metadat_d14 %>% 
  group_by(lineage_fate) %>% 
  summarise(n_cells = n())

lineage_fate_d14 <- metadat_d14 %>% 
  group_by(lineage_barcode, lineage_fate, majority_fate) %>% 
  summarise(n_cells = n())
lineage_size_d14 <- metadat_d14 %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_total = n())

lineage_fate_d14 <- merge(lineage_fate_d14, lineage_size_d14, by = 'lineage_barcode')
lineage_fate_d14$freq <- lineage_fate_d14$n_cells / lineage_fate_d14$n_total


# ==============================================================================
# D14 status (predicted)
# ==============================================================================
fate_potential_df <- merge(metadat[, c('cell_barcode', 'lineage_barcode')], fate_potential_df, by = 'cell_barcode')
predicted_lineage_size <- fate_potential_df %>% 
  group_by(lineage_barcode) %>% 
  summarise(n_cycling = sum(cycling_fate),
            n_non_cycling = sum(non_cycling_fate))

predicted_lineage_size$n_total <- predicted_lineage_size$n_cycling + predicted_lineage_size$n_non_cycling
predicted_lineage_size$freq_cycling <- predicted_lineage_size$n_cycling / predicted_lineage_size$n_total
predicted_lineage_size$freq_non_cycling <- predicted_lineage_size$n_non_cycling / predicted_lineage_size$n_total

length(intersect(unique(lineage_fate_d14$lineage_barcode), unique(predicted_lineage_size$lineage_barcode)))

lineage_fate_d14_cycling <- lineage_fate_d14[lineage_fate_d14$lineage_fate == 'cycling', ]
# lineage_fate_d14_cycling <- lineage_fate_d14_cycling[lineage_fate_d14_cycling$n_total > 10, ]
comp <- merge(lineage_fate_d14_cycling[, c('lineage_barcode', 'freq')], predicted_lineage_size[, c('lineage_barcode', 'freq_cycling')], by = 'lineage_barcode')

ggplot(comp, aes(x = freq, y = freq_cycling)) +
  geom_point(size = 0.5) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab('Observed cycling fate percentage') +
  ylab('Predicted cycling fate percentage') +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()

lineage_fate_d14_non_cycling <- lineage_fate_d14[lineage_fate_d14$lineage_fate == 'non-cycling', ]
# lineage_fate_d14_non_cycling <- lineage_fate_d14_non_cycling[lineage_fate_d14_non_cycling$n_total > 10, ]
comp <- merge(lineage_fate_d14_non_cycling[, c('lineage_barcode', 'freq')], predicted_lineage_size[, c('lineage_barcode', 'freq_non_cycling')], by = 'lineage_barcode')

ggplot(comp, aes(x = freq, y = freq_non_cycling)) +
  geom_point(size = 0.5) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab('Observed non-cycling fate percentage') +
  ylab('Predicted non-cycling fate percentage') +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw()
