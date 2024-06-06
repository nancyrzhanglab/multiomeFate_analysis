library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================

cospar_out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/scripts/task5_cospar/'
cospar_day4to6 <- read.csv(paste0(cospar_out_dir, 'Hematopoiesis_day4.csv'))


in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Kevin_Fate_Potentials/'
load(paste0(in_dir, 'Writeup8_Neutrophil-6_from_day4_postprocess.RData'))
day4_to_Neu <- cell_imputed_score

load(paste0(in_dir, 'Writeup8_Monocyte-6_from_day4_postprocess.RData'))
day4_to_Mo <- cell_imputed_score

load(paste0(in_dir, 'Writeup8_Undifferentiated-6_from_day4_postprocess.RData'))
day4_to_Undiff <- cell_imputed_score

# ==============================================================================
# Sanity check
# ==============================================================================
ggplot(cospar_day4to6, aes(x = `SPRING.x`, y = `SPRING.y`)) +
  geom_point(aes(color = fate_bias_intraclone_transition_map_Neutrophil.Monocyte))

# ==============================================================================
# Compare
# ==============================================================================
cospar_day4to6_output_use <- cospar_day4to6[, c('SPRING.x', 'SPRING.y', 'time_info', 'state_info', 'fate_bias_intraclone_transition_map_Neutrophil.Monocyte',
                                                     'fate_bias_transition_map_Neutrophil.Monocyte', "Cell.ID", "assigned_lineage")]

day4_to_Neu <- as.data.frame(day4_to_Neu)
day4_to_Neu$Cell.ID <- rownames(day4_to_Neu)

day4_to_Mo <- as.data.frame(day4_to_Mo)
day4_to_Mo$Cell.ID <- rownames(day4_to_Mo)

day4_to_Undiff <- as.data.frame(day4_to_Undiff)
day4_to_Undiff$Cell.ID <- rownames(day4_to_Undiff)

cospar_day4to6_output_use <- merge(cospar_day4to6_output_use, day4_to_Neu, by = 'Cell.ID')
cospar_day4to6_output_use <- merge(cospar_day4to6_output_use, day4_to_Mo, by = 'Cell.ID')
cospar_day4to6_output_use <- merge(cospar_day4to6_output_use, day4_to_Undiff, by = 'Cell.ID')

cospar_day4to6_output_use$day4_to_Neu_Num <- 10^cospar_day4to6_output_use$day4_to_Neu
cospar_day4to6_output_use$day4_to_Mo_Num <- 10^cospar_day4to6_output_use$day4_to_Mo
cospar_day4to6_output_use$day4_to_Undiff_Num <- 10^cospar_day4to6_output_use$day4_to_Undiff

cospar_day4to6_output_use$Neu_Mo_bias <- cospar_day4to6_output_use$day4_to_Neu_Num / (cospar_day4to6_output_use$day4_to_Neu_Num + cospar_day4to6_output_use$day4_to_Mo_Num)
hist(cospar_day4to6_output_use$Neu_Mo_bias)

ggplot(cospar_day4to6_output_use, aes(x = `Neu_Mo_bias`, y = `fate_bias_transition_map_Neutrophil.Monocyte`)) +
  geom_point(alpha = 0.2)

cor.test(cospar_day4to6_output_use$Neu_Mo_bias, 
         cospar_day4to6_output_use$fate_bias_transition_map_Neutrophil.Monocyte,
         method = 'spearman')


cospar_day4to6_output_use$diff <- cospar_day4to6_output_use$fate_bias_transition_map_Neutrophil.Monocyte - cospar_day4to6_output_use$Neu_Mo_bias

to_plot <- cospar_day4to6_output_use[cospar_day4to6_output_use$diff > 0.3 | cospar_day4to6_output_use$diff < -0.3, ]
ggplot(cospar_day4to6_output_use, aes(x = `SPRING.x`, y = `SPRING.y`, color = diff)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red") +
  theme_bw()

nrow(to_plot) / nrow(cospar_day4to6_output_use)
