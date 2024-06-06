library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================

cospar_out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/scripts/task5_cospar/'
cospar_day2to4 <- read.csv(paste0(cospar_out_dir, 'Hematopoiesis_day2.csv'))
cospar_all <- read.csv(paste0(cospar_out_dir, 'Hematopoiesis_all.csv'))

day4_lineage <- cospar_all %>% 
  filter(Time.point == 4) %>%
  group_by(assigned_lineage) %>%
  summarise(count_Mo = sum(Cell.type.annotation == 'Monocyte'),
            count_Neu = sum(Cell.type.annotation == 'Neutrophil'),
            count_Undiff = sum(Cell.type.annotation == 'Undifferentiated'))
day4_lineage$day4_Neu_bias <- day4_lineage$count_Neu / (day4_lineage$count_Mo + day4_lineage$count_Neu)

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Kevin_Fate_Potentials/'
load(paste0(in_dir, 'Writeup8_Neutrophil-4_from_day2_postprocess.RData'))
day2_to_Neu <- cell_imputed_score

load(paste0(in_dir, 'Writeup8_Monocyte-4_from_day2_postprocess.RData'))
day2_to_Mo <- cell_imputed_score

load(paste0(in_dir, 'Writeup8_Undifferentiated-4_from_day2_postprocess.RData'))
day2_to_Undiff <- cell_imputed_score

# ==============================================================================
# Compare
# ==============================================================================
cospar_day2to4_output_use <- cospar_day2to4[, c('SPRING.x', 'SPRING.y', 'time_info', 'state_info', 'fate_bias_intraclone_transition_map_Neutrophil.Monocyte',
                                                'fate_bias_transition_map_Neutrophil.Monocyte', "Cell.ID", "assigned_lineage")]


# colnames(cospar_day2to4_output_use) <- c('SPRING.x', 'SPRING.y', 'time_info', 'state_info', 'fate_bias_intraclone_transition_map_Neutrophil.Monocyte',
#                                          'fate_bias_transition_map_Neutrophil.Monocyte', 'Cell.ID', "Cell.barcode", "assigned_lineage")

day2_to_Neu <- as.data.frame(day2_to_Neu)
day2_to_Neu$Cell.ID <- rownames(day2_to_Neu)

day2_to_Mo <- as.data.frame(day2_to_Mo)
day2_to_Mo$Cell.ID <- rownames(day2_to_Mo)

day2_to_Undiff <- as.data.frame(day2_to_Undiff)
day2_to_Undiff$Cell.ID <- rownames(day2_to_Undiff)

cospar_day2to4_output_use <- merge(cospar_day2to4_output_use, day2_to_Neu, by = 'Cell.ID')
cospar_day2to4_output_use <- merge(cospar_day2to4_output_use, day2_to_Mo, by = 'Cell.ID')
cospar_day2to4_output_use <- merge(cospar_day2to4_output_use, day2_to_Undiff, by = 'Cell.ID')

cospar_day2to4_output_use$day2_to_Neu_Num <- 10^cospar_day2to4_output_use$day2_to_Neu
cospar_day2to4_output_use$day2_to_Mo_Num <- 10^cospar_day2to4_output_use$day2_to_Mo
cospar_day2to4_output_use$day2_to_Undiff_Num <- 10^cospar_day2to4_output_use$day2_to_Undiff

cospar_day2to4_output_use$Neu_Mo_bias <- cospar_day2to4_output_use$day2_to_Neu_Num / (cospar_day2to4_output_use$day2_to_Neu_Num + cospar_day2to4_output_use$day2_to_Mo_Num)

cospar_day2to4_output_use <- merge(cospar_day2to4_output_use, day4_lineage, by = 'assigned_lineage')
# cospar_day2to4_output_use$day4_Neu_bias <- cospar_day2to4_output_use$count_Neu / (cospar_day2to4_output_use$count_Mo + cospar_day2to4_output_use$count_Neu)



p1 <- ggplot(cospar_day2to4_output_use) +
  geom_point(aes(x = day4_Neu_bias, y = Neu_Mo_bias), alpha=0.5) +
  ggtitle('Fate Potential fate bias') +
  theme_bw()
p2 <- ggplot(cospar_day2to4_output_use) +
  geom_point(aes(x = day4_Neu_bias, y = fate_bias_intraclone_transition_map_Neutrophil.Monocyte), alpha=0.5) +
  ggtitle('CoSpar intraclone fate bias') +
  theme_bw()
p3 <- ggplot(cospar_day2to4_output_use) +
  geom_point(aes(x = day4_Neu_bias, y = fate_bias_transition_map_Neutrophil.Monocyte), alpha=0.5) +
  ggtitle('CoSpar state fate bias') +
  theme_bw()
p1 / p2 / p3
ggsave('~/Downloads/comp.png', width = 5, height = 12, dpi = 300)

par(mfrow = c(2, 2))
hist(cospar_day2to4_output_use$day4_Neu_bias, breaks = 100, main = 'Day 4 lineage bias')
hist(cospar_day2to4_output_use$Neu_Mo_bias, breaks = 100, main = 'Fate Potential fate bias')
hist(cospar_day2to4_output_use$fate_bias_intraclone_transition_map_Neutrophil.Monocyte, breaks = 100, main = 'CoSpar Intraclone fate bias')
hist(cospar_day2to4_output_use$fate_bias_transition_map_Neutrophil.Monocyte, breaks = 100, main = 'CoSpar state fate bias')

cor.test(cospar_day2to4_output_use$day4_Neu_bias, cospar_day2to4_output_use$fate_bias_transition_map_Neutrophil.Monocyte, method = 'spearman')


fate_potential_var <- cospar_day2to4_output_use %>% 
  group_by(assigned_lineage) %>% 
  summarise(var_FatePotential = var(Neu_Mo_bias),
            var_CoSpar_Intraclone = var(fate_bias_intraclone_transition_map_Neutrophil.Monocyte),
            var_CoSpar_State = var(fate_bias_transition_map_Neutrophil.Monocyte))

par(mfrow = c(2, 2))
hist(fate_potential_var$var_FatePotential, breaks = 100, main = 'Fate Potential fate bias')
hist(fate_potential_var$var_CoSpar_Intraclone, breaks = 100, main = 'CoSpar intraclone fate bias')
hist(fate_potential_var$var_CoSpar_State, breaks = 100, main = 'CoSpar fate bias')

p1 <- ggplot(fate_potential_var) +
  geom_point(aes(x = var_FatePotential, y = var_CoSpar_Intraclone), alpha=0.5) +
  ylim(0, 0.5) +
  xlim(0, 0.5) +
  theme_bw()
p2 <- ggplot(fate_potential_var) +
  geom_point(aes(x = var_FatePotential, y = var_CoSpar_State), alpha=0.5) +
  ylim(0, 0.5) +
  xlim(0, 0.5) +
  theme_bw()
p3 <- ggplot(fate_potential_var) +
  geom_point(aes(x = var_CoSpar_State, y = var_CoSpar_Intraclone), alpha=0.5) +
  ylim(0, 0.5) +
  xlim(0, 0.5) +
  theme_bw()
(p1 + p2) / (p3 + p3)


