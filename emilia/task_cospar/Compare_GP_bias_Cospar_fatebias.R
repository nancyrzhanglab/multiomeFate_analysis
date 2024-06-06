library(ggplot2)
library(Seurat)
# ==============================================================================
# Read data
# ==============================================================================

cospar_out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/scripts/task5_cospar/'
cospar_day2to4 <- read.csv(paste0(cospar_out_dir, 'Hematopoiesis_day2.csv'))
kevin_metadat <- read.csv(paste0(cospar_out_dir, 'Hematopoiesis_all_metadat_kevin_processed.csv'))

# cospar_day2to4_full <- merge(cospar_day2to4, kevin_metadat, by = c('SPRING.x', 'SPRING.y'))
cospar_day2to4_full <- cospar_day2to4

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Kevin_Fate_Potentials/'
load(paste0(in_dir, 'Writeup8_Neutrophil-4_from_day2_postprocess.RData'))
day2_to_Neu <- cell_imputed_score

load(paste0(in_dir, 'Writeup8_Monocyte-4_from_day2_postprocess.RData'))
day2_to_Mo <- cell_imputed_score

load(paste0(in_dir, 'Writeup8_Undifferentiated-4_from_day2_postprocess.RData'))
day2_to_Undiff <- cell_imputed_score

# ==============================================================================
# Sanity check
# ==============================================================================
ggplot(cospar_day2to4_full, aes(x = `SPRING.x`, y = `SPRING.y`)) +
  geom_point(aes(color = fate_bias_intraclone_transition_map_Neutrophil.Monocyte))

hist(cospar_day2to4_full$`fate_bias_intraclone_transition_map_Neutrophil.Monocyte`)

# ==============================================================================
# Compare
# ==============================================================================
cospar_day2to4_output_use <- cospar_day2to4_full[, c('SPRING.x', 'SPRING.y', 'time_info', 'state_info', 'fate_bias_intraclone_transition_map_Neutrophil.Monocyte',
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
hist(cospar_day2to4_output_use$Neu_Mo_bias, breaks = 100)
hist(cospar_day2to4_output_use$fate_bias_transition_map_Neutrophil.Monocyte)

ggplot(cospar_day2to4_output_use, aes(x = `SPRING.x`, y = `SPRING.y`)) +
  geom_point(aes(color = Neu_Mo_bias), size = 0.5) +
  scale_color_gradient2(low = "#3e76fa",midpoint = 0.5, mid = 'gray', high = "#e02626") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank())

# ggsave('~/Downloads/FatePotential_fate_bias_day2_to_day4.png', width = 6.5, height = 4.5, dpi = 300)

summary(cospar_day2to4_output_use)
ggplot(cospar_day2to4_output_use, aes(x = `Neu_Mo_bias`, y = `fate_bias_transition_map_Neutrophil.Monocyte`)) +
  geom_point(alpha = 0.5) +
  theme_bw()

cor.test(cospar_day2to4_output_use$Neu_Mo_bias,
         cospar_day2to4_output_use$fate_bias_intraclone_transition_map_Neutrophil.Monocyte)

cospar_day2to4_output_use$diff <- cospar_day2to4_output_use$fate_bias_transition_map_Neutrophil.Monocyte - cospar_day2to4_output_use$Neu_Mo_bias
# cospar_day2to4_output_use$diff <- cospar_day2to4_output_use$Neu_Mo_bias - cospar_day2to4_output_use$fate_bias_transition_map_Neutrophil.Monocyte


to_plot <- cospar_day2to4_output_use[cospar_day2to4_output_use$diff > 0.3 | cospar_day2to4_output_use$diff < -0.3, ]
ggplot(to_plot, aes(x = `SPRING.x`, y = `SPRING.y`, color = diff)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red") +
  theme_bw()+
  theme(panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank())
ggsave('~/Downloads/CoSpar_fate_bias_MINUS_FatePotential_fate_bias_day2_to_day4_ALL.png', width = 6.5, height = 4.5, dpi = 300)

par(mfrow = c(1, 1))
hist(cospar_day2to4_output_use$diff, breaks = 100 )

lineage_size <- cospar_day2to4_output_use %>% 
  group_by(assigned_lineage) %>% 
  summarise(lineage_size = dplyr::n())

cospar_day2to4_output_use <- merge(cospar_day2to4_output_use, lineage_size, by = 'assigned_lineage')
ggplot(cospar_day2to4_output_use, aes(x = lineage_size, y = diff)) +
  geom_point() +
  xlim(0, 10) +
  ylab('CoSpar_fate_bias - Fate_potential_fate_bias') +
  theme_bw()

# ==============================================================================
# Correlation with gene expression
# ==============================================================================

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/LARRY_hematopoiesis/'
load(paste0(in_dir,'/Writeup8_larry-dataset_step2_lineage-plotting.RData'))

seurat_object <- FindVariableFeatures(seurat_object)
day2 <- subset(seurat_object, Time.point == 2)
day2.RNA.meta <- day2@meta.data
var_features <- day2@assays[["RNA"]]@meta.data[["var.features"]]
var_features <- var_features[!is.na(var_features)]
day2.RNA <- day2@assays[["RNA"]]@layers[["scale.data"]]
day2.RNA <- as.data.frame(day2.RNA)
rownames(day2.RNA) <- var_features
colnames(day2.RNA) <- rownames(day2.RNA.meta)
day2.RNA <- as.data.frame(t(day2.RNA))
day2.RNA$Cell.ID <- rownames(day2.RNA)
day2.RNA <- merge(day2.RNA, cospar_day2to4_output_use, by = 'Cell.ID')

cor.df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(cor.df) <- c('gene', 'cor.CoSpar.state', 'p.CoSpar.state', 'cor.CoSpar.intraclone', 'p.CoSpar.intraclone', 'cor.FP', 'p.FP')
for (g in var_features) {
  df <- day2.RNA[, c('Cell.ID', g, 'fate_bias_transition_map_Neutrophil.Monocyte', 'fate_bias_intraclone_transition_map_Neutrophil.Monocyte', 'Neu_Mo_bias')]
  res.CoSpar.state <- cor.test(df[[g]], df$`fate_bias_transition_map_Neutrophil.Monocyte`)
  cor.CoSpar.state <- res.CoSpar.state[["estimate"]][["cor"]]
  p.CoSpar.state <- res.CoSpar.state[["p.value"]]
  
  res.CoSpar.intraclone <- cor.test(df[[g]], df$`fate_bias_intraclone_transition_map_Neutrophil.Monocyte`)
  cor.CoSpar.intraclone <- res.CoSpar.intraclone[["estimate"]][["cor"]]
  p.CoSpar.intraclone <- res.CoSpar.intraclone[["p.value"]]
  
  res.FP <- cor.test(df[[g]], df$`Neu_Mo_bias`)
  cor.FP <- res.FP[["estimate"]][["cor"]]
  p.FP <- res.FP[["p.value"]]
  
  cor.df[nrow(cor.df) + 1, ] <- c(g, cor.CoSpar.state, p.CoSpar.state, cor.CoSpar.intraclone, p.CoSpar.intraclone, cor.FP, p.FP)
}


neu.genes <- c('Camp', 'Itgb2l', 'Ngp', 'Ltf', 'Ceacam10', 'G0s2', 'Cd177', 'Chil3', 'S100a9', 'Elane', 'Mpo', 'Lipg', 'Gfi1', 'Syne1', 'Lcn2', 'Gstm1', 'Ctsg')
mo.genes <- c('S100a4', 'Ctss', 'Fabp5', 'Psap', 'Ctsc', 'Ahnak', 'Spp1', 'Mmp8', 'Gpnmb', 'Saa3', 'Mrc1', 'Lpl', 'Mmp12', 'Tgfbi', 'Ms4a6d', 'Wfdc17', 'Clec4n')

neu.genes <- c('Cd15', 'Cd11b', 'Cd16', 'Cd10')
cor.df$cor.CoSpar.state <- as.numeric(cor.df$cor.CoSpar.state)
cor.df$cor.CoSpar.intraclone <- as.numeric(cor.df$cor.CoSpar.intraclone)
cor.df$cor.FP <- as.numeric(cor.df$cor.FP)
cor.df.neu.subset <- cor.df[cor.df$gene %in% neu.genes, ]
cor.df.mo.subset <- cor.df[cor.df$gene %in% mo.genes, ]

res <- cor.test(cor.df$cor.CoSpar.state, cor.df$cor.FP)
p1 <- ggplot(cor.df, aes(x = cor.CoSpar.state, y = cor.FP)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = cor.df.neu.subset, color = 'red') +
  geom_point(data = cor.df.mo.subset, color = 'blue') +
  ggtitle(paste0('r = ', res[['estimate']][['cor']])) +
  xlim(-0.7, 0.7) +
  ylim(-0.7, 0.7) +
  theme_bw()

res <- cor.test(cor.df$cor.CoSpar.intraclone, cor.df$cor.FP)
p2 <- ggplot(cor.df, aes(x = cor.CoSpar.intraclone, y = cor.FP)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = cor.df.neu.subset, color = 'red') +
  geom_point(data = cor.df.mo.subset, color = 'blue') +
  ggtitle(paste0('r = ', res[['estimate']][['cor']])) +
  xlim(-0.7, 0.7) +
  ylim(-0.7, 0.7) +
  theme_bw()

res <- cor.test(cor.df$cor.CoSpar.intraclone, cor.df$cor.CoSpar.state)
p3 <- ggplot(cor.df, aes(x = cor.CoSpar.state, y = cor.CoSpar.intraclone)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = cor.df.neu.subset, color = 'red') +
  geom_point(data = cor.df.mo.subset, color = 'blue') +
  ggtitle(paste0('r = ', res[['estimate']][['cor']])) +
  xlim(-0.7, 0.7) +
  ylim(-0.7, 0.7) +
  theme_bw()

p1 + p2 + p3


