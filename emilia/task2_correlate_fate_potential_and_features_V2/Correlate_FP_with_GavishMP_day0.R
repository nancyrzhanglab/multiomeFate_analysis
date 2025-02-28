rm(list = ls())
library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggdensity)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '~/Downloads/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
# load('~/Downloads/Writeup10a_data_wnn.RData')

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
# all_data[["wnn.umap"]] <- all_data_wnn

gavish.mp <- read_csv(paste0(ref_dir, 'Gavish_Malignant_Meta_Programs.csv'))
isg_rs <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/ISG.RS.txt', col_names = F)
isg_mem <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Persistent IFN ISG Groups/Memory ISGs Human.csv', col_names = F)

colnames(isg_rs) <- 'Gene'
colnames(isg_mem) <- 'Gene'

mitf_axl <- read.csv('~/Downloads/Triosh_2016_MITF_AXL_Programs.csv')

# ==============================================================================
# Wrangle data
# ==============================================================================
gavish.mp.list <- list()
for (mp in colnames(gavish.mp)) {
  gavish.mp.list[[mp]] <- c(gavish.mp[[mp]] %>% as.character())
}

gavish.mp.list <- lapply(1:length(gavish.mp.list), function(x) {
  gs <- gavish.mp.list[[x]]
  gs <- unique(na.omit(gs[gs %in% all_data_saver@data@Dimnames[[1]]]))
  return(gs)
})
names(gavish.mp.list) <- colnames(gavish.mp)

isg_rs <- unique(na.omit(isg_rs$Gene[isg_rs$Gene %in% all_data_saver@data@Dimnames[[1]]]))
isg_mem <- unique(na.omit(isg_mem$Gene[isg_mem$Gene %in% all_data_saver@data@Dimnames[[1]]]))

mitf_program <- unique(na.omit(mitf_axl$MITF_Program[mitf_axl$MITF_Program %in% all_data_saver@data@Dimnames[[1]]]))
axl_program <- unique(na.omit(mitf_axl$AXL_Program[mitf_axl$AXL_Program %in% all_data_saver@data@Dimnames[[1]]]))

# ==============================================================================
# Calcualte module scores 
# ==============================================================================

metadat <- all_data@meta.data
scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=c(gavish.mp.list, 
                                           list(isg_rs), 
                                           list(isg_mem), 
                                           list(mitf_program), 
                                           list(axl_program)))

colnames(scores) <- c(names(gavish.mp.list), 'ISG.RS', 'ISG.Mem', 'MITF_Program', 'AXL_Program')
scores.df <- as.data.frame(scores) 

# ==============================================================================
# Compare with fate potentials
# ==============================================================================
fatepot.DABTRAM.d0.d10 <- all_data_fatepotential[['fatepotential_DABTRAM_d0_d10']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.DABTRAM.d0.d10) <- 'fatepotential_DABTRAM_d0_d10'

fatepot.DABTRAM.d10.w5 <- all_data_fatepotential[['fatepotential_DABTRAM_d10_w5']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.DABTRAM.d10.w5) <- 'fatepotential_DABTRAM_d10_w5'
fatepot.DABTRAM.d10.w5$cell_id <- rownames(fatepot.DABTRAM.d10.w5)

fatepot.COCL2.d0.d10 <- all_data_fatepotential[['fatepotential_COCL2_d0_d10']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.COCL2.d0.d10) <- 'fatepotential_COCL2_d0_d10'
fatepot.COCL2.d0.d10$cell_id <- rownames(fatepot.COCL2.d0.d10)

fatepot.COCL2.d10.w5 <- all_data_fatepotential[['fatepotential_COCL2_d10_w5']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.COCL2.d10.w5) <- 'fatepotential_COCL2_d10_w5'
fatepot.COCL2.d10.w5$cell_id <- rownames(fatepot.COCL2.d10.w5)

fatepot.CIS.d0.d10 <- all_data_fatepotential[['fatepotential_CIS_d0_d10']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.CIS.d0.d10) <- 'fatepotential_CIS_d0_d10'
fatepot.CIS.d0.d10$cell_id <- rownames(fatepot.CIS.d0.d10)

fatepot.CIS.d10.w5 <- all_data_fatepotential[['fatepotential_CIS_d10_w5']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.CIS.d10.w5) <- 'fatepotential_CIS_d10_w5'
fatepot.CIS.d10.w5$cell_id <- rownames(fatepot.CIS.d10.w5)

# add fate potential
scores.df <- merge(scores.df, fatepot.DABTRAM.d0.d10, by='row.names', all = T)
colnames(scores.df)[1] <- 'cell_id'
scores.df <- merge(scores.df, fatepot.DABTRAM.d10.w5, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.COCL2.d0.d10, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.COCL2.d10.w5, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.CIS.d0.d10, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.CIS.d10.w5, by='cell_id', all = T)

# ==============================================================================
# Subset to day0
# ==============================================================================
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
scores.df.day0 <- scores.df[scores.df$cell_id %in% rownames(metadat.day0), ]

fatepot.d10.w5 <- c('fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d10_w5')
scores.df.day0 <- scores.df.day0[, -which(names(scores.df.day0) %in% fatepot.d10.w5)]

# ==============================================================================
# Correlation test
# ==============================================================================

# DABTRAM
cor.vec.dabtram <- sapply(colnames(scores), function(x) {
  res <- cor.test(scores.df.day0[, x], scores.df.day0$fatepotential_DABTRAM_d0_d10, method = 'spearman')
  return(res$estimate)
})
cor.vec.dabtram <- as.data.frame(cor.vec.dabtram)

# COCL2
cor.vec.cocl2 <- sapply(colnames(scores), function(x) {
  res <- cor.test(scores.df.day0[, x], scores.df.day0$fatepotential_COCL2_d0_d10, method = 'spearman')
  return(res$estimate)
})
cor.vec.cocl2 <- as.data.frame(cor.vec.cocl2)

# CIS
cor.vec.cis <- sapply(colnames(scores), function(x) {
  res <- cor.test(scores.df.day0[, x], scores.df.day0$fatepotential_CIS_d0_d10, method = 'spearman')
  return(res$estimate)
})
cor.vec.cis <- as.data.frame(cor.vec.cis)

cor.df <- merge(cor.vec.dabtram, cor.vec.cocl2, by='row.names')
colnames(cor.df) <- c('MetaProgram', 'DABTRAM.D0.D10', 'COCL2.D0.D10')

cor.vec.cis$MetaProgram <- rownames(cor.vec.cis)
cor.df <- merge(cor.df, cor.vec.cis, by='MetaProgram')
colnames(cor.df) <- c('MetaProgram', 'DABTRAM.D0.D10', 'COCL2.D0.D10', 'CIS.D0.D10')

ggpairs(cor.df[, -1]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = 'Correlation of MetaPrograms with Fate Potential (Day0 to Day10)')

cor.df.melt <- melt(cor.df, id.vars = 'MetaProgram')
cor.df.melt$MetaProgram <- as.character(cor.df.melt$MetaProgram)

ggplot(cor.df.melt, aes(y = reorder(MetaProgram, value), x =  value)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~variable, scales = 'free_x') +
  xlim(-0.5, 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = 'Correlation of MetaPrograms with Fate Potential (Day0 to Day10)',
       x = 'Correlation',
       y = 'MetaProgram')

write.csv(cor.df, paste0(result_dir, 'GavishMP_UCell_cor_d0_d10.csv'), row.names = F)


# ==============================================================================
# Plot feature
# ==============================================================================
scores.df.day0.plot <- scores.df.day0[, c('cell_id', 'Chromatin', 'fatepotential_DABTRAM_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_CIS_d0_d10')]
scores.df.day0.plot <- melt(scores.df.day0.plot, id.vars = c('cell_id', 'Chromatin'))
ggplot(scores.df.day0.plot, aes(x = Chromatin, y = value)) + 
  geom_point(alpha = 1) +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  facet_wrap(~variable, ncol = 3) +
  stat_cor() +
  # geom_smooth(method = 'lm', alpha = 0.1) +
  ylab('Fate potential') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12))

umap <- all_data_wnn@cell.embeddings
umap <- as.data.frame(umap)
umap$cell_id <- rownames(umap)

scores.df.day0.plot <- merge(scores.df.day0.plot, umap, by='cell_id')
ggplot(scores.df.day0.plot, aes(x = wnnUMAP_1, y = wnnUMAP_2)) + 
  geom_point(aes(color = value), alpha = 0.8) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  facet_wrap(~variable, ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12))

scores.df.chromatin <- scores.df[, c('cell_id', 'Chromatin')]
scores.df.chromatin$Chromatin_Scaled <- scale(scores.df.chromatin$Chromatin)
colnames(scores.df.chromatin) <- c('cell_id', 'Chromatin', 'Chromatin_Scaled')

scores.df.chromatin <- merge(scores.df.chromatin, umap, by='cell_id')
metadat$cell_id <- rownames(metadat)
scores.df.chromatin <- merge(scores.df.chromatin, metadat[, c("cell_id", "dataset", "Phase", "assigned_lineage")], by='cell_id')
ggplot(scores.df.chromatin, aes(x = wnnUMAP_1, y = wnnUMAP_2)) + 
  geom_point(aes(color = Chromatin_Scaled), alpha = 0.8) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  facet_wrap(~dataset, ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12))

ggplot(scores.df.chromatin, aes(x = Phase, y = Chromatin_Scaled)) + 
  geom_violin() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12))

ggplot(scores.df.chromatin, aes(x = dataset, y = Chromatin_Scaled)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  theme_bw()

# ==============================================================================
# Low FP cells
# ==============================================================================
Lin.day0 <- metadat[metadat$dataset == 'day0', 'assigned_lineage']
Lin.day10.dabtram <- metadat[metadat$dataset == 'day10_DABTRAM', 'assigned_lineage']
Lin.week5.dabtram <- metadat[metadat$dataset == 'week5_DABTRAM', 'assigned_lineage']
in.d0.w5 <- intersect(Lin.day0, Lin.week5.dabtram)
in.d0.d10.w5 <- intersect(in.d0.w5, Lin.day10.dabtram)

metadat.day0.fatepot <- metadat[metadat$dataset == 'day0', c('cell_id', 'fatepotential_DABTRAM_d0_d10', 'assigned_lineage')]
metadat.day0.fatepot$inWeek5 <- ifelse(metadat.day0.fatepot$assigned_lineage %in%Lin.week5.dabtram, 'Yes', 'No')
metadat.day0.fatepot <- metadat.day0.fatepot %>% 
  dplyr::group_by(assigned_lineage, inWeek5) %>%
  dplyr::summarise(mean_fatepotential = mean(fatepotential_DABTRAM_d0_d10),
                   max_fatepotential = max(fatepotential_DABTRAM_d0_d10),
                   var_fatepotential = var(fatepotential_DABTRAM_d0_d10),
                   n = n())
metadat.day0.fatepot <- metadat.day0.fatepot[metadat.day0.fatepot$n > 3,]

ggplot(metadat.day0.fatepot, aes(x = inWeek5, y = n)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.1) +
  geom_point() +
  theme_bw()

lin.size <- metadat %>% 
  dplyr::group_by(assigned_lineage, dataset) %>%
  dplyr::summarise(n = n()) %>% 
  spread(key = dataset, value = n)
lin.size.dabtram <- lin.size[, c('assigned_lineage', 'day0', 'week5_DABTRAM')]
lin.size.dabtram <- lin.size.dabtram %>% drop_na()
ggplot(lin.size.dabtram, aes(x = day0, y  = week5_DABTRAM)) +
  geom_point() +
  theme_bw()

metadat.day0$cell_id <- rownames(metadat.day0)
metadat.day0.umap <- merge(metadat.day0, umap, by='cell_id')
metadat.day0.umap$inWeek5 <- ifelse(metadat.day0.umap$assigned_lineage %in% Lin.week5.dabtram, 'Yes', 'No')
metadat.day0.umap$inDay10Week5 <- ifelse(metadat.day0.umap$assigned_lineage %in% in.d0.d10.w5, 'Yes', 'No')

ggplot(metadat.day0.umap, aes(x = wnnUMAP_1, y = wnnUMAP_2)) + 
  geom_point(aes(color = inDay10Week5), alpha = 0.8) +
  scale_color_manual(values = c('No' = 'blue', 'Yes' = 'red')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12))
