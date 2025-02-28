rm(list = rm())
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
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver

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

# write.csv(scores.df, paste0(out_dir, 'saver.GavishMP.UCellScores.csv'))

# ==============================================================================
# Subset to day10
# ==============================================================================
metadat.day10.dabtram <- metadat[metadat$dataset == 'day10_DABTRAM', ]
scores.df.day10.dabtram <- scores.df[scores.df$cell_id %in% rownames(metadat.day10.dabtram), ]

fatepot.d10.w5 <- c('fatepotential_DABTRAM_d0_d10', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d10_w5')
scores.df.day10.dabtram <- scores.df.day10.dabtram[, -which(names(scores.df.day10.dabtram) %in% fatepot.d10.w5)]

metadat.day10.cocl2 <- metadat[metadat$dataset == 'day10_COCL2', ]
scores.df.day10.cocl2 <- scores.df[scores.df$cell_id %in% rownames(metadat.day10.cocl2), ]

fatepot.d10.w5 <- c('fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d0_d10', 'fatepotential_CIS_d10_w5')
scores.df.day10.cocl2 <- scores.df.day10.cocl2[, -which(names(scores.df.day10.cocl2) %in% fatepot.d10.w5)]

metadat.day10.cis <- metadat[metadat$dataset == 'day10_CIS', ]
scores.df.day10.cis <- scores.df[scores.df$cell_id %in% rownames(metadat.day10.cis), ]

fatepot.d10.w5 <- c('fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d0_d10')
scores.df.day10.cis <- scores.df.day10.cis[, -which(names(scores.df.day10.cis) %in% fatepot.d10.w5)]

# ==============================================================================
# Correlation test
# ==============================================================================

# DABTRAM
metadat.day10.dabtram <- metadat.day10.dabtram[scores.df.day10.dabtram$cell_id, ]
cor.vec.dabtram <- sapply(colnames(scores), function(x) {
  res <- cor.test(scores.df.day10.dabtram[, x], metadat.day10.dabtram$fatepotential_DABTRAM_d10_w5, method = 'spearman')
  return(res$estimate)
})
cor.vec.dabtram <- as.data.frame(cor.vec.dabtram)

# COCL2
metadat.day10.cocl2 <- metadat.day10.cocl2[scores.df.day10.cocl2$cell_id, ]
cor.vec.cocl2 <- sapply(colnames(scores), function(x) {
  res <- cor.test(scores.df.day10.cocl2[, x], metadat.day10.cocl2$fatepotential_COCL2_d10_w5, method = 'spearman')
  return(res$estimate)
})
cor.vec.cocl2 <- as.data.frame(cor.vec.cocl2)

# CIS
metadat.day10.cis <- metadat.day10.cis[scores.df.day10.cis$cell_id, ]
cor.vec.cis <- sapply(colnames(scores), function(x) {
  res <- cor.test(scores.df.day10.cis[, x], metadat.day10.cis$fatepotential_CIS_d10_w5, method = 'spearman')
  return(res$estimate)
})
cor.vec.cis <- as.data.frame(cor.vec.cis)

cor.df <- merge(cor.vec.dabtram, cor.vec.cocl2, by='row.names')
colnames(cor.df) <- c('MetaProgram', 'DABTRAM.D10.W5', 'COCL2.D10.W5')

cor.vec.cis$MetaProgram <- rownames(cor.vec.cis)
cor.df <- merge(cor.df, cor.vec.cis, by='MetaProgram')
colnames(cor.df) <- c('MetaProgram', 'DABTRAM.D10.W5', 'COCL2.D10.W5', 'CIS.D10.W5')

ggpairs(cor.df[, -1]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = 'Correlation of MetaPrograms with Fate Potential (Day10 to Week5)')

cor.df.melt <- melt(cor.df, id.vars = 'MetaProgram')
cor.df.melt$MetaProgram <- as.character(cor.df.melt$MetaProgram)

ggplot(cor.df.melt, aes(y = reorder(MetaProgram, value), x =  value)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~variable, scales = 'free_x') +
  # xlim(-0.2, 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = 'Correlation of MetaPrograms with Fate Potential (Day10 to Week5)',
       x = 'Correlation',
       y = 'MetaProgram')

write.csv(cor.df, paste0(result_dir, 'GavishMP_UCell_cor_d10_w5.csv'), row.names = F)

