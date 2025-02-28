rm(list=ls())
library(Seurat)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data[['saver']] <- all_data_saver

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# ==============================================================================
# Get high fp cells from d10 and cells in w5
# ==============================================================================

# d0
metadat.d0 <- metadat %>% filter(dataset == 'day0')

# high fp d10
fp.dabtram.d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(fp.dabtram.d10_w5) <- c('fatepotential_DABTRAM_d10_w5')
fp.dabtram.d10_w5$cell_id <- rownames(fp.dabtram.d10_w5)
high_fp_d10 <- fp.dabtram.d10_w5 %>% filter(fatepotential_DABTRAM_d10_w5 > 0)

high_fp_d10 <- merge(high_fp_d10, metadat, by=c('fatepotential_DABTRAM_d10_w5', 'cell_id'))

low_fp_d10 <- fp.dabtram.d10_w5 %>% filter(fatepotential_DABTRAM_d10_w5 < 0)
low_fp_d10 <- merge(low_fp_d10, metadat, by=c('fatepotential_DABTRAM_d10_w5', 'cell_id'))

# week5 
metadat.week5_DABTRAM <- metadat %>% filter(dataset == 'week5_DABTRAM')

metadat.DABTRAM <- read.csv(paste0(data_dir, 'metadat.all_data_dabtram_recluster.csv'), row.names = 1)
metadat.week5_DABTRAM <- metadat.DABTRAM[metadat.DABTRAM$dataset == 'week5_DABTRAM', ]
metadat.week5_DABTRAM.clust0 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 0)
metadat.week5_DABTRAM.clust3 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 3)

# ==============================================================================
# Signatures
# ==============================================================================

gavish.mp <- read_csv(paste0(ref_dir, 'Gavish_Malignant_Meta_Programs.csv'))
isg_rs <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/ISG.RS.txt', col_names = F)
isg_mem <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Persistent IFN ISG Groups/Memory ISGs Human.csv', col_names = F)

colnames(isg_rs) <- 'Gene'
colnames(isg_mem) <- 'Gene'

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

# ==============================================================================
# Calculate module scores 
# ==============================================================================

metadat <- all_data@meta.data

scores <- ScoreSignatures_UCell(all_data@assays[["saver"]]@data, 
                                features=c(gavish.mp.list, 
                                           list(isg_rs), 
                                           list(isg_mem)))
colnames(scores) <- c(names(gavish.mp.list), 'ISG.RS', 'ISG.Mem')
scores.df <- as.data.frame(scores) 


# ==============================================================================
# Plot module scores
# ==============================================================================
# scores.df.scaled <- scale(scores.df)
# scores.df.scaled <- as.data.frame(scores.df.scaled)

scores.df.d0 <- scores.df[metadat.d0$cell_id, ]
scores.df.high.fp.d10 <- scores.df[high_fp_d10$cell_id, ]
scores.df.low.fp.d10 <- scores.df[low_fp_d10$cell_id, ]
scores.df.week5 <- scores.df[metadat.week5_DABTRAM$cell_id, ]
scores.df.week5.cl0 <- scores.df[metadat.week5_DABTRAM.clust0$cell_id, ]
scores.df.week5.cl3 <- scores.df[metadat.week5_DABTRAM.clust3$cell_id, ]

df2 <- data.frame(x = scores.df.week5.cl0$`Interferon/MHC-II (I)`)
df2$group <- 'clust0'
df3 <- data.frame(x = scores.df.week5.cl3$`Interferon/MHC-II (I)`)
df3$group <- 'clust3'
df4 <- rbind(df2, df3)

ggplot(df4, aes(x = x, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('clust0' = '#F8766D', 'clust3' = '#00B0F6')) +
  xlab('MP17: Interferon/MHC-II (I) expression') +
  theme_Publication()

df2 <- data.frame(x = scores.df.week5.cl0$Stress)
df2$group <- 'clust0'
df3 <- data.frame(x = scores.df.week5.cl3$Stress)
df3$group <- 'clust3'
df4 <- rbind(df2, df3)

ggplot(df4, aes(x = x, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('clust0' = '#F8766D', 'clust3' = '#00B0F6')) +
  xlab('MP5: Stress expression') +
  theme_Publication()


scores.df.d0.summarized <- scores.df.d0 %>% 
  summarise_all(mean)
scores.df.d0.summarized <- t(scores.df.d0.summarized)
colnames(scores.df.d0.summarized) <- 'd0'

scores.df.high.fp.d10.summarized <- scores.df.high.fp.d10 %>% 
  summarise_all(mean)
scores.df.high.fp.d10.summarized <- t(scores.df.high.fp.d10.summarized)
colnames(scores.df.high.fp.d10.summarized) <- 'high.fp.d10'

scores.df.low.fp.d10.summarized <- scores.df.low.fp.d10 %>% 
  summarise_all(mean)
scores.df.low.fp.d10.summarized <- t(scores.df.low.fp.d10.summarized)
colnames(scores.df.low.fp.d10.summarized) <- 'low.fp.d10'

scores.df.week5.summarized <- scores.df.week5 %>% 
  summarise_all(mean)
scores.df.week5.summarized <- t(scores.df.week5.summarized)
colnames(scores.df.week5.summarized) <- 'week5'

scores.df.week5.cl0.summarized <- scores.df.week5.cl0 %>% 
  summarise_all(mean)
scores.df.week5.cl0.summarized <- t(scores.df.week5.cl0.summarized)
colnames(scores.df.week5.cl0.summarized) <- 'week5.cl0'

scores.df.week5.cl3.summarized <- scores.df.week5.cl3 %>% 
  summarise_all(mean)
scores.df.week5.cl3.summarized <- t(scores.df.week5.cl3.summarized)
colnames(scores.df.week5.cl3.summarized) <- 'week5.cl3'

to_plot <- cbind(scores.df.d0.summarized,
                 scores.df.low.fp.d10.summarized, scores.df.high.fp.d10.summarized, 
                 scores.df.week5.cl0.summarized, scores.df.week5.cl3.summarized)

to_plot <- cbind(t(scores.df.low.fp.d10), t(scores.df.high.fp.d10), 
                 t(scores.df.week5.cl0), t(scores.df.week5.cl3))

to_plot <- cbind(t(scores.df.week5.cl0), t(scores.df.week5.cl3))

to_plot <- cbind(scores.df.week5.cl0.summarized, scores.df.week5.cl3.summarized)


mps.to_plot <- c('EMT-I', 'Stress', 'Unfolded protein response', 'Interferon/MHC-II (I)',  'Hypoxia')
mps.to_plot <- c("Cell Cycle - G2/M" , "Cell Cycle - G1/S",  "Stress",
                 "Proteasomal degradation", "Unfolded protein response", 
                 "EMT-I", "Interferon/MHC-II (I)", "Interferon/MHC-II (II)",
                 "Epithelial Senescence",  "MYC", "Respiration", "Secreted I", "Secreted II", "Skin-pigmentation",
                 'ISG.RS')

pheatmap::pheatmap(to_plot[mps.to_plot, ], show_rownames=TRUE, show_colnames=TRUE, scale = 'row',
                   cluster_rows=FALSE, cluster_cols=FALSE, fontsize=8, border_color='black',
                   breaks = seq(-1.2, 1.2, length.out = 100),
                   color = colorRampPalette(c("#2192FF", "#FFFBCA", "#EB5A3C"))(100))

# ==============================================================================
# DE in week5
# ==============================================================================

de.res <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(de.res) <- c('MP', 'diff', 'pval')

for(mp in colnames(scores.df)){
  
  res <- t.test(scores.df.week5.cl0[, mp], 
                scores.df.week5.cl3[, mp])
  diff <- mean(scores.df.week5.cl0[, mp]) - mean(scores.df.week5.cl3[, mp])
  de.res <- rbind(de.res, c(mp, diff, res$p.value))
}

colnames(de.res) <- c('MP', 'diff', 'pval')
de.res$diff <- as.numeric(de.res$diff)
de.res$pval <- as.numeric(de.res$pval)
de.res$pval.adj <- p.adjust(de.res$pval, method = 'BH')

de.res$neg.log10.pval.adj <- -log10(de.res$pval.adj)
de.res$neg.log10.pval.adj <- ifelse(de.res$neg.log10.pval.adj > 100, 100, de.res$neg.log10.pval.adj)

ggplot(de.res, aes(x=diff, y=neg.log10.pval.adj)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x='Difference in mean module score', y='-log10(p-value)')

scores.df.week5.cl0.new <- scores.df.week5.cl0
scores.df.week5.cl3.new <- scores.df.week5.cl3

scores.df.week5.cl0.new$group <- 'week5.cl0'
scores.df.week5.cl3.new$group <- 'week5.cl3'

scores.df.week5.new <- rbind(scores.df.week5.cl0.new, scores.df.week5.cl3.new)
scores.df.week5.new <- melt(scores.df.week5.new, id.vars = 'group')

signatures <- c('Cell Cycle - G2/M', "Cell Cycle - G1/S", "Cell Cycle HMG-rich",      
                'Chromatin',  'Stress', 'Hypoxia', 'Stress (in vitro)', 'Proteasomal degradation',
                'Unfolded protein response', 'Protein maturation', 'Translation initiation', 'EMT-I',                   
                'EMT-II', 'EMT-III', 'EMT-IV', 'Interferon/MHC-II (I)', 'Interferon/MHC-II (II)',
                'Epithelial Senescence', 'MYC', 'Respiration', 'Secreted I', 'Secreted II',
                'Cilia', 'Skin-pigmentation', 'Unassigned', 'ISG.RS', 'ISG.Mem' )
scores.df.week5.new <- scores.df.week5.new[scores.df.week5.new$variable %in% signatures, ]
ggplot(scores.df.week5.new, aes(x=group, y=value)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  facet_wrap(~variable, scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x='Group', y='Module score')
