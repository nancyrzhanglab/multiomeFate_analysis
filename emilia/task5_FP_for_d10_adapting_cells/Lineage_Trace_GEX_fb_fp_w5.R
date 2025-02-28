rm(list=ls())
library(Seurat)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE
treatment <- 'COCL2'

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

# fate bias from d0 to week5 adapting
df.bias <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_', treatment, '.csv'))

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
# ==============================================================================
# Check if high bias cells are in day10 or week5
# ==============================================================================

# get high bias cells
high_bias_cells <- df.bias %>% filter(bias > 0.5)
high_bias_cells <- merge(high_bias_cells, metadat, by='cell_id')


metadat.high_bias_cells <- metadat %>% filter(assigned_lineage %in% high_bias_cells$assigned_lineage)
check <- as.data.frame(table(metadat.high_bias_cells$assigned_lineage, metadat.high_bias_cells$dataset))
check <- spread(check, Var2, Freq)
check <- check[, c('Var1', 'day0', 'day10_DABTRAM', 'week5_DABTRAM')]

# high fp d10
fp.dabtram.d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(fp.dabtram.d10_w5) <- c('fatepotential_DABTRAM_d10_w5')
fp.dabtram.d10_w5$cell_id <- rownames(fp.dabtram.d10_w5)
high_fp_d10 <- fp.dabtram.d10_w5 %>% filter(fatepotential_DABTRAM_d10_w5 > 0)

high_fp_d10 <- merge(high_fp_d10, metadat, by=c('fatepotential_DABTRAM_d10_w5', 'cell_id'))
high_fp_d10.sub <- high_fp_d10 %>% filter(assigned_lineage %in% high_bias_cells$assigned_lineage)

# week5 
metadat.week5_DABTRAM <- metadat %>% filter(dataset == 'week5_DABTRAM')

metadat.DABTRAM <- read.csv(paste0(data_dir, 'metadat.all_data_dabtram_recluster.csv'), row.names = 1)
metadat.week5_DABTRAM <- metadat.DABTRAM[metadat.DABTRAM$dataset == 'week5_DABTRAM', ]
metadat.week5_DABTRAM.clust0 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 0)
metadat.week5_DABTRAM.clust3 <- metadat.week5_DABTRAM %>% filter(seurat_clusters == 3)

metadat.week5_DABTRAM.sub <- metadat.week5_DABTRAM %>% filter(assigned_lineage %in% high_bias_cells$assigned_lineage)

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
scores.df.scaled <- scale(scores.df)
scores.df.scaled <- as.data.frame(scores.df.scaled)

scores.df.high.bias.cells <- scores.df.scaled[high_bias_cells$cell_id, ]
scores.df.high.fp.d10 <- scores.df.scaled[high_fp_d10.sub$cell_id, ]
scores.df.week5 <- scores.df.scaled[metadat.week5_DABTRAM.sub$cell_id, ]

scores.df.week5.cl0 <- scores.df.scaled[metadat.week5_DABTRAM.clust0$cell_id, ]
scores.df.week5.cl3 <- scores.df.scaled[metadat.week5_DABTRAM.clust3$cell_id, ]

scores.df.high.bias.cells.summarized <- scores.df.high.bias.cells %>% 
  summarise_all(mean)
scores.df.high.bias.cells.summarized <- t(scores.df.high.bias.cells.summarized)
colnames(scores.df.high.bias.cells.summarized) <- 'high.bias.cells.d0'

scores.df.high.fp.d10.summarized <- scores.df.high.fp.d10 %>% 
  summarise_all(mean)
scores.df.high.fp.d10.summarized <- t(scores.df.high.fp.d10.summarized)
colnames(scores.df.high.fp.d10.summarized) <- 'high.fp.d10'

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

to_plot <- cbind(scores.df.high.bias.cells.summarized, scores.df.high.fp.d10.summarized, scores.df.week5.summarized)

to_plot <- cbind(scores.df.high.bias.cells.summarized, scores.df.high.fp.d10.summarized,
                 scores.df.week5.cl0.summarized, scores.df.week5.cl3.summarized)

to_plot <- rbind(scores.df.week5.cl0, scores.df.week5.cl3)

to_plot <- rbind(scores.df.high.bias.cells, scores.df.high.fp.d10, scores.df.week5)
to_plot <- t(to_plot)

# to_plot <- scale(to_plot)

pheatmap::pheatmap(to_plot, show_rownames=FALSE, show_colnames=TRUE, 
                  cluster_rows=FALSE, cluster_cols=FALSE, fontsize=8, border_color='black',
                  gaps_row = c(4173, 6824),
                  breaks=seq(-2, 2, length.out = 101), color = colorRampPalette(c("#2192FF", "#FFFBCA", "#EB5A3C"))(100))

row.to.keep <- c("Cell Cycle - G2/M", "Cell Cycle - G1/S", "Cell Cycle HMG-rich", "Chromatin",
                 "Stress", "Hypoxia", "Stress (in vitro)", "Proteasomal degradation", 
                 "Unfolded protein response", "Protein maturation", "Translation initiation",
                 "EMT-I", "EMT-II", "EMT-III", "EMT-IV","Interferon/MHC-II (I)", "Interferon/MHC-II (II)",
                 "Epithelial Senescence", "MYC", "Respiration", "Secreted I", "Secreted II", "Cilia",                  
                 "Skin-pigmentation", "Metal-response", "ISG.RS" )

to_plot <- to_plot[row.to.keep, ]

pheatmap::pheatmap(to_plot, show_rownames=TRUE, show_colnames=TRUE,
                   cluster_rows=FALSE, cluster_cols=FALSE, fontsize=8, border_color='black',
                   gaps_col = c(44, 54, 84),
                   breaks=seq(-2, 2, length.out = 101), color = colorRampPalette(c("#2192FF", "#FFFBCA", "#EB5A3C"))(100))
