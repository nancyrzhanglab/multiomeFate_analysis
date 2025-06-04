rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'


remove_unassigned_cells <- TRUE

theme_Publication<- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

theme_Clean<- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            axis.text = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour=NA,fill=NA),
            panel.background = element_blank(),
            strip.text = element_text(face="bold")
    ))
}
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

# =============================================================================
# Wrangle
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# DABTRAM
dabtram_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(dabtram_d10_w5) <- 'DABTRAM_d10_w5'

dabtram_d10_w5$cell_id <- rownames(dabtram_d10_w5)

dabtram_d10_w5 <- merge(dabtram_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
dabtram_d10_w5$progeny_size <- 10**dabtram_d10_w5$DABTRAM_d10_w5

gini.df.dabtram <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(gini.df.dabtram) <- c('lineage', 'gini_index', 'n_cells')
lineages <- unique(dabtram_d10_w5$assigned_lineage)
for(lin in lineages){
  df <- dabtram_d10_w5[dabtram_d10_w5$assigned_lineage == lin, ]
  

  progeny_sizes <- df$progeny_size
  
  # calculate gini index
  progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  increment = 1/length(progeny_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  
  gini.df.dabtram <- rbind(gini.df.dabtram, data.frame(lineage = lin, gini_index = gini, n_cells = nrow(df)))
  
}
gini.df.dabtram$dataset <-'DABTRAM_d10_w5'


# COCL2
cocl2_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_COCL2_d10_w5"]][["cell_imputed_score"]])
colnames(cocl2_d10_w5) <- 'COCL2_d10_w5'

cocl2_d10_w5$cell_id <- rownames(cocl2_d10_w5)

cocl2_d10_w5 <- merge(cocl2_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cocl2_d10_w5$progeny_size <- 10**cocl2_d10_w5$COCL2_d10_w5

gini.df.cocl2 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(gini.df.cocl2) <- c('lineage', 'gini_index', 'n_cells')
lineages <- unique(cocl2_d10_w5$assigned_lineage)
for(lin in lineages){
  df <- cocl2_d10_w5[cocl2_d10_w5$assigned_lineage == lin, ]
  
  
  progeny_sizes <- df$progeny_size
  
  # calculate gini index
  progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  increment = 1/length(progeny_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  
  gini.df.cocl2 <- rbind(gini.df.cocl2, data.frame(lineage = lin, gini_index = gini, n_cells = nrow(df)))
  
}
gini.df.cocl2$dataset <-'COCL2_d10_w5'

gini.df <- rbind(gini.df.dabtram, gini.df.cocl2)
gini.df <- gini.df[gini.df$n_cells > 10, ] # 5

gini.df$dataset <- factor(gini.df$dataset, levels = c('COCL2_d10_w5', 'DABTRAM_d10_w5'))
ggplot(gini.df, aes(y = dataset, x = gini_index)) +
  geom_violin(scale = 'width') +
  geom_jitter(aes(size = n_cells, fill = dataset), shape = 21, color = 'black') +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("#6DC49C", "#9D85BE")) +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('') +
  xlab('Gini index') +
  xlim(0, 1) +
  theme_Publication() +
  theme(legend.position = 'none')

# ggsave(paste0(figure_dir, 'fate_potential_gini_index.pdf'), width = 5.5, height = 2.3)
write.csv(gini.df, paste0(output_dir, 'fate_potential_gini_index.csv'), row.names = FALSE)


# CIS
cis_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_CIS_d10_w5"]][["cell_imputed_score"]])
colnames(cis_d10_w5) <- 'CIS_d10_w5'

cis_d10_w5$cell_id <- rownames(cis_d10_w5)

cis_d10_w5 <- merge(cis_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cis_d10_w5$progeny_size <- 10**cis_d10_w5$CIS_d10_w5

gini.df.cis <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(gini.df.cis) <- c('lineage', 'gini_index', 'n_cells')
lineages <- unique(cis_d10_w5$assigned_lineage)
for(lin in lineages){
  df <- cis_d10_w5[cis_d10_w5$assigned_lineage == lin, ]
  
  
  progeny_sizes <- df$progeny_size
  
  # calculate gini index
  progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  increment = 1/length(progeny_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  
  gini.df.cis <- rbind(gini.df.cis, data.frame(lineage = lin, gini_index = gini, n_cells = nrow(df)))
  
}
gini.df.cis$dataset <-'CIS_d10_w5'
write.csv(gini.df.cis, paste0(output_dir, 'fate_potential_gini_index_CIS.csv'), row.names = FALSE)


# DABTRAM
progeny_sizes <- dabtram_d10_w5$progeny_size

# calculate gini index
progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
increment = 1/length(progeny_sizes)
B=sum(increment*sizeprop)-increment/2
gini = 1-2*B
gini

pdf(paste0(figure_dir, 'fate_potential_gini_index_DABTRAM.pdf'), width = 5, height = 5)
plot(NA,
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "",
       ylab = "",
       asp = TRUE,
       xaxt = "n",
       yaxt = "n",
       bty = "n")
axis(1); axis(2)
lines(c(0,1), c(0,1), lty = 2, col = "coral", lwd = 3)
lines(x = seq(increment, 1, length.out = length(progeny_sizes)), 
      y = sizeprop,
      col = '#9D85BE', lwd = 8)
dev.off()

# COCL2

# calculate gini index
progeny_sizes <- cocl2_d10_w5$progeny_size
progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
increment = 1/length(progeny_sizes)
B=sum(increment*sizeprop)-increment/2
gini = 1-2*B
gini


pdf(paste0(figure_dir, 'fate_potential_gini_index_COCL2.pdf'), width = 5, height = 5)
plot(NA,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "",
     ylab = "",
     asp = TRUE,
     xaxt = "n",
     yaxt = "n",
     bty = "n")
axis(1); axis(2)
lines(c(0,1), c(0,1), lty = 2, col = "coral", lwd = 3)
lines(x = seq(increment, 1, length.out = length(progeny_sizes)),
      y = sizeprop,
      col = '#6DC49C', lwd = 8)
dev.off()


# CIS
cis_d10_w5 <- as.data.frame(all_data_fatepotential[["fatepotential_CIS_d10_w5"]][["cell_imputed_score"]])
colnames(cis_d10_w5) <- 'CIS_d10_w5'
cis_d10_w5$cell_id <- rownames(cis_d10_w5)
cis_d10_w5 <- merge(cis_d10_w5, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cis_d10_w5$progeny_size <- 10**cis_d10_w5$CIS_d10_w5


progeny_sizes <- cis_d10_w5$progeny_size

# calculate gini index
progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
increment = 1/length(progeny_sizes)
B=sum(increment*sizeprop)-increment/2
gini = 1-2*B
gini
pdf(paste0(figure_dir, 'Supp_fate_potential_gini_index_CIS.pdf'), width = 5, height = 5)
plot(NA,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "",
     ylab = "",
     asp = TRUE,
     xaxt = "n",
     yaxt = "n",
     bty = "n")
axis(1); axis(2)
lines(c(0,1), c(0,1), lty = 2, col = "coral", lwd = 3)
lines(x = seq(increment, 1, length.out = length(progeny_sizes)), 
      y = sizeprop,
      col = '#C96D29', lwd = 8)
dev.off()


# ==============================================================================
# Calculate the fate potential gini for all datasets
# ==============================================================================

# DABTRAM
dabtram_d0_d10 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d0_d10"]][["cell_imputed_score"]])
colnames(dabtram_d0_d10) <- 'DABTRAM_d0_d10'

dabtram_d0_d10$cell_id <- rownames(dabtram_d0_d10)

dabtram_d0_d10 <- merge(dabtram_d0_d10, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
dabtram_d0_d10$progeny_size <- 10**dabtram_d0_d10$DABTRAM_d0_d10

gini.df.dabtram.d0 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(gini.df.dabtram.d0) <- c('lineage', 'gini_index', 'n_cells')
lineages <- unique(dabtram_d0_d10$assigned_lineage)
for(lin in lineages){
  df <- dabtram_d0_d10[dabtram_d0_d10$assigned_lineage == lin, ]
  
  
  progeny_sizes <- df$progeny_size
  
  # calculate gini index
  progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  increment = 1/length(progeny_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  
  gini.df.dabtram.d0 <- rbind(gini.df.dabtram.d0, data.frame(lineage = lin, gini_index = gini, n_cells = nrow(df)))
  
}
gini.df.dabtram.d0$dataset <-'DABTRAM_d0_d10'


# COCL2
cocl2_d0_d10 <- as.data.frame(all_data_fatepotential[["fatepotential_COCL2_d0_d10"]][["cell_imputed_score"]])
colnames(cocl2_d0_d10) <- 'COCL2_d0_d10'

cocl2_d0_d10$cell_id <- rownames(cocl2_d0_d10)

cocl2_d0_d10 <- merge(cocl2_d0_d10, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cocl2_d0_d10$progeny_size <- 10**cocl2_d0_d10$COCL2_d0_d10

gini.df.cocl2.d0 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(gini.df.cocl2.d0) <- c('lineage', 'gini_index', 'n_cells')
lineages <- unique(cocl2_d0_d10$assigned_lineage)
for(lin in lineages){
  df <- cocl2_d0_d10[cocl2_d0_d10$assigned_lineage == lin, ]
  
  
  progeny_sizes <- df$progeny_size
  
  # calculate gini index
  progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  increment = 1/length(progeny_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  
  gini.df.cocl2.d0 <- rbind(gini.df.cocl2.d0, data.frame(lineage = lin, gini_index = gini, n_cells = nrow(df)))
  
}
gini.df.cocl2.d0$dataset <-'COCL2_d0_d10'


# CIS
cis_d0_d10 <- as.data.frame(all_data_fatepotential[["fatepotential_CIS_d0_d10"]][["cell_imputed_score"]])
colnames(cis_d0_d10) <- 'CIS_d0_d10'

cis_d0_d10$cell_id <- rownames(cis_d0_d10)

cis_d0_d10 <- merge(cis_d0_d10, metadat[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
cis_d0_d10$progeny_size <- 10**cis_d0_d10$CIS_d0_d10

gini.df.cis.d0 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(gini.df.cis.d0) <- c('lineage', 'gini_index', 'n_cells')
lineages <- unique(cis_d0_d10$assigned_lineage)
for(lin in lineages){
  df <- cis_d0_d10[cis_d0_d10$assigned_lineage == lin, ]
  
  
  progeny_sizes <- df$progeny_size
  
  # calculate gini index
  progeny_sizes=sort(progeny_sizes, decreasing=FALSE)
  sizeprop = cumsum(progeny_sizes/sum(progeny_sizes))
  increment = 1/length(progeny_sizes)
  B=sum(increment*sizeprop)-increment/2
  gini = 1-2*B
  
  gini.df.cis.d0 <- rbind(gini.df.cis.d0, data.frame(lineage = lin, gini_index = gini, n_cells = nrow(df)))
  
}
gini.df.cis.d0$dataset <-'CIS_d0_d10'


gini.df.all <- rbind(gini.df.dabtram, gini.df.cocl2)
gini.df.all <- rbind(gini.df.all, gini.df.cis)
gini.df.all <- rbind(gini.df.all, gini.df.dabtram.d0)
gini.df.all <- rbind(gini.df.all, gini.df.cocl2.d0)
gini.df.all <- rbind(gini.df.all, gini.df.cis.d0)

gini.df.all <- gini.df.all[gini.df.all$n_cells > 5, ]
ggplot(gini.df.all, aes(x = dataset, y = gini_index)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.1, outlier.shape = NA)

write.csv(gini.df.all, paste0(output_dir, 'fate_potential_gini_index_all.csv'), row.names = FALSE)
