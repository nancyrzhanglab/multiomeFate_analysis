rm(list = ls())

library(tidyverse)
library(ggplot2)
library(data.table)
library(GGally)

remove_unassigned_cells <- TRUE
dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")
theme_Publication<- function(base_size=12, base_family="sans") {
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
            axis.text = element_text(),
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
# =============================================================================
# reading data
# =============================================================================
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig5/'

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_CIS.RData'))


all_data[['FastTopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[["ft.DABTRAM.umap"]] <- all_data_ft_DABTRAM_umap
all_data[['FastTopic_COCL2']] <- all_data_fasttopic_COCL2
all_data[["ft.COCL2.umap"]] <- all_data_ft_COCL2_umap
all_data[['FastTopic_CIS']] <- all_data_fasttopic_CIS
all_data[["ft.CIS.umap"]] <- all_data_ft_CIS_umap


# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

rna.gavish.mp <- read.csv(paste0(results_dir, 'GavishMP_UCell_scores.csv'), row.names = 1)
rna.gavish.mp$cell_id <- rownames(rna.gavish.mp)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.DABTRAM <- metadat[metadat$dataset %in% c('day0', 'day10_DABTRAM', 'week5_DABTRAM'),]
metadat.COCL2 <- metadat[metadat$dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2'),]
metadat.CIS <- metadat[metadat$dataset %in% c('day0', 'day10_CIS', 'week5_CIS'),]

umap.DABTRAM <- as.data.frame(all_data@reductions[["ft.DABTRAM.umap"]]@cell.embeddings) %>% drop_na()
umap.COCL2 <- as.data.frame(all_data@reductions[["ft.COCL2.umap"]]@cell.embeddings) %>% drop_na()
umap.CIS <- as.data.frame(all_data@reductions[["ft.CIS.umap"]]@cell.embeddings) %>% drop_na()

umap.DABTRAM$cell_id <- rownames(umap.DABTRAM)
umap.COCL2$cell_id <- rownames(umap.COCL2)
umap.CIS$cell_id <- rownames(umap.CIS)

mps <- c('Cell.Cycle...G2.M_UCell', 'Cell.Cycle...G1.S_UCell', 'Cell.Cycle.HMG.rich_UCell', 'Chromatin_UCell',
         'Stress_UCell', 'Hypoxia_UCell', 'Stress..in.vitro._UCell', 'Proteasomal.degradation_UCell', 'Unfolded.protein.response_UCell', 
         'Protein.maturation_UCell', 'Translation.initiation_UCell', 'EMT.I_UCell', 'EMT.II_UCell', 
         'EMT.III_UCell', 'EMT.IV_UCell', 'Interferon.MHC.II..I._UCell', 'Interferon.MHC.II..II._UCell',
         'Epithelial.Senescence_UCell', 'MYC_UCell', 'Respiration_UCell', 'Secreted.I_UCell', 'Secreted.II_UCell',
         'Skin.pigmentation_UCell')

# =============================================================================
# Plot violins and umap
# =============================================================================

to_plot.COCL2 <- merge(umap.COCL2, rna.gavish.mp, by = 'cell_id')
to_plot.COCL2 <- merge(to_plot.COCL2, metadat.COCL2[, c('cell_id', 'dataset')], by = 'cell_id')

ggplot(to_plot.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = dataset)) +
  geom_point() +
  theme_minimal() +
  ggtitle("COCL2")

to_plot.COCL2 <- to_plot.COCL2[order(to_plot.COCL2$`EMT.I_UCell`), ]
ggplot(to_plot.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2, color = `EMT.I_UCell`)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.3) +
  theme_minimal() +
  ggtitle("COCL2")

to_plot.COCL2_m <- melt(to_plot.COCL2, id.vars = c('cell_id', 'dataset', 'ftCOCL2umap_1', 'ftCOCL2umap_2'))
to_plot.COCL2_m$dataset <- factor(to_plot.COCL2_m$dataset, levels = c('day0', 'day10_COCL2', 'week5_COCL2'))
ggplot(to_plot.COCL2_m, aes(x = dataset, y = value)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.1) +
  facet_wrap(~variable, scales = 'free_y')


to_plot.CIS <- merge(umap.CIS, rna.gavish.mp, by = 'cell_id')
to_plot.CIS <- merge(to_plot.CIS, metadat.CIS[, c('cell_id', 'dataset')], by = 'cell_id')
to_plot.CIS_m <- melt(to_plot.CIS, id.vars = c('cell_id', 'dataset', 'ftCISumap_1', 'ftCISumap_2'))
to_plot.CIS_m$dataset <- factor(to_plot.CIS_m$dataset, levels = c('day0', 'day10_CIS', 'week5_CIS'))
ggplot(to_plot.CIS_m, aes(x = dataset, y = value)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.1) +
  facet_wrap(~variable, scales = 'free_y')

# =============================================================================
# Plot heatmap
# =============================================================================

mps.to.plot <- c('Cell.Cycle...G2.M_UCell', 'Cell.Cycle...G1.S_UCell', 'Stress_UCell', 'Hypoxia_UCell',
                  'EMT.I_UCell', 'EMT.II_UCell', 'EMT.III_UCell', 'EMT.IV_UCell',
                 'Interferon.MHC.II..I._UCell', 'Interferon.MHC.II..II._UCell')

mat.DABTRAM <- merge(rna.gavish.mp, metadat.DABTRAM[, c('cell_id', 'dataset')], by = 'cell_id')
mat.DABTRAM <- mat.DABTRAM[, -c(1)]
mat.DABTRAM_m <- melt(mat.DABTRAM, id.vars = c('dataset'))
mat.DABTRAM_m <- mat.DABTRAM_m[mat.DABTRAM_m$variable %in% mps.to.plot,]

mat.COCL2 <- merge(rna.gavish.mp, metadat.COCL2[, c('cell_id', 'dataset')], by = 'cell_id')
mat.COCL2 <- mat.COCL2[, -c(1)]
mat.COCL2_m <- melt(mat.COCL2, id.vars = c('dataset'))
mat.COCL2_m <- mat.COCL2_m[mat.COCL2_m$variable %in% mps.to.plot,]

mat.CIS <- merge(rna.gavish.mp, metadat.CIS[, c('cell_id', 'dataset')], by = 'cell_id')
mat.CIS <- mat.CIS[, -c(1)]
mat.CIS_m <- melt(mat.CIS, id.vars = c('dataset'))
mat.CIS_m <- mat.CIS_m[mat.CIS_m$variable %in% mps.to.plot,]

ggplot(mat.DABTRAM_m, aes(x = dataset, y = value)) +
  geom_violin(aes(fill = dataset), scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = 'black', lwd = 0.5) +
  scale_fill_manual(values = dataset_colors) +
  facet_wrap(~variable, scales = 'free_y', ncol = 5) +
  ggtitle("DABTRAM") +
  theme_Publication(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(figure_dir, 'DABTRAM_violin.pdf'), width = 6, height = 4.5)

ggplot(mat.COCL2_m, aes(x = dataset, y = value)) +
  geom_violin(aes(fill = dataset), scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = 'black', lwd = 0.5) +
  scale_fill_manual(values = dataset_colors) +
  facet_wrap(~variable, scales = 'free_y', ncol = 5) +
  ggtitle("COCL2") +
  theme_Publication(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(figure_dir, 'COCL2_violin.pdf'), width = 6, height = 4.5)

ggplot(mat.CIS_m, aes(x = dataset, y = value)) +
  geom_violin(aes(fill = dataset), scale = 'width') +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = 'black', lwd = 0.5) +
  scale_fill_manual(values = dataset_colors) +
  facet_wrap(~variable, scales = 'free_y', ncol = 5) +
  ggtitle("CIS") +
  theme_Publication(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(figure_dir, 'CIS_violin.pdf'), width = 6, height = 4.5)


# =============================================================================
# Plot correlations ggpairs
# =============================================================================

mat.DABTRAM.day10 <- mat.DABTRAM[mat.DABTRAM$dataset == 'day10_DABTRAM', ]
mat.DABTRAM.week5 <- mat.DABTRAM[mat.DABTRAM$dataset == 'week5_DABTRAM', ]

ggpairs(mat.DABTRAM.day10[, c('EMT.II_UCell', 'Epithelial.Senescence_UCell', 'Interferon.MHC.II..II._UCell', 'Stress_UCell')])
ggpairs(mat.DABTRAM.week5[, c('EMT.II_UCell', 'Epithelial.Senescence_UCell', 'Interferon.MHC.II..II._UCell', 'Stress_UCell')])


mat.COCL2.day10 <- mat.COCL2[mat.COCL2$dataset == 'day10_COCL2', ]
mat.COCL2.week5 <- mat.COCL2[mat.COCL2$dataset == 'week5_COCL2', ]

ggpairs(mat.COCL2.week5[, c('EMT.IV_UCell', 'Hypoxia_UCell', 'Interferon.MHC.II..II._UCell', 'Secreted.II_UCell')])
ggpairs(mat.COCL2.day10[, c('EMT.IV_UCell', 'Hypoxia_UCell', 'Interferon.MHC.II..II._UCell', 'Secreted.II_UCell')])


mat.CIS.day10 <- mat.CIS[mat.CIS$dataset == 'day10_CIS', ]
mat.CIS.week5 <- mat.CIS[mat.CIS$dataset == 'week5_CIS', ]
ggpairs(mat.CIS.day10[, c('Cell.Cycle...G1.S_UCell', 'EMT.I_UCell', 'MYC_UCell', 'Proteasomal.degradation_UCell')])
ggpairs(mat.CIS.week5[, c('Cell.Cycle...G1.S_UCell', 'EMT.I_UCell', 'MYC_UCell', 'Proteasomal.degradation_UCell')])



# DABTRAM day10
p1 <- ggpairs(mat.DABTRAM.day10[, c('EMT.II_UCell', 'Epithelial.Senescence_UCell', 'Interferon.MHC.II..II._UCell', 'Stress_UCell')], 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(mat.DABTRAM.day10[, c('EMT.II_UCell', 'Epithelial.Senescence_UCell', 'Interferon.MHC.II..II._UCell', 'Stress_UCell')], 
             label = TRUE, 
             limits = c(-0.2, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.4,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- 4 # length(mps.to_plot)

for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p1,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p1,k1,k2) +
      theme_bw()
    p1 <- putPlot(p1,plt,k2,k1)
    p1 <- putPlot(p1,plt2,k1,k2)
    idx <- idx+1
    
  }
}
print(p1)

# DABTRAM week5
p1 <- ggpairs(mat.DABTRAM.week5[, c('EMT.II_UCell', 'Epithelial.Senescence_UCell', 'Interferon.MHC.II..II._UCell', 'Stress_UCell')], 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(mat.DABTRAM.week5[, c('EMT.II_UCell', 'Epithelial.Senescence_UCell', 'Interferon.MHC.II..II._UCell', 'Stress_UCell')], 
             label = TRUE, 
             limits = c(-0.2, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.4,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- 4 # length(mps.to_plot)

for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p1,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p1,k1,k2) +
      theme_bw()
    p1 <- putPlot(p1,plt,k2,k1)
    p1 <- putPlot(p1,plt2,k1,k2)
    idx <- idx+1
    
  }
}
print(p1)


# COCL2 day10
# Matrix of plots
p1 <- ggpairs(mat.COCL2.day10[, c('EMT.IV_UCell', 'Hypoxia_UCell', 'Interferon.MHC.II..II._UCell', 'Secreted.II_UCell')], 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(mat.COCL2.day10[, c('EMT.IV_UCell', 'Hypoxia_UCell', 'Interferon.MHC.II..II._UCell', 'Secreted.II_UCell')], 
             label = TRUE, 
             limits = c(-0.2, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.4,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- 4 # length(mps.to_plot)

for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p1,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p1,k1,k2) +
      theme_bw()
    p1 <- putPlot(p1,plt,k2,k1)
    p1 <- putPlot(p1,plt2,k1,k2)
    idx <- idx+1
    
  }
}
print(p1)

# COCL2 week5
# Matrix of plots
p1 <- ggpairs(mat.COCL2.week5[, c('EMT.IV_UCell', 'Hypoxia_UCell', 'Interferon.MHC.II..II._UCell', 'Secreted.II_UCell')], 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(mat.COCL2.week5[, c('EMT.IV_UCell', 'Hypoxia_UCell', 'Interferon.MHC.II..II._UCell', 'Secreted.II_UCell')], 
             label = TRUE, 
             limits = c(-0.2, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.4,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- 4 # length(mps.to_plot)

for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p1,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p1,k1,k2) +
      theme_bw()
    p1 <- putPlot(p1,plt,k2,k1)
    p1 <- putPlot(p1,plt2,k1,k2)
    idx <- idx+1
    
  }
}
print(p1)


# CIS day10
# Matrix of plots
p1 <- ggpairs(mat.CIS.day10[, c('Cell.Cycle...G1.S_UCell', 'EMT.I_UCell', 'MYC_UCell', 'Proteasomal.degradation_UCell')], 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(mat.CIS.day10[, c('Cell.Cycle...G1.S_UCell', 'EMT.I_UCell', 'MYC_UCell', 'Proteasomal.degradation_UCell')],
             label = TRUE, 
             limits = c(-0.2, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.4,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- 4 # length(mps.to_plot)

for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p1,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p1,k1,k2) +
      theme_bw()
    p1 <- putPlot(p1,plt,k2,k1)
    p1 <- putPlot(p1,plt2,k1,k2)
    idx <- idx+1
    
  }
}
print(p1)

# CIS week5
# Matrix of plots
p1 <- ggpairs(mat.CIS.week5[, c('Cell.Cycle...G1.S_UCell', 'EMT.I_UCell', 'MYC_UCell', 'Proteasomal.degradation_UCell')], 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(mat.CIS.week5[, c('Cell.Cycle...G1.S_UCell', 'EMT.I_UCell', 'MYC_UCell', 'Proteasomal.degradation_UCell')],
             label = TRUE, 
             limits = c(-0.2, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.4,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- 4 # length(mps.to_plot)

for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p1,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p1,k1,k2) +
      theme_bw()
    p1 <- putPlot(p1,plt,k2,k1)
    p1 <- putPlot(p1,plt2,k1,k2)
    idx <- idx+1
    
  }
}

print(p1)

