rm(list = ls())

library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggdensity)
library(gridExtra)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task6_examine_subpop/'
ref_dir <- '~/Downloads/'

remove_unassigned_cells <- TRUE
# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))


all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[['FastTopic_DABTRAM']] <- all_data_fasttopic_DABTRAM

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

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
metadat.d10_dabtram <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.week5_dabtram <- metadat[metadat$dataset == 'week5_DABTRAM', ]

scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=c(gavish.mp.list, 
                                           list(isg_rs), 
                                           list(isg_mem)))

colnames(scores) <- c(names(gavish.mp.list), 'ISG.RS', 'ISG.Mem')
scores.df <- as.data.frame(scores) 


mps.to_plot <- c('EMT-I', 'Stress', 'Unfolded protein response', 'Interferon/MHC-II (I)',  'Hypoxia')


scores.df.d10_dabtram <- scores.df[rownames(metadat.d10_dabtram), mps.to_plot]
scores.df.week5_dabtram <- scores.df[rownames(metadat.week5_dabtram), mps.to_plot]

# Matrix of plots
p1 <- ggpairs(scores.df.d10_dabtram, upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p2 <- ggcorr(scores.df.d10_dabtram, label = TRUE, 
             limits = c(-0.1, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.45,
             label_round = 2)
p2

g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- length(mps.to_plot)

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


# Matrix of plots
p3 <- ggpairs(scores.df.week5_dabtram, 
              upper = list(continuous = wrap("cor", size = 5, color = "black")))
# Correlation matrix plot
p4 <- ggcorr(scores.df.week5_dabtram, label = TRUE, limits = c(-0.1, 1), 
             low = "#3EDBF0",
             mid = "#F8E7F6",
             high = "#9400FF",
             midpoint = 0.45,
             label_round = 2)


g4 <- ggplotGrob(p4)
colors <- g4$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p <- length(mps.to_plot)
for (k1 in 2: p) {
  for (k2 in 1:(k1 -1)) {
    plt <- getPlot(p3,k2,k1) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    plt2 <- getPlot(p3,k1,k2) +
      theme_bw()
    p3 <- putPlot(p3,plt,k2,k1)
    p3 <- putPlot(p3,plt2,k1,k2)
    idx <- idx+1
  }
}
print(p3)

ggsave(p1, filename = paste0(out_dir, 'CorrelationMatrix_D10_DABTRAM.png'), width = 6, height = 6)
ggsave(p3, filename = paste0(out_dir, 'CorrelationMatrix_W5_DABTRAM.png'), width = 6, height = 6)


