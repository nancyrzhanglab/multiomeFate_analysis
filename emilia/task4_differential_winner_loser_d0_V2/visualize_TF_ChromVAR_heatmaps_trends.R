library(tidyverse)
library(data.table)
library(ggplot2)


out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_differential_winner_loser_d0_V2/'

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

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Downloads/'
treatment <- 'DABTRAM'
# ==============================================================================
# read data
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))

load(paste0(result_dir, 'differential_winner_loser_chromVAR_lineage_specific_adapation_TF_', treatment, '_t_test_results.RData'))
# ==============================================================================
# wrangle data
# ==============================================================================

fp_dabtram_d0_d10 <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d0_d10"]][["cell_imputed_score"]])
colnames(fp_dabtram_d0_d10) <- c('fate_potential_d0_d10_dabtram')
fp_dabtram_d0_d10$cell_id <- rownames(fp_dabtram_d0_d10)

# order by fate potential
fp_dabtram_d0_d10 <- fp_dabtram_d0_d10[order(fp_dabtram_d0_d10$fate_potential_d0_d10_dabtram, decreasing = F),]
fp_dabtram_d0_d10$order <- seq(1: nrow(fp_dabtram_d0_d10))

# subset to day0 data
all_data_chromVar_day0 <- all_data_chromVar_day0@data
all_data_chromVar_day0 <- all_data_chromVar_day0[, fp_dabtram_d0_d10$cell_id]
all_data_chromVar_day0 <- as.matrix(all_data_chromVar_day0)

all_data_chromVar_day0 <- as.data.frame(t(all_data_chromVar_day0))

# Order chromaVAR data by FP
all_data_chromVar_day0 <- all_data_chromVar_day0[fp_dabtram_d0_d10$cell_id, ]

tf_of_interest <- c('JUN(var.2)', 'FOSL1::JUND', 
                    'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 
                    'SNAI1', 'SNAI2', 'SNAI3', 'MITF', 'SOX10',
                    'ASCL1(var.2)', 'USF1', 'MAFK', 'RUNX2', 'NFE2L1', 'POU5F1')
tf_of_interest <- t_test_results$feature

all_data_chromVar_day0_subset <- all_data_chromVar_day0[, tf_of_interest]
all_data_chromVar_day0_subset <- all_data_chromVar_day0 %>% drop_na()
# all_data_chromVar_day0_subset <- scale(all_data_chromVar_day0_subset)
pheatmap::pheatmap(t(all_data_chromVar_day0_subset), 
                  cluster_rows = T, cluster_cols = F, 
                  show_rownames = F, show_colnames = F,
                  breaks = seq(-2, 2, length.out = 101))
