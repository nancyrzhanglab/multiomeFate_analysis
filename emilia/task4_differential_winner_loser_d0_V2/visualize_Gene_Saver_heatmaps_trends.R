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

############################################################################################################
# Read data general
############################################################################################################

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver

genes_of_interest <- c('MITF', 'FN1', 'TRPM1', 'MYO1D', 'SLCO5A1', 'ACTB')

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
all_data_d0 <- subset(all_data, dataset == 'day0')
all_data_d0_saver <- all_data_d0@assays[["Saver"]]@data
all_data_d0_saver <- all_data_d0_saver[, fp_dabtram_d0_d10$cell_id]
all_data_d0_saver <- as.matrix(all_data_d0_saver)


to_plot <- as.data.frame(all_data_d0_saver[genes_of_interest, fp_dabtram_d0_d10$cell_id])
to_plot$feature <- rownames(to_plot)

to_plot_m <- melt(to_plot)
colnames(to_plot_m) <- c('Gene', 'cell_id', 'Saver')
to_plot_m <- merge(to_plot_m, fp_dabtram_d0_d10, by = 'cell_id')

ggplot(to_plot_m, aes(x = order, y = Saver, color = Gene)) +
  # geom_point(alpha = 0.02) +
  geom_smooth(method = 'loess', alpha = 0.05) +
  xlab('Fate potential order (low -> high)') +
  scale_color_manual(values = c("MITF" = "#D31816", 'FN1' = '#B8E186')) +
  facet_wrap(~Gene) +
  # coord_cartesian(ylim = c(-2, 2)) +
  theme_Publication()
