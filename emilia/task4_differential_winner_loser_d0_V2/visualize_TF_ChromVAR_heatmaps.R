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

# ==============================================================================
# read data
# ==============================================================================
load('~/Downloads/Writeup10a_data_fatepotential.RData')
load('~/Downloads/Writeup10a_data_chromVar_day0.RData')

tf_of_interest <- c('JUN(var.2)', 'FOSL1::JUND', 
                    'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 
                    'SNAI1', 'SNAI2', 'SNAI3', 'MITF', 'SOX10',
                    'ASCL1(var.2)', 'USF1', 'MAFK', 'RUNX2', 'NFE2L1', 'POU5F1')
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

diff_res_dabtram <- diff_res_dabtram[diff_res_dabtram$neg_log10_p_val > 8, ]
to_plot <- as.data.frame(all_data_chromVar_day0[tf_of_interest, fp_dabtram_d0_d10$cell_id])
to_plot <- as.data.frame(all_data_chromVar_day0[c(diff_res_dabtram$feature, 'MITF'), fp_dabtram_d0_d10$cell_id])
to_plot$feature <- rownames(to_plot)
to_plot$TF_label <- ifelse(grepl('JUN', to_plot$feature), 'AP1', 'Other')
to_plot$TF_label <- ifelse(grepl('FOS', to_plot$feature), 'AP1', to_plot$TF_label)
to_plot$TF_label <- ifelse(grepl('TEAD', to_plot$feature), 'TEAD', to_plot$TF_label)
to_plot$TF_label <- ifelse(grepl('SNAI', to_plot$feature), 'SNAI', to_plot$TF_label)
to_plot$TF_label <- ifelse(grepl('MITF', to_plot$feature), 'MITF', to_plot$TF_label)
to_plot$TF_label <- ifelse(grepl('SOX10', to_plot$feature), 'SOX10', to_plot$TF_label)
to_plot$TF_label <- factor(to_plot$TF_label, levels = c('AP1', 'TEAD', 'SNAI', 'MITF', 'SOX10', 'Other'))

to_plot <- to_plot[, c(fp_dabtram_d0_d10$cell_id[1:200],
                       fp_dabtram_d0_d10$cell_id[1001:1201],
                       fp_dabtram_d0_d10$cell_id[2001:2201],
                       fp_dabtram_d0_d10$cell_id[3001:3201])]
pheatmap(to_plot, color = colorRampPalette(c("blue", "white", "red"))(12),
         breaks = seq(-2, 2, length.out = 12), 
         cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F)

to_plot_m <- melt(to_plot)
colnames(to_plot_m) <- c('TF', 'TF_label', 'cell_id', 'chromVAR')
to_plot_m <- merge(to_plot_m, fp_dabtram_d0_d10, by = 'cell_id')

ggplot(to_plot_m, aes(x = order, y = chromVAR, color = TF)) +
  # geom_point(alpha = 0.02) +
  geom_smooth(method = 'loess', alpha = 0.05) +
  xlab('Fate potential order (low -> high)') +
  scale_color_manual(values = c("JUN(var.2)" = "#D31816", 
                                "TEAD1" = "#A269A2FF","TEAD2" = "#A269A2FF","TEAD3" = "#A269A2FF","TEAD4" = "#A269A2FF",
                                "SNAI1" =  '#33608CFF', "SNAI2" =  '#33608CFF', "SNAI3" =  '#33608CFF',
                                "MITF" =  '#276419',
                                "SOX10" = '#66C2A5',
                                "BACH1" = '#DE708DFF', 'BATF' = '#FEE08B', 'NFE2' = "#99540FFF", 'USF1' = '#B8E186')) +
  facet_wrap(~TF_label) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_Publication()

ggsave('~/Downloads/TF_day0_trends.pdf', width = 6, height = 5.5)


colors_epi_Clusters<-c("#D31816", "#99540FFF","#F2A36BFF","#A269A2FF","#33608CFF","#DE708DFF","#276419","#B8E186","#66C2A5","#FEE08B","#F1B6DA","#87A3B8")
