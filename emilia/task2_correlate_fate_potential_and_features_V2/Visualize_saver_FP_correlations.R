library(tidyverse)
library(ggplot2)
library(GGally)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

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
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
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
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/saver_cor_vec.RData')

dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']] 

dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']] 

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_saver_cor_vec) <- paste0(colnames(dabtram_d0_saver_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_saver_cor_vec) <- paste0(colnames(cocl2_d0_saver_cor_vec), '.COCL2_d0')
colnames(cis_d0_saver_cor_vec) <- paste0(colnames(cis_d0_saver_cor_vec), '.CIS_d0')

colnames(dabtram_d10_saver_cor_vec) <- paste0(colnames(dabtram_d10_saver_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_saver_cor_vec) <- paste0(colnames(cocl2_d10_saver_cor_vec), '.COCL2_d10')
colnames(cis_d10_saver_cor_vec) <- paste0(colnames(cis_d10_saver_cor_vec), '.CIS_d10')

# Day0
d0_cor <- merge(dabtram_d0_saver_cor_vec, cocl2_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

d0_cor <- merge(d0_cor, cis_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

# Day10
d10_cor <- merge(dabtram_d10_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

# ==============================================================================
# Plotting
# ==============================================================================

ggpairs(d0_cor[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0')], mapping = aes(alpha = 0.7, color = 'pink')) +
  theme_Publication()

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_RNA_saver_FP_day0.png', 
       width = 8, height = 8, dpi = 300)

ggpairs(d10_cor[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')], mapping = aes(alpha = 0.7, color = 'pink')) +
  theme_Publication()

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_RNA_saver_FP_day10.png', 
       width = 8, height = 8, dpi = 300)


