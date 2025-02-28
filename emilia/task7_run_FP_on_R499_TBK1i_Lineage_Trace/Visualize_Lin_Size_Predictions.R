rm(list = ls())

library(tidyverse)
library(ggplot2)


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
            # axis.ticks = element_blank(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

use.pal_SCP_Signature1 <- c("#6BAED6", "#08519C", "#74C476", "#FCBBA1", "#EF3B2C", "#67000D", "#A63603", "grey75", "grey50")
names(use.pal_SCP_Signature1) <- c("SCP-Sensitive", "SCP-Sensitive2", "SCP 11", "SCP-Resistant1", "SCP-Resistant2", "SCP 11/SCP-Resistant2", "SCP-Resistant2;IFNmem-hi", "Plastic", "Missing")

use.pal_SCP_Signature2 <- c("#6BAED6", "#74C476", "#FCBBA1", "#EF3B2C", "#67000D", "grey75", "grey50")
names(use.pal_SCP_Signature2) <- c("SCP-Sen1/2", "SCP 11", "SCP-Res1", "SCP-Res2", "SCP 11/SCP-Res2", "Plastic", "Missing")

# ==============================================================================
# read data
# ==============================================================================
data_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Lineage Trace/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task7_run_FP_on_R499_TBK1i_Lineage_Trace/'

cell_anno <- readRDS(paste0(data_dir, 'Lineage Barcode Annotations 061624.rds'))
metadat.LB_ALL_FILT_Annotated <- readRDS(paste0(data_dir, 'metadat.LB_ALL_FILT_Annotated.rds'))
metadat.LB_Anno <- readRDS(paste0(data_dir, 'Lineage Barcode Annotations 061624.rds'))

# lin.size.comp.df <- read.csv(paste0(out_dir, 'NT_d14_d19_w_d14_size_lineage_size_comparison.csv')) 
lin.size.comp.df <- read.csv(paste0(out_dir, 'JAKi_TBK1i_d14_d19_lineage_size_comparison.csv')) 



# ==============================================================================
# wrangle
# ==============================================================================

# metadat.LB_Anno <- metadat.LB_Anno[metadat.LB_Anno$Label == 'Naive NT', ]
metadat.LB_Anno <- metadat.LB_Anno[metadat.LB_Anno$Label == 'JAKi TBK1i NT', ]

lin.size.comp.df <- merge(lin.size.comp.df, 
                          metadat.LB_Anno[, c('lineage_barcode_assignment', "Annotation", "LB_Annotation", "LB_Annotation_V2", 
                                              "LB_Annotation_V3", "InVitro_NT_Annotation", "D14_NT_Annotation")], 
                          by = 'lineage_barcode_assignment', all.x = T)


# ==============================================================================
# plot
# ==============================================================================
lin.size.comp.df.na <- lin.size.comp.df[is.na(lin.size.comp.df$LB_Annotation_V3), ]
lin.size.comp.df.anno <- lin.size.comp.df[!is.na(lin.size.comp.df$LB_Annotation_V3), ]

ggplot(lin.size.comp.df, aes(x = lineage_future_count, y = lineage_imputed_count)) +
  geom_jitter(width = 2, color = '#DCDCDC', size = 2) +
  geom_jitter(data = lin.size.comp.df.anno, width = 0.1, shape = 21, stroke = 2, fill = 'white', aes(color = LB_Annotation_V3, size = lineage_current_count)) +
  scale_color_manual(values = use.pal_SCP_Signature2) +
  stat_cor(method = 'spearman') +
  labs(title = 'Lineage Size Predictions',
       x = 'Future Lineage Size (D19)',
       y = 'Imputed Lineage Size (FatePotential)',
       color = 'LB Annotation',
       size = 'Current Lineage Size (D14)') +
  theme_Publication()
ggsave(paste0(out_dir, 'Lineage_Size_Predictions_JAKi_TBK1i.png'), width = 6, height = 4, dpi = 300)

ggplot(lin.size.comp.df, aes(x = lineage_current_count, y = lineage_future_count)) +
  geom_jitter(width = 2, color = '#DCDCDC', size = 2) +
  geom_jitter(data = lin.size.comp.df.anno, width = 0.1, shape = 21, stroke = 2, fill = 'white', aes(color = LB_Annotation_V3, size = lineage_current_count)) +
  scale_color_manual(values = use.pal_SCP_Signature2) +
  stat_cor() +
  labs(title = 'Lineage Size Comparison',
       x = 'Current Lineage Size (D14)',
       y = 'Future Lineage Size (D19)',
       color = 'LB Annotation',
       size = 'Current Lineage Size (D14)') +
  theme_Publication()

ggsave(paste0(out_dir, 'Lineage_Size_Current_JAKi_TBK1i.png'), width = 6, height = 4, dpi = 300)






  
  
