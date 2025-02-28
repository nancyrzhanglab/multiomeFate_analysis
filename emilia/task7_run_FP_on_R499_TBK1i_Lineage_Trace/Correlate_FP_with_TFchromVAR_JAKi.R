rm(list = ls())
library(tidyverse)
library(ggplot2)
library(matrixStats)
library(Seurat)
library(gplots)
library(ggpubr)

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
# Read data
# ==============================================================================


# For in vivo data
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Minn_lab/Final/'
base_name <- 'TBK1i_Multiome_StringentMT_InVivo_Run3_ATAC_Cancer_Harmony_eigs18'

chromVar_signature <- readRDS(paste0(data_dir, base_name, '/', base_name, '.scMultiome_combined_object_CHROMVAR_OTHER_SIGNATURES_022324.rds'))
chromVar_motif <- readRDS(paste0(data_dir, base_name, '/',  base_name, '.scMultiome_combined_object_CHROMVAR_MOTIFS.rds'))
chromVar_motif <- t(chromVar_motif)

metadat <- readRDS(paste0(data_dir, base_name, '/', base_name, '.scMultiome_combined_object_METADAT.rds'))

fp <- as.data.frame(readRDS('~/Downloads/final_fit_d14_d19_JAKi_TBK1i_w_d14_size_scaled.fatepotential.rds'))
# fp <- readRDS('~/Downloads/final_fit_d14_d19_scaled.rds')
# fp <- as.data.frame(fp[["cell_imputed_score"]])

# =============================================================================
# Prepare the response variables
# =============================================================================
metadat <- metadat[!is.na(metadat$lineage_barcode_assignment), ]
metadat <- metadat[metadat$lineage_barcode_assignment != 'Too many barcodes', ]
metadat.JAKi_TBK1i.d14 <- metadat[metadat$Condition == 'JAKi TBK1i D14 InVivo', ]

chromVAR.signature.JAKi_TBK1i.d14 <- chromVar_signature[rownames(metadat.JAKi_TBK1i.d14), ]

chromVAR.motif.JAKi_TBK1i.d14 <- chromVar_motif[rownames(metadat.JAKi_TBK1i.d14), ]
# chromVAR.JAKi_TBK1i.d14 <- t(chromVAR.JAKi_TBK1i.d14)

colnames(fp) <- 'FatePotential'
fp$cell_id <- rownames(fp)
fp <- fp[rownames(metadat.JAKi_TBK1i.d14), ]

chromVAR_cor_vec_motif <- sapply(1:ncol(chromVAR.motif.JAKi_TBK1i.d14), function(j){
  # skip if there are too many NAs
  if(sum(is.na(chromVAR.motif.JAKi_TBK1i.d14[,j])) > 0.5 * length(chromVAR.motif.JAKi_TBK1i.d14[,j])){
    return(c(NA, NA))
  }
  res <- stats::cor.test(fp[['FatePotential']], chromVAR.motif.JAKi_TBK1i.d14[,j],
                         alternative = "two.sided",
                         method = "pearson",
                         rm.na = T)
  c(res$estimate, res$p.value)
})

chromVAR_cor_vec_motif <- as.data.frame(t(chromVAR_cor_vec_motif))
colnames(chromVAR_cor_vec_motif) <- c("correlation", "p.value")
chromVAR_cor_vec_motif$Signature <- colnames(chromVAR.motif.JAKi_TBK1i.d14)
# rownames(chromVAR_cor_vec) <- colnames(chromVAR.NT.d14)

chromVAR_cor_vec_signature <- sapply(1:ncol(chromVAR.signature.JAKi_TBK1i.d14), function(j){
  # skip if there are too many NAs
  if(sum(is.na(chromVAR.signature.JAKi_TBK1i.d14[,j])) > 0.5 * length(chromVAR.signature.JAKi_TBK1i.d14[,j])){
    return(c(NA, NA))
  }
  res <- stats::cor.test(fp[['FatePotential']], chromVAR.signature.JAKi_TBK1i.d14[,j],
                         alternative = "two.sided",
                         method = "pearson",
                         rm.na = T)
  c(res$estimate, res$p.value)
})

chromVAR_cor_vec_signature <- as.data.frame(t(chromVAR_cor_vec_signature))
colnames(chromVAR_cor_vec_signature) <- c("correlation", "p.value")
chromVAR_cor_vec_signature$Signature <- colnames(chromVAR.signature.JAKi_TBK1i.d14)

chromVAR_cor_vec <- rbind(chromVAR_cor_vec_motif, chromVAR_cor_vec_signature)

# hist(chromVAR_cor_vec$correlation, breaks = 50)

chromVAR_cor_vec <- chromVAR_cor_vec[order(chromVAR_cor_vec$correlation, decreasing = T), ]
chromVAR_cor_vec <- chromVAR_cor_vec %>% drop_na()
chromVAR_cor_vec$order <- seq(1: nrow(chromVAR_cor_vec))

# =============================================================================
# plot
# =============================================================================
signatures <- c('AC0672:IRF/STAT:IRF', 'AC0691:RELA/NFKB:Rel', 'AC0051:BATF/IRF(bZIP)', 'AC0154:ZNF354A/ZNF:C2H2-ZF',
                'AC0652:TEAD:TEA', 'C0121:ESRRA/NR5A:Nuclear-receptor', 'AC0681:NFAT:Rel',
                'AC0049:JUND/JUN:bZIP')
signatures2 <- c('gained_other_promoters_ATAC', 'activated_enhancer_ATAC_inclusive', 'DAR_Resistant2_UP', 'DAR_Resistant3_UP', 'Chronic_IFNG_Resistance_ATAC_B16y',
                'Res499_SCP_11',  'Res499_SCP_Resistant1', 'deactivated_enhancer_ATAC_inclusive', 'lost_other_promoters', 'DAR_Sensitive1_UP', 'DAR_Sensitive2_UP', 'lost_other_promoters_ATAC')


top.signatures <- c('gained_other_promoters_ATAC', 'gained_TSS_promoters', 'DAR_Resistant2_UP', 'DAR_Resistant3_UP', 
                    'DAR_Resistant1_UP', 'DAR_Sensitive3_UP', 'Res499_SCP_Resistant1')
bottom.signatures <- c('deactivated_enhancer_ATAC_inclusive', 'lost_other_promoters', 'DAR_Sensitive1_UP', 'DAR_Sensitive2_UP', 'lost_other_promoters_ATAC')

# top.TFs <- c('AC0672:IRF/STAT:IRF','AC0673:IRF:IRF', 'AC0154:ZNF354A/ZNF:C2H2-ZF','AC0681:NFAT:Rel')
# bottom.TFs <- c('AC0691:RELA/NFKB:Rel', 'AC0051:BATF/IRF:bZIP', 'AC0652:TEAD:TEA', 'AC0668:STAT/ETV:STAT')


ggplot(chromVAR_cor_vec, aes(x = order, y = correlation)) +
  geom_point() +
  geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures, ], color = 'red', size = 3) +
  ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures, ], 
                           aes(label = Signature), size = 3, color = 'red', nudge_x = 10, nudge_y = 0.2) +
  geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures2, ], color = 'blue', size = 3) +
  ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures2, ], 
                           aes(label = Signature), size = 3, color = 'blue', nudge_x = -10, nudge_y = -0.2) +
  labs(title = 'Correlation between Fate Potential and TF activity (JAKi TBK1i)',
       x = 'Signature',
       y = 'Correlation') +
  xlim(-60, 1300) +
  ylim(-0.6, 0.4) +
  theme_bw()

chromVAR_cor_vec$p.adj <- p.adjust(chromVAR_cor_vec$p.value, method = 'BH')

thres.top <- chromVAR_cor_vec %>% 
  filter(p.adj < 0.05) %>% 
  filter(correlation > 0) %>% 
  arrange(correlation) %>%
  head(1) %>% 
  pull(correlation)

thres.bottom <- chromVAR_cor_vec %>% 
  filter(p.adj < 0.05) %>% 
  filter(correlation < 0) %>% 
  arrange(correlation) %>%
  tail(1) %>% 
  pull(correlation)

ggplot(chromVAR_cor_vec, aes(x = order, y = correlation)) +
  geom_point() +
  geom_hline(yintercept = thres.top, linetype = 'dashed', color = 'gray') +
  geom_hline(yintercept = thres.bottom, linetype = 'dashed', color = 'gray') +
  geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% top.signatures, ], color = 'red', size = 3) +
  ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% top.signatures, ],
                           aes(label = Signature), size = 3, color = 'red', nudge_x = 15, nudge_y = 0.2) +
  geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% bottom.signatures, ], color = 'blue', size = 3) +
  ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% bottom.signatures, ],
                           aes(label = Signature), size = 3, color = 'blue', nudge_x = -20, nudge_y = -0.15) +
  # geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% top.TFs, ], color = 'red', size = 3) +
  # ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% top.TFs, ],
  #                          aes(label = Signature), size = 3, color = 'red', nudge_x = 10, nudge_y = 0.2) +
  # geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% bottom.TFs, ], color = 'blue', size = 3) +
  # ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% bottom.TFs, ],
  #                          aes(label = Signature), size = 3, color = 'blue', nudge_x = -20, nudge_y = -0.15) +
  # geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures2, ], color = 'blue', size = 3) +
  # ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures2, ], 
  #                          aes(label = Signature), size = 3, color = 'blue', nudge_x = -10, nudge_y = -0.2) +
  labs(title = 'Correlation between Fate Potential and chromVAR (JAKi TBK1i)',
       x = 'Signature Order',
       y = 'Correlation') +
  xlim(-60, 1300) +
  ylim(-0.6, 0.5) +
  theme_Publication()



ggplot(chromVAR_cor_vec, aes(x = order, y = correlation)) +
  geom_point() +
  geom_hline(yintercept = thres.top, linetype = 'dashed', color = 'gray') +
  geom_hline(yintercept = thres.bottom, linetype = 'dashed', color = 'gray') +
  geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature == 'Res 499 Memory', ], color = 'red', size = 3) +
  ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature == 'Res 499 Memory', ],
                           aes(label = Signature), size = 3, color = 'red', nudge_x = 10, nudge_y = 0.2) +
  # geom_point(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures2, ], color = 'blue', size = 3) +
  # ggrepel::geom_text_repel(data = chromVAR_cor_vec[chromVAR_cor_vec$Signature %in% signatures2, ], 
  #                          aes(label = Signature), size = 3, color = 'blue', nudge_x = -10, nudge_y = -0.2) +
  labs(title = 'Correlation between Fate Potential and chromVAR (JAKi TBK1i)',
       x = 'Signature Order',
       y = 'Correlation') +
  xlim(-60, 1300) +
  ylim(-0.6, 0.5) +
  theme_Publication()

