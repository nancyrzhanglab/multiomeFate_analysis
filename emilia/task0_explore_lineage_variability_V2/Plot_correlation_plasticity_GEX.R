library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler) 
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
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/Archive_Sydney/task0_explore_lineage_variability/features_associated_with_plasticity/'
sample_name <- 'day10_DABTRAM'
load(paste0(results_dir, sample_name, '_gene_exp_mean_day10_plasticity_correlation.RData'))


singature_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
hallmark <- read.gmt(paste0(singature_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(singature_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
threeCA <- read.gmt(paste0(singature_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))

# ==============================================================================
# Enrichment analysis
# ==============================================================================
getAndSortTable <- function(df) {
  cor_vec <- df %>% drop_na()
  cor_vec$p.value <- as.numeric(cor_vec$p.value)
  cor_vec$p_adj <- p.adjust(cor_vec$p.value, method = 'BH')
  cor_vec <- cor_vec[cor_vec$p_adj < 0.05, ]
  
  cor_vec$gene <- rownames(cor_vec)
  cor_vec$correlation <- as.numeric(cor_vec$correlation)
  cor_vec <- cor_vec[order(cor_vec$correlation, decreasing = TRUE),]
  
  return(cor_vec)
}

day10_dabtram <- getAndSortTable(cor_vec)

# GSEA DABTRAM
gsea_input_dabtram <- day10_dabtram$correlation
names(gsea_input_dabtram) <- day10_dabtram$gene

set.seed(123)
# Hallmark
GSEA_res_dabtram <- GSEA(geneList = gsea_input_dabtram, 
                         TERM2GENE = hallmark, 
                         pvalueCutoff = 0.2,
                         seed = T,
                         verbose = F)
GSEA_res_dabtram.df <- as_tibble(GSEA_res_dabtram@result)
GSEA_res_dabtram.df <- GSEA_res_dabtram.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')


# Reactome
GSEA_res_dabtram.reactome <- GSEA(geneList = gsea_input_dabtram, 
                         TERM2GENE = reactome, 
                         pvalueCutoff = 0.2,
                         seed = T,
                         verbose = F)
GSEA_res_dabtram.reactome.df <- as_tibble(GSEA_res_dabtram.reactome@result)
GSEA_res_dabtram.reactome.df <- GSEA_res_dabtram.reactome.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.reactome.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')


# 3CA
GSEA_res_dabtram.3CA <- GSEA(geneList = gsea_input_dabtram, 
                                  TERM2GENE = threeCA, 
                                  pvalueCutoff = 0.2,
                                  seed = T,
                                  verbose = F)
GSEA_res_dabtram.3CA.df <- as_tibble(GSEA_res_dabtram.3CA@result)
GSEA_res_dabtram.3CA.df <- GSEA_res_dabtram.3CA.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.3CA.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')


GSEA_res <- rbind(GSEA_res_dabtram.df, GSEA_res_dabtram.3CA.df)
GSEA_res <- GSEA_res[order(GSEA_res$NES.DABTRAM, decreasing = TRUE),]
GSEA_res$order <- as.numeric(rownames(GSEA_res))

features_to_label_top <- head(GSEA_res, 2) %>% select(ID) %>% pull()
emt_related <- c(grep('EMT', GSEA_res$ID, value = T) , 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')

features_to_label_top <- c(features_to_label_top, 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
features_to_label_top <- c(features_to_label_top, 'GAVISH_3CA_MALIGNANT_METAPROGRAM_32_SKIN_PIGMENTATION')

ggplot(GSEA_res, aes(x = order, y = `NES.DABTRAM`)) +
  geom_point(size=0.5) +
  geom_point(data = subset(GSEA_res, ID %in% emt_related), size=2, color = '#0079FF') +
  # geom_point(data = subset(cor_vec, feature %in% sox_related), size=2, color = '#5DB996') +
  # geom_point(data = subset(cor_vec, feature %in% tead_related), size=2, color = '#F0BB78') +
  # geom_point(data = subset(cor_vec, feature %in% mitf), size=2, color = '#F72C5B') +
  # geom_point(data = subset(cor_vec, feature %in% snai_related), size=2, color = '#FFA07A') +
  ggrepel::geom_text_repel(data = subset(GSEA_res, ID %in% features_to_label_top),
                           ggplot2::aes(label = ID),
                           box.padding = ggplot2::unit(0.4, 'lines'),
                           point.padding = ggplot2::unit(0.2, 'lines'),
                           nudge_y = 0.3,
                           nudge_x = 0.1,
                           color = 'black',
                           size=4,
                           max.overlaps = 80) +
  # ggrepel::geom_text_repel(data = subset(cor_vec, feature %in% features_to_label_bottom),
  #                          ggplot2::aes(label = Num),
  #                          box.padding = ggplot2::unit(0.4, 'lines'),
  #                          point.padding = ggplot2::unit(0.2, 'lines'),
  #                          nudge_y = 0.5,
  #                          nudge_x = -2,
  #                          color = 'black',
  #                          size=4,
  #                          max.overlaps = 80) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylab('Normalized Enrichment Score') +
  ggtitle(sample_name) +
  ylim(-3, 3.5) +
  # xlim(-0.5, 645) +
  theme_Publication()
ggsave(paste0(results_dir, sample_name, '_gene_exp_mean_day10_plasticity_byRNA_GSEA_NES.png'), width = 4, height = 5, units = 'in', dpi = 300)
