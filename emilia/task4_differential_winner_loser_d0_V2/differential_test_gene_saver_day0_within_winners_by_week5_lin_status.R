library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
remove_unassigned_cells <- TRUE

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
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_wnn.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[["wnn.umap"]] <- all_data_wnn

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

## Read features to test
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task3_identify_lineage_specific_adapation_features_V2/lineage_specific_adapation_genes_saver.RData')

keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"))
)

# DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
#                  "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
#                  "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
#                  "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
#                  "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44"))

# COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                "CAV1"))
keygenes <- unlist(keygenes)
# ==============================================================================
# Get day0 winner cells
# ==============================================================================

# get fate potential
cur_time = 'd0'
fut_time = 'd10'
treatment = 'CIS'
fp_name = paste0("fatepotential_", treatment, '_', cur_time, '_', fut_time)

fp <- as.data.frame(all_data@misc[[fp_name]][["cell_imputed_score"]])
fp$cell_id <- rownames(fp)
colnames(fp) <- c(fp_name, 'cell_id')

# get winner cells
fp_median <- median(fp[[fp_name]])
fp_winner <- fp[fp[[fp_name]] > fp_median, ]

# get saver
metadat <- all_data@meta.data
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
metadat.day0$cell_id <- rownames(metadat.day0)

saver_day0 <- t(all_data[["Saver"]]@data)
saver_day0 <- saver_day0[metadat.day0$cell_id, ]

saver_day0_winner <- saver_day0[fp_winner$cell_id, ]

# ==============================================================================
# Get lineages survived to Week5 DABTRAM
# ==============================================================================
lin.week5 <- metadat[metadat$dataset == 'week5_CIS', ]$assigned_lineage
lin.day0 <- metadat[metadat$dataset == 'day0', ]$assigned_lineage

lin.day0.week5 <- intersect(lin.day0, lin.week5)

fp_winner <- merge(fp_winner, metadat.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')
fp_winner$LinInWeek5 <- ifelse(fp_winner$assigned_lineage %in% lin.day0.week5, 'Yes', 'No')

table(fp_winner$LinInWeek5)

# ==============================================================================
# Differential tests
# ==============================================================================
columns <- c('feature', 'mean_winner_lin_in_week5', 'mean_winner_lin_notin_week5', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

# features <- lineage_specific_adapation_TFs[[paste0(tolower(treatment), "_", fut_time, "_saver_cor_vec_top25")]][[1]]
features <- colnames(saver_day0_winner)

day0_winner_lin_in_week5 <- fp_winner[fp_winner$LinInWeek5 == 'Yes', 'cell_id']
day0_winner_lin_Notin_week5 <- fp_winner[fp_winner$LinInWeek5 == 'No', 'cell_id']

for (f in features) {
  
  feature_day0_winner_lin_in_week5 <- saver_day0_winner[day0_winner_lin_in_week5, f]
  feature_day0_winner_lin_Notin_week5 <- saver_day0_winner[day0_winner_lin_Notin_week5, f]
  
  feature_day0_winner_lin_in_week5 <- feature_day0_winner_lin_in_week5[!is.na(feature_day0_winner_lin_in_week5)]
  feature_day0_winner_lin_Notin_week5 <- feature_day0_winner_lin_Notin_week5[!is.na(feature_day0_winner_lin_Notin_week5)]
  
  variance <- var(feature_day0_winner_lin_in_week5) + var(feature_day0_winner_lin_Notin_week5)
  
  if(variance == 0) {
    next
  }
  res <- t.test(feature_day0_winner_lin_in_week5,
                feature_day0_winner_lin_Notin_week5,
                alternative = 'two.sided')
  
  t_statistics <- res[["statistic"]][["t"]]
  t_test_p_val <- res[["p.value"]] 
  
  t_test_results[nrow(t_test_results) + 1, ] <- c(
    f, 
    mean(feature_day0_winner_lin_in_week5), 
    mean(feature_day0_winner_lin_Notin_week5), 
    t_statistics, 
    t_test_p_val
  )
}

t_test_results$p_val <- as.numeric(t_test_results$p_val)
t_test_results$p_val_adj <- p.adjust(t_test_results$p_val, method = 'BH')
t_test_results$neg_log10_p_val <- -log10(t_test_results$p_val)

t_test_results$mean_winner_lin_in_week5 <- as.numeric(t_test_results$mean_winner_lin_in_week5)
t_test_results$mean_winner_lin_notin_week5 <- as.numeric(t_test_results$mean_winner_lin_notin_week5)
t_test_results$mean_diff <- log2(t_test_results$mean_winner_lin_in_week5) - log2(t_test_results$mean_winner_lin_notin_week5)
t_test_results$t_statistic <- as.numeric(t_test_results$t_statistic)
# write.csv(t_test_results, paste0('~/Downloads/t_test.csv'), row.names = FALSE)

t_test_results.key.gene <- t_test_results[t_test_results$feature %in% keygenes, ]
p_val_thres <- t_test_results[t_test_results$p_val_adj < 0.05, ]
p_val_thres <- min(p_val_thres$neg_log10_p_val)

ggplot(t_test_results, aes(x = mean_diff, 
                           y = neg_log10_p_val)) +
  geom_point(color = 'gray') +
  geom_point(data = t_test_results.key.gene, aes(x = mean_diff, y = neg_log10_p_val), color = 'red') +
  ggrepel::geom_text_repel(data = t_test_results.key.gene, 
                           aes(label = feature), 
                           box.padding = 0.5,
                           nudge_x = 0.1,
                           nudge_y = 0.1) +
  geom_hline(yintercept = p_val_thres, linetype = 'dashed') +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = 'log2(Fold Change)', y = '-log10(p_val)', title = 'Day 0 (CIS)') +
  theme_bw() +
  theme_Publication() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))
ggsave('~/Downloads/volcano_d0_CIS.png', width = 5, height = 5, dpi = 300)

t_test_results.sig <- t_test_results[t_test_results$p_val_adj < 0.1, ]
t_test_results.sig <- t_test_results.sig[t_test_results.sig$mean_diff > 0, ]
write.csv(t_test_results.sig, '~/Downloads/t_test_cis.csv', row.names = F)


t_test_results <- t_test_results[order(t_test_results$mean_diff, decreasing = T), ]
rownames(t_test_results) <- seq(1, nrow(t_test_results))
t_test_gsea <- t_test_results[, c('feature', 'mean_diff')]
t_test_gsea <- t_test_gsea$mean_diff
names(t_test_gsea) <- t_test_results$feature

# run GSEA using the 'GSEA' function from clusterProfiler
hallmark <- read.gmt("/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/h.all.v2024.1.Hs.symbols.gmt")
t_testGSEA.res <- GSEA(t_test_gsea, TERM2GENE=hallmark, verbose=FALSE)
t_testGSEA.df <- as_tibble(t_testGSEA.res@result)
t_testGSEA.df$treatment <- 'CIS'

# create 'bubble plot' to summarize y signatures across x phenotypes
t_testGSEA.df <- t_testGSEA.df[order(t_testGSEA.df$NES, decreasing = T), ]

write.csv(t_testGSEA.df, '~/Downloads/day0_CIS_GSEA.csv', row.names = F)

ggplot(t_testGSEA.df, aes(x = treatment,  y= ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

# create enrichment plots using the enrichplot package
gseaplot2(t_testGSEA.res, 
          geneSetID = c(3, 8, 11), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = "Top 3 Significant Pathways") #can also turn off this title
