library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# in_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task1_ANOVA_lineage_specific_features/'
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_ANOVA_lineage_specific_features_V2/'

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
anova.day0 <- read.csv(paste0(in_dir, 'day0_SAVER_RNA_processed_ANOVA_pvals.csv'))
anova.day10CIS <- read.csv(paste0(in_dir, 'day10_CIS_SAVER_RNA_processed_ANOVA_pvals.csv'))
anova.day10COCL2 <- read.csv(paste0(in_dir, 'day10_COCL2_SAVER_RNA_processed_ANOVA_pvals.csv'))
anova.day10DABTRAM <- read.csv(paste0(in_dir, 'day10_DABTRAM_SAVER_RNA_processed_ANOVA_pvals.csv'))
anova.week5CIS <- read.csv(paste0(in_dir, 'week5_CIS_SAVER_RNA_processed_ANOVA_pvals.csv'))
anova.week5COCL2 <- read.csv(paste0(in_dir, 'week5_COCL2_SAVER_RNA_processed_ANOVA_pvals.csv'))
anova.week5DABTRAM <- read.csv(paste0(in_dir, 'week5_DABTRAM_SAVER_RNA_processed_ANOVA_pvals.csv'))

# ==============================================================================
# Genes of interest
# ==============================================================================
keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
  DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
  COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
                 "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
                 "CAV1"))
)

memGene <- list('GAPDH', 'LMNA', 'SOX10', 'BABAM1', 'MITF', 'KDM5A', 'FGFR1', 'KDM5B', 'LOXL2', 
                'CCNA2', 'SERPINE1', 'NRG1', 'VEGFC', 'EGFR', 'JUN', 'WNT5A', 'AXL', 'NGFR', 'PDGFRB',
                'PDGFC', 'RUNX2', 'VGF', 'FOSL1')

keygenes <- unlist(keygenes)
keygenes <- unlist(memGene)
# ==============================================================================
# Wrangle data
# ==============================================================================

anova.day0$timepoint <- 'day0'
anova.day0$treatment <- 'baseline'
anova.day0$p_val <- as.numeric(anova.day0$p_val)
anova.day0$p_val_adj <- p.adjust(anova.day0$p_val, method = 'BH')

anova.day10CIS$timepoint <- 'day10'
anova.day10CIS$treatment <- 'CIS'
anova.day10CIS$p_val <- as.numeric(anova.day10CIS$p_val)
anova.day10CIS$p_val_adj <- p.adjust(anova.day10CIS$p_val, method = 'BH')

anova.day10COCL2$timepoint <- 'day10'
anova.day10COCL2$treatment <- 'COCL2'
anova.day10COCL2$p_val <- as.numeric(anova.day10COCL2$p_val)
anova.day10COCL2$p_val_adj <- p.adjust(anova.day10COCL2$p_val, method = 'BH')

anova.day10DABTRAM$timepoint <- 'day10'
anova.day10DABTRAM$treatment <- 'DABTRAM'
anova.day10DABTRAM$p_val <- as.numeric(anova.day10DABTRAM$p_val)
anova.day10DABTRAM$p_val_adj <- p.adjust(anova.day10DABTRAM$p_val, method = 'BH')

anova.week5CIS$timepoint <- 'week5'
anova.week5CIS$treatment <- 'CIS'
anova.week5CIS$p_val <- as.numeric(anova.week5CIS$p_val)
anova.week5CIS$p_val_adj <- p.adjust(anova.week5CIS$p_val, method = 'BH')

anova.week5COCL2$timepoint <- 'week5'
anova.week5COCL2$treatment <- 'COCL2'
anova.week5COCL2$p_val <- as.numeric(anova.week5COCL2$p_val)
anova.week5COCL2$p_val_adj <- p.adjust(anova.week5COCL2$p_val, method = 'BH')

anova.week5DABTRAM$timepoint <- 'week5'
anova.week5DABTRAM$treatment <- 'DABTRAM'
anova.week5DABTRAM$p_val <- as.numeric(anova.week5DABTRAM$p_val)
anova.week5DABTRAM$p_val_adj <- p.adjust(anova.week5DABTRAM$p_val, method = 'BH')


to_plot <- rbind(anova.day0, anova.day10CIS)
to_plot <- rbind(to_plot, anova.day10COCL2)
to_plot <- rbind(to_plot, anova.day10DABTRAM)
to_plot <- rbind(to_plot, anova.week5CIS)
to_plot <- rbind(to_plot, anova.week5COCL2)
to_plot <- rbind(to_plot, anova.week5DABTRAM)

to_plot$dataset <- paste0(to_plot$timepoint, '_', to_plot$treatment)
to_plot$neg_log10_pval <- (-1) * log10(to_plot$p_val)
to_plot$neg_log10_pval <- ifelse(to_plot$neg_log10_pval > 300, 300, to_plot$neg_log10_pval)

to_plot <- to_plot[to_plot$feature %in% keygenes,]
to_plot$isSig <- ifelse(to_plot$p_val_adj < 0.05, "sig", "not_sig")
to_plot
# ==============================================================================
# Plotting
# ==============================================================================
order <-  c("day0_baseline", "day10_DABTRAM", "day10_COCL2", "day10_CIS", "week5_DABTRAM", "week5_COCL2", "week5_CIS")
dataset_colors <- c(day0_baseline = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

# dot plot 
to_plot.sig <- to_plot[to_plot$isSig == "sig",]
to_plot.not_sig <- to_plot[to_plot$isSig == "not_sig",]

p <- ggplot() +
  geom_point(data = to_plot.sig, aes(x = dataset, y = feature, fill = neg_log10_pval, size = neg_log10_pval), shape = 21) +
  geom_point(data = to_plot.not_sig, aes(x = dataset, y = feature), shape = 21, color = 'gray') +
  # scale_fill_manual(values = dataset_colors) +
  scale_fill_gradient(low = "gray", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "right") +
  # facet_wrap(~dataset, scales = "free_x", ncol = 1) +
  # coord_flip() +
  # scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) +
  labs(subtitle = "ANOVA p-values",
       x = "Gene",
       y = "-log10(p-value)")
p
ggsave('~/Downloads/ANOVA_pvalues.png', p, width = 4, height = 5, units = 'in', dpi = 300)


# ==============================================================================
# Plot ANOVA illustration violin plot
# ==============================================================================
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver

remove_unassigned_cells <- T
# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

gene <- 'AXL'
metadat <- all_data@meta.data
metadat.day0 <- metadat[metadat$dataset == 'day0',]
metadat.day0$cell_id <- rownames(metadat.day0)

metadat.day10Dabtram <- metadat[metadat$dataset == 'day10_DABTRAM',]
metadat.day10Dabtram$cell_id <- rownames(metadat.day10Dabtram)

metadat.week5Dabtram <- metadat[metadat$dataset == 'week5_DABTRAM',]
metadat.week5Dabtram$cell_id <- rownames(metadat.week5Dabtram)

# day 0
saver.gene <- as.data.frame(all_data[["Saver"]]@data[gene, rownames(metadat.day0)])
colnames(saver.gene) <- gene
saver.gene$cell_id <- rownames(saver.gene)
saver.gene <- merge(saver.gene, metadat.day0[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

saver.gene.meta <- saver.gene %>% 
  group_by(assigned_lineage) %>% 
  summarise(n_cells = n(),
            mean_value = mean(AXL),
            var_value = var(AXL)) %>% 
  filter(n_cells > 3) %>% 
  arrange(desc(n_cells))

saver.gene.fil <- saver.gene[saver.gene$assigned_lineage %in% saver.gene.meta$assigned_lineage, ]
summary(aov(AXL ~ assigned_lineage, saver.gene.fil))

saver.gene.toplot <- saver.gene[saver.gene$assigned_lineage %in% saver.gene.meta$assigned_lineage[1:10], ]
# length(unique(saver.gene.toplot$assigned_lineage))

ggplot(saver.gene.toplot, aes(x = reorder(assigned_lineage, -AXL, median), y = AXL)) +
  geom_violin(scale = 'width', fill = 'gray') +
  # geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave('~/downloads/ANOVA_ns.png', width = 4, height = 3, dpi = 300)


# day 10 / week 5 DABTRAM
saver.gene <- as.data.frame(all_data[["Saver"]]@data[gene, rownames(metadat.week5Dabtram)])
colnames(saver.gene) <- gene
saver.gene$cell_id <- rownames(saver.gene)
saver.gene <- merge(saver.gene, metadat.week5Dabtram[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

saver.gene.meta <- saver.gene %>% 
  group_by(assigned_lineage) %>% 
  summarise(n_cells = n(),
            mean_value = mean(AXL),
            var_value = var(AXL)) %>% 
  filter(n_cells > 3) %>% 
  arrange(desc(var_value))

saver.gene.fil <- saver.gene[saver.gene$assigned_lineage %in% saver.gene.meta$assigned_lineage, ]
summary(aov(AXL ~ assigned_lineage, saver.gene.fil))

saver.gene.toplot <- saver.gene[saver.gene$assigned_lineage %in% saver.gene.meta$assigned_lineage[1:10], ]
# length(unique(saver.gene.toplot$assigned_lineage))

ggplot(saver.gene.toplot, aes(x = reorder(assigned_lineage, -AXL, median), y = AXL)) +
  geom_violin(scale = 'width', fill = '#623594', alpha = 0.5) +
  # geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave('~/downloads/ANOVA_sig.png', width = 4, height = 3, dpi = 300)
