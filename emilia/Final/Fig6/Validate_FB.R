rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'


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

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data@misc <- all_data_fatepotential

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

df.bias.DABTRAM <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
df.bias.COCL2 <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_COCL2.csv'))
df.bias.CIS <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_CIS.csv'))

# =============================================================================
# QC
# =============================================================================

nrow(df.bias.DABTRAM[df.bias.DABTRAM$bias > 0.5, ]) / 3326
nrow(df.bias.COCL2[df.bias.COCL2$bias > 0.5, ]) / 3326
nrow(df.bias.CIS[df.bias.CIS$bias > 0.5, ]) / 3326

p1 <- ggplot(df.bias.DABTRAM, aes(x = bias)) +
  geom_histogram(bins = 100, color = 'black', fill = 'lightgray') +
  theme_Publication() +
  xlab('Adapting fate bias') +
  ylab('Number of cells') +
  ggtitle('DABTRAM') +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(df.bias.COCL2, aes(x = bias)) +
  geom_histogram(bins = 100, color = 'black', fill = 'lightgray') +
  theme_Publication() +
  xlab('Adapting fate bias') +
  ylab('Number of cells') +
  ggtitle('COCL2') +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(df.bias.CIS, aes(x = bias)) +
  geom_histogram(bins = 100, color = 'black', fill = 'lightgray') +
  theme_Publication() +
  xlab('Adapting fate bias') +
  ylab('Number of cells') +
  ggtitle('CIS') +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggarrange(p1, p2, p3,
          ncol = 3,
          label.x = 0.5,
          label.y = 0.95,
          font.label = list(size = 10, face = 'bold'))

ggsave(paste0(figure_dir, 'Supp_adapting_bias_histograms.pdf'), width = 9, height = 3)

# =============================================================================
# Wrangle
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.day0 <- metadat[metadat$dataset == 'day0', ]


metadata.week5.COCL2 <- subset(metadat, dataset == 'week5_COCL2')
metadata.week5.CIS <- subset(metadat, dataset == 'week5_CIS')
metadata.week5.DABTRAM <- subset(metadat, dataset == 'week5_DABTRAM')

df.bias.DABTRAM <- merge(df.bias.DABTRAM, metadat.day0[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
df.bias.COCL2 <- merge(df.bias.COCL2, metadat.day0[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
df.bias.CIS <- merge(df.bias.CIS, metadat.day0[, c("cell_id", 'assigned_lineage')], by = 'cell_id')
# =============================================================================
# Plot violin
# =============================================================================
mean.bias.DABTRAM <- df.bias.DABTRAM %>% 
  group_by(assigned_lineage) %>% 
  summarise(bias.mean = mean(bias), 
            n.cells = n()) %>% 
  arrange(desc(n.cells))

lin.to.plot <- head(mean.bias.DABTRAM, 10) %>% pull(assigned_lineage)

df.bias.DABTRAM$assigned_lineage_label <- ifelse(df.bias.DABTRAM$assigned_lineage %in% lin.to.plot, df.bias.DABTRAM$assigned_lineage, 'Other')

lin.to.anova <- mean.bias.DABTRAM %>% filter(n.cells > 3) %>% pull(assigned_lineage)
anova_res <- stats::oneway.test(bias ~ assigned_lineage, data = df.bias.DABTRAM[df.bias.DABTRAM$assigned_lineage %in% lin.to.anova, ])


ggplot(df.bias.DABTRAM, aes(x = assigned_lineage_label, y = bias)) +
  geom_violin(trim=T, scale = "width", fill = '#A9A9A9') +
  geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.9, size = 1) +
  scale_x_discrete(limits = c(lin.to.plot, "Other"),
                   guide = ggplot2::guide_axis(angle = 45)) +
  ylab("Day 0 to Day 10 adapting fate bias") +
  xlab("") +
  stat_summary(fun=mean, geom="point", shape=16, size=3, color="red") +
  stat_summary(fun=max, geom="point", shape=10, size=5, color="blue") +
  theme_Publication()

ggsave(paste0(figure_dir, 'FB_violin.pdf'), width = 7, height = 3)

# =============================================================================
# Compare in and out of W5
# =============================================================================

# DABTRAM 
fp.d10.w5.DABTRAM <- all_data@misc[["fatepotential_DABTRAM_d10_w5"]][["lineage_imputed_count"]]
fp.d10.w5.DABTRAM <- as.data.frame(fp.d10.w5.DABTRAM)
fp.d10.w5.DABTRAM$assigned_lineage <- rownames(fp.d10.w5.DABTRAM)
fp.d10.w5.DABTRAM.winner <- fp.d10.w5.DABTRAM %>% 
  filter(fp.d10.w5.DABTRAM > 1)

df.bias.DABTRAM$adaptingCount <- 10**df.bias.DABTRAM$adaptingFP
df.bias.DABTRAM$is.lin.in.w5 <- df.bias.DABTRAM$assigned_lineage %in% metadata.week5.DABTRAM$assigned_lineage
df.bias.DABTRAM$is.lin.in.w5.imputed <- df.bias.DABTRAM$assigned_lineage %in% fp.d10.w5.DABTRAM.winner$assigned_lineage

df.DABTRAM.summary <- df.bias.DABTRAM %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(adaptingCount = sum(adaptingCount))

p1 <- ggplot(df.DABTRAM.summary, aes(x = is.lin.in.w5, y = adaptingCount)) +
  geom_violin(scale = 'width', aes(fill = is.lin.in.w5)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  scale_fill_manual(values = c("#3EDBF0", "#FFB433")) +
  ylab('Imputed Adapting Progeny in Lineage') +
  ggtitle('DABTRAM') +
  theme_Publication() 

# COCL2
fp.d10.w5.COCL2 <- all_data@misc[["fatepotential_COCL2_d10_w5"]][["lineage_imputed_count"]]
fp.d10.w5.COCL2 <- as.data.frame(fp.d10.w5.COCL2)
fp.d10.w5.COCL2$assigned_lineage <- rownames(fp.d10.w5.COCL2)
fp.d10.w5.COCL2.winner <- fp.d10.w5.COCL2 %>% 
  filter(fp.d10.w5.COCL2 > 1)

df.bias.COCL2$adaptingCount <- 10**df.bias.COCL2$adaptingFP
df.bias.COCL2$is.lin.in.w5 <- df.bias.COCL2$assigned_lineage %in% metadata.week5.COCL2$assigned_lineage
df.bias.COCL2$is.lin.in.w5.imputed <- df.bias.COCL2$assigned_lineage %in% fp.d10.w5.COCL2.winner$assigned_lineage

df.COCL2.summary <- df.bias.COCL2 %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(adaptingCount = sum(adaptingCount))

p2 <- ggplot(df.COCL2.summary, aes(x = is.lin.in.w5, y = adaptingCount)) +
  geom_violin(scale = 'width', aes(fill = is.lin.in.w5)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  scale_fill_manual(values = c("#3EDBF0", "#FFB433")) +
  ylab('Imputed Adapting Progeny in Lineage') +
  ggtitle('COCL2') +
  theme_Publication() 

# CIS
fp.d10.w5.CIS <- all_data@misc[["fatepotential_CIS_d10_w5"]][["lineage_imputed_count"]]
fp.d10.w5.CIS <- as.data.frame(fp.d10.w5.CIS)
fp.d10.w5.CIS$assigned_lineage <- rownames(fp.d10.w5.CIS)
fp.d10.w5.CIS.winner <- fp.d10.w5.CIS %>% 
  filter(fp.d10.w5.CIS > 1)

df.bias.CIS$adaptingCount <- 10**df.bias.CIS$adaptingFP
df.bias.CIS$is.lin.in.w5 <- df.bias.CIS$assigned_lineage %in% metadata.week5.CIS$assigned_lineage
df.bias.CIS$is.lin.in.w5.imputed <- df.bias.CIS$assigned_lineage %in% fp.d10.w5.CIS.winner$assigned_lineage

df.CIS.summary <- df.bias.CIS %>% 
  group_by(assigned_lineage, is.lin.in.w5, is.lin.in.w5.imputed) %>% 
  summarize(adaptingCount = sum(adaptingCount))

p3 <- ggplot(df.CIS.summary, aes(x = is.lin.in.w5, y = adaptingCount)) +
  geom_violin(scale = 'width', aes(fill = is.lin.in.w5)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE')),
                     paired = FALSE, 
                     method = "wilcox", aes(label =  after_stat(p.signif))) +
  scale_fill_manual(values = c("#3EDBF0", "#FFB433")) +
  ylab('Imputed Adapting Progeny in Lineage') +
  ggtitle('CIS') +
  theme_Publication() 

p4 <- ggarrange(p1, p2, p3, ncol = 3, legend = F)

ggsave(paste0(figure_dir, 'Supp_Violin_Sum_AdaptingCount_panel.pdf'), width = 9, height = 3)
