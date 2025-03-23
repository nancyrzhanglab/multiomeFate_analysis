rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

# data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
# results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/'

data_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup10a/'
results_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig2/'

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

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594",
                    shuffled.day0 = '#F0F0F0',
                    shuffled.day10_DABTRAM = '#F0F0F0',
                    shuffled.day10_COCL2 = '#F0F0F0',
                    shuffled.day10_CIS = '#F0F0F0',
                    shuffled.week5_DABTRAM = '#F0F0F0',
                    shuffled.week5_COCL2 = '#F0F0F0',
                    shuffled.week5_CIS = '#F0F0F0')

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# read lineage variability scores
lin_var.day0 <- read.csv(paste0(results_dir, 'day0/lineage_variability_day0_saver_sample.csv'))
lin_var.day0_Shuffled <- read.csv(paste0(results_dir, 'day0/lineage_variability_shuffledday0_saver_sample.csv'))

lin_var.day10_DABTRAM <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'))
lin_var.day10_DABTRAM_Shuffled <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_shuffledday10_DABTRAM_saver_sample.csv'))

lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_COCL2_Shuffled <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_shuffledday10_COCL2_saver_sample.csv'))
lin_var.day10_COCL2_atac <- read.csv(paste0(results_dir, 'day10_COCL2/embeding_treatment/lineage_variability_day10_COCL2_peakvi.csv'))

lin_var.day10_CIS <- read.csv(paste0(results_dir, 'day10_CIS/lineage_variability_day10_CIS_saver_sample.csv'))
lin_var.day10_CIS_Shuffled <- read.csv(paste0(results_dir, 'day10_CIS/lineage_variability_shuffledday10_CIS_saver_sample.csv'))

lin_var.week5_DABTRAM <- read.csv(paste0(results_dir, 'week5_DABTRAM/lineage_variability_week5_DABTRAM_saver_sample.csv'))
lin_var.week5_DABTRAM_Shuffled <- read.csv(paste0(results_dir, 'week5_DABTRAM/lineage_variability_shuffledweek5_DABTRAM_saver_sample.csv'))

lin_var.week5_COCL2 <- read.csv(paste0(results_dir, 'week5_COCL2/lineage_variability_week5_COCL2_saver_sample.csv'))
lin_var.week5_COCL2_Shuffled <- read.csv(paste0(results_dir, 'week5_COCL2/lineage_variability_shuffledweek5_COCL2_saver_sample.csv'))

lin_var.week5_CIS <- read.csv(paste0(results_dir, 'week5_CIS/lineage_variability_week5_CIS_saver_sample.csv'))
lin_var.week5_CIS_Shuffled <- read.csv(paste0(results_dir, 'week5_CIS/lineage_variability_shuffledweek5_CIS_saver_sample.csv'))

# =============================================================================
# assembliing data
# =============================================================================
lin_var.day0$dataset <- 'day0'
lin_var.day0$category <- 'day0'
lin_var.day0_Shuffled$dataset <- 'day0'
lin_var.day0_Shuffled$category <- 'shuffled.day0'

lin_var.day10_COCL2$dataset <- 'day10_COCL2'
lin_var.day10_COCL2$category <- 'day10_COCL2'
lin_var.day10_COCL2_Shuffled$dataset <- 'day10_COCL2'
lin_var.day10_COCL2_Shuffled$category <- 'shuffled.day10_COCL2'

lin_var.day10_DABTRAM$dataset <- 'day10_DABTRAM'
lin_var.day10_DABTRAM$category <- 'day10_DABTRAM'
lin_var.day10_DABTRAM_Shuffled$dataset <- 'day10_DABTRAM'
lin_var.day10_DABTRAM_Shuffled$category <- 'shuffled.day10_DABTRAM'

lin_var.day10_CIS$dataset <- 'day10_CIS'
lin_var.day10_CIS$category <- 'day10_CIS'
lin_var.day10_CIS_Shuffled$dataset <- 'day10_CIS'
lin_var.day10_CIS_Shuffled$category <- 'shuffled.day10_CIS'

lin_var.week5_DABTRAM$dataset <- 'week5_DABTRAM'
lin_var.week5_DABTRAM$category <- 'week5_DABTRAM'
lin_var.week5_DABTRAM_Shuffled$dataset <- 'week5_DABTRAM'
lin_var.week5_DABTRAM_Shuffled$category <- 'shuffled.week5_DABTRAM'

lin_var.week5_COCL2$dataset <- 'week5_COCL2'
lin_var.week5_COCL2$category <- 'week5_COCL2'
lin_var.week5_COCL2_Shuffled$dataset <- 'week5_COCL2'
lin_var.week5_COCL2_Shuffled$category <- 'shuffled.week5_COCL2'

lin_var.week5_CIS$dataset <- 'week5_CIS'
lin_var.week5_CIS$category <- 'week5_CIS'
lin_var.week5_CIS_Shuffled$dataset <- 'week5_CIS'
lin_var.week5_CIS_Shuffled$category <- 'shuffled.week5_CIS'


df.day10_COCL2 <- rbind(lin_var.day10_COCL2, lin_var.day10_COCL2_Shuffled)
df.day10_COCL2$category <- factor(df.day10_COCL2$category, levels = c('shuffled', 'day10_COCL2'))
# plot density plot in df.day10_COCL2, filled by category
ggplot(df.day10_COCL2, aes(x = normalized_avg_eud_dist_by_shuffle)) +
  geom_density(aes(fill = category), alpha = 0.5) +
  scale_fill_manual(values = c('day10_COCL2' = '#6DC49C', 'shuffled' = '#F0F0F0')) +
  theme_Publication()

df <- rbind(lin_var.day0, lin_var.day0_Shuffled)
df <- rbind(df, lin_var.day10_COCL2)
df <- rbind(df, lin_var.day10_COCL2_Shuffled)
df <- rbind(df, lin_var.day10_DABTRAM)
df <- rbind(df, lin_var.day10_DABTRAM_Shuffled)
df <- rbind(df, lin_var.day10_CIS)
df <- rbind(df, lin_var.day10_CIS_Shuffled)
df <- rbind(df, lin_var.week5_DABTRAM)
df <- rbind(df, lin_var.week5_DABTRAM_Shuffled)
df <- rbind(df, lin_var.week5_COCL2)
df <- rbind(df, lin_var.week5_COCL2_Shuffled)
df <- rbind(df, lin_var.week5_CIS)
df <- rbind(df, lin_var.week5_CIS_Shuffled)

df$category <- factor(df$category, levels = c('shuffled.day0', 'day0', 
                                              'shuffled.day10_DABTRAM', 'day10_DABTRAM', 
                                              'shuffled.day10_COCL2', 'day10_COCL2',
                                              'shuffled.day10_CIS', 'day10_CIS', 
                                              'shuffled.week5_DABTRAM', 'week5_DABTRAM',
                                              'shuffled.week5_COCL2', 'week5_COCL2',
                                              'shuffled.week5_CIS', 'week5_CIS'))


# =============================================================================
# plot overview of lineage variability of all data
# =============================================================================
ggplot(df, aes(x = dataset, y= normalized_avg_eud_dist_by_shuffle)) +
  geom_violin(aes(fill = category), scale = 'width', width = 0.8) +
  geom_boxplot(aes(group = category), width = 0.2, position = position_dodge(0.8), fill = 'white', outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  theme_Publication()


# =============================================================================
# Plot day10 COCL2
# =============================================================================

# plot COCL2 violin plot
df1 <- rbind(lin_var.day10_COCL2, lin_var.day10_COCL2_Shuffled)
ggplot(df1, aes(x = dataset, y= normalized_avg_eud_dist_by_shuffle)) +
  geom_violin(aes(fill = category), scale = 'width', width = 0.8) +
  geom_boxplot(aes(group = category), width = 0.2, position = position_dodge(0.8), fill = 'white', outlier.shape = NA) +
  scale_fill_manual(values = dataset_colors) +
  ylab('Lineage variability (RNA)') +
  theme_Publication()
ggsave(paste0(figure_dir, 'day10_COCL2_lin_var_RNA.pdf'), width = 2.5, height = 3)


# get lineage variability by RNA
lin_var.rna <- lin_var.day10_COCL2[, c('assigned_lineage', 'normalized_avg_eud_dist_by_shuffle')]
lin_var.rna$modality <- 'RNA'
colnames(lin_var.rna) <- c('assigned_lineage', 'lin_var.RNA', 'modality')

# get lineage variability by ATAC
lin_var.atac <- lin_var.day10_COCL2_atac[, c('assigned_lineage', 'normalized_avg_eud_dist_by_shuffle', 'n_cells')]
lin_var.atac$modality <- 'ATAC'
colnames(lin_var.atac) <- c('assigned_lineage', 'lin_var.ATAC', 'n_cells', 'modality')

# merge
rna.atac.comp <- merge(lin_var.rna, lin_var.atac, by = 'assigned_lineage')

p1 <- ggplot(rna.atac.comp, aes(x = lin_var.RNA, y = lin_var.ATAC)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  stat_cor() +
  xlab('Lineage variability (RNA)') +
  ylab('Lineage variability (ATAC)') +
  theme_Publication()
ggsave(paste0(figure_dir, '/day10_COCL2_RNA_ATAC_comparison.pdf'), width = 3.5, height = 4)


# overlay with UMAP
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))
umap <- as.data.frame(all_data_ft_COCL2_umap@cell.embeddings)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.day10_COCL2 <- metadat[metadat$dataset == 'day10_COCL2', ]
metadat.COCL2 <- metadat[metadat$dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2'), ]

umap.COCL2 <- umap[metadat.COCL2$cell_id, ]
umap.day10_COCL2 <- umap[metadat.day10_COCL2$cell_id, ]

lin1.cells <- metadat.COCL2[metadat.COCL2$assigned_lineage == 'Lin70618', 'cell_id']
umap.day10_COCL2.lin1 <- umap.day10_COCL2[lin1.cells, ] # high variance

lin2.cells <- metadat.COCL2[metadat.COCL2$assigned_lineage == 'Lin34625', 'cell_id']
umap.day10_COCL2.lin2 <- umap.day10_COCL2[lin2.cells, ] # low variance

p2 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#F0F0F0') +
  geom_point(data = umap.day10_COCL2.lin1,  size = 2, color = '#6DC49C') +
  ggtitle('High var.') +
  theme_Publication()

p3 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#F0F0F0') +
  geom_point(data = umap.day10_COCL2.lin2,  size = 2, color = '#6DC49C') +
  ggtitle('Low var.') +
  theme_Publication()

p4 <- ggarrange(p1, 
          ggarrange(p2, p3, ncol = 1, nrow = 2),
          widths = c(1, 0.6))

ggsave(paste0(figure_dir, '/Fig2D_day10_COCL2_lineage_variability_examples.pdf'), p4, width = 5, height = 4)
