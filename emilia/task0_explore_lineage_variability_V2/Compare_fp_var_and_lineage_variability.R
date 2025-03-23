rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs//'

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
                    week5_DABTRAM = "#623594")

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
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

lin_var.day10_DABTRAM <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'))
lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_CIS <- read.csv(paste0(results_dir, 'day10_CIS/lineage_variability_day10_CIS_saver_sample.csv'))

# =============================================================================
# Week5 metadata
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

metadat.week5 <- metadat[metadat$dataset %in% c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS'),]
lin.size.w5 <- table(metadat.week5$assigned_lineage, metadat.week5$dataset)
lin.size.w5 <- as.data.frame(lin.size.w5)
lin.size.w5 <- spread(lin.size.w5, key = 'Var2', value = 'Freq')
colnames(lin.size.w5)[1] <- 'assigned_lineage'


# ==============================================================================
# Calculate fate potential variance
# ==============================================================================

# DABTRAM
fp.d10_w5.DABTRAM <- as.data.frame(all_data@misc[["fatepotential_DABTRAM_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.DABTRAM) <- 'fatepotential_DABTRAM_d10_w5'
fp.d10_w5.DABTRAM$cell_id <- rownames(fp.d10_w5.DABTRAM)
fp.d10_w5.DABTRAM <- merge(fp.d10_w5.DABTRAM, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.DABTRAM <- fp.d10_w5.DABTRAM %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
            n_cells = n())

# COCL2
fp.d10_w5.COCL2 <- as.data.frame(all_data@misc[["fatepotential_COCL2_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.COCL2) <- 'fatepotential_COCL2_d10_w5'
fp.d10_w5.COCL2$cell_id <- rownames(fp.d10_w5.COCL2)
fp.d10_w5.COCL2 <- merge(fp.d10_w5.COCL2, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.COCL2 <- fp.d10_w5.COCL2 %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            n_cells = n())

# CIS
fp.d10_w5.CIS <- as.data.frame(all_data@misc[["fatepotential_CIS_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.CIS) <- 'fatepotential_CIS_d10_w5'
fp.d10_w5.CIS$cell_id <- rownames(fp.d10_w5.CIS)
fp.d10_w5.CIS <- merge(fp.d10_w5.CIS, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.CIS <- fp.d10_w5.CIS %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_CIS_d10_w5, na.rm = TRUE),
            n_cells = n())

# ==============================================================================
# Assemble data
# ==============================================================================

# DABTRAM
df.day10_DABTRAM <- merge(lin_var.day10_DABTRAM, var.fp.d10_w5.DABTRAM, by = 'assigned_lineage')

p1 <- ggplot(df.day10_DABTRAM, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#9D85BE') +
  ylab('Lineage variability') +
  stat_cor(method = 'spearman') +
  theme_Publication()

# COCL2
df.day10_COCL2 <- merge(lin_var.day10_COCL2, var.fp.d10_w5.COCL2, by = 'assigned_lineage')

p2 <- ggplot(df.day10_COCL2, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#6DC49C') +
  ylab('Lineage variability') +
  stat_cor(method = 'spearman') +
  theme_Publication()

# CIS
df.day10_CIS <- merge(lin_var.day10_CIS, var.fp.d10_w5.CIS, by = 'assigned_lineage')

p3 <- ggplot(df.day10_CIS, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#FBD08C') +
  ylab('Lineage variability') +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggarrange(p1, p2, p3, ncol = 3)

# ==============================================================================
# Assemble data 2
# ==============================================================================
# DABTRAM
df.day10_DABTRAM2 <- merge(lin_var.day10_DABTRAM, lin.size.w5[, c('assigned_lineage', 'week5_DABTRAM')], by = 'assigned_lineage', all.x = T)
df.day10_DABTRAM2$week5_DABTRAM[is.na(df.day10_DABTRAM2$week5_DABTRAM)] <- 0
df.day10_DABTRAM2$log10_week5_DABTRAM <- log10(df.day10_DABTRAM2$week5_DABTRAM + 1)

p1 <- ggplot(df.day10_DABTRAM2, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_DABTRAM)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#9D85BE') +
  stat_cor(method = 'spearman') +
  theme_Publication()

# COCL2
df.day10_COCL22 <- merge(lin_var.day10_COCL2, lin.size.w5[, c('assigned_lineage', 'week5_COCL2')], by = 'assigned_lineage', all.x = T)
df.day10_COCL22$week5_COCL2[is.na(df.day10_COCL22$week5_COCL2)] <- 0
df.day10_COCL22$log10_week5_COCL2 <- log10(df.day10_COCL22$week5_COCL2 + 1)

p2 <- ggplot(df.day10_COCL22, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_COCL2)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  stat_cor(method = 'spearman') +
  theme_Publication()

# CIS
df.day10_CIS2 <- merge(lin_var.day10_CIS, lin.size.w5[, c('assigned_lineage', 'week5_CIS')], by = 'assigned_lineage', all.x = T)
df.day10_CIS2$week5_CIS[is.na(df.day10_CIS2$week5_CIS)] <- 0
df.day10_CIS2$log10_week5_CIS <- log10(df.day10_CIS2$week5_CIS + 1)

p3 <- ggplot(df.day10_CIS2, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_CIS)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#FBD08C') +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggarrange(p1, p2, p3, ncol = 3)
