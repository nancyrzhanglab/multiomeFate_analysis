rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig4/'

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

theme_Clean<- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            axis.text = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour=NA,fill=NA),
            panel.background = element_blank(),
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

var.fp.d10_w5.DABTRAM <- read.csv(paste0(out_dir, 'var.fp.d10_w5.DABTRAM.csv'))
var.fp.d10_w5.COCL2 <- read.csv(paste0(out_dir, 'var.fp.d10_w5.COCL2.csv'))
var.fp.d10_w5.CIS <- read.csv(paste0(out_dir, 'var.fp.d10_w5.CIS.csv'))

lin_var.day10_DABTRAM <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'))
lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_CIS <- read.csv(paste0(results_dir, 'day10_CIS/lineage_variability_day10_CIS_saver_sample.csv'))

# ==============================================================================
# Calculate future lineage size
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

metadat.week5 <- metadat[metadat$dataset %in% c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS'),]
lin.size.w5 <- table(metadat.week5$assigned_lineage, metadat.week5$dataset)
lin.size.w5 <- as.data.frame(lin.size.w5)
lin.size.w5 <- spread(lin.size.w5, key = 'Var2', value = 'Freq')
colnames(lin.size.w5)[1] <- 'assigned_lineage'

# =============================================================================
# Assemble data
# =============================================================================
var.fp.d10_w5.DABTRAM <- merge(var.fp.d10_w5.DABTRAM, lin.size.w5[, c('assigned_lineage', 'week5_DABTRAM')], by = 'assigned_lineage', all.x = T)
var.fp.d10_w5.COCL2 <- merge(var.fp.d10_w5.COCL2, lin.size.w5[, c('assigned_lineage', 'week5_COCL2')], by = 'assigned_lineage', all.x = T)
var.fp.d10_w5.CIS <- merge(var.fp.d10_w5.CIS, lin.size.w5[, c('assigned_lineage', 'week5_CIS')], by = 'assigned_lineage', all.x = T)


lin.size.thres <- 5
var.fp.d10_w5.DABTRAM <- var.fp.d10_w5.DABTRAM %>% filter(n_cells > lin.size.thres)
var.fp.d10_w5.DABTRAM$week5_DABTRAM[is.na(var.fp.d10_w5.DABTRAM$week5_DABTRAM)] <- 0
var.fp.d10_w5.DABTRAM$log10_week5_DABTRAM <- log10(var.fp.d10_w5.DABTRAM$week5_DABTRAM + 1)
var.fp.d10_w5.DABTRAM$category <- ifelse(var.fp.d10_w5.DABTRAM$variance > 0.07693912, 'High variance', 'Low variance')

var.fp.d10_w5.COCL2 <- var.fp.d10_w5.COCL2 %>% filter(n_cells > lin.size.thres)
var.fp.d10_w5.COCL2$week5_COCL2[is.na(var.fp.d10_w5.COCL2$week5_COCL2)] <- 0
var.fp.d10_w5.COCL2$log10_week5_COCL2 <- log10(var.fp.d10_w5.COCL2$week5_COCL2 + 1)
var.fp.d10_w5.COCL2$category <- ifelse(var.fp.d10_w5.COCL2$variance > 3.9212267, 'High variance', NA)
var.fp.d10_w5.COCL2$category <- ifelse(var.fp.d10_w5.COCL2$variance < 1.8213771, 'Low variance', var.fp.d10_w5.COCL2$category)

var.fp.d10_w5.CIS <- var.fp.d10_w5.CIS %>% filter(n_cells > lin.size.thres)
var.fp.d10_w5.CIS$week5_CIS[is.na(var.fp.d10_w5.CIS$week5_CIS)] <- 0
var.fp.d10_w5.CIS$log10_week5_CIS <- log10(var.fp.d10_w5.CIS$week5_CIS + 1)

# ==============================================================================
# Assemble data
# ==============================================================================

# DABTRAM
df.day10_DABTRAM <- merge(lin_var.day10_DABTRAM, var.fp.d10_w5.DABTRAM, by = 'assigned_lineage')

p1 <- ggplot(df.day10_DABTRAM, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#9D85BE') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('Lineage variability') +
  xlab('Growth potential standard deviation') +
  stat_cor(method = 'spearman') +
  theme_Publication()

p1

# COCL2
df.day10_COCL2 <- merge(lin_var.day10_COCL2, var.fp.d10_w5.COCL2, by = 'assigned_lineage')

p2 <- ggplot(df.day10_COCL2, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#6DC49C') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('Lineage variability') +
  xlab('Growth potential standard deviation') +
  stat_cor(method = 'spearman') +
  theme_Publication()
p2

p3 <- ggarrange(p2, p1, ncol = 1)

ggsave(paste0(figure_dir, 'scatter_plot_lin_var_growth_potential_sd_DABTRAM_COCL2.png'), p3, width = 5, height = 8.5)

p1.cl <- ggplot(df.day10_DABTRAM, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#9D85BE') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('') +
  xlab('') +
  theme_Clean()


# COCL2
p2.cl <- ggplot(df.day10_COCL2, aes(x = sqrt(variance), y = normalized_avg_eud_dist_by_shuffle)) +
  geom_point(aes(size = n_cells.x), shape = 21, fill = '#6DC49C') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('') +
  xlab('') +
  theme_Clean()

p3.cl <- ggarrange(p2.cl, p1.cl, ncol = 1)
ggsave(paste0(figure_dir, 'scatter_plot_lin_var_growth_potential_sd_DABTRAM_COCL2_clean.pdf'), p3.cl, width = 3, height = 3.8)

# =============================================================================
# Plot
# =============================================================================

p.fp.var.dabtram <- ggplot(var.fp.d10_w5.DABTRAM, aes(x = sqrt(variance), y = log10_week5_DABTRAM)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#9D85BE') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  stat_cor(method = 'spearman') +
  theme_Publication()
ggsave(paste0(figure_dir, 'Supp_scatter_plot_growth_potential_sd_w5_size_DABTRAM.pdf'), p.fp.var.dabtram, width = 3, height = 2.4)

ggplot(var.fp.d10_w5.COCL2, aes(x = sqrt(variance), y = log10_week5_COCL2)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggplot(var.fp.d10_w5.CIS, aes(x = sqrt(variance), y = log10_week5_CIS)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#FBD08C') +
  stat_cor(method = 'spearman') +
  theme_Publication()



# comp.COCL2 <- merge(fp.d10_w5.COCL2, cell_imputed_score_df, by = 'cell_id', all = T)
# comp.DABTRAM <- merge(fp.d10_w5.DABTRAM, cell_imputed_score_df, by = 'cell_id', all = T)
# comp.CIS <- merge(fp.d10_w5.CIS, cell_imputed_score_df, by = 'cell_id', all = T)
# 
# 
# ggplot(comp.COCL2, aes(x = fatepotential_COCL2_d10_w5, y = cell_imputed_score)) +
#   geom_point(shape = 21, fill = '#6DC49C') +
#   stat_cor(method = 'spearman') 
# 
# ggplot(comp.DABTRAM, aes(x = fatepotential_DABTRAM_d10_w5, y = cell_imputed_score)) +
#   geom_point(shape = 21, fill = '#9D85BE') +
#   stat_cor(method = 'spearman') 
# 
# ggplot(comp.CIS, aes(x = fatepotential_CIS_d10_w5, y = cell_imputed_score)) +
#   geom_point(shape = 21, fill = '#FBD08C') +
#   stat_cor(method = 'spearman') 

# ==============================================================================
# Assemble data 2
# ==============================================================================
# DABTRAM
df.day10_DABTRAM2 <- merge(lin_var.day10_DABTRAM, lin.size.w5[, c('assigned_lineage', 'week5_DABTRAM')], by = 'assigned_lineage', all.x = T)
df.day10_DABTRAM2$week5_DABTRAM[is.na(df.day10_DABTRAM2$week5_DABTRAM)] <- 0
df.day10_DABTRAM2$log10_week5_DABTRAM <- log10(df.day10_DABTRAM2$week5_DABTRAM + 1)

p4 <- ggplot(df.day10_DABTRAM2, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_DABTRAM)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#9D85BE') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  stat_cor(method = 'spearman') +
  ylab('log10(week5 + 1)') +
  xlab('Lineage variability') +
  theme_Publication()

# COCL2
df.day10_COCL22 <- merge(lin_var.day10_COCL2, lin.size.w5[, c('assigned_lineage', 'week5_COCL2')], by = 'assigned_lineage', all.x = T)
df.day10_COCL22$week5_COCL2[is.na(df.day10_COCL22$week5_COCL2)] <- 0
df.day10_COCL22$log10_week5_COCL2 <- log10(df.day10_COCL22$week5_COCL2 + 1)

p5 <- ggplot(df.day10_COCL22, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_COCL2)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  stat_cor(method = 'spearman') +
  ylab('log10(week5 + 1)') +
  xlab('Lineage variability') +
  theme_Publication()


p6 <- ggarrange(p5, p4, ncol = 1)

ggsave(paste0(figure_dir, 'scatter_plot_lin_var_w5_size_DABTRAM_COCL2.png'), p6, width = 5, height = 8.5)


p4.cl <- ggplot(df.day10_DABTRAM2, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_DABTRAM)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#9D85BE') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('') +
  xlab('') +
  theme_Clean()

p5.cl <- ggplot(df.day10_COCL22, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_COCL2)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  scale_size_continuous(range = c(1, 8), breaks = c(20, 40, 60, 80)) +
  ylab('') +
  xlab('') +
  theme_Clean()


p6.cl <- ggarrange(p5.cl, p4.cl, ncol = 1)

ggsave(paste0(figure_dir, 'scatter_plot_lin_var_w5_size_DABTRAM_COCL2_clean.pdf'), p6.cl, width = 3, height = 3.8)

# CIS
df.day10_CIS2 <- merge(lin_var.day10_CIS, lin.size.w5[, c('assigned_lineage', 'week5_CIS')], by = 'assigned_lineage', all.x = T)
df.day10_CIS2$week5_CIS[is.na(df.day10_CIS2$week5_CIS)] <- 0
df.day10_CIS2$log10_week5_CIS <- log10(df.day10_CIS2$week5_CIS + 1)

p3 <- ggplot(df.day10_CIS2, aes(x = normalized_avg_eud_dist_by_shuffle, y = log10_week5_CIS)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#FBD08C') +
  stat_cor(method = 'spearman') +
  ylab('log10(week5 + 1)') +
  xlab('Lineage variability') +
  theme_Publication()

ggarrange(p1, p2, p3, ncol = 3)
