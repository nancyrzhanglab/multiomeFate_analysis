rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs//'
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

# ==============================================================================
# Calculate fate potential variance
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)


# DABTRAM
fp.d10_w5.DABTRAM <- as.data.frame(all_data@misc[["fatepotential_DABTRAM_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.DABTRAM) <- 'fatepotential_DABTRAM_d10_w5'
fp.d10_w5.DABTRAM$cell_id <- rownames(fp.d10_w5.DABTRAM)
fp.d10_w5.DABTRAM <- merge(fp.d10_w5.DABTRAM, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.DABTRAM <- fp.d10_w5.DABTRAM %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
            max = max(fatepotential_DABTRAM_d10_w5, na.rm = TRUE),
            n_cells = n())
# write.csv(var.fp.d10_w5.DABTRAM, paste0(out_dir, 'var.fp.d10_w5.DABTRAM.csv'), row.names = F)


# COCL2
fp.d10_w5.COCL2 <- as.data.frame(all_data@misc[["fatepotential_COCL2_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.COCL2) <- 'fatepotential_COCL2_d10_w5'
fp.d10_w5.COCL2$cell_id <- rownames(fp.d10_w5.COCL2)
fp.d10_w5.COCL2 <- merge(fp.d10_w5.COCL2, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.COCL2 <- fp.d10_w5.COCL2 %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            max = max(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            n_cells = n())

range.fp.d10_w5.COCL2 <- fp.d10_w5.COCL2 %>% 
  group_by(assigned_lineage) %>%
  summarise(range = max(fatepotential_COCL2_d10_w5, na.rm = TRUE) - min(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            max = max(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            variance = var(fatepotential_COCL2_d10_w5, na.rm = TRUE),
            n_cells = n())
# write.csv(var.fp.d10_w5.COCL2, paste0(out_dir, 'var.fp.d10_w5.COCL2.csv'), row.names = F)


# CIS
fp.d10_w5.CIS <- as.data.frame(all_data@misc[["fatepotential_CIS_d10_w5"]][['cell_imputed_score']])
colnames(fp.d10_w5.CIS) <- 'fatepotential_CIS_d10_w5'
fp.d10_w5.CIS$cell_id <- rownames(fp.d10_w5.CIS)
fp.d10_w5.CIS <- merge(fp.d10_w5.CIS, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

var.fp.d10_w5.CIS <- fp.d10_w5.CIS %>% 
  group_by(assigned_lineage) %>%
  summarise(variance = var(fatepotential_CIS_d10_w5, na.rm = TRUE),
            n_cells = n())
# write.csv(var.fp.d10_w5.CIS, paste0(out_dir, 'var.fp.d10_w5.CIS.csv'), row.names = F)

# ==============================================================================
# Calculate future lineage size
# ==============================================================================
metadat.week5 <- metadat[metadat$dataset %in% c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS'),]
lin.size.w5 <- table(metadat.week5$assigned_lineage, metadat.week5$dataset)
lin.size.w5 <- as.data.frame(lin.size.w5)
lin.size.w5 <- spread(lin.size.w5, key = 'Var2', value = 'Freq')
colnames(lin.size.w5) <- c('assigned_lineage','week5_CIS.size','week5_COCL2.size', 'week5_DABTRAM.size')

# =============================================================================
# Assemble data
# =============================================================================
var.fp.d10_w5.DABTRAM <- merge(var.fp.d10_w5.DABTRAM, lin.size.w5[, c('assigned_lineage', 'week5_DABTRAM.size')], by = 'assigned_lineage')
var.fp.d10_w5.COCL2 <- merge(var.fp.d10_w5.COCL2, lin.size.w5[, c('assigned_lineage', 'week5_COCL2.size')], by = 'assigned_lineage')
var.fp.d10_w5.CIS <- merge(var.fp.d10_w5.CIS, lin.size.w5[, c('assigned_lineage', 'week5_CIS.size')], by = 'assigned_lineage')

range.fp.d10_w5.COCL2 <- merge(range.fp.d10_w5.COCL2, lin.size.w5[, c('assigned_lineage', 'week5_COCL2.size')], by = 'assigned_lineage')
# =============================================================================
# plot
# =============================================================================

var.fp.d10_w5.DABTRAM <- var.fp.d10_w5.DABTRAM[var.fp.d10_w5.DABTRAM$n_cells > 5, ]
var.fp.d10_w5.DABTRAM <- var.fp.d10_w5.DABTRAM[var.fp.d10_w5.DABTRAM$week5_DABTRAM.size > 0, ]
var.fp.d10_w5.DABTRAM$winner_week5 <- ifelse(var.fp.d10_w5.DABTRAM$week5_DABTRAM.size > 10, 'yes', 'no')
var.fp.d10_w5.DABTRAM$large.max <- ifelse(var.fp.d10_w5.DABTRAM$max > 0.35916223, 'yes', 'no')
p1 <- ggplot(var.fp.d10_w5.DABTRAM, aes(x = sqrt(variance), y = log10(week5_DABTRAM.size+ 1))) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#9D85BE') +
  stat_cor(method = 'spearman') +
  theme_Publication()
p1

ggplot(var.fp.d10_w5.DABTRAM, aes(x = sqrt(variance), y = log10(week5_DABTRAM.size+ 1))) +
  geom_point(aes(size = n_cells, fill = large.max), shape = 21) +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggplot(var.fp.d10_w5.DABTRAM, aes(x = log10(n_cells + 1), y = sqrt(variance))) +
  geom_point(aes(size = log10(week5_DABTRAM.size + 1)), shape = 21) +
  stat_cor(method = 'spearman') +
  ggtitle('DABTRAM Day10') +
  # xlim(c(0, 2.5)) +
  # ylim(c(0, 2)) +
  theme_Publication()

var.fp.d10_w5.COCL2 <- var.fp.d10_w5.COCL2[var.fp.d10_w5.COCL2$n_cells > 2, ]
var.fp.d10_w5.COCL2 <- var.fp.d10_w5.COCL2[var.fp.d10_w5.COCL2$week5_COCL2.size > 0, ]
var.fp.d10_w5.COCL2$winner_week5 <- ifelse(var.fp.d10_w5.COCL2$week5_COCL2.size > 10, 'yes', 'no')
var.fp.d10_w5.COCL2$large.max <- ifelse(var.fp.d10_w5.COCL2$max > -0.3649756, 'yes', 'no')
p2 <- ggplot(var.fp.d10_w5.COCL2, aes(x = sqrt(variance), y = log10(week5_COCL2.size+ 1))) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggplot(var.fp.d10_w5.COCL2, aes(x = sqrt(variance), y = log10(week5_COCL2.size+ 1))) +
  geom_point(aes(size = n_cells, fill = large.max), shape = 21) +
  facet_wrap(~large.max) +
  stat_cor(method = 'spearman') +
  theme_Publication()

range.fp.d10_w5.COCL2 <- range.fp.d10_w5.COCL2[range.fp.d10_w5.COCL2$n_cells > 1, ]
ggplot(range.fp.d10_w5.COCL2, aes(x = range, y = log10(week5_COCL2.size+ 1))) +
  geom_point(aes(size = n_cells), shape = 21) +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggplot(range.fp.d10_w5.COCL2, aes(x = range, y = max)) +
  geom_point(aes(size = n_cells), shape = 21) +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggplot(var.fp.d10_w5.COCL2, aes(y = sqrt(variance), x = log10(n_cells + 1))) +
  geom_point(aes(size = log10(week5_COCL2.size + 1)), shape = 21) +
  stat_cor(method = 'spearman') +
  ggtitle('COCL2 Day10') +
  ylim(c(0, 4)) +
  theme_Publication()

var.fp.d10_w5.CIS <- var.fp.d10_w5.CIS[var.fp.d10_w5.CIS$n_cells > 5, ]
var.fp.d10_w5.CIS <- var.fp.d10_w5.CIS[var.fp.d10_w5.CIS$week5_CIS.size > 0, ]
ggplot(var.fp.d10_w5.CIS, aes(x = sqrt(variance), y = log10(week5_CIS.size+ 1))) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#FBD08C') +
  stat_cor(method = 'spearman') +
  theme_Publication()

ggplot(var.fp.d10_w5.CIS, aes(x = log10(n_cells + 1), y = sqrt(variance))) +
  geom_point(aes(size = log10(week5_CIS.size + 1)), shape = 21) +
  stat_cor(method = 'spearman') +
  ggtitle('CIS Day10') +
  # xlim(c(0, 2.5)) +
  # ylim(c(0, 2)) +
  theme_Publication()

p3 <- ggarrange(p1, p2, p3, ncol = 2, nrow = 1)

p3
