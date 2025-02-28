library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggthemes)

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

treatment <- 'DABTRAM'
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

# ==============================================================================
# Get imputed day 10 status per cell
# ==============================================================================
# imputed progeny size day 0 to day10 
fp_name <- paste0('fatepotential_', treatment, '_d0_d10')
imputedCell.d0_d10 <- as.data.frame(all_data_fatepotential[[fp_name]][["cell_imputed_score"]])
colnames(imputedCell.d0_d10) <- paste0('cell_imputed_score_', treatment, '_d0_d10')
imputedCell.d0_d10$cell_id <- rownames(imputedCell.d0_d10)

imputedCell.d0_d10 <- merge(imputedCell.d0_d10, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')


# imputed lineage size day 10 to week5 
fp_name_2 <- paste0('fatepotential_', treatment, '_d10_w5')
imputedLineage.d10_w5 <- as.data.frame(all_data_fatepotential[[fp_name_2]][["lineage_imputed_count"]])
colnames(imputedLineage.d10_w5) <- paste0('lineage_imputed_count_', treatment, '_d10_w5')
imputedLineage.d10_w5$assigned_lineage <- rownames(imputedLineage.d10_w5)


# combine day0 and week5 through day 10
combine.d0.d10.w5 <- merge(imputedCell.d0_d10, imputedLineage.d10_w5, by = 'assigned_lineage', all.x = T)
# combine.d0.d10.w5 <- combine.d0.d10.w5 %>% drop_na()


fp_median <- median(combine.d0.d10.w5[[paste0('cell_imputed_score_', treatment, '_d0_d10')]])
combine.d0.d10.w5$winner_day0 <- ifelse(combine.d0.d10.w5[[paste0('cell_imputed_score_', treatment, '_d0_d10')]] > fp_median, 'Yes', 'No')
combine.d0.d10.w5$winnerLin_week5 <- ifelse(combine.d0.d10.w5[[paste0('lineage_imputed_count_', treatment, '_d10_w5')]] > 1, 'Yes', 'No')
combine.d0.d10.w5$winnerLin_week5 <- ifelse(is.na(combine.d0.d10.w5$winnerLin_week5), 'No', combine.d0.d10.w5$winnerLin_week5)
combine.d0.d10.w5 %>% 
  filter(winner_day0 == 'Yes',
         winnerLin_week5 == 'Yes') %>% 
  nrow()

combine.d0.d10.w5 %>% 
  filter(winner_day0 == 'Yes',
         winnerLin_week5 != 'Yes') %>% 
  nrow()

combine.d0.d10.w5.plot <- combine.d0.d10.w5[sample(nrow(combine.d0.d10.w5), 100), ]
ggplot(combine.d0.d10.w5) +
  geom_point(aes(x = cell_imputed_score_DABTRAM_d0_d10, 
                 y = log10(lineage_imputed_count_DABTRAM_d10_w5)), 
             shape = 21,
             size = 3, color = 'gray') +
  geom_vline(xintercept = -0, color = 'red', linewidth = 2) +
  geom_hline(yintercept = 0, color = 'red', linewidth = 2) +
  xlim(-1, 1) +
  ylim(-2, 2) +
  labs(x = 'Imputed Cell Fate Potential from Day 0 to Day 10',
       y = 'Imputed Lineage Size at Week5 (log10)') +
  theme_Publication()

ggsave('~/Downloads/filter_scheme_diagram.png', width = 3, height = 3, dpi = 300)




# ==============================================================================
# Calcualte lineage size
# ==============================================================================

metadat <- all_data@meta.data
metadat.day10DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.day10DABTRAM$cell_id <- rownames(metadat.day10DABTRAM)

lineage_size.day10DABTRAM <- metadat.day10DABTRAM %>% 
  group_by(assigned_lineage) %>% 
  summarise(n_cells = n())

# imputed lineage size
imputedLinSize.day10DABTRAM <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["lineage_imputed_count"]])
colnames(imputedLinSize.day10DABTRAM) <- 'lineage_imputed_count'
imputedLinSize.day10DABTRAM$assigned_lineage <- rownames(imputedLinSize.day10DABTRAM)

fp_vs_day10Size <- merge(lineage_size.day10DABTRAM, imputedLinSize.day10DABTRAM, by = 'assigned_lineage')


# imputed cell fate potential
fp.day10DABTRAM <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(fp.day10DABTRAM) <- 'cell_imputed_score'
fp.day10DABTRAM$cell_id <- rownames(fp.day10DABTRAM)
fp.day10DABTRAM <- merge(fp.day10DABTRAM, metadat.day10DABTRAM[, c('assigned_lineage', 'cell_id')], by = 'cell_id')
fp.day10DABTRAM.summary <- fp.day10DABTRAM %>% 
  group_by(assigned_lineage) %>% 
  summarise(mean_fp = mean(cell_imputed_score))

fp_vs_day10Size <- merge(fp_vs_day10Size, fp.day10DABTRAM.summary, by = 'assigned_lineage')

ggplot(fp_vs_day10Size, aes(x = log10(n_cells), y = log10(lineage_imputed_count))) +
  geom_jitter(aes(fill = mean_fp), shape = 21, size = 3, width = 0.1) +
  scale_fill_gradient2(low = 'red', mid = 'gray', high = 'blue') +
  theme_classic()


# imputed lineage size day 0 to day10 DABTRAM
imputedLinSize.day0DABTRAM <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d0_d10"]][["lineage_imputed_count"]])
colnames(imputedLinSize.day0DABTRAM) <- 'lineage_imputed_count_DABTRAM_d0_d10'
imputedLinSize.day0DABTRAM$assigned_lineage <- rownames(imputedLinSize.day0DABTRAM)

fp_vs_day0 <- merge(imputedLinSize.day0DABTRAM, imputedLinSize.day10DABTRAM, by = 'assigned_lineage')

ggplot(fp_vs_day0, aes(x = log10(lineage_imputed_count_DABTRAM_d0_d10), 
                       y = log10(lineage_imputed_count))) +
  geom_point() +
  stat_cor()
