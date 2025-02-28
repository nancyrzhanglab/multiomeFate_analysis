rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

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

# =============================================================================
# FP
# =============================================================================
d10_w5 <- all_data@misc[["fatepotential_DABTRAM_d10_w5"]]
d10_w5 <- as.data.frame(d10_w5[["cell_imputed_score"]])
colnames(d10_w5) <- c('cell_imputed_score_d10_w5')
d10_w5$cell_id <- rownames(d10_w5)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

d10_w5 <- merge(d10_w5, metadat[, c('assigned_lineage', 'cell_id')], by = 'cell_id')

candidates <- seq(min(d10_w5$cell_imputed_score_d10_w5), max(d10_w5$cell_imputed_score_d10_w5), length.out = 100)

thres_cor <- data.frame(thres = numeric(), rho = numeric())
for(thres in candidates) {
  adapting_thres <- thres
  fp_vs_lineage_size <- d10_w5 %>% 
    group_by(assigned_lineage) %>%
    summarize(n_cells = n(), 
              mean_fp = mean(cell_imputed_score_d10_w5),
              num_greater_than_thres = sum(cell_imputed_score_d10_w5 > adapting_thres))
  
  res <- cor.test(fp_vs_lineage_size$num_greater_than_thres, fp_vs_lineage_size$n_cells, method = 'spearman')
  rho <- res$estimate
  
  thres_cor <- rbind(thres_cor, data.frame(thres = thres, rho = rho))
}

ggplot(thres_cor) +
  geom_point(aes(x = thres, y = rho)) +
  ylim(c(0, 1))

adapting_thres <- 0.5
fp_vs_lineage_size <- d10_w5 %>% 
  group_by(assigned_lineage) %>%
  summarize(n_cells = n(), 
            mean_fp = mean(cell_imputed_score_d10_w5),
            num_greater_than_thres = sum(cell_imputed_score_d10_w5 > adapting_thres))

ggplot(fp_vs_lineage_size, aes(x = n_cells, y = num_greater_than_thres)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  geom_smooth(method = 'lm') +
  theme_minimal()
  

# =============================================================================
# Number of adapting cell fate
# =============================================================================
thres_adapt_cell_num <- data.frame()
for(thres in candidates) {
  adapting_thres <- thres
  fp_vs_lineage_size <- d10_w5 %>% 
    group_by(assigned_lineage) %>%
    summarize(n_cells = n(), 
              mean_fp = mean(cell_imputed_score_d10_w5),
              num_greater_than_thres = sum(cell_imputed_score_d10_w5 > adapting_thres))
  fp_vs_lineage_size$adapting_thres <- adapting_thres
  
  # thres_adapt_cell_num <- rbind(thres_adapt_cell_num, 
  #                               data.frame(thres = thres, 
  #                                          num_min = min(fp_vs_lineage_size$num_greater_than_thres),
  #                                          num_max = max(fp_vs_lineage_size$num_greater_than_thres)))
  thres_adapt_cell_num <- rbind(thres_adapt_cell_num,
                                fp_vs_lineage_size)
}

# thres_adapt_cell_num$adapting_thres <- as.factor(thres_adapt_cell_num$adapting_thres)

ggplot(thres_adapt_cell_num, aes(x = adapting_thres, y = num_greater_than_thres)) +
  geom_jitter(width = 0.01, size = 0.5) +
  theme_minimal()
